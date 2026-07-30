#!/usr/bin/env bash
# scan_full_coverage.sh — one command that closes every coverage gap a discovery run leaves,
# so the package can say what WAS and WAS NOT searched with nothing left implicit.
#
#   bash scan_full_coverage.sh <run_dir> [out_dir] [--email you@inst.edu] [--seeds seeds.faa]
#       [--skip-catalogues] [--skip-hosts] [--skip-seeds] [--max-contigs N] [--cpu N]
#
# --seeds is needed only if the run has no seed_qc/seed_status.csv yet (the family census
# writes it). Pass it explicitly: the recorded path in run_manifest.json is redacted by the
# privacy scrubber, so it cannot always be recovered.
#
# A default phage run searches 6 databases, of which only the two genome databases can contain
# an unannotated gene. THREE gaps remain, and this closes all three:
#
#   1. SEED SOURCES — every seed whose own genome was never searched is fetched and scanned, so
#      "the family is accounted for" stops resting on sibling matches.  (scan_missed_seeds.sh)
#   2. METAGENOMES  — GPD (~142k gut phage genomes) and GVD-AVrC (~300k viral representatives)
#      are in the catalog but in no default set, so a homolog living only there is invisible.
#      Streamed, batched, discarded.                                (stream_scan_catalogue.py)
#   3. PROPHAGES    — the host genera's RefSeq representative genomes, where a phage gene would
#      sit as a prophage.                                             (scan_host_genera.sh)
#
# NOT covered even by this, and deliberately so: the FULL RefSeq bacterial genome set (~600 GB
# compressed, ~2 TB six-frame — days-to-weeks, a server job, and the host-genera scan is the
# bounded answer to the same question), and Pfam/PHROGs, which are annotation databases that
# cannot contain an unannotated gene by construction.
#
# Nothing is added to the database cache. The catalogue scans use a TEMPORARY compressed file
# that is deleted when they finish.
#
# Hours to run. Combine with `run.sh --detach`, or drive it via run.sh --scan-full-coverage.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
         "/opt/homebrew/Caskroom/miniforge/base" /opt/conda; do
  [ -f "$b/etc/profile.d/conda.sh" ] && { . "$b/etc/profile.d/conda.sh"; break; }
done
conda activate hmm-discovery 2>/dev/null || true

usage() {
  sed -n '2,31p' "$0" | sed 's/^# \{0,1\}//'
  echo
  echo "Example:"
  echo "  bash scan_full_coverage.sh /path/to/discovery_run --email you@inst.edu --cpu 4"
}
# `--help` used to be taken as the RUN DIRECTORY: mkdir/tee failed on "--help/full_coverage",
# `find --help` printed its usage to STDOUT so the "no profile HMM" guard saw a non-empty
# string, and the script went on to start a 1.5 GB catalogue download with that usage text as
# the model path. Handle the flag, and require a real directory.
case "${1:-}" in
  -h|--help|"") usage; [ -n "${1:-}" ] && exit 0; echo; echo "ERROR: need <run_dir>" >&2; exit 2;;
esac
RUN_DIR="$1"; shift
[ -d "$RUN_DIR" ] || { echo "ERROR: run_dir is not a directory: $RUN_DIR" >&2; exit 2; }
OUT=""; EMAIL=""; CPU=4; MAXC=0; DO_CAT=1; DO_HOST=1; DO_SEED=1; SEEDS=""
while [ $# -gt 0 ]; do
  case "$1" in
    --seeds) SEEDS="${2:?--seeds needs a value}"; shift 2;;
    --email) EMAIL="${2:?--email needs a value}"; shift 2;;
    --cpu) CPU="${2:?--cpu needs a value}"; shift 2;;
    --max-contigs) MAXC="${2:?--max-contigs needs a value}"; shift 2;;
    --skip-catalogues) DO_CAT=0; shift;;
    --skip-hosts) DO_HOST=0; shift;;
    --skip-seeds) DO_SEED=0; shift;;
    -h|--help) usage; exit 0;;
    # An unknown --flag used to fall through and silently become the OUTPUT DIRECTORY.
    -*) echo "ERROR: unknown option: $1" >&2; echo; usage >&2; exit 2;;
    *) [ -z "$OUT" ] && OUT="$1"; shift;;
  esac
done
[ -z "$OUT" ] && OUT="$RUN_DIR/full_coverage"
mkdir -p "$OUT" || { echo "ERROR: cannot create output directory: $OUT" >&2; exit 2; }
LOG="$OUT/full_coverage.log"

HMM="$(find "$RUN_DIR" -path '*03_hmm_profile/profile.hmm' 2>/dev/null | head -1)"
[ -z "$HMM" ] && HMM="$(find "$RUN_DIR" -name 'benchmark_profile.hmm' 2>/dev/null | sort | tail -1)"
[ -z "$HMM" ] && { echo "ERROR: no profile HMM found under $RUN_DIR"; exit 2; }
# Existence is not usability: a 0-byte / truncated profile.hmm (what an interrupted run leaves)
# passed the check above, and every downstream scan then returned 0 hits, which this pipeline
# writes into coverage_summary.csv as "the family is absent from this catalogue".
[ -s "$HMM" ] || { echo "ERROR: profile HMM is EMPTY (0 bytes): $HMM"; exit 2; }
head -c 6 "$HMM" 2>/dev/null | grep -q '^HMMER' || {
  echo "ERROR: not a HMMER profile: $HMM"; exit 2; }

say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
say "FULL COVERAGE SCAN — model: $HMM"
say "  run dir: $RUN_DIR"
say "  output : $OUT"
[ "$MAXC" -gt 0 ] && say "  BOUNDED TEST: --max-contigs $MAXC per catalogue"

# ---- 1. the seeds' own source genomes -------------------------------------------------
if [ "$DO_SEED" = 1 ]; then
  say "[1/3] seed source genomes (every seed whose own genome was never searched)"
  # This step needs seed_qc/seed_status.csv, which the family census writes. Produce it here if
  # it is absent, so a run that predates the census is not silently skipped.
  if [ ! -f "$RUN_DIR/seed_qc/seed_status.csv" ]; then
    if [ -z "$SEEDS" ]; then
      # Best effort from the manifest. Note the run's own privacy scrubber rewrites the home
      # directory and username in shared provenance, so the recorded path may be redacted
      # (".../<user>/...") and unusable — hence --seeds being the reliable route.
      SEEDS="$(python3 - "$RUN_DIR" <<'PY'
import json, os, sys, re
from pathlib import Path
try:
    m = json.loads((Path(sys.argv[1]) / "run_manifest.json").read_text(encoding="utf-8"))
except Exception:
    sys.exit(0)
p = str(((m.get("parameters") or {}).get("fasta")) or m.get("fasta") or "")
if not p:
    sys.exit(0)
p = p.replace("<user>", os.environ.get("USER", "")).replace("<host>", "")
p = os.path.expanduser(p.replace("~", str(Path.home()), 1) if p.startswith("~") else p)
print(p if Path(p).is_file() else "")
PY
)"
      [ -n "$SEEDS" ] && say "  seed FASTA recovered from the manifest: $SEEDS"
    fi
    if [ -n "$SEEDS" ] && [ -f "$SEEDS" ]; then
      say "  seed_status.csv absent — running the family census first"
      python3 "$REPO/scripts/family_census.py" --discovery-dir "$RUN_DIR" \
              --seeds "$SEEDS" --cpu "$CPU" 2>&1 | tee -a "$LOG" | tail -4
    else
      say "  (cannot run the census: pass --seeds <seeds.faa>; the manifest path is redacted"
      say "   by the privacy scrubber, so it cannot always be recovered automatically)"
    fi
  fi
  if [ -n "$EMAIL" ]; then
    bash "$REPO/scripts/scan_missed_seeds.sh" "$RUN_DIR" "$OUT/seed_sources" \
         --email "$EMAIL" 2>&1 | tee -a "$LOG" | tail -6
  else
    say "  (skipped: needs --email for NCBI)"
  fi
fi

# ---- 2. the metagenome catalogues -----------------------------------------------------
if [ "$DO_CAT" = 1 ]; then
  for c in gpd gvd; do
    say "[2/3] metagenome catalogue: $c"
    EXTRA=""; [ "$MAXC" -gt 0 ] && EXTRA="--max-contigs $MAXC"
    python3 "$REPO/scripts/stream_scan_catalogue.py" --hmm "$HMM" --catalogue "$c" \
        --out "$OUT/catalogue_$c" --cpu "$CPU" --resume $EXTRA 2>&1 | tee -a "$LOG" | tail -4
  done
fi

# ---- 3. host genera (prophage) --------------------------------------------------------
if [ "$DO_HOST" = 1 ]; then
  say "[3/3] host genera (prophage copies in the hosts' RefSeq representative genomes)"
  bash "$REPO/scripts/scan_host_genera.sh" "$RUN_DIR" "$OUT/host_genera" 2>&1 \
      | tee -a "$LOG" | tail -5
fi

# ---- consolidated coverage statement --------------------------------------------------
python3 "$REPO/scripts/coverage_report.py" --run-dir "$RUN_DIR" --coverage-dir "$OUT" 2>&1 \
    | tee -a "$LOG"

say "FULL COVERAGE SCAN COMPLETE -> $OUT/coverage_summary.csv"

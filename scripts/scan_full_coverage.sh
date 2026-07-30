#!/usr/bin/env bash
# scan_full_coverage.sh — one command that closes every coverage gap a discovery run leaves,
# so the package can say what WAS and WAS NOT searched with nothing left implicit.
#
#   bash scan_full_coverage.sh <run_dir> [out_dir] [--email you@inst.edu]
#       [--skip-catalogues] [--skip-hosts] [--skip-seeds] [--max-contigs N] [--cpu N]
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

RUN_DIR="${1:?need <run_dir>}"; shift || true
OUT=""; EMAIL=""; CPU=4; MAXC=0; DO_CAT=1; DO_HOST=1; DO_SEED=1
while [ $# -gt 0 ]; do
  case "$1" in
    --email) EMAIL="$2"; shift 2;;
    --cpu) CPU="$2"; shift 2;;
    --max-contigs) MAXC="$2"; shift 2;;
    --skip-catalogues) DO_CAT=0; shift;;
    --skip-hosts) DO_HOST=0; shift;;
    --skip-seeds) DO_SEED=0; shift;;
    *) [ -z "$OUT" ] && OUT="$1"; shift;;
  esac
done
[ -z "$OUT" ] && OUT="$RUN_DIR/full_coverage"
mkdir -p "$OUT"
LOG="$OUT/full_coverage.log"

HMM="$(find "$RUN_DIR" -path '*03_hmm_profile/profile.hmm' 2>/dev/null | head -1)"
[ -z "$HMM" ] && HMM="$(find "$RUN_DIR" -name 'benchmark_profile.hmm' 2>/dev/null | sort | tail -1)"
[ -z "$HMM" ] && { echo "ERROR: no profile HMM found under $RUN_DIR"; exit 2; }

say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
say "FULL COVERAGE SCAN — model: $HMM"
say "  run dir: $RUN_DIR"
say "  output : $OUT"
[ "$MAXC" -gt 0 ] && say "  BOUNDED TEST: --max-contigs $MAXC per catalogue"

# ---- 1. the seeds' own source genomes -------------------------------------------------
if [ "$DO_SEED" = 1 ]; then
  say "[1/3] seed source genomes (every seed whose own genome was never searched)"
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

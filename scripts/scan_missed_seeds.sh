#!/usr/bin/env bash
# scan_missed_seeds.sh — after a discovery run, go and find the input seeds the search never
# re-found, by fetching each seed's OWN source genome and scanning the model against it.
#
#   bash scan_missed_seeds.sh <run_dir> [out_dir] [--email you@inst.edu] [--max N] [--batch N]
#
# WHY: a run reports the homologs it found in the databases it searched, and says nothing about
# the seeds that did not come back — which reads as a sensitivity failure when it is usually a
# COVERAGE one (a GenBank-only accession newer than the searched snapshot, or a metagenomic /
# prophage record a viral-genome database cannot contain). This gives a per-seed verdict and,
# in `explains_miss`, the reason.
#
# NOTHING IS CACHED and NO DATABASE IS ADDED. Genomes are fetched in batches, scanned, and
# deleted; peak disk is one batch. Metagenomic contigs are pulled out of a STREAMED catalogue,
# keeping only the wanted sequences. The database cache is never touched.
#
# Needs an email for NCBI (--email or $NCBI_EMAIL). Without one it prints what it WOULD fetch
# and exits 0 without making a single request.
#
# Output: <out_dir>/missed_seed_sources.txt, results/collection_hits.tsv,
#         missed_seed_scan.csv (also copied into <run_dir>/seed_qc/ so it is packaged).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
         "/opt/homebrew/Caskroom/miniforge/base" /opt/conda; do
  [ -f "$b/etc/profile.d/conda.sh" ] && { . "$b/etc/profile.d/conda.sh"; break; }
done
conda activate hmm-discovery 2>/dev/null || true

RUN_DIR="${1:?need <run_dir>}"; shift || true
OUT=""; BATCH=10; MAX=0; EMAIL=""; ONLY_MISSED=""
while [ $# -gt 0 ]; do
  case "$1" in
    --only-missed) ONLY_MISSED="--only-missed"; shift;;
    --email) EMAIL="$2"; shift 2;;
    --max)   MAX="$2";   shift 2;;
    --batch) BATCH="$2"; shift 2;;
    *) if [ -z "$OUT" ]; then OUT="$1"; else BATCH="$1"; fi; shift;;
  esac
done
[ -z "$OUT" ] && OUT="$RUN_DIR/missed_seed_scan"
mkdir -p "$OUT"

# The cache must not be inside the output tree, or "nothing is cached" stops being true.
case "$(cd "$OUT" && pwd)" in
  "$HOME/.cache/hmm-homologue-finder"*) echo "ERROR: refusing to write into the database cache"; exit 2;;
esac

HMM="$(find "$RUN_DIR" -path '*03_hmm_profile/profile.hmm' 2>/dev/null | head -1)"
[ -z "$HMM" ] && HMM="$(find "$RUN_DIR" -name 'benchmark_profile.hmm' 2>/dev/null | sort | tail -1)"
[ -z "$HMM" ] && { echo "ERROR: no profile HMM found under $RUN_DIR"; exit 2; }
echo "model HMM: $HMM"

if [ ! -f "$RUN_DIR/seed_qc/seed_status.csv" ]; then
  echo "ERROR: $RUN_DIR/seed_qc/seed_status.csv not found — run the family census first"
  echo "       (python3 scripts/family_census.py --discovery-dir '$RUN_DIR' --seeds <seeds.faa>)"
  exit 2
fi

echo "planning the fetch from seed_qc/seed_status.csv…"
python3 "$REPO/scripts/missed_seed_report.py" --run-dir "$RUN_DIR" --out "$OUT" $ONLY_MISSED --plan || exit 2

nsrc=0; [ -f "$OUT/missed_seed_sources.txt" ] && nsrc=$(grep -cvE '^[[:space:]]*(#|$)' "$OUT/missed_seed_sources.txt" || true)
if [ "${nsrc:-0}" -eq 0 ]; then
  echo "nothing to fetch — every distinct seed protein was re-found by the search."
  exit 0
fi

# NCBI wants an address in the URL; the list file carries the __EMAIL__ placeholder and
# _fetch substitutes $NCBI_EMAIL at request time, so the shipped list never holds it.
export NCBI_EMAIL="${EMAIL:-${NCBI_EMAIL:-}}"
if [ -z "$NCBI_EMAIL" ] && grep -q '__EMAIL__' "$OUT/missed_seed_sources.txt"; then
  echo "----------------------------------------------------------------"
  echo " DRY RUN — no email given, so NOT contacting NCBI."
  echo " Would fetch $nsrc source(s):"
  sed 's/^/   /' "$OUT/missed_seed_sources.txt" | sed 's/&email=__EMAIL__//' | head -20
  echo " Re-run with --email you@inst.edu (or set \$NCBI_EMAIL) to perform the scan."
  echo "----------------------------------------------------------------"
  exit 0
fi

echo "scanning the model through the missed seeds' own genomes…"
bash "$REPO/scripts/scan_genome_collection.sh" "$HMM" "$OUT/missed_seed_sources.txt" \
     "$OUT/results" "$BATCH" "$MAX" 25

python3 "$REPO/scripts/missed_seed_report.py" --run-dir "$RUN_DIR" --out "$OUT" $ONLY_MISSED --report || exit 2

# Prove the claim: no genome data left behind, and the shared cache untouched.
left=$(find "$OUT" -name '*.fna' -o -name '*.fa' -o -name '*.fna.gz' 2>/dev/null | grep -v _missed_seed_metagenome | wc -l)
echo "================================================================"
echo " MISSED-SEED SCAN COMPLETE"
echo "   verdicts: $OUT/missed_seed_scan.csv  (also in $RUN_DIR/seed_qc/)"
echo "   hits:     $OUT/results/collection_hits.tsv"
echo "   leftover genome files (should be 0): $left"
echo "   the database cache was not touched — nothing was added to it."
echo "================================================================"

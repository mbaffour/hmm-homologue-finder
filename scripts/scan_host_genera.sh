#!/usr/bin/env bash
# scan_host_genera.sh — after a discovery run, scan its model six-frame through the RefSeq
# representative genomes of the phages' HOST GENERA (where a phage gene would sit as a
# prophage). End-to-end: detect the host genera from the run, assemble the host-genome URL
# list from NCBI assembly_summary, then run the batched six-frame collection scan.
#
#   bash scan_host_genera.sh <run_dir> [out_dir] [batch=25]
#       [--genera "Escherichia,Klebsiella,..."]   # override the auto-detected genera
#       [--max N]                                  # only the first N genomes (testing)
#
# Output: <out_dir>/host_genome_urls.txt, results/collection_hits.tsv, results/collection_hits_aa.faa.
# A model with 0 hits here is phage-specific (no prophage homolog in the host genera).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
         "/opt/homebrew/Caskroom/miniforge/base" /opt/conda; do
  [ -f "$b/etc/profile.d/conda.sh" ] && { . "$b/etc/profile.d/conda.sh"; break; }
done
conda activate hmm-discovery 2>/dev/null || true

RUN_DIR="${1:?need <run_dir>}"; shift || true
OUT=""; BATCH=25; GENERA=""; MAX=0
while [ $# -gt 0 ]; do
  case "$1" in
    --genera) GENERA="$2"; shift 2;;
    --max)    MAX="$2";    shift 2;;
    --batch)  BATCH="$2";  shift 2;;
    *) if [ -z "$OUT" ]; then OUT="$1"; else BATCH="$1"; fi; shift;;
  esac
done
[ -z "$OUT" ] && OUT="$RUN_DIR/host_genera_scan"
mkdir -p "$OUT"

HMM="$(find "$RUN_DIR" -path '*03_hmm_profile/profile.hmm' 2>/dev/null | head -1)"
[ -z "$HMM" ] && HMM="$(find "$RUN_DIR" -name 'benchmark_profile.hmm' 2>/dev/null | sort | tail -1)"
[ -z "$HMM" ] && { echo "ERROR: no profile HMM found under $RUN_DIR"; exit 2; }
echo "model HMM: $HMM"

if [ -n "$GENERA" ]; then
  GRE=$(echo "$GENERA" | tr ',' '|' | tr -d ' ')
else
  GRE=$(python3 - "$RUN_DIR" <<'PY'
import csv, re, sys, glob
run = sys.argv[1]; orgs = set()
for fn in set(glob.glob(run + "/all_runs_hits.csv") + glob.glob(run + "/**/all_runs_hits.csv", recursive=True)):
    try:
        for r in csv.DictReader(open(fn, encoding="utf-8")):
            o = (r.get("organism") or "").strip()
            if o: orgs.add(o)
    except Exception: pass
drop = {"Gamaleyavirus", "MAG", "UNVERIFIED", "Enterobacteria", "Uncultured", "Caudoviricetes"}
gen = set()
for o in orgs:
    w = re.split(r"[ _]", o)[0]; w = {"Echerichia": "Escherichia"}.get(w, w)
    if w[:1].isupper() and w.isalpha() and len(w) > 3 and w not in drop \
       and "phage" not in w.lower() and "virus" not in w.lower():
        gen.add(w)
print("|".join(sorted(gen)))
PY
)
fi
[ -z "$GRE" ] && { echo "ERROR: no host genera detected (pass --genera \"A,B,...\")"; exit 2; }
echo "host genera: $GRE"

echo "fetching NCBI assembly_summary + filtering to reference/representative genomes…"
curl -sS "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" \
| awk -F'\t' -v re="^($GRE)$" '
    /^#/ {next}
    ($5=="reference genome" || $5=="representative genome") {
      split($8,a," ")
      if (a[1] ~ re) { ftp=$20; sub(/\/$/,"",ftp); n=split(ftp,b,"/"); print ftp"/"b[n]"_genomic.fna.gz"; c++ }
    }
    END { printf "  %d host genome(s)\n", c > "/dev/stderr" }' > "$OUT/host_genome_urls.txt"
ng=$(wc -l < "$OUT/host_genome_urls.txt")
echo "  -> $OUT/host_genome_urls.txt  ($ng genomes)"
[ "$ng" -eq 0 ] && { echo "ERROR: no genomes matched the host genera"; exit 2; }

bash "$REPO/scripts/scan_genome_collection.sh" "$HMM" "$OUT/host_genome_urls.txt" "$OUT/results" "$BATCH" "$MAX"

AGG="$OUT/results/collection_hits.tsv"
nh=0; [ -f "$AGG" ] && nh=$(($(wc -l < "$AGG") - 1)); [ "$nh" -lt 0 ] && nh=0
echo "================================================================"
echo " HOST-GENERA SCAN COMPLETE: $nh hit(s) of the model across the host genomes"
echo "   genera: $GRE"
echo "   hits: $AGG ; proteins: $OUT/results/collection_hits_aa.faa"
echo "   (0 hits = the family is phage-specific — no prophage homolog in these hosts.)"
echo "================================================================"

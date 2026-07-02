#!/usr/bin/env bash
# scan_genome_collection.sh — scan a profile HMM through a COLLECTION of genomes,
# six-frame read-through, batch by batch so memory + disk stay bounded, aggregating
# every hit into one table + FASTA. This is how you search a model against a large,
# bounded genome set (e.g. the host genera of your phages, or a custom genome list)
# without loading it all at once — the complement to the streamed catalog databases.
#
# Usage:
#   bash scan_genome_collection.sh <hmm> <list.txt> <out_dir> [batch=25] [max=0]
#
# <list.txt>: one genome source per line — an http(s)/ftp URL to a (gzipped) nucleotide
#   FASTA, or a local .fna[.gz] path. Lines starting with '#' are skipped.
# Reuses scan_genome.py for the per-batch six-frame + read-through (overprinting) scan.
# Output: <out_dir>/collection_hits.tsv, collection_hits_aa.faa, scan_console.log.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
         "/opt/homebrew/Caskroom/miniforge/base" /opt/conda; do
  [ -f "$b/etc/profile.d/conda.sh" ] && { . "$b/etc/profile.d/conda.sh"; break; }
done
conda activate hmm-discovery 2>/dev/null || true

HMM="${1:?need <hmm>}"; LIST="${2:?need <list.txt>}"; OUT="${3:?need <out_dir>}"
BATCH="${4:-25}"; MAX="${5:-0}"
cd "$REPO" || exit 2
mkdir -p "$OUT"
AGG="$OUT/collection_hits.tsv"; AGGAA="$OUT/collection_hits_aa.faa"; LOG="$OUT/scan_console.log"

# one fetch helper: URL -> curl|gunzip ; local .gz -> zcat ; local plain -> cat
_fetch() {
  case "$1" in
    http*|ftp*) curl -sSL --retry 2 "$1" 2>/dev/null | gunzip -c 2>/dev/null ;;
    *.gz)       gunzip -c "$1" 2>/dev/null ;;
    *)          cat "$1" 2>/dev/null ;;
  esac
}

# read the source list into an array (portable to macOS bash 3.2, which lacks mapfile)
SRC=(); while IFS= read -r _line; do [ -n "$_line" ] && SRC+=("$_line"); done < <(grep -vE '^[[:space:]]*(#|$)' "$LIST")
total=${#SRC[@]}; [ "$MAX" -gt 0 ] && [ "$MAX" -lt "$total" ] && total=$MAX
nb=$(( (total + BATCH - 1) / BATCH ))
: > "$AGG"; : > "$AGGAA"; hdr=0; hits=0
echo "$(date '+%F %T') START: $total genome source(s), $nb batch(es) of $BATCH; HMM=$(basename "$HMM")" | tee "$LOG"
i=0; b=0
while [ "$i" -lt "$total" ]; do
  b=$((b+1)); end=$((i+BATCH)); [ "$end" -gt "$total" ] && end=$total
  : > "$OUT/_batch.fna"
  j=$i; while [ "$j" -lt "$end" ]; do _fetch "${SRC[$j]}" >> "$OUT/_batch.fna"; j=$((j+1)); done
  nseq=$(grep -c '^>' "$OUT/_batch.fna" 2>/dev/null); nseq=${nseq:-0}   # grep -c prints 0 itself; no '|| echo 0' (double-counts on empty)
  if [ "${nseq:-0}" -gt 0 ]; then
    python3 scripts/scan_genome.py --hmm "$HMM" --genome "$OUT/_batch.fna" \
        --out "$OUT/_bout" --find-interrupted >/dev/null 2>&1 || true
    if [ -f "$OUT/_bout/scan_hits.tsv" ]; then
      [ "$hdr" = 0 ] && { head -1 "$OUT/_bout/scan_hits.tsv" > "$AGG"; hdr=1; }
      tail -n +2 "$OUT/_bout/scan_hits.tsv" >> "$AGG"
      [ -f "$OUT/_bout/scan_hits_aa.faa" ] && cat "$OUT/_bout/scan_hits_aa.faa" >> "$AGGAA"
    fi
  fi
  hits=$(( $(wc -l < "$AGG" 2>/dev/null || echo 1) - 1 )); [ "$hits" -lt 0 ] && hits=0
  echo "$(date '+%T') batch $b/$nb ($nseq contigs)  cumulative hits: $hits" | tee -a "$LOG"
  rm -f "$OUT/_batch.fna"; rm -rf "$OUT/_bout"
  i=$end
done
echo "$(date '+%F %T') DONE. total hits: $hits  ->  $AGG" | tee -a "$LOG"

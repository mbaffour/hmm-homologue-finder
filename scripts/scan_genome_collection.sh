#!/usr/bin/env bash
# scan_genome_collection.sh — scan a profile HMM through a COLLECTION of genomes,
# six-frame read-through, batch by batch so memory + disk stay bounded, aggregating
# every hit into one table + FASTA. This is how you search a model against a large,
# bounded genome set (e.g. the host genera of your phages, or a custom genome list)
# without loading it all at once — the complement to the streamed catalog databases.
#
# Usage:
#   bash scan_genome_collection.sh <hmm> <list.txt> <out_dir> [batch=25] [max=0] [min_bit=25]
#
# <list.txt>: one genome source per line — an http(s)/ftp URL to a nucleotide FASTA,
#   gzipped OR plain text (an NCBI efetch URL is plain text), or a local .fna[.gz] path.
#   Lines starting with '#' are skipped.
# Reuses scan_genome.py for the per-batch six-frame + read-through (overprinting) scan.
# Output: <out_dir>/collection_hits.tsv, collection_hits_aa.faa, scan_console.log.
#
# NOTHING IS CACHED: each batch is fetched, scanned, and deleted, so peak disk is one
# batch (_batch.fna) no matter how large the collection. The database cache is never
# touched — this is the streaming complement to the cached catalog databases.
#
# If the list contains an address placeholder __EMAIL__, it is substituted from
# $NCBI_EMAIL at request time so the list file itself never carries the address.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
         "/opt/homebrew/Caskroom/miniforge/base" /opt/conda; do
  [ -f "$b/etc/profile.d/conda.sh" ] && { . "$b/etc/profile.d/conda.sh"; break; }
done
conda activate hmm-discovery 2>/dev/null || true

HMM="${1:?need <hmm>}"; LIST="${2:?need <list.txt>}"; OUT="${3:?need <out_dir>}"
BATCH="${4:-25}"; MAX="${5:-0}"; MINBIT="${6:-25}"
cd "$REPO" || exit 2

# REFUSE TO REPORT A NEGATIVE FROM A SCAN THAT NEVER RAN.
# The per-batch scan below is `>/dev/null 2>&1 || true`, so a missing, empty or corrupt HMM
# used to produce "DONE. total hits: 0" and exit 0 — indistinguishable from a real absence,
# and the callers then print "0 hits = the family is phage-specific". Check the model, and
# the genome list, up front instead.
if [ ! -f "$HMM" ]; then
  echo "ERROR: profile HMM not found: $HMM" >&2; exit 2
fi
if [ ! -s "$HMM" ]; then
  echo "ERROR: profile HMM is EMPTY (0 bytes): $HMM" >&2
  echo "       A 0-byte model would silently return 0 hits for every genome." >&2; exit 2
fi
if ! head -c 6 "$HMM" 2>/dev/null | grep -q '^HMMER'; then
  echo "ERROR: not a HMMER profile (no HMMER magic in the first line): $HMM" >&2; exit 2
fi
if [ ! -f "$LIST" ]; then
  echo "ERROR: genome list not found: $LIST" >&2; exit 2
fi
mkdir -p "$OUT"
AGG="$OUT/collection_hits.tsv"; AGGAA="$OUT/collection_hits_aa.faa"; LOG="$OUT/scan_console.log"

# one fetch helper: URL -> curl|decompress ; local .gz -> decompress ; local plain -> cat
#
# `gzip -dcf` not `gunzip -c`: -f passes NON-gzip input through verbatim. With `gunzip -c`
# a plain-text URL (an NCBI efetch response is plain text) produced ZERO BYTES and no
# error, so the batch reported "0 contigs" and the run concluded "gene absent" from a
# fetch that had in fact succeeded. A silent empty fetch must never read as a real negative.
_fetch() {
  local src="$1"
  # substitute the address only at request time; the list file keeps the placeholder
  case "$src" in *__EMAIL__*) src="${src//__EMAIL__/${NCBI_EMAIL:-}}" ;; esac
  case "$src" in
    http*|ftp*) curl -sSL --retry 2 "$src" 2>/dev/null | gzip -dcf 2>/dev/null; sleep 0.4 ;;
    *.gz)       gzip -dcf "$src" 2>/dev/null ;;
    *)          cat "$src" 2>/dev/null ;;
  esac
}

# read the source list into an array (portable to macOS bash 3.2, which lacks mapfile)
SRC=(); while IFS= read -r _line; do [ -n "$_line" ] && SRC+=("$_line"); done < <(grep -vE '^[[:space:]]*(#|$)' "$LIST")
total=${#SRC[@]}; [ "$MAX" -gt 0 ] && [ "$MAX" -lt "$total" ] && total=$MAX
nb=$(( (total + BATCH - 1) / BATCH ))
: > "$AGG"; : > "$AGGAA"; hdr=0; hits=0
echo "$(date '+%F %T') START: $total genome source(s), $nb batch(es) of $BATCH; HMM=$(basename "$HMM")" | tee "$LOG"
i=0; b=0; nscanned=0; nfailed=0; nempty=0
while [ "$i" -lt "$total" ]; do
  b=$((b+1)); end=$((i+BATCH)); [ "$end" -gt "$total" ] && end=$total
  : > "$OUT/_batch.fna"
  j=$i; while [ "$j" -lt "$end" ]; do _fetch "${SRC[$j]}" >> "$OUT/_batch.fna"; j=$((j+1)); done
  nseq=$(grep -c '^>' "$OUT/_batch.fna" 2>/dev/null); nseq=${nseq:-0}   # grep -c prints 0 itself; no '|| echo 0' (double-counts on empty)
  [ "${nseq:-0}" -eq 0 ] && nempty=$((nempty+1))
  if [ "${nseq:-0}" -gt 0 ]; then
    # --no-neighbours: the per-batch Prodigal gene-calling and genome-map rendering were
    # being deleted with _bout a moment later, so they were pure wasted time.
    nscanned=$((nscanned+1))
    # Keep the scan's own diagnostics: a batch that CRASHED must not be indistinguishable
    # from a batch that legitimately found nothing.
    if ! python3 scripts/scan_genome.py --hmm "$HMM" --genome "$OUT/_batch.fna" \
        --out "$OUT/_bout" --find-interrupted --no-neighbours \
        --min-bit "$MINBIT" >"$OUT/_batch_scan.log" 2>&1; then
      nfailed=$((nfailed+1))
      echo "  WARNING: the scan FAILED on batch $b (this batch contributes no evidence either way)" | tee -a "$LOG"
      tail -3 "$OUT/_batch_scan.log" | sed 's/^/    /' | tee -a "$LOG"
    fi
    rm -f "$OUT/_batch_scan.log"
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
echo "  batches with sequence: $nscanned ; batches that fetched NOTHING: $nempty ; scans that FAILED: $nfailed" | tee -a "$LOG"
# A run in which nothing was ever successfully scanned must not be read as "the gene is absent".
if [ "$nscanned" -eq 0 ]; then
  echo "  ERROR: no batch ever contained a sequence — every fetch came back empty." | tee -a "$LOG"
  echo "         'total hits: 0' here is NOT evidence of absence; the search did not happen." | tee -a "$LOG"
  exit 3
fi
if [ "$nfailed" -eq "$nscanned" ]; then
  echo "  ERROR: the scan failed on EVERY batch ($nfailed/$nscanned)." | tee -a "$LOG"
  echo "         'total hits: 0' here is NOT evidence of absence; the search did not happen." | tee -a "$LOG"
  exit 3
fi
[ "$nfailed" -gt 0 ] && echo "  WARNING: $nfailed of $nscanned batch scan(s) failed — this table is INCOMPLETE." | tee -a "$LOG"
exit 0

#!/bin/bash
# start.sh — friendly interactive launcher for the HMM Homologue Finder.
# Walks you through choosing a seed FASTA, databases, and options, then runs the
# pipeline via run.sh (which handles conda activation). Works on macOS, Linux, WSL2.
#
#   bash start.sh
#
set -u -o pipefail   # NOT -e: interactive prompts / probes return non-zero by design
ORIG_PWD="$(pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

bold() { printf '\033[1m%s\033[0m\n' "$1"; }
# ask "prompt" "default"  -> prints the answer on stdout (prompt goes to stderr)
ask() {
  local p="$1" d="${2:-}" a
  if [ -n "$d" ]; then read -r -p "$p [$d]: " a; printf '%s' "${a:-$d}"
  else read -r -p "$p: " a; printf '%s' "$a"; fi
}
# resolve a user-typed path (absolute / ~ / relative-to-original-cwd) to absolute
resolve_path() {
  local p="$1"
  p="${p%\"}"; p="${p#\"}"; p="${p%\'}"; p="${p#\'}"   # strip drag-and-drop quotes
  case "$p" in
    /*) printf '%s' "$p" ;;
    "~"/*|"~") printf '%s' "${p/#\~/$HOME}" ;;
    *) printf '%s/%s' "$ORIG_PWD" "$p" ;;
  esac
}

echo
bold "=== HMM Homologue Finder — guided run ==="
case "$(uname -s)" in
  Darwin*) echo "Platform: macOS" ;;
  Linux*) if grep -qi microsoft /proc/version 2>/dev/null; then echo "Platform: Windows (WSL2)"; else echo "Platform: Linux"; fi ;;
  *) echo "Platform: $(uname -s)" ;;
esac

# default CPU = core count
if command -v nproc >/dev/null 2>&1; then DEF_CPU="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then DEF_CPU="$(sysctl -n hw.ncpu 2>/dev/null || echo 8)"
else DEF_CPU=8; fi

# 1. seed FASTA
echo
echo "Step 1 — seed FASTA (protein, OR nucleotide CDS — auto-detected & translated)"
echo "  Drag the file into this window, or type its path. Enter = bundled example."
FASTA="$(resolve_path "$(ask '  seed FASTA' "$HERE/examples/example_seeds.fasta")")"
if [ ! -f "$FASTA" ]; then echo "  !! File not found: $FASTA"; exit 1; fi
echo "  A nucleotide CDS seed is auto-detected and translated (default genetic code 11,"
echo "  bacterial/phage). Press Enter to accept, or set a different NCBI translation table."
TT="$(ask '  genetic code table for a nucleotide seed' '')"
[ -n "$TT" ] && ARGS_TT=(--trans-table "$TT") || ARGS_TT=()

# 2. mode
echo
echo "Step 2 — run mode"
echo "  1) Smoke test       — fast (1 iteration, 1 small DB); confirms everything works"
echo "  2) Full run         — discovery across your databases, multiple iterations"
echo "  3) Scan ONE genome  — build the HMM from your seeds and check a single genome for the gene"
MODE="$(ask '  choose 1, 2 or 3' '1')"

# Mode 3 — single-genome scan: build the HMM from the Step-1 seeds and scan one genome.
if [ "$MODE" = "3" ]; then
  echo
  echo "Step 3 — the genome to scan"
  echo "  Give a LOCAL genome FASTA (drag it in / type its path), OR an NCBI nucleotide"
  echo "  accession to fetch from NCBI (e.g. KX098390, NC_031062 — comma-separate several)."
  GIN="$(ask '  genome FASTA path OR NCBI accession' '')"
  OUT="$(resolve_path "$(ask '  output directory' "$ORIG_PWD/genome_scan")")"
  FI="$(ask '  also report stop-interrupted / overprinted copies? (y/N)' 'N')"
  SCAN=(--scan --seeds "$FASTA" --out "$OUT")
  [ "${#ARGS_TT[@]}" -gt 0 ] && SCAN+=("${ARGS_TT[@]}")   # nucleotide-seed genetic code, if set
  case "$FI" in y|Y|yes|YES) SCAN+=(--find-interrupted) ;; esac
  GPATH="$(resolve_path "$GIN")"
  if [ -f "$GPATH" ]; then
    SCAN+=(--genome "$GPATH")
  else
    echo "  '$GIN' is not a local file — treating it as an NCBI accession (will fetch it)."
    EMAIL="$(ask '  NCBI email (required to fetch)' '')"
    [ -n "$EMAIL" ] && SCAN+=(--email "$EMAIL")
    SCAN+=(--accession "$GIN")
  fi
  echo
  bold "About to run:"; echo "  bash run.sh ${SCAN[*]}"; echo
  C="$(ask '  press Enter to run (or type n to cancel)' '')"
  case "$C" in n|N|no|NO) echo "Cancelled."; exit 0 ;; esac
  echo
  exec bash run.sh "${SCAN[@]}"
fi

ARGS=(--fasta "$FASTA")
[ "${#ARGS_TT[@]}" -gt 0 ] && ARGS+=("${ARGS_TT[@]}")   # nucleotide-seed genetic code, if set
NAME="$(ask '  output label (folder name)' "$(basename "${FASTA%.*}")")"; ARGS+=(--name "$NAME")
EMAIL="$(ask '  NCBI email for genome/protein fetch (blank = run offline)' '')"
[ -n "$EMAIL" ] && ARGS+=(--email "$EMAIL")   # never send a placeholder address to NCBI

# where the results folder should be created (Enter = next to the seed FASTA)
echo
echo "  Output location — where should the results folder go?"
OUTLOC="$(ask '  output directory (Enter = next to your FASTA)' '')"
if [ -n "$OUTLOC" ]; then
  OUTLOC="$(resolve_path "$OUTLOC")"
  ARGS+=(--out-dir "$OUTLOC/${NAME}_discovery")
  echo "  -> results will go to: $OUTLOC/${NAME}_discovery"
fi

if [ "$MODE" = "2" ]; then
  ARGS+=(--iterations "$(ask '  iterations' '3')")
  ARGS+=(--cpu "$(ask '  CPU threads' "$DEF_CPU")")
  echo
  echo "Step 3 — databases"
  echo "  1) Default set (10 phage/viral DBs)   2) Choose interactively   3) List them first"
  DBMODE="$(ask '  choose 1/2/3' '1')"
  if [ "$DBMODE" = "3" ]; then bash run.sh --list-databases; DBMODE="$(ask '  now 1 (default) or 2 (pick)' '2')"; fi
  [ "$DBMODE" = "2" ] && ARGS+=(--pick-databases)
  echo
  PG="$(ask '  strict mode — also require Prodigal coding-locus overlap? (y/N)' 'N')"
  case "$PG" in y|Y|yes|YES) ARGS+=(--prodigal-gate) ;; esac
  echo
  echo "  Stop-interrupted / overprinted homologs: a read-through scan (stops kept) that"
  echo "  finds homologs broken by a premature stop the normal search can't see. Slower"
  echo "  (re-scans each nucleotide DB). Writes interrupted_homologs.tsv."
  FI="$(ask '  also scan for stop-interrupted / overprinted homologs? (y/N)' 'N')"
  case "$FI" in y|Y|yes|YES) ARGS+=(--find-interrupted) ;; esac

  echo
  echo "Step 4 — figures"
  ST="$(ask '  include the pre-run seed QC tree + alignment? (Y/n)' 'Y')"
  case "$ST" in n|N|no|NO) ARGS+=(--no-seed-tree) ;; esac
  GL="$(ask '  label neighbour genes in the synteny figures? (y/N)' 'N')"
  case "$GL" in y|Y|yes|YES) ARGS+=(--synteny-gene-labels) ;; esac
  echo "  colour synteny genes by:   1) function   2) conservation   3) both"
  CB="$(ask '  choose 1/2/3' '3')"
  case "$CB" in
    1) ARGS+=(--color-by function) ;;
    2) ARGS+=(--color-by conservation) ;;
    *) ARGS+=(--color-by both) ;;
  esac
else
  ARGS+=(--smoke --cpu "$DEF_CPU")
fi

echo
bold "About to run:"
echo "  bash run.sh ${ARGS[*]}"
echo
echo "  Tip: inspect or free cached databases later with:  python3 scripts/manage_cache.py  (add --clear-all to wipe)"
echo
C="$(ask '  press Enter to run (or type n to cancel)' '')"
case "$C" in n|N|no|NO) echo "Cancelled."; exit 0 ;; esac
echo
exec bash run.sh "${ARGS[@]}"

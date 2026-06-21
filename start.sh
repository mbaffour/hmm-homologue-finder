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
echo "Step 1 — seed protein FASTA"
echo "  Drag the file into this window, or type its path. Enter = bundled example."
FASTA="$(resolve_path "$(ask '  seed FASTA' "$HERE/examples/example_seeds.fasta")")"
if [ ! -f "$FASTA" ]; then echo "  !! File not found: $FASTA"; exit 1; fi

# 2. mode
echo
echo "Step 2 — run mode"
echo "  1) Smoke test  — fast (1 iteration, 1 small DB); confirms everything works"
echo "  2) Full run    — your databases, multiple iterations"
MODE="$(ask '  choose 1 or 2' '1')"

ARGS=(--fasta "$FASTA")
NAME="$(ask '  output label (folder name)' "$(basename "${FASTA%.*}")")"; ARGS+=(--name "$NAME")
EMAIL="$(ask '  NCBI email for genome/protein fetch (blank = run offline)' '')"
[ -n "$EMAIL" ] && ARGS+=(--email "$EMAIL")   # never send a placeholder address to NCBI

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
else
  ARGS+=(--smoke --cpu "$DEF_CPU")
fi

echo
bold "About to run:"
echo "  bash run.sh ${ARGS[*]}"
echo
C="$(ask '  press Enter to run (or type n to cancel)' '')"
case "$C" in n|N|no|NO) echo "Cancelled."; exit 0 ;; esac
echo
exec bash run.sh "${ARGS[@]}"

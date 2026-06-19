#!/bin/bash
# Double-click this file in Finder to run the HMM Homologue Finder.
# It verifies the required software (installing anything missing on first run),
# then launches the interactive pipeline, which asks for your seed FASTA.
cd "$(dirname "$0")"
TOOL_DIR="$(cd .. && pwd)"

# activate the conda environment (create/install it if needed)
if ! command -v conda >/dev/null 2>&1; then
  for base in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge"; do
    [ -f "$base/bin/activate" ] && source "$base/bin/activate" && break
  done
fi
conda activate hmm-discovery 2>/dev/null

# verify tools; if any are missing, run the one-time setup
if ! python3 scripts/check_tools.py >/dev/null 2>&1; then
  echo "Some required software is missing — running one-time setup…"
  bash "$TOOL_DIR/setup.sh"
  conda activate hmm-discovery 2>/dev/null
fi

# caffeinate keeps the Mac awake for the (long) run; the pipeline prompts for the FASTA.
caffeinate -i python3 scripts/hmm_finder.py
echo ""
echo "Finished. Press Enter to close this window."
read

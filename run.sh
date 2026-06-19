#!/bin/bash
# HMM Homologue Finder — terminal launcher for macOS and Linux (incl. WSL2).
# Verifies the software (installs on first run), then starts the interactive
# pipeline, which asks for your seed FASTA.
set -u
cd "$(dirname "$0")"

# Find conda (or tell the user to install Miniforge).
if ! command -v conda >/dev/null 2>&1; then
  for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
    [ -f "$base/bin/activate" ] && . "$base/bin/activate" && break
  done
fi
conda activate hmm-discovery 2>/dev/null || true

# Install/verify tools on first run.
if ! python3 scripts/check_tools.py >/dev/null 2>&1; then
  echo "Setting up the software (one-time)…"
  bash setup.sh || { echo "Setup failed — see messages above."; exit 1; }
  conda activate hmm-discovery 2>/dev/null || true
fi

# Keep the machine awake during the (long) run if a tool is available.
KEEP_AWAKE=""
command -v caffeinate >/dev/null 2>&1 && KEEP_AWAKE="caffeinate -i"          # macOS
command -v systemd-inhibit >/dev/null 2>&1 && KEEP_AWAKE="systemd-inhibit"   # Linux

$KEEP_AWAKE python3 scripts/hmm_finder.py "$@"

#!/bin/bash
# HMM Homologue Finder — terminal launcher for macOS and Linux (incl. WSL2).
# Verifies the software (installs on first run), then starts the interactive
# pipeline, which asks for your seed FASTA.
set -u
cd "$(dirname "$0")"

ENV_NAME="hmm-discovery"

# Load a conda *shell hook* so `conda activate` works in a non-interactive shell.
# NOTE: we deliberately source etc/profile.d/conda.sh and NOT $base/bin/activate —
# sourcing bin/activate inherits THIS launcher's positional parameters ("$@") and
# feeds them to `conda activate`, which then errors with
# "activate does not accept more than one argument".
load_conda() {
  local base
  for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
              "$HOME/opt/miniconda3" "$HOME/opt/anaconda3" \
              "/opt/homebrew/Caskroom/miniforge/base" "/usr/local/Caskroom/miniforge/base" \
              "/opt/miniconda3" "/opt/anaconda3" "/opt/conda"; do
    if [ -f "$base/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1091
      . "$base/etc/profile.d/conda.sh"
      return 0
    fi
  done
  # conda already on PATH (e.g. via condabin) but hook not yet sourced: derive base.
  if command -v conda >/dev/null 2>&1; then
    local b
    b="$(conda info --base 2>/dev/null || true)"
    if [ -n "$b" ] && [ -f "$b/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1091
      . "$b/etc/profile.d/conda.sh"
      return 0
    fi
  fi
  return 1
}

if ! load_conda; then
  echo "conda/mamba was not found. Install Miniforge first, then re-run ./run.sh"
  echo "  (see README / docs/INSTALL.md)"
  exit 1
fi

# Activate the env; create it on first run if activation fails.
if ! conda activate "$ENV_NAME" 2>/dev/null; then
  echo "Setting up the software (one-time)…"
  bash setup.sh || { echo "Setup failed — see messages above."; exit 1; }
  conda activate "$ENV_NAME" || { echo "Could not activate '$ENV_NAME' after setup."; exit 1; }
fi

# Verify required tools; install on first run if anything is missing.
if ! python3 scripts/check_tools.py >/dev/null 2>&1; then
  echo "Setting up the software (one-time)…"
  bash setup.sh || { echo "Setup failed — see messages above."; exit 1; }
  conda activate "$ENV_NAME" || { echo "Could not activate '$ENV_NAME' after setup."; exit 1; }
fi

# Keep the machine awake during the (long) run if a working tool is available.
# systemd-inhibit is probed because it exists but is denied under WSL2.
KEEP_AWAKE=""
if command -v caffeinate >/dev/null 2>&1; then
  KEEP_AWAKE="caffeinate -i"                                   # macOS
elif command -v systemd-inhibit >/dev/null 2>&1 && systemd-inhibit --what=idle true >/dev/null 2>&1; then
  KEEP_AWAKE="systemd-inhibit --what=idle"                     # Linux (skip if denied, e.g. WSL)
fi

# exec so the pipeline's exit code becomes this script's exit code (no masking).
exec $KEEP_AWAKE python3 scripts/hmm_finder.py "$@"

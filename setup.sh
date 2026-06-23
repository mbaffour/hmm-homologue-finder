#!/bin/bash
# setup.sh — one-time environment setup for the HMM Homologue Finder.
# Checks for conda, creates the `hmm-discovery` environment with all required
# tools, and verifies everything is installed. Safe to re-run.
#
#   bash setup.sh            # create/verify the environment
#   bash setup.sh --check    # dry run: print the platform plan, install nothing
#
# NOTE: we use `set -u -o pipefail` but deliberately NOT `set -e` — this script
# sources the conda shell hook and runs `command -v` / `grep -q` probes that
# legitimately exit non-zero; `-e` would abort on those. Errors are handled
# explicitly below.
set -u -o pipefail
TOOL_DIR="$(cd "$(dirname "$0")" && pwd)"
# Use the engine bundled with the tool; fall back to the dev repo if absent.
if [ -f "$TOOL_DIR/engine/environment.yml" ]; then
  DEPLOY="$TOOL_DIR/engine"
else
  DEPLOY="$HOME/Documents/HMM-Discovery-Deployable-20260602"
fi
ENV_NAME="hmm-discovery"

CHECK=0
case "${1:-}" in --check|--dry-run|-n) CHECK=1 ;; esac

# OS check: bioconda tools are Linux/macOS only. On Windows this must run in WSL2.
case "$(uname -s)" in
  Linux*|Darwin*) : ;;
  MINGW*|MSYS*|CYGWIN*)
    echo "Native Windows is not supported (the bioinformatics tools have no Windows builds)."
    echo "Please run this inside WSL2 (Ubuntu). See README / run.bat."; exit 1 ;;
esac

# macOS Apple Silicon (M1/M2/M3): several bioconda tools (HMMER, MAFFT, Prodigal,
# IQ-TREE, MEME, …) lack reliable osx-arm64 builds. Without this, conda either
# fails to solve or silently mixes channels. Force the x86-64 subdir so the env
# resolves consistently (runs under Rosetta 2). Respect a user-set CONDA_SUBDIR
# and only do this on arm64 Darwin; Linux and Intel macOS are untouched.
if [ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ] && [ -z "${CONDA_SUBDIR:-}" ]; then
  export CONDA_SUBDIR=osx-64
  echo "(Apple Silicon detected: setting CONDA_SUBDIR=osx-64 so bioconda tools resolve "
  echo " via x86-64/Rosetta 2 — override by exporting CONDA_SUBDIR before running.)"
fi

echo "=== HMM Homologue Finder — environment setup ($(uname -s) $(uname -m)) ==="

# --check / dry run: report the platform plan and exit before touching anything.
if [ "$CHECK" = "1" ]; then
  echo ""
  echo "--- dry run (--check): nothing will be created or installed ---"
  echo "  Platform:      $(uname -s) $(uname -m)"
  echo "  CONDA_SUBDIR:  ${CONDA_SUBDIR:-(platform default)}"
  echo "  Target env:    $ENV_NAME"
  if [ -f "$DEPLOY/environment.yml" ]; then
    echo "  Env recipe:    $DEPLOY/environment.yml"
  else
    echo "  Env recipe:    (none found at $DEPLOY — would use a minimal inline recipe)"
  fi
  if command -v conda >/dev/null 2>&1; then
    echo "  conda:         $(command -v conda)"
    if conda env list 2>/dev/null | grep -qE "^[[:space:]]*${ENV_NAME}[[:space:]]"; then
      echo "  env status:    exists — setup would verify it and install any missing tools"
    else
      echo "  env status:    missing — setup would create it from the recipe"
    fi
  else
    echo "  conda:         NOT FOUND — setup would print Miniforge install instructions"
  fi
  echo "Plan looks right? Re-run without --check to apply."
  exit 0
fi

# 1. conda / mamba present?  Load a conda *shell hook* (etc/profile.d/conda.sh)
# rather than bin/activate, so that `conda activate` works later even when conda
# is only reachable via condabin (which otherwise errors "Run 'conda init' first").
for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
            "$HOME/opt/miniconda3" "$HOME/opt/anaconda3" \
            "/opt/homebrew/Caskroom/miniforge/base" "/usr/local/Caskroom/miniforge/base" \
            "/opt/miniconda3" "/opt/anaconda3" "/opt/conda"; do
  if [ -f "$base/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "$base/etc/profile.d/conda.sh"
    break
  fi
done
# Fall back: conda on PATH but hook not found above — derive its base and source it.
if ! type conda 2>/dev/null | grep -q 'function'; then
  if command -v conda >/dev/null 2>&1; then
    _cbase="$(conda info --base 2>/dev/null || true)"
    [ -n "$_cbase" ] && [ -f "$_cbase/etc/profile.d/conda.sh" ] && source "$_cbase/etc/profile.d/conda.sh"
  fi
fi
if ! command -v conda >/dev/null 2>&1; then
  case "$(uname -s)" in
    Darwin*) _mf="Miniforge3-MacOSX-$(uname -m).sh" ;;
    *)       _mf="Miniforge3-Linux-$(uname -m).sh" ;;
  esac
  echo "conda/mamba was not found. Install Miniforge first (one time):"
  echo
  echo "  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/$_mf"
  echo "  bash $_mf"
  echo "  # then close & reopen the terminal and re-run:  bash setup.sh"
  echo
  exit 1
fi
echo "conda found: $(command -v conda)"

# 2. create the environment if missing (reuse the deployable repo's recipe)
if ! conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
  echo "Creating the '${ENV_NAME}' environment (this can take several minutes)…"
  if [ -f "$DEPLOY/environment.yml" ]; then
    conda env create -f "$DEPLOY/environment.yml" || mamba env create -f "$DEPLOY/environment.yml"
  else
    echo "environment.yml not found at $DEPLOY; creating a minimal env."
    conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
      python=3.11 hmmer mafft trimal prodigal seqkit cd-hit iqtree meme curl
  fi
fi

# 3. activate and verify; install any stragglers
conda activate "$ENV_NAME"
if [ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]; then
  echo "ERROR: could not activate the '$ENV_NAME' environment (CONDA_DEFAULT_ENV='${CONDA_DEFAULT_ENV:-}')." >&2
  echo "Try: conda activate $ENV_NAME   — or re-run after 'conda init bash' && restarting the shell." >&2
  exit 1
fi
echo "Active env: ${CONDA_DEFAULT_ENV:-?}"

need=""
for t in hmmsearch hmmbuild mafft trimal prodigal seqkit cd-hit iqtree meme fimo curl; do
  command -v "$t" >/dev/null 2>&1 || need="$need $t"
done
if [ -n "$need" ]; then
  echo "Installing missing tools:$need"
  conda install -y -c conda-forge -c bioconda $need || mamba install -y -c conda-forge -c bioconda $need
fi
# clinker (+ biopython/pandas) via pip if absent
command -v clinker >/dev/null 2>&1 || pip install clinker
python -c "import Bio, pandas" 2>/dev/null || pip install biopython pandas

# Optional: headless browser for STATIC clinker figures (PNG export of clinker's
# JS-rendered plot). Entirely optional and non-fatal — the pipeline runs fine
# without it (static SVG/PNG/PDF synteny panels are produced regardless).
if python -c "import playwright" 2>/dev/null; then
  echo "Setting up the headless browser for static clinker figures (optional)…"
  python -m playwright install chromium >/dev/null 2>&1 || echo "  (chromium download skipped)"
  # chromium needs some OS libraries, installed via apt/dnf (root). Try unattended;
  # if sudo needs a password, print the one-time command instead of blocking setup.
  if sudo -n "$(command -v python)" -m playwright install-deps chromium >/dev/null 2>&1; then
    echo "  static clinker enabled (browser + OS libs installed)."
  else
    echo "  Static clinker PNGs need chromium OS libraries (one-time, needs your password):"
    echo "      sudo \$(command -v python) -m playwright install-deps chromium"
    echo "  (Skipping for now — the static synteny panels in downstream/synteny/ are produced regardless.)"
  fi
fi

# 4. final report
echo ""
python3 "$(dirname "$0")/scripts/check_tools.py"
echo ""
echo "Setup complete. To run the pipeline:"
echo "  conda activate ${ENV_NAME}"
echo "  python3 scripts/hmm_finder.py        # prompts for your seed FASTA"

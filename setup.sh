#!/bin/bash
# setup.sh — one-time environment setup for the HMM Homologue Finder.
# Checks for conda, creates the `hmm-discovery` environment with all required
# tools, and verifies everything is installed. Safe to re-run.
#
#   bash setup.sh
#
set -u
TOOL_DIR="$(cd "$(dirname "$0")" && pwd)"
# Use the engine bundled with the tool; fall back to the dev repo if absent.
if [ -f "$TOOL_DIR/engine/environment.yml" ]; then
  DEPLOY="$TOOL_DIR/engine"
else
  DEPLOY="$HOME/Documents/HMM-Discovery-Deployable-20260602"
fi
ENV_NAME="hmm-discovery"

# OS check: bioconda tools are Linux/macOS only. On Windows this must run in WSL2.
case "$(uname -s)" in
  Linux*|Darwin*) : ;;
  MINGW*|MSYS*|CYGWIN*)
    echo "Native Windows is not supported (the bioinformatics tools have no Windows builds)."
    echo "Please run this inside WSL2 (Ubuntu). See README / run.bat."; exit 1 ;;
esac

echo "=== HMM Homologue Finder — environment setup ($(uname -s)) ==="

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

# 4. final report
echo ""
python3 "$(dirname "$0")/scripts/check_tools.py"
echo ""
echo "Setup complete. To run the pipeline:"
echo "  conda activate ${ENV_NAME}"
echo "  python3 scripts/hmm_finder.py        # prompts for your seed FASTA"

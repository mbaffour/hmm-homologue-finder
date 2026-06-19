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

# 1. conda / mamba present?
if ! command -v conda >/dev/null 2>&1; then
  # try common install locations
  for base in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge"; do
    [ -f "$base/bin/activate" ] && source "$base/bin/activate" && break
  done
fi
if ! command -v conda >/dev/null 2>&1; then
  cat <<'MSG'
conda/mamba was not found. Install Miniforge first (one time):

  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh
  bash Miniforge3-MacOSX-$(uname -m).sh
  # then close & reopen the terminal and re-run:  bash setup.sh

MSG
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
source activate "$ENV_NAME" 2>/dev/null || conda activate "$ENV_NAME"
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

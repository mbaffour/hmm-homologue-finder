#!/bin/bash
# HMM Homologue Finder — terminal launcher for macOS and Linux (incl. WSL2).
# Verifies the software (installs on first run), then starts the interactive
# pipeline, which asks for your seed FASTA.
set -u -o pipefail   # NOT -e: conda hook + `command -v` probes return non-zero by design
cd "$(dirname "$0")"
HERE="$(pwd)"

# --detach: re-launch this run in its OWN session (setsid) so closing the terminal
# window — which sends SIGHUP to the foreground job and would otherwise kill a long
# run — cannot stop it. The detached copy does the normal setup + search; its console
# output is redirected to a log file (next to your seed, else your home folder).
# --detach is stripped before relaunch so the run happens exactly once. Detached runs
# have no terminal, so they never prompt: pass --email for online lookups (else they
# run offline) and --databases to choose DBs (else the defaults are used).
case " $* " in
  *" --detach "*)
    detach_args=()
    for a in "$@"; do [ "$a" = "--detach" ] || detach_args+=("$a"); done
    seed=""; want=0
    for a in "$@"; do
      if [ "$want" = 1 ]; then seed="$a"; want=0; fi
      [ "$a" = "--fasta" ] && want=1
    done
    case "$seed" in [A-Za-z]:[\\/]*) seed="$(wslpath -u "$seed" 2>/dev/null || printf '%s' "$seed")";; esac
    logdir="$HOME"
    if [ -n "$seed" ]; then sd="$(dirname "$seed" 2>/dev/null || true)"; [ -n "$sd" ] && [ -d "$sd" ] && logdir="$sd"; fi
    log="$logdir/hmm_run_console_$(date +%Y%m%d_%H%M%S).log"
    setsid nohup bash "$HERE/run.sh" ${detach_args[@]+"${detach_args[@]}"} </dev/null >"$log" 2>&1 &
    bgpid=$!
    printf '\n============================================================\n'
    printf ' Your search is now running IN THE BACKGROUND (pid %s).\n' "$bgpid"
    printf ' You can safely CLOSE THIS WINDOW — the run keeps going.\n\n'
    printf '   progress log : %s\n' "$log"
    printf '   when it finishes, open  report.html  in your results folder.\n'
    printf '   (double-click CHECK_RUN.bat any time to see how it is going)\n'
    printf '============================================================\n\n'
    exit 0
    ;;
esac

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

# exec so the tool's exit code becomes this script's exit code (no masking).
# `run.sh --scan ...` routes to the single-genome scan mode; otherwise the
# discovery pipeline. Both share the conda activation + tool check above.
if [ "${1:-}" = "--scan" ]; then
  shift
  exec $KEEP_AWAKE python3 scripts/scan_genome.py "$@"
fi
# `run.sh --scan-host-genera <run_dir> …` searches the run's model six-frame through the
# RefSeq representative genomes of the phages' host genera (prophage search). Combine with
# the top-of-file --detach to run it in the background.
if [ "${1:-}" = "--scan-host-genera" ]; then
  shift
  exec $KEEP_AWAKE bash scripts/scan_host_genera.sh "$@"
fi
exec $KEEP_AWAKE python3 scripts/hmm_finder.py "$@"

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
    # Do not announce a background run without checking one started. A bad flag or an
    # unreadable seed file kills the detached copy — but only AFTER conda activation and the
    # tool check, i.e. ~20 s in — and printing "your search is now running, you can close this
    # window" for a run that is already dead is the worst possible message: the user waits for
    # a result that will never come. Watch it until it is past setup (or gone).
    printf '\nStarting the background run — checking that it comes up'
    _alive=0
    for _i in $(seq 1 45); do
      if ! kill -0 "$bgpid" 2>/dev/null; then _alive=0; break; fi
      _alive=1
      # hmm_finder prints this banner once every argument has been validated and the real
      # work has begun — past every fast-failure path, so seeing it means the run is real.
      if grep -q '=== family pipeline:' "$log" 2>/dev/null; then break; fi
      printf '.'
      sleep 1
    done
    printf '\n'
    if [ "$_alive" != 1 ]; then
      printf '\n============================================================\n'
      printf ' THE BACKGROUND RUN IS NOT RUNNING. Nothing was started.\n\n'
      printf ' Last lines of %s:\n\n' "$log"
      tail -12 "$log" 2>/dev/null | sed 's/^/   /'
      printf '\n Fix the problem above and re-run.\n'
      printf '============================================================\n\n'
      exit 1
    fi
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
# `run.sh --scan-missed-seeds <run_dir> …` chases the input seeds the search never re-found by
# fetching each seed's OWN source genome and scanning the model against it — answering whether
# the miss was database coverage or a genuine absence. Streams: nothing is cached, no database
# is added. Combine with the top-of-file --detach to run it in the background.
if [ "${1:-}" = "--scan-missed-seeds" ]; then
  shift
  exec $KEEP_AWAKE bash scripts/scan_missed_seeds.sh "$@"
fi
# `run.sh --scan-catalogue gpd|gvd --hmm H --out D` streams one large metagenome catalogue
# (GPD / GVD-AVrC) through the model, batch by batch, discarding as it goes.
if [ "${1:-}" = "--scan-catalogue" ]; then
  shift
  exec $KEEP_AWAKE python3 scripts/stream_scan_catalogue.py "$@"
fi
# `run.sh --scan-full-coverage <run_dir> --email E` closes EVERY coverage gap in one command:
# the seeds' own source genomes, the GPD + GVD metagenome catalogues, and the host genera
# (prophage) — then writes a coverage_summary.csv that also NAMES what is still not searched.
# Hours to run; combine with the top-of-file --detach.
if [ "${1:-}" = "--scan-full-coverage" ]; then
  shift
  exec $KEEP_AWAKE bash scripts/scan_full_coverage.sh "$@"
fi
exec $KEEP_AWAKE python3 scripts/hmm_finder.py "$@"

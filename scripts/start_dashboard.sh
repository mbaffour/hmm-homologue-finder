#!/usr/bin/env bash
# start_dashboard.sh — start the run dashboard so it STAYS up, and print the URL.
#
#   bash scripts/start_dashboard.sh [port] [--stop]
#
# WHY A SCRIPT RATHER THAN "just run dashboard.py &": under WSL a background job dies when the
# shell that launched it goes away, and WSL tears the whole VM down once its last process
# exits. A dashboard started with a bare `&` from a one-shot `wsl bash …` command is therefore
# gone by the time you open a browser — which looks exactly like "the page can't be reached".
# This detaches with setsid and holds a keepalive so the VM does not shut underneath it.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
PORT="${1:-8765}"
case "${1:-}" in --stop|stop) PORT=""; ;; esac
cd "$REPO" || exit 2

if [ "${1:-}" = "--stop" ] || [ "${1:-}" = "stop" ]; then
  pkill -f 'scripts/dashboard.py' 2>/dev/null && echo "dashboard stopped" || echo "not running"
  exit 0
fi

for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
         "/opt/homebrew/Caskroom/miniforge/base" /opt/conda; do
  [ -f "$b/etc/profile.d/conda.sh" ] && { . "$b/etc/profile.d/conda.sh"; break; }
done
conda activate hmm-discovery 2>/dev/null || true
PY="$(command -v python3)"

# already up? then just report it rather than starting a second one on a busy port
if pgrep -f 'scripts/dashboard.py' >/dev/null 2>&1; then
  echo "dashboard already running:"
  pgrep -fa 'scripts/dashboard.py' | sed 's/^/  /'
  echo "  ->  http://127.0.0.1:${PORT}"
  exit 0
fi

# Keepalive: WSL shuts the VM down when nothing is left running, which would take the
# dashboard with it. Only start one if none is already holding the door open.
pgrep -f 'sleep 604800' >/dev/null 2>&1 || \
  setsid nohup sleep 604800 </dev/null >/dev/null 2>&1 &

LOG="$HOME/hmm_dashboard.log"
WINUSER="$(ls /mnt/c/Users 2>/dev/null | grep -viE '^(public|default|all users|defaultuser0)$' | head -1)"
export WINUSER
setsid nohup "$PY" "$REPO/scripts/dashboard.py" --port "$PORT" --interval 12 \
    --root "$HOME/hmm_runs" \
    ${WINUSER:+--root "/mnt/c/Users/$WINUSER/Downloads"} \
    </dev/null >"$LOG" 2>&1 &
sleep 3

if pgrep -f 'scripts/dashboard.py' >/dev/null 2>&1; then
  echo "================================================================"
  echo " DASHBOARD RUNNING — open this in your browser:"
  echo ""
  echo "     http://127.0.0.1:${PORT}"
  echo ""
  echo "   It keeps running after you close this window."
  echo "   First page load takes a few seconds while it scans the run folders."
  echo "   log:  $LOG"
  echo "   stop: bash scripts/start_dashboard.sh --stop"
  echo "================================================================"
else
  echo "FAILED to start. Last lines of $LOG:"; tail -15 "$LOG" 2>/dev/null | sed 's/^/  /'
  exit 1
fi

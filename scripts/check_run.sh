#!/usr/bin/env bash
# Show whether a background (detached) search is still running and its latest progress.
# Invoked by final/CHECK_RUN.bat; safe to run any time.
echo "=== Is a search running right now? ==="
if ps -eo etime,args 2>/dev/null | grep -E '[h]mm_finder\.py|[s]can_genome\.py|[s]can_genome_collection\.sh|[s]can_host_genera\.sh|[r]un_pipeline' | grep -v grep; then
  echo "  ^ YES — a search is running (elapsed time shown on the left)."
else
  echo "  No search process is currently running (it may be finished, or not started)."
fi
echo
echo "=== Latest background progress log ==="
log=$(ls -t /mnt/c/Users/*/Downloads/hmm_run_console_*.log "$HOME"/hmm_run_console_*.log 2>/dev/null | head -1)
if [ -n "$log" ]; then
  echo "  $log"
  echo "  --- last 18 lines ---"
  tail -18 "$log" 2>/dev/null | sed 's/^/    /'
else
  echo "  (no detached console log found yet)"
fi
echo
echo "=== Most recent finished report(s) ==="
found=$(ls -t /mnt/c/Users/*/Downloads/*discovery*/report.html "$HOME"/hmm_runs/*/report.html 2>/dev/null | head -3)
if [ -n "$found" ]; then
  printf '%s\n' "$found" | sed 's/^/  /'
  echo "  (open one of these in your browser for the results)"
else
  echo "  (no finished report.html yet)"
fi

#!/usr/bin/env bash
# Pauses the pipeline's heavy processes (SIGSTOP) when swap is nearly exhausted and resumes them (SIGCONT) when it recovers,
# so that the laptop does not freeze. Only touches miniprot, minimap2 and diamond started from this folder.
set -uo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
log="$root/logs/watchdog.log"
paused=0
while true; do
  free=$(sysctl -n vm.swapusage | sed -E 's/.*free = ([0-9]+)\.[0-9]+M.*/\1/')
  pids=$(pgrep -f "$root" | xargs -I{} sh -c 'ps -o pid=,comm= -p {} | grep -E "miniprot|minimap2|diamond" | awk "{print \$1}"' 2>/dev/null | tr '\n' ' ')
  if (( paused == 0 )) && (( free < 150 )) && [[ -n "$pids" ]]; then
    kill -STOP $pids 2>/dev/null && paused=1 && echo "$(date '+%F %T') PAUSED swap_free=${free}M pids=$pids" >> "$log"
  elif (( paused == 1 )) && (( free > 500 )); then
    kill -CONT $pids 2>/dev/null; paused=0; echo "$(date '+%F %T') RESUMED swap_free=${free}M pids=$pids" >> "$log"
  fi
  # exit when the run it guards is finished: WATCHDOG_DONE_FILE / WATCHDOG_DONE_PAT / WATCHDOG_DONE_N override the pilot defaults
  n_done=$(grep -cE "${WATCHDOG_DONE_PAT:-ALL_DONE}" "${WATCHDOG_DONE_FILE:-$root/logs/driver.log}" 2>/dev/null || true)
  (( ${n_done:-0} >= ${WATCHDOG_DONE_N:-1} )) && (( paused == 0 )) && exit 0
  sleep 5
done

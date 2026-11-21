#!/usr/bin/env bash
# One-off for the night of 2026-11-12/13 (laptop is shut down at 08:00).
# At CUT_AT: stop the pair3 HiFi download (pstream) so that hifi_stream.sh finishes with what was downloaded
# and the driver runs HIFI_CHECK pair3. Then stop the HiFi queue so pair1/pair2 do not start before shutdown
# (they are rerun in the evening), and write a state snapshot for the restart.
set -euo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
log="$root/logs/verify.log"
mylog="$root/logs/cut_hifi.log"
say() { echo "$(date '+%F %T') $*" | tee -a "$mylog" >> "$log"; }
CUT_AT="${CUT_AT:-2026-11-13 06:40}"
DRIVER_PID="${DRIVER_PID:-17853}"

target=$(date -j -f "%Y-%m-%d %H:%M" "$CUT_AT" +%s)
now=$(date +%s)
echo "$(date '+%F %T') armed: cut pair3 HiFi at $CUT_AT, driver pid $DRIVER_PID" >> "$mylog"
(( target > now )) && sleep $(( target - now ))

if pgrep -f "pstream.py.*SRR39556278" >/dev/null; then
  got=$(tail -n 1 "$root/verify/pair3/hifi/SRR39556278.pstream.log" 2>/dev/null || true)
  say "HIFI_CUT pair3 by plan (shutdown at 08:00): $got"
  pkill -f "pstream.py.*SRR39556278" || true
else
  say "HIFI_CUT pair3: download already finished, nothing to cut"
fi

# wait for the driver to finish HIFI_CHECK pair3 (minimap2 drains its last batch first, up to ~15 min)
for _ in $(seq 1 180); do
  grep -qE "HIFI_CHECK_(OK|FAILED) pair3" "$log" && break
  sleep 30
done
grep -qE "HIFI_CHECK_(OK|FAILED) pair3" "$log" || say "HIFI_CHECK pair3 did not finish within 90 min; stopping queue anyway"

# stop the queue: driver first, then anything already started for pair1
kill "$DRIVER_PID" 2>/dev/null || true
pkill -f "hifi_stream.sh .*verify/pair1" 2>/dev/null || true
pkill -f "pstream.py.*SRR38296513" 2>/dev/null || true
pkill -f "minimap2 .*verify/pair1/hifi/ref_both.fa" 2>/dev/null || true
sleep 5
rm -f "$root/verify/pair1/hifi/SRR38296513.paf.part"
say "HIFI_QUEUE_STOPPED pair1 pair2 deferred to the evening run (full volume, faster stream)"
say "VERIFY_DONE partial=pair3_only"

# state snapshot for the evening restart
snap="$root/logs/STATE_before_shutdown_$(date '+%Y-%m-%d_%H%M').txt"
{
  echo "state at $(date '+%F %T'); laptop to be shut down at 08:00"
  echo; echo "== verify.log (tail)"; tail -n 25 "$log"
  echo; echo "== verify2.log (tail)"; tail -n 15 "$root/logs/verify2.log" 2>/dev/null
  echo; echo "== summaries present"; ls -la "$root"/verify/pair*/*summary*.txt "$root"/verify/pair*/rna/*summary*.txt 2>/dev/null
  echo; echo "== hifi files"; ls -la "$root"/verify/pair*/hifi/ 2>/dev/null
  echo; echo "== rna files"; ls -la "$root"/verify/pair*/rna/ 2>/dev/null
  echo; echo "== still running"; pgrep -fl "minimap2|rna_check|hifi_check|verify_driver" || true
  echo; echo "== TODO evening: HiFi pair1 (18 GB) and pair2 (20 GB) full volume; fix pstream speed (8 parallel = 0.33 MB/s vs single curl 0.97 MB/s); RNA pair3 if not finished"
} > "$snap"
say "STATE_SNAPSHOT $snap"

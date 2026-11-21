#!/usr/bin/env bash
# RNA-seq part of the night run: wait for the downloads, map each pair's RNA reads to its windows and run rna_check.py.
# Runs alongside verify_driver.sh (HiFi stream is network-bound). Logged to logs/verify.log.
set -uo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
log="$root/logs/verify.log"
say() { echo "$(date '+%F %T') $*" | tee -a "$log"; }
say "RNA_CHECK_START"
while ! grep -q "RNA_DONE verify/pair3/rna" "$root/logs/get_rna.log" 2>/dev/null; do
  if ! pgrep -f get_rna.sh >/dev/null; then say "RNA chain not running and not finished; continuing with what is there"; break; fi
  sleep 60
done
for pair in pair2 pair1 pair3; do
  if ! ls "$root/verify/$pair/rna/"*.fastq.gz >/dev/null 2>&1; then say "RNA_CHECK_SKIP $pair (no reads)"; continue; fi
  say "RNA_CHECK $pair"
  if nice -n 10 python3 -I "$root/scripts/rna_check.py" "$root/$pair" "$root/verify/$pair" --map > "$root/logs/${pair}.rna_check.log" 2>&1; then
    say "RNA_CHECK_OK $pair"
  else
    say "RNA_CHECK_FAILED $pair (see logs/${pair}.rna_check.log)"
  fi
done
say "RNA_CHECK_DONE"

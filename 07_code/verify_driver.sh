#!/usr/bin/env bash
# Night run: wait for the RNA-seq downloads, then stream HiFi read prefixes of each assembly's own individual
# through minimap2 (pilot windows of both haplotypes) and run the read-level checks. Everything logged to logs/verify.log.
set -uo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
log="$root/logs/verify.log"
say() { echo "$(date '+%F %T') $*" | tee -a "$log"; }
say "VERIFY_START"
# 1. wait for RNA downloads (get_rna.sh chain writes RNA_DONE for the last directory)
while ! grep -q "RNA_DONE verify/pair3/rna" "$root/logs/get_rna.log" 2>/dev/null; do
  if ! pgrep -f get_rna.sh >/dev/null; then say "RNA chain not running and not finished; continuing anyway"; break; fi
  sleep 60
done
say "RNA_READY"
# 2. HiFi prefixes: pair3 Eidolon (Revio, 2.7 Gbp/GB), pair1 Myotis (Revio, 2.3 Gbp/GB), pair2 Rhinolophus (Sequel II, 1.1 Gbp/GB)
declare -a JOBS=(
  "pair3|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR395/078/SRR39556278/SRR39556278_subreads.fastq.gz|15000000000"
  "pair1|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR382/013/SRR38296513/SRR38296513_subreads.fastq.gz|18000000000"
  "pair2|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR369/062/SRR36951662/SRR36951662_subreads.fastq.gz|20000000000"
)
for job in "${JOBS[@]}"; do
  IFS='|' read -r pair url limit <<< "$job"
  # SKIP_PAIRS="pair3" re-runs the queue without the pairs that are already done
  if [[ " ${SKIP_PAIRS:-} " == *" $pair "* ]]; then say "HIFI_SKIP $pair (already done)"; continue; fi
  free_gb=$(df -g "$root" | awk 'NR==2{print $4}')
  if [[ "$free_gb" -lt 5 ]]; then say "DISK_LOW ${free_gb}G, stop"; exit 1; fi
  say "HIFI_STREAM $pair $(basename "$url") limit=$limit"
  "$root/scripts/hifi_stream.sh" "$root/verify/$pair" "$url" "$limit" >> "$log" 2>&1
  say "HIFI_CHECK $pair"
  if nice -n 10 python3 -I "$root/scripts/hifi_check.py" "$root/$pair" "$root/verify/$pair" > "$root/logs/${pair}.hifi_check.log" 2>&1; then
    say "HIFI_CHECK_OK $pair"
  else
    say "HIFI_CHECK_FAILED $pair (see logs/${pair}.hifi_check.log)"
  fi
done
say "VERIFY_DONE"

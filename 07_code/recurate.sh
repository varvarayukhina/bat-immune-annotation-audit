#!/usr/bin/env bash
# Re-download the two genomes of a pair that has already been through the pipeline, run curate.py again, delete the genomes.
# usage: recurate.sh PAIRTAG ACC_HAP1 ACC_HAP2
set -euo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
tag="$1"; acc1="$2"; acc2="$3"
w="$root/$tag"
for hap in hap1 hap2; do
  acc="$acc1"; [[ "$hap" == hap2 ]] && acc="$acc2"
  if ! ls "$w"/genomes/"$hap"/*_genomic.fna >/dev/null 2>&1; then
    "$root/scripts/dl_genome.sh" "$acc" "$w/genomes/$hap" > "$w/logs/dl_${hap}_recurate.log" 2>&1
    rm -f "$w"/genomes/"$hap"/*.zip
  fi
  "$root/bin/samtools" faidx "$(ls "$w"/genomes/"$hap"/*_genomic.fna)"
done
if [[ -d "$w/curation" ]]; then
  mv "$w/curation" "$w/curation_first_run_$(date +%H%M)"
fi
nice -n 10 python3 -I "$root/scripts/curate.py" "$w" "$tag"
rm -f "$w"/genomes/hap1/*_genomic.fna "$w"/genomes/hap2/*_genomic.fna "$w"/genomes/hap*/*.fai
echo "$(date '+%F %T') RECURATION_OK $tag" >> "$root/logs/driver.log"

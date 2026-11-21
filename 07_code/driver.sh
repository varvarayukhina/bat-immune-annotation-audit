#!/usr/bin/env bash
# Runs the three pairs one after another: pipeline -> curation -> delete genome FASTA.
# Everything is nice'd; a failure in one pair keeps its genomes and moves on to the next pair.
set -uo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
log="$root/logs/driver.log"
mkdir -p "$root/logs"
ts() { date '+%F %T'; }
free_gb() { df -g "$HOME" | awk 'NR==2{print $4}'; }

run_one() { # tag acc1 acc2
  local tag="$1" a1="$2" a2="$3"
  if (( $(free_gb) < 9 )); then
    echo "$(ts) SKIP $tag: only $(free_gb) GB free on disk" >> "$log"
    return
  fi
  echo "$(ts) START $tag" >> "$log"
  if nice -n 10 "$root/scripts/run_pair.sh" "$tag" "$a1" "$a2" >> "$root/logs/$tag.run.log" 2>&1; then
    echo "$(ts) PIPELINE_OK $tag" >> "$log"
  else
    echo "$(ts) PIPELINE_FAILED $tag (genomes kept; see logs/$tag.run.log)" >> "$log"
    return
  fi
  if nice -n 10 python3 -I "$root/scripts/curate.py" "$root/$tag" "$tag" >> "$root/logs/$tag.curate.log" 2>&1; then
    echo "$(ts) CURATION_OK $tag" >> "$log"
    rm -f "$root/$tag"/genomes/hap1/*_genomic.fna "$root/$tag"/genomes/hap2/*_genomic.fna "$root/$tag"/genomes/hap*/*.fai
    echo "$(ts) GENOMES_DELETED $tag" >> "$log"
  else
    echo "$(ts) CURATION_FAILED $tag (genomes kept; see logs/$tag.curate.log)" >> "$log"
  fi
}

echo "$(ts) DRIVER_START host=$(hostname) free=$(free_gb)GB" >> "$log"
run_one pair1 GCA_057945255.1 GCA_057945245.1
run_one pair2 GCA_054645015.1 GCA_054645025.1
run_one pair3 GCA_059749095.1 GCA_059749105.1
echo "$(ts) ALL_DONE" >> "$log"

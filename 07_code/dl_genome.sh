#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
acc="$1"; dir="$2"
mkdir -p "$dir"
curl -fsS --retry 5 --retry-delay 5 --retry-all-errors -o "$dir/$acc.zip" \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=GENOME_FASTA&filename=${acc}.zip"
unzip -l "$dir/$acc.zip"
unzip -o -j "$dir/$acc.zip" '*_genomic.fna' -d "$dir"
ls -l "$dir"

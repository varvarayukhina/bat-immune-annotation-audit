#!/usr/bin/env bash
# Download (prefixes of) paired RNA-seq runs from ENA with parallel range requests.
# usage: get_rna.sh OUTDIR LIMIT_BYTES RUN [RUN ...]   (LIMIT_BYTES per mate; 0 = whole file)
set -euo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
out="$1"; limit="$2"; shift 2
mkdir -p "$out"
tmp="${TMPDIR:-/tmp}/pstream_rna"
for run in "$@"; do
  rep="$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp,fastq_bytes,read_count,base_count,instrument_model,library_layout" | tail -n +2)"
  [[ -n "$rep" ]] || { echo "no ENA report for $run" >&2; exit 1; }
  echo -e "${run}\t${rep}" >> "$out/manifest.tsv"
  # ENA prepends run_accession: f1=run, f2=fastq_ftp, f3=fastq_bytes
  IFS=';' read -r -a urls <<< "$(echo "$rep" | cut -f2)"
  IFS=';' read -r -a sizes <<< "$(echo "$rep" | cut -f3)"
  for i in "${!urls[@]}"; do
    u="https://${urls[$i]}"; sz="${sizes[$i]}"
    n="$sz"; if [[ "$limit" -gt 0 && "$limit" -lt "$sz" ]]; then n="$limit"; fi
    f="$out/$(basename "${urls[$i]}")"
    if [[ -s "$f" ]]; then echo "$(date '+%T') have $f"; continue; fi
    echo "$(date '+%T') $run mate $((i+1)): $n of $sz bytes -> $f"
    python3 -I "$root/scripts/pstream.py" "$u" "$n" 8 33554432 "$tmp" > "$f.part"
    mv "$f.part" "$f"
  done
done
echo "$(date '+%T') RNA_DONE $out"

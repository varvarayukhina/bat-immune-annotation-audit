#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'
# usage: run_miniprot.sh genome.fna hapLabel
g="$1"; lab="$2"; base="$(cd "$(dirname "$0")/.." && pwd)"
/usr/bin/time -l miniprot -t4 -j2 -p 0.3 -N 100 --outs=0.3 --outc=0.3 --outn=200 --gff -P "$lab" \
  "$g" "$base/queries/queries_final.faa" > "$base/mapping/$lab.miniprot.gff" 2> "$base/logs/miniprot_$lab.log"

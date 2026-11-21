#!/usr/bin/env bash
# Reference proteomes for the back-search (family verification): UniProt Swiss-Prot, reviewed entries, human + mouse.
set -euo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
d="$root/ref"
mkdir -p "$d"
for tax in 9606 10090; do
  url="https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A${tax}%29%20AND%20%28reviewed%3Atrue%29&format=fasta"
  curl -fsS --retry 5 --retry-all-errors -D "$d/sp_$tax.headers" -o "$d/sp_$tax.fasta" "$url"
  printf '%s entries: %s  release: %s\n' "$tax" "$(grep -c '>' "$d/sp_$tax.fasta")" "$(grep -i '^x-uniprot-release:' "$d/sp_$tax.headers" | tr -d '\r' | cut -d' ' -f2-)"
done
cat "$d/sp_9606.fasta" "$d/sp_10090.fasta" > "$d/swissprot_hm.fasta"
diamond makedb --in "$d/swissprot_hm.fasta" -d "$d/swissprot_hm" --quiet
ls -l "$d"

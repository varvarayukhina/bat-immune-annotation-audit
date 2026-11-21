#!/usr/bin/env bash
# Stream a prefix of an ENA HiFi fastq.gz through minimap2 against both haplotypes' pilot windows; keep only the PAF.
# usage: hifi_stream.sh VERIFY_DIR URL LIMIT_BYTES
# Reads are never stored. A gzip prefix decompresses up to the cut; the last (truncated) FASTQ record is dropped by the awk filter.
set -uo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
v="$1"; url="$2"; limit="$3"
run="$(basename "$url" | cut -d_ -f1)"
mkdir -p "$v/hifi"
ref="$v/hifi/ref_both.fa"
[[ -s "$ref" ]] || cat "$v/windows_hap1.fa" "$v/windows_hap2.fa" > "$ref"
paf="$v/hifi/${run}.paf"
echo "$(date '+%F %T') HIFI_START $run limit=$limit -> $paf"
python3 -I "$root/scripts/pstream.py" "$url" "$limit" 8 33554432 "${TMPDIR:-/tmp}/pstream_hifi" 2>> "$v/hifi/${run}.pstream.log" \
  | gzip -dc 2>/dev/null \
  | paste - - - - | awk -F'\t' 'NF==4 && length($2)==length($4) && length($2)>0 {print $1"\n"$2"\n"$3"\n"$4}' \
  | nice -n 10 minimap2 -x map-hifi -t 3 --secondary=no -c --cs "$ref" - > "$paf.part" 2> "$v/hifi/${run}.minimap2.log"
st=("${PIPESTATUS[@]}")
mv "$paf.part" "$paf"
echo "$(date '+%F %T') HIFI_DONE $run pipestatus=${st[*]} paf_lines=$(wc -l < "$paf")"

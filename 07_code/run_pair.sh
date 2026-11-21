#!/usr/bin/env bash
# usage: run_pair.sh PAIRTAG ACC_HAP1 ACC_HAP2   (work dir: <root>/PAIRTAG, shares queries/ and bin/)
# Rerun of 2026-11-12: identical to the 2026-11-09 pipeline except that the genome FASTA files are NOT deleted here
# (driver.sh deletes them after the curation step has finished).
set -euo pipefail
IFS=$'\n\t'
root="$(cd "$(dirname "$0")/.." && pwd)"
tag="$1"; acc1="$2"; acc2="$3"
w="$root/$tag"
mkdir -p "$w"/{mapping,out,logs,genomes}
ln -sfn "$root/queries" "$w/queries"
ln -sfn "$root/bin" "$w/bin"
api="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
for hap in hap1 hap2; do
  acc="$acc1"; [[ "$hap" == hap2 ]] && acc="$acc2"
  curl -fsS --retry 5 --retry-all-errors "$api/$acc/dataset_report" -o "$w/out/report_$acc.json"
  curl -fsS --retry 5 --retry-all-errors "$api/$acc/sequence_reports?page_size=1000" -o "$w/out/seqreport_$acc.json"
  if ! ls "$w"/genomes/"$hap"/*_genomic.fna >/dev/null 2>&1; then
    "$root/scripts/dl_genome.sh" "$acc" "$w/genomes/$hap" > "$w/logs/dl_$hap.log" 2>&1
    rm -f "$w"/genomes/"$hap"/*.zip
  fi
  fna="$(ls "$w"/genomes/"$hap"/*_genomic.fna)"
  "$root/bin/samtools" faidx "$fna"
  python3 -I "$root/scripts/run_miniprot_chunked.py" "$fna" "$hap" "$w" > "$w/logs/chunked_$hap.out" 2>&1
  python3 -I "$root/scripts/analyze_genes.py" "$w" "$tag" "$hap" "$acc" "$w/out/seqreport_$acc.json"
done
f1="$(ls "$w"/genomes/hap1/*_genomic.fna)"; f2="$(ls "$w"/genomes/hap2/*_genomic.fna)"
python3 -I "$root/scripts/compare_regions.py" "$w" "$tag" "$f1" "$f2" "$acc1" "$acc2" "$w/out/seqreport_$acc1.json" "$w/out/seqreport_$acc2.json" > "$w/logs/compare.log" 2>&1
python3 -I "$root/scripts/finalize.py" "$w" "$tag"
# higher-sensitivity pass restricted to the region windows
s="$w/sens"
python3 -I "$root/scripts/run_sens_windows.py" "$w" "$s" hap1 "$f1" "$tag"
python3 -I "$root/scripts/run_sens_windows.py" "$w" "$s" hap2 "$f2" "$tag"
python3 -I "$root/scripts/analyze_genes.py" "$s" "${tag}sens" hap1 "$acc1" "$w/out/seqreport_$acc1.json"
python3 -I "$root/scripts/analyze_genes.py" "$s" "${tag}sens" hap2 "$acc2" "$w/out/seqreport_$acc2.json"
python3 -I "$root/scripts/compare_regions.py" "$s" "${tag}sens" "$f1" "$f2" "$acc1" "$acc2" "$w/out/seqreport_$acc1.json" "$w/out/seqreport_$acc2.json" > "$s/logs/compare.log" 2>&1
python3 -I "$root/scripts/finalize.py" "$s" "${tag}sens"
echo "pair $tag pipeline finished (genomes kept for curation)"

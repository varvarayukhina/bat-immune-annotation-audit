#!/usr/bin/env bash
# Rebuild every figure and derived statistic from the deposited tables. No raw data.
# Needs matplotlib and scipy, or:
#   PYTHON="uv run --with matplotlib --with scipy python3" ./reproduce.sh
set -euo pipefail
IFS=$'\n\t'

# PYTHON may carry arguments, so split it on spaces into an array.
PY=()
IFS=' ' read -r -a PY <<< "${PYTHON:-python3}"
IFS=$'\n\t'

cd "$(dirname "$0")"

# The analysis code reads the layout it ran in; recreate it as links, tables stored once.
for p in pair1 pair2 pair3; do
  mkdir -p "verify/$p"
  for t in hifi_depth refseq_start rna_loci windows; do
    ln -sf "../../05_tables/readcheck/${p}_${t}.tsv" "verify/$p/${t}.tsv"
  done
done

mkdir -p summary
for p in pair1 pair2 pair3; do
  ln -sf "../05_tables/summary/${p}.refined_differences.tsv" "summary/${p}.refined_differences.tsv"
  mkdir -p "$p/out" "$p/curation" "$p/sens/out"
  ln -sf "../../05_tables/summary/${p}.family_summary.tsv" "$p/out/${p}.family_summary.tsv"
  ln -sf "../../05_tables/curation/${p}.locus_curation.tsv" "$p/curation/${p}.locus_curation.tsv"
  ln -sf "../../../05_tables/summary/${p}sens.family_summary.tsv" "$p/sens/out/${p}sens.family_summary.tsv"
done

# analysis_extra.py writes 07_code/derived/; make_figures.py reads manuscript/derived/.
mkdir -p 07_code/derived manuscript/derived
# The KLRD1 junction check is a separate step, deposited as a table.
ln -sf ../../05_tables/derived/klrd1_junction_mapq.tsv manuscript/derived/klrd1_junction_mapq.tsv

# Section 4 needs verify/<pair>/windows_<hap>.fa, not deposited; see 07_code/fetch_windows.py.
# Without them the script stops after section 3. Expected, not a failure.
have_fa=yes
for p in pair1 pair2 pair3; do
  for h in hap1 hap2; do
    [[ -f "verify/$p/windows_$h.fa" ]] || have_fa=no
  done
done

echo "== per-family counts, summed depth, depth screen, confounders =="
if [[ "$have_fa" == yes ]]; then
  "${PY[@]}" 07_code/analysis_extra.py
else
  "${PY[@]}" 07_code/analysis_extra.py || {
    echo
    echo "note: sections 1-3 above are complete. Section 4, the GC and length confounder" >&2
    echo "      check, needs verify/<pair>/windows_<hap>.fa; its deposited result is" >&2
    echo "      05_tables/derived/confounders_gc_length.tsv." >&2
  }
fi

for t in copies_implied_vs_assembled depth_null_fdr; do
  cp "07_code/derived/$t.tsv" "manuscript/derived/$t.tsv"
done

echo
echo "== phasing checks: counterpart identity, KLRD1 paralogue identity =="
"${PY[@]}" 07_code/phase_checks.py

echo "== figures 1-5 =="
"${PY[@]}" 07_code/make_figures.py

echo
echo "derived tables: 07_code/derived/"
echo "figures:        07_code/figures/"
echo "Compare them against the deposited 05_tables/derived/ and 06_figures/."

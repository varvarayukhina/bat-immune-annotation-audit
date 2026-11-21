#!/usr/bin/env bash
# Assemble everything a reviewer needs into one folder, one single-file digest and one zip.
set -euo pipefail
IFS=$'\n\t'

ROOT="/Users/varvarauhina/bat-nk-pilot/rerun_2026-11-12"
OUT="/Users/varvarauhina/bat-nk-pilot/REVIEW_BUNDLE_2026-11-20"

mkdir -p "$OUT"/{04_numbers,05_tables/{derived,summary,curation,readcheck,queries},06_figures,07_code}

# 1. texts
cp "$ROOT/manuscript/submission/manuscript_submission.md" "$OUT/01_manuscript.md"
cp "$ROOT/METHODS_CURATION.txt"                           "$OUT/02_methods_full.txt"
cp "$ROOT/REPORT_2026-11-12.md"                           "$OUT/03_internal_report_ru.md"
cp "$ROOT/manuscript/submission/WHERE_TO_SUBMIT.md"       "$OUT/08_where_to_submit.md"

# 2. machine-written number files
cp "$ROOT/manuscript/analysis_extra.out" "$ROOT/manuscript/make_figures.out" "$OUT/04_numbers/"
cp "$ROOT/summary/overview.tsv" "$ROOT/summary/digest.txt" "$OUT/04_numbers/"
cp "$ROOT/results_final_2026-11-15/verify.log" "$OUT/04_numbers/"
for p in pair1 pair2 pair3; do
  cp "$ROOT/$p/curation/$p.curation_summary.txt" "$OUT/04_numbers/"
  cp "$ROOT/results_final_2026-11-15/${p}_hifi_summary.txt" \
     "$ROOT/results_final_2026-11-15/${p}_rna_summary.txt" \
     "$ROOT/results_final_2026-11-15/${p}_refseq_start_summary.txt" "$OUT/04_numbers/"
done

# 3. tables
cp "$ROOT"/manuscript/derived/*.tsv            "$OUT/05_tables/derived/"
cp "$ROOT"/summary/pair*.refined_differences.tsv "$OUT/05_tables/summary/"
for p in pair1 pair2 pair3; do
  cp "$ROOT/$p/curation/$p.locus_curation.tsv" \
     "$ROOT/$p/curation/$p.difference_curation.tsv" \
     "$ROOT/$p/curation/$p.gap_curation.tsv" \
     "$ROOT/$p/curation/$p.backsearch.tsv" "$OUT/05_tables/curation/"
done
cp "$ROOT"/results_final_2026-11-15/pair*_{hifi_depth,hifi_sites,rna_loci,rna_sites,refseq_start}.tsv \
   "$OUT/05_tables/readcheck/"
cp "$ROOT"/queries/queries_final_table.tsv "$ROOT"/queries/queries_missing.tsv \
   "$ROOT"/queries/queries_alias_excluded.tsv "$OUT/05_tables/queries/"

# 4. figures
cp "$ROOT"/manuscript/figures/*.png "$ROOT"/manuscript/figures/*.pdf "$OUT/06_figures/"

# 5. code
cp "$ROOT"/manuscript/analysis_extra.py "$ROOT"/manuscript/make_figures.py "$OUT/07_code/"
cp "$ROOT"/scripts/*.py "$ROOT"/scripts/*.sh "$OUT/07_code/"

# 6. single-file digest for pasting into a chat: texts + every number file, no big tables
ONE="$OUT/ONEFILE_FOR_GPT.md"
{
  cat "$OUT/00_START_HERE.md"
  printf '\n\n\n%s\n\n' '==================== 01_manuscript.md ===================='
  cat "$OUT/01_manuscript.md"
  printf '\n\n\n%s\n\n```text\n' '==================== 02_methods_full.txt ===================='
  cat "$OUT/02_methods_full.txt"
  printf '```\n\n\n%s\n\n' '==================== 03_internal_report_ru.md ===================='
  cat "$OUT/03_internal_report_ru.md"
  printf '\n\n\n%s\n' '==================== 04_numbers/ ===================='
  for f in "$OUT"/04_numbers/*; do
    printf '\n\n### %s\n\n```text\n' "$(basename "$f")"
    cat "$f"
    printf '```\n'
  done
} > "$ONE"

# 7. zip for upload (figures included)
cd "$(dirname "$OUT")"
staging="$(mktemp -d "$PWD/.bat-nk-bundle.XXXXXX")"
zip -qr "$staging/review_bundle_2026-11-20.zip" "$(basename "$OUT")" -x '*.DS_Store'
mv -f "$staging/review_bundle_2026-11-20.zip" "$PWD/review_bundle_2026-11-20.zip"
rmdir "$staging"

printf 'bundle: %s\n' "$OUT"
printf 'onefile: %s (%s)\n' "$ONE" "$(du -h "$ONE" | cut -f1)"
printf 'zip:     %s (%s)\n' "$PWD/review_bundle_2026-11-20.zip" "$(du -h review_bundle_2026-11-20.zip | cut -f1)"
find "$OUT" -type f | wc -l | xargs printf 'files:   %s\n'

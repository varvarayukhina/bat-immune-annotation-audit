#!/usr/bin/env python3
"""Combine per-haplotype outputs into the three deliverable TSVs. usage: finalize.py base_dir pair_tag"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

base, pair = Path(sys.argv[1]), sys.argv[2]
out = base / "out"


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


# 1. genes per haplotype
loci = rd(out / f"{pair}.hap1.loci.tsv") + rd(out / f"{pair}.hap2.loci.tsv")
keys = list(loci[0].keys())
with (out / f"{pair}.genes_per_haplotype.tsv").open("w") as f:
    f.write("\t".join(keys) + "\n")
    for r in loci:
        f.write("\t".join(r[k] for k in keys) + "\n")

# 1b. family x haplotype
fam = defaultdict(lambda: defaultdict(Counter))
for r in loci:
    fam[r["family"]][r["haplotype"]][r["status"]] += 1
    fam[r["family"]][r["haplotype"]]["total"] += 1
with (out / f"{pair}.family_summary.tsv").open("w") as f:
    f.write("family\thap\tloci\tintact\tpartial\tpseudogene_like\tuncertain\n")
    for fm in sorted(fam):
        for h in ("hap1", "hap2"):
            c = fam[fm][h]
            f.write(f"{fm}\t{h}\t{c['total']}\t{c['intact']}\t{c['partial']}\t{c['pseudogene-like']}\t{c['uncertain']}\n")

# 2. region summary (one row per hap1/hap2 window pair)
regs = {r["region_id"]: r for h in ("hap1", "hap2") for r in rd(out / f"{pair}.{h}.regions.tsv")}
pairs = rd(out / f"{pair}.region_pairs.tsv")
diffs = rd(out / f"{pair}.haplotype_differences.tsv")


def side(rid, h):
    if not rid:
        return [""] * 11
    r = regs[rid]
    return [r["chromosome"].split(":")[0], r["contig"], r["start"], r["end"], r["span_bp"], r["n_loci"], r["intact"], r["partial"], r["pseudogene_like"], r["uncertain"], r["families"]]


hdr = ["region_hap1", "region_hap2", "type"]
for h in ("hap1", "hap2"):
    hdr += [f"{h}_chr", f"{h}_contig", f"{h}_window_start", f"{h}_window_end", f"{h}_window_bp", f"{h}_loci", f"{h}_intact", f"{h}_partial", f"{h}_pseudogene_like", f"{h}_uncertain", f"{h}_families", f"{h}_extract_len_pm100kb", f"{h}_frac_aligned"]
hdr += ["gaps_ge5kb_hap1only", "gaps_ge5kb_hap2only", "internal_indels_ge5kb", "inverted_blocks", "strand_switches", "loci_matched_identical_call", "n_diff_rows", "flags_hap1_(within200kb_of_end_or_gap)", "flags_hap2_(within200kb_of_end_or_gap)", "note"]
rows = []
for p in pairs:
    r1, r2 = p["region_hap1"], p["region_hap2"]
    nd = sum(1 for d in diffs if d["region_hap1"] == r1 and d["region_hap2"] == r2)
    rows.append([r1, r2, p["type"]] + side(r1, "hap1")[:] + [p["hap1_len"], p["hap1_frac_aligned"]] + side(r2, "hap2") + [p["hap2_len"], p["hap2_frac_aligned"]] +
                [p["gaps_ge5kb_hap1"], p["gaps_ge5kb_hap2"], p["internal_indels_ge5kb"], p["inverted_blocks"], p["strand_switches"], p["loci_matched_identical_call"], nd,
                 p["flags_hap1"], p["flags_hap2"], p["method_note"]])
# fix column order: side() has 11 entries then extract_len and frac -> 13 per hap, matches hdr
with (out / f"{pair}.region_summary.tsv").open("w") as f:
    f.write("\t".join(hdr) + "\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print("genes", len(loci), "families", len(fam), "region rows", len(rows), "diff rows", len(diffs))

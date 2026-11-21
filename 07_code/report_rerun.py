#!/usr/bin/env python3
"""Digest of the rerun + curation for all pairs that have a curation/ folder.
usage: report_rerun.py rerun_root original_root
Writes rerun_root/summary/overview.tsv, <pair>.refined_differences.tsv, digest.txt"""
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

root, orig = Path(sys.argv[1]), Path(sys.argv[2])
outd = root / "summary"
outd.mkdir(exist_ok=True)
NAMES = {"pair1": "Myotis sodalis", "pair2": "Rhinolophus clivosus", "pair3": "Eidolon dupreanum"}
BIG = 50  # indel longer than this = structural change, not a point frameshift


def rd(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def identical(a, b):
    return a.exists() and b.exists() and a.read_bytes() == b.read_bytes()


def refine(r):
    """final_class of difference_curation -> short class for the report"""
    fc = r["final_class"]
    base = fc.split(";")[0]
    if base.startswith("same_status"):
        return "artefact_sensitivity_or_threshold"
    if base.startswith("status_differs"):
        items = r["variants_in_cds_of_better_copy"]
        big = [m for m in re.findall(r"(?:ins|del)(\d+)@", items) if int(m) > BIG]
        if big:
            return f"status_differs_large_indel_in_gene({max(map(int, big))}bp)"
        if int(r["n_frameshift_indels"] or 0) or int(r["n_stop_gains"] or 0) or "splice_site" in items:
            return "status_differs_with_disruptive_variant"
        return "status_differs_no_disruptive_variant_found(threshold_effect_likely)"
    if "absent" in base:
        return "locus_absent_from_other_haplotype"
    if base.startswith("sequence_aligned_but_no_gene"):
        return "sequence_present_no_gene_hit"
    if base.startswith("gene_hit_missing"):
        return "gene_hit_missing_one_side_at_sensitive"
    return base


ov_rows, digest = [], []
for pair in ("pair1", "pair2", "pair3"):
    cur = root / pair / "curation"
    if not (cur / f"{pair}.curation_summary.txt").exists():
        continue
    L = rd(cur / f"{pair}.locus_curation.tsv")
    D = rd(cur / f"{pair}.difference_curation.tsv")
    G = rd(cur / f"{pair}.gap_curation.tsv")
    rep = {t: identical(orig / pair / f"{pair}.{t}.tsv", root / pair / "out" / f"{pair}.{t}.tsv") for t in ("family_summary", "genes_per_haplotype", "haplotype_differences", "region_summary")}
    reps = {t: identical(orig / pair / "sensitive_windows" / f"{pair}sens.{t}.tsv", root / pair / "sens" / "out" / f"{pair}sens.{t}.tsv") for t in ("family_summary", "genes_per_haplotype", "haplotype_differences")}
    o = {"pair": pair, "species": NAMES[pair], "reproduced_main_tables": "yes" if all(rep.values()) else "NO:" + ",".join(k for k, v in rep.items() if not v),
         "reproduced_sensitive_tables": "yes" if all(reps.values()) else "NO:" + ",".join(k for k, v in reps.items() if not v)}
    for h in ("hap1", "hap2"):
        rows = [r for r in L if r["haplotype"] == h]
        intact = [r for r in rows if r["pilot_status"] == "intact"]
        o[f"{h}_loci"] = len(rows)
        o[f"{h}_realign_same_status"] = sum(r["realign_consistent"] == "yes" for r in rows)
        o[f"{h}_intact_pilot"] = len(intact)
        o[f"{h}_intact_full_gene_model"] = sum(r["curated_intact"] == "yes" for r in intact)
        o[f"{h}_intact_fail_start_only"] = sum(r["curated_intact"] == "no" and r["curation_note"].startswith("start:") and "," not in r["curation_note"] for r in intact)
        o[f"{h}_noncanonical_splice_loci"] = sum(r["noncanonical_splice"] not in ("-", "") for r in rows)
        o[f"{h}_backsearch_confirmed"] = sum(r["family_check"].startswith("confirmed") for r in rows)
        o[f"{h}_backsearch_relabel"] = sum(r["family_check"].startswith("relabel") for r in rows)
        o[f"{h}_backsearch_ambiguous"] = sum(r["family_check"].startswith("ambiguous") for r in rows)
        o[f"{h}_no_protein"] = sum(r["family_check"] in ("no_protein", "no_swissprot_hit") for r in rows)
        o[f"{h}_loci_with_ge98pct_paralog"] = sum(1 for r in rows if r["paralog_pident"] not in ("-", "") and float(r["paralog_pident"]) >= 98)
    ref = Counter(refine(r) for r in D)
    o["diff_rows_pilot"] = None
    o["locus_differences_examined"] = len(D)
    for k in ("artefact_sensitivity_or_threshold", "status_differs_with_disruptive_variant", "status_differs_no_disruptive_variant_found(threshold_effect_likely)",
              "locus_absent_from_other_haplotype", "sequence_present_no_gene_hit", "gene_hit_missing_one_side_at_sensitive"):
        o[k] = ref.get(k, 0)
    o["status_differs_large_indel"] = sum(v for k, v in ref.items() if k.startswith("status_differs_large_indel"))
    o["homopolymer_frameshifts"] = sum("homopolymer" in r["variants_in_cds_of_better_copy"] for r in D)
    o["multi_copy_landings"] = sum("multi_copy_landing" in r["final_class"] for r in D)
    gc = Counter(r["gap_class"].split(";")[0].split("(")[0] for r in G)
    o["gaps_unique"] = len(G)
    o["gaps_present_in_paired_window"] = gc.get("present_in_paired_window", 0)
    o["gaps_present_elsewhere"] = gc.get("present_elsewhere_in_hap1", 0) + gc.get("present_elsewhere_in_hap2", 0)
    o["gaps_absent"] = gc.get("absent_from_hap1", 0) + gc.get("absent_from_hap2", 0)
    o["gaps_N_run"] = gc.get("N_run", 0)
    o["gaps_bp_absent"] = sum(int(r["length"]) for r in G if r["gap_class"].startswith("absent"))
    summ = (cur / f"{pair}.curation_summary.txt").read_text()
    m = re.search(r"pilot difference rows: (\d+); of these (\d+)", summ)
    o["diff_rows_pilot"] = m.group(1) if m else ""
    o["diff_rows_duplicated_across_windows"] = m.group(2) if m else ""
    ov_rows.append(o)
    # refined table
    with open(outd / f"{pair}.refined_differences.tsv", "w") as f:
        cols = ["refined_class", "diff_index", "feature", "haplotype", "locus_id", "family", "pilot_status", "src_status_sensitive", "counterpart_status_sensitive",
                "locus_span_aligned_frac", "flank_aligned_frac_left_right", "n_alt_landings", "variants_in_cds_of_better_copy", "flags", "final_class"]
        f.write("\t".join(cols) + "\n")
        for r in sorted(D, key=lambda r: refine(r)):
            f.write("\t".join([refine(r)] + [r[c] for c in cols[1:]]) + "\n")
    # digest text
    digest.append(f"=== {pair} {NAMES[pair]}")
    digest.append(f"reproduced: main {o['reproduced_main_tables']}, sensitive {o['reproduced_sensitive_tables']}")
    for h in ("hap1", "hap2"):
        digest.append(f"{h}: loci {o[h + '_loci']}, same status on +-5kb realignment {o[h + '_realign_same_status']}; intact {o[h + '_intact_pilot']} -> full gene model {o[h + '_intact_full_gene_model']} "
                      f"(start codon only reason: {o[h + '_intact_fail_start_only']}); back-search confirmed {o[h + '_backsearch_confirmed']}, relabel {o[h + '_backsearch_relabel']}, "
                      f"ambiguous {o[h + '_backsearch_ambiguous']}, no protein {o[h + '_no_protein']}; >=98% paralog {o[h + '_loci_with_ge98pct_paralog']}; non-canonical splice {o[h + '_noncanonical_splice_loci']}")
    rel = Counter(f'{r["haplotype"]} {r["family"]}->{r["back_symbol"]}' for r in L if r["family_check"].startswith("relabel"))
    digest.append("relabels: " + "; ".join(f"{k} x{v}" for k, v in rel.most_common()) if rel else "relabels: none")
    amb = Counter(f'{r["family"]}:{r["family_check"]}' for r in L if r["family_check"].startswith("ambiguous"))
    if amb:
        digest.append("ambiguous: " + "; ".join(f"{k} x{v}" for k, v in amb.most_common(10)))
    digest.append(f"differences: pilot rows {o['diff_rows_pilot']} (dup {o['diff_rows_duplicated_across_windows']}), locus-level examined {len(D)}:")
    for k, v in ref.most_common():
        digest.append(f"   {v:3d} {k}")
    digest.append("   by family (locus-level): " + ", ".join(f"{k}={v}" for k, v in Counter(r["family"].split("|")[0] for r in D).most_common()))
    digest.append(f"   artefact share by family: " + ", ".join(f"{fam}={n_art}/{n}" for fam, n, n_art in sorted(
        ((fam, sum(1 for r in D if r["family"].split("|")[0] == fam), sum(1 for r in D if r["family"].split("|")[0] == fam and refine(r) == "artefact_sensitivity_or_threshold"))
         for fam in {r["family"].split("|")[0] for r in D}), key=lambda x: -x[1])))
    digest.append(f"gaps >=5kb unique {len(G)}: " + ", ".join(f"{k}={v}" for k, v in gc.most_common()) + f"; bp truly absent {o['gaps_bp_absent']}")
    # per-family paralog clusters
    par = defaultdict(list)
    for r in L:
        if r["paralog_pident"] not in ("-", "") and float(r["paralog_pident"]) >= 98:
            par[(r["haplotype"], r["family"])].append(f'{r["locus_id"]}~{r["closest_paralog"]}({r["paralog_pident"]})')
    for k, v in sorted(par.items()):
        digest.append(f"   >=98% paralogs {k[0]} {k[1]}: " + ", ".join(v))
    nc = [r for r in L if r["noncanonical_splice"] not in ("-", "")]
    if nc:
        digest.append("non-canonical splice sites: " + "; ".join(f'{r["haplotype"]} {r["locus_id"]} {r["family"]} {r["pilot_status"]} [{r["noncanonical_splice"]}]' for r in nc))
    # KLRD1 detail (query-bias suspicion in Eidolon)
    kl = [r for r in L if r["family"] == "KLRD1"]
    if len(kl) > 4:
        digest.append("KLRD1 loci: " + "; ".join(f'{r["haplotype"]} {r["locus_id"]} {r["pilot_status"]} back={r["back_symbol"]}({r["back_pident"]}) second={r["back_second"]} paralog={r["paralog_pident"]} start={r["start_codon"]}' for r in kl))
    digest.append("")

cols = list(ov_rows[0].keys())
with open(outd / "overview.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for o in ov_rows:
        f.write("\t".join(str(o[c]) for c in cols) + "\n")
(outd / "digest.txt").write_text("\n".join(digest) + "\n")
print("\n".join(digest))

#!/usr/bin/env python3
"""Read-level checks of the curation calls with HiFi reads of the SAME individual (SRA runs of the assembly BioSample),
mapped to the pilot windows of both haplotypes (verify/<pair>/windows_hap1.fa + windows_hap2.fa).

usage: hifi_check.py pair_dir verify_dir
inputs : verify_dir/hifi/*.paf   (minimap2 -x map-hifi -c --cs --secondary=no, both haplotypes' windows as one reference)
         verify_dir/windows.tsv, pair_dir/curation/<pair>.{difference,locus,gap}_curation.tsv
outputs: verify_dir/hifi_sites.tsv   read support at every disruptive variant (better copy and the copy that carries it)
         verify_dir/hifi_depth.tsv   depth of every locus / hap-only locus / absent gap relative to the window median
         verify_dir/hifi_summary.txt
Filters: primary alignments, aligned block >= 2 kb, divergence (de) <= 0.02."""
import sys
from array import array
from collections import Counter
from itertools import accumulate
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readcheck_lib import DIFFFLANK, VAR, events, fmt_counts, lift_table, load_paf, load_windows, rd, site_counts, variant_sides, window_of  # noqa: E402

pair, vdir = Path(sys.argv[1]), Path(sys.argv[2])
tag = pair.name
cur = pair / "curation"
tmp = vdir / "tmp"
tmp.mkdir(exist_ok=True)
MIN_BLOCK, MAX_DE, MAPQ = 2_000, 0.02, 20

W, WSEQ = load_windows(vdir)
recs, n_in, n_keep = load_paf(sorted((vdir / "hifi").glob("*.paf")), MIN_BLOCK, MAX_DE)

# ---------- depth per window ----------
DEPTH = {}  # window name -> (prefix_all, prefix_mq, median_all)
for w in W:
    L = w["end"] - w["start"] + 1
    d_all, d_mq = array("i", bytes(4 * (L + 1))), array("i", bytes(4 * (L + 1)))
    for r in recs.get(w["name"], []):
        d_all[r["ts"]] += 1
        d_all[min(L, r["te"])] -= 1
        if r["mq"] >= MAPQ:
            d_mq[r["ts"]] += 1
            d_mq[min(L, r["te"])] -= 1
    dep_all = list(accumulate(d_all))[:L]
    dep_mq = list(accumulate(d_mq))[:L]
    med = median(dep_all[::100]) if L > 100 else 0
    DEPTH[w["name"]] = ([0] + list(accumulate(dep_all)), [0] + list(accumulate(dep_mq)), med)


def mean_depth(wname, off_s, off_e):
    """mean depth over window offsets [off_s, off_e] (1-based inclusive) -> (all, mapq20)"""
    P, Q, _ = DEPTH[wname]
    s, e = max(0, off_s - 1), min(len(P) - 1, off_e)
    if e <= s:
        return 0.0, 0.0
    return (P[e] - P[s]) / (e - s), (Q[e] - Q[s]) / (e - s)


def spanning(wname, off_s, off_e, pad=1000, mq=MAPQ):
    return sum(1 for r in recs.get(wname, []) if r["mq"] >= mq and r["ts"] + 1 <= off_s - pad and r["te"] >= off_e + pad)


def depth_row(hap, contig, s, e, pad=3000):
    """mean depth over the region +-pad (a locus of 2-5 kb alone gives a noisy estimate at 10x), spanning reads over the region +-1 kb"""
    loc = window_of(W, hap, contig, s, e)
    if not loc:
        return None
    w, a, b = loc
    L = w["end"] - w["start"] + 1
    da, dm = mean_depth(w["name"], max(1, a - pad), min(L, b + pad))
    med = DEPTH[w["name"]][2]
    return {"window": w["region"], "depth_all": f"{da:.1f}", "depth_mapq20": f"{dm:.1f}", "window_median": f"{med:.1f}",
            "ratio_all": f"{da / med:.2f}" if med else "-", "ratio_mapq20": f"{dm / med:.2f}" if med else "-", "spanning_reads_mapq20": spanning(w["name"], a, b)}


def ratio_verdict(ratio, med, n_span):
    if not med or med < 4:
        return "low_coverage"
    if ratio >= 1.7:
        return "about_2x:sequence_present_in_both_haplotypes(other_assembly_lacks_it)"
    if ratio >= 0.55:
        return "about_1x:supported_as_assembled" + ("" if n_span >= 2 else "(no_spanning_read)")
    if ratio >= 0.2:
        return "below_1x:reads_split_over_copies(possible_assembly_duplication)"
    return "near_0:no_read_support(assembly_artefact?)"


def add_depth(rows, kind, hap, ident, family, status, contig, s, e, paralog="-", closest="-"):
    d = depth_row(hap, contig, s, e)
    row = {"kind": kind, "hap": hap, "id": ident, "family": family, "status": status, "contig": contig, "start": s, "end": e, "length": e - s + 1,
           "paralog_pident": paralog, "closest_paralog": closest}
    if d:
        row.update(d)
        row["verdict"] = ratio_verdict(float(d["ratio_all"]) if d["ratio_all"] != "-" else 0, float(d["window_median"]), d["spanning_reads_mapq20"])
    else:
        row["verdict"] = "outside_windows"
    rows.append(row)


# ---------- A. depth of loci, hap-only loci, absent gaps ----------
drows = []
for lc in rd(cur / f"{tag}.locus_curation.tsv"):
    add_depth(drows, "locus", lc["haplotype"], lc["locus_id"], lc["family"], lc["pilot_status"], lc["contig"], int(lc["start"]), int(lc["end"]), lc["paralog_pident"], lc["closest_paralog"])
D = rd(cur / f"{tag}.difference_curation.tsv")
for r in D:
    if "absent" in r["final_class"] and r["haplotype"] in ("hap1", "hap2"):
        add_depth(drows, "hap_only_locus", r["haplotype"], r["locus_id"], r["family"], r["pilot_status"], r["contig"], int(r["start"]), int(r["end"]))
for g in rd(cur / f"{tag}.gap_curation.tsv"):
    if g["gap_class"].startswith("absent"):
        add_depth(drows, "absent_gap", g["haplotype"], f'gap{g["diff_index"]}', "-", g["gap_class"].split(";")[0], g["contig"], int(g["start"]), int(g["end"]))
DCOLS = ["kind", "hap", "id", "family", "status", "contig", "start", "end", "length", "window", "depth_all", "depth_mapq20", "window_median", "ratio_all", "ratio_mapq20",
         "spanning_reads_mapq20", "paralog_pident", "closest_paralog", "verdict"]
with open(vdir / "hifi_depth.tsv", "w") as f:
    f.write("\t".join(DCOLS) + "\n")
    for r in drows:
        f.write("\t".join(str(r.get(c, "-")) for c in DCOLS) + "\n")

# ---------- B. variant sites ----------
srows = []
for r in D:
    if not r["final_class"].startswith("status_differs_confirmed") or r["variants_in_cds_of_better_copy"] in ("-", ""):
        continue
    ref_h, ref_loc, ref_st, alt_h, alt_loc, alt_st, classes = variant_sides(r)
    a_ref, a_alt = max(1, ref_loc[1] - DIFFFLANK), max(1, alt_loc[1] - DIFFFLANK)
    wr = window_of(W, ref_h, ref_loc[0], a_ref, ref_loc[2] + DIFFFLANK)
    wa = window_of(W, alt_h, alt_loc[0], a_alt, alt_loc[2] + DIFFFLANK)
    base = {"diff_index": r["diff_index"], "locus": r["locus_id"], "family": r["family"], "classes": classes, "better_hap": ref_h, "alt_hap": alt_h}
    if not wr or not wa:
        srows.append({**base, "variant": r["variants_in_cds_of_better_copy"], "verdict": "outside_windows(" + ("ref" if not wr else "alt") + ")"})
        continue
    ref_seq = WSEQ[wr[0]["name"]][wr[1] - 1:wr[2]]
    alt_seq = WSEQ[wa[0]["name"]][wa[1] - 1:wa[2]]
    lift, strand, alen = lift_table(ref_seq, alt_seq, tmp, f"lift{r['diff_index']}")
    for item in r["variants_in_cds_of_better_copy"].split(";"):
        mv = VAR.match(item)
        if not mv:
            continue
        exon, rb, ab, spos, _stop, ikind, ilen, ipos, hp = mv.groups()
        pos1 = int(spos or ipos)
        n = int(ilen) if ilen else 1
        hp = int(hp) if hp else 0
        ref_off = wr[1] + pos1 - 1
        alt_p = lift(pos1)
        row = {**base, "variant": item, "ref_contig_pos": f"{ref_loc[0]}:{a_ref + pos1 - 1}", "alt_contig_pos": "-", "alt_strand": strand}
        slop = max(3, hp + 2) if ikind else 2
        if ikind and n > 50:  # large indel: depth over the segment on the copy that has it + junction reads on the copy that lacks it
            if ikind.endswith("del"):
                da, dm = mean_depth(wr[0]["name"], ref_off, ref_off + n - 1)
                med = DEPTH[wr[0]["name"]][2]
                ratio = da / med if med else 0
                row["ref_site"] = f"segment_depth={da:.1f}(mapq20 {dm:.1f}) window_median={med:.1f} ratio={ratio:.2f}"
                if alt_p:
                    alt_off = wa[1] + alt_p - 1
                    row["alt_contig_pos"] = f"{alt_loc[0]}:{a_alt + alt_p - 1}"
                    span = [rr for rr in recs.get(wa[0]["name"], []) if rr["ts"] + 1 <= alt_off - 500 and rr["te"] >= alt_off + 500]
                    clean = sum(1 for rr in span if not any(e.startswith(("ins", "del")) and int(e[3:]) > 50 for e in (events(rr, alt_off, 500) or [])))
                    row["alt_site"] = f"junction_spanning_reads={len(span)} clean(no_indel>50bp within 500bp)={clean}"
                    row["verdict"] = (f"deletion_real_in_{alt_h}" if (ratio < 1.5 and clean >= 2) else "segment_at_2x:present_in_both_haplotypes(alt_assembly_lacks_it)" if ratio >= 1.7 else "ambiguous")
                else:
                    row["alt_site"], row["verdict"] = "lift_failed", "ambiguous(lift_failed)"
            else:
                row["ref_site"] = row["alt_site"] = row["verdict"] = "large_insertion_not_evaluated"
            srows.append(row)
            continue
        cr = site_counts(recs.get(wr[0]["name"], []), ref_off, slop)
        row["ref_site"] = fmt_counts(cr)
        if alt_p is None:
            row["alt_site"], row["verdict"] = "lift_failed", "ambiguous(lift_failed)"
            srows.append(row)
            continue
        alt_off = wa[1] + alt_p - 1
        row["alt_contig_pos"] = f"{alt_loc[0]}:{a_alt + alt_p - 1}"
        ca = site_counts(recs.get(wa[0]["name"], []), alt_off, slop)
        row["alt_site"] = fmt_counts(ca)
        # reverse event at the alt copy = reads carry the better copy's allele -> the alt assembly is wrong there
        if ikind:
            rev_alt = ca["del"] if ikind.endswith("ins") else ca["ins"]
            same_ref = cr["ins"] if ikind.endswith("ins") else cr["del"]
        else:
            rev_alt, same_ref = ca["sub"], cr["sub"]
        row["alt_reads_with_better_allele"], row["ref_reads_with_alt_allele"] = rev_alt, same_ref
        if ca["n"] < 3:
            v = "low_coverage_at_alt"
        elif rev_alt <= max(1, 0.2 * ca["n"]):  # reads agree with the assembled alt copy (one stray read = HiFi error)
            v = f"variant_supported_by_reads_of_{alt_h}"
        elif rev_alt >= 2 and rev_alt / ca["n"] >= 0.6:
            v = f"variant_NOT_in_reads:assembly_error_in_{alt_h}"
        else:
            v = "mixed_reads_at_alt(heterozygous_or_collapsed_copies)"
        if cr["n"] >= 3 and cr["match"] / cr["n"] < 0.6:
            v += f";better_copy_{ref_h}_poorly_supported"
        row["verdict"] = v
        srows.append(row)
SCOLS = ["diff_index", "locus", "family", "classes", "better_hap", "alt_hap", "variant", "ref_contig_pos", "alt_contig_pos", "alt_strand", "ref_site", "alt_site",
         "alt_reads_with_better_allele", "ref_reads_with_alt_allele", "verdict"]
with open(vdir / "hifi_sites.tsv", "w") as f:
    f.write("\t".join(SCOLS) + "\n")
    for r in srows:
        f.write("\t".join(str(r.get(c, "-")) for c in SCOLS) + "\n")

# ---------- summary ----------
out = [f"{tag}: PAF records {n_in}, kept (primary, block>={MIN_BLOCK}, de<={MAX_DE}) {n_keep}"]
for w in W:
    out.append(f"  {w['name']} median depth {DEPTH[w['name']][2]:.1f} reads {len(recs.get(w['name'], []))}")
for kind in ("locus", "hap_only_locus", "absent_gap"):
    c = Counter(r["verdict"].split("(")[0].split(":")[0] for r in drows if r["kind"] == kind)
    out.append(f"{kind}: " + ", ".join(f"{k}={v}" for k, v in c.most_common()))
par = [r for r in drows if r["kind"] == "locus" and r["paralog_pident"] not in ("-", "") and float(r["paralog_pident"]) >= 98]
out.append(">=98% paralog loci: " + "; ".join(f'{r["hap"]} {r["id"]} {r["family"]} ratio={r.get("ratio_all")} mq20={r.get("ratio_mapq20")} span={r.get("spanning_reads_mapq20")} -> {r["verdict"].split("(")[0]}' for r in par))
for r in drows:
    if r["kind"] != "locus":
        out.append(f'  {r["kind"]} {r["hap"]} {r["id"]} {r["family"]} {r["status"]} len={r["length"]} ratio={r.get("ratio_all")} span={r.get("spanning_reads_mapq20")} -> {r["verdict"]}')
out.append("variant sites: " + ", ".join(f"{k}={v}" for k, v in Counter(r["verdict"].split(";")[0] for r in srows).most_common()))
for r in srows:
    out.append(f'  {r["locus"]} {r["family"]} {r["classes"]} {r["variant"]} | better {r["better_hap"]}: {r.get("ref_site")} | alt {r["alt_hap"]}: {r.get("alt_site")} -> {r["verdict"]}')
(vdir / "hifi_summary.txt").write_text("\n".join(out) + "\n")
print("\n".join(out))

#!/usr/bin/env python3
"""Additional analyses for the short report, all from tables already on disk.

Run: uv run --with scipy python manuscript/analysis_extra.py
Writes manuscript/derived/*.tsv and prints every number for the text.

1. Per-family |hap1 - hap2| locus-count differences in the two annotation passes.
2. Summed normalized depth and annotated locus counts, using per-locus eligibility.
3. Exploratory depth flags from a robust fit to log2(ratio) and nominal BH q values.
4. Confounder check: does the depth ratio track GC content or locus length?
5. Near-identical paralogs vs the rest: are the depth anomalies where copies are similar?
6. Locus-difference counts, pairing flags and descriptive Wilson intervals.
"""
from __future__ import annotations

import csv
import math
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import mannwhitneyu, norm, pearsonr

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "derived"
OUT.mkdir(exist_ok=True)
PAIRS = [("pair1", "Myotis sodalis"), ("pair2", "Rhinolophus clivosus"), ("pair3", "Eidolon dupreanum")]


def rd(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def wr(name: str, rows: list[dict], cols: list[str]) -> None:
    with open(OUT / name, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a proportion; stable at the extremes, unlike normal approximation."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def bh(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, in the input order."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = min(prev, pvals[i] * m / k)
        adj[i] = prev
    return adj


def read_fa(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


print("=" * 70)
print("1. OBSERVED HAPLOTYPE DIFFERENCES PER FAMILY (annotated loci)")
print("=" * 70)
floor_rows: list[dict] = []
for pair, sp in PAIRS:
    for mode, path in (("default", ROOT / pair / "out" / f"{pair}.family_summary.tsv"),
                       ("sensitive", ROOT / pair / "sens" / "out" / f"{pair}sens.family_summary.tsv")):
        if not path.exists():
            print(f"  missing: {path}")
            continue
        by_fam: dict[str, dict[str, int]] = defaultdict(dict)
        for r in rd(path):
            by_fam[r["family"]][r["hap"]] = int(r["loci"])
        for fam, h in sorted(by_fam.items()):
            a, b = h.get("hap1", 0), h.get("hap2", 0)
            floor_rows.append(dict(pair=pair, species=sp, mode=mode, family=fam, hap1=a, hap2=b,
                                   diff=abs(a - b), total=a + b))
wr("haplotype_count_differences.tsv", floor_rows, ["pair", "species", "mode", "family", "hap1", "hap2", "diff", "total"])
for pair, sp in PAIRS:
    for mode in ("default", "sensitive"):
        rows = [r for r in floor_rows if r["pair"] == pair and r["mode"] == mode]
        if not rows:
            continue
        diffs = [r["diff"] for r in rows]
        nz = [r for r in rows if r["diff"]]
        worst = sorted(nz, key=lambda r: -r["diff"])[:3]
        print(f"{sp:24s} {mode:10s} families={len(rows):3d} with a difference={len(nz):3d} "
              f"sum|hap1-hap2|={sum(diffs):3d} max={max(diffs)} "
              f"worst={', '.join(f'{r['family']} {r['hap1']}v{r['hap2']}' for r in worst)}")
print("\nFamilies detected in both annotation passes:")
for pair, sp in PAIRS:
    d = {r["family"]: r for r in floor_rows if r["pair"] == pair and r["mode"] == "default"}
    s = {r["family"]: r for r in floor_rows if r["pair"] == pair and r["mode"] == "sensitive"}
    flips = [(f, d[f]["diff"], s[f]["diff"]) for f in d if f in s and d[f]["diff"] != s[f]["diff"]]
    sign = [(f, d[f]["hap1"] - d[f]["hap2"], s[f]["hap1"] - s[f]["hap2"]) for f in d if f in s
            and (d[f]["hap1"] - d[f]["hap2"]) * (s[f]["hap1"] - s[f]["hap2"]) < 0]
    print(f"  {sp:24s} families whose difference changes between passes: {len(flips)}; "
          f"families where the SIGN flips: {len(sign)} {sign if sign else ''}")

print()
print("=" * 70)
print("2. SUMMED NORMALIZED DEPTH PER FAMILY vs ANNOTATED LOCUS COUNT")
print("=" * 70)
cn_rows: list[dict] = []
for pair, sp in PAIRS:
    rows = [r for r in rd(ROOT / "verify" / pair / "hifi_depth.tsv")
            if r["kind"] == "locus" and r["ratio_mapq20"] not in ("-", "")
            and r["window_median"] not in ("-", "") and float(r["window_median"]) >= 4]
    by = defaultdict(list)
    for r in rows:
        by[(r["family"], r["hap"])].append(r)
    for (fam, hap), rs in sorted(by.items()):
        assembled = len(rs)
        implied = sum(float(r["ratio_mapq20"]) for r in rs)
        cn_rows.append(dict(pair=pair, species=sp, family=fam, haplotype=hap, assembled_loci=assembled,
                            implied_copies=round(implied, 1), ratio=round(implied / assembled, 2)))
wr("copies_implied_vs_assembled.tsv", cn_rows,
   ["pair", "species", "family", "haplotype", "assembled_loci", "implied_copies", "ratio"])
for r in sorted(cn_rows, key=lambda r: -(r["implied_copies"] - r["assembled_loci"]))[:12]:
    print(f"  {r['species']:22s} {r['family']:8s} {r['haplotype']} assembled={r['assembled_loci']:3d} "
          f"implied={r['implied_copies']:6.1f} (x{r['ratio']})")
above = sum(r["implied_copies"] - r["assembled_loci"] > 2 for r in cn_rows)
below = sum(r["assembled_loci"] - r["implied_copies"] > 2 for r in cn_rows)
print(f"  summed normalized depth within two units: {len(cn_rows) - above - below} of {len(cn_rows)}")
print(f"  absolute deviation above two units: {above} excesses and {below} deficits")

print()
print("=" * 70)
print("3. EXPLORATORY DEPTH SCREEN (nominal BH q < 0.05)")
print("=" * 70)
null_rows: list[dict] = []
for pair, sp in PAIRS:
    rows = [r for r in rd(ROOT / "verify" / pair / "hifi_depth.tsv")
            if r["kind"] == "locus" and r["ratio_mapq20"] not in ("-", "")
            and r["window_median"] not in ("-", "") and float(r["window_median"]) >= 4]
    vals = [(r, float(r["ratio_mapq20"])) for r in rows]
    logs = [math.log2(v) for _, v in vals if v > 0]
    centre = st.median(logs)
    mad = st.median([abs(x - centre) for x in logs]) or 1e-9
    sigma = 1.4826 * mad
    ps, keep = [], []
    for r, v in vals:
        z = (math.log2(v) - centre) / sigma if v > 0 else -8.0
        ps.append(2 * norm.sf(abs(z)))
        keep.append((r, v, z))
    adj = bh(ps)
    sig = 0
    for (r, v, z), p, q in zip(keep, ps, adj):
        flag = q < 0.05
        sig += flag
        null_rows.append(dict(pair=pair, species=sp, hap=r["hap"], locus=r["id"], family=r["family"],
                              ratio_mapq20=v, log2=round(math.log2(v) if v > 0 else -8, 3), z=round(z, 2),
                              p=f"{p:.3g}", q=f"{q:.3g}", significant="yes" if flag else "no",
                              old_call=r["verdict"].split("(")[0].split(":")[0]))
    print(f"{sp:24s} loci tested={len(vals):3d} centre(log2)={centre:+.3f} robust sd={sigma:.3f} "
          f"screen flags (nominal q < 0.05)={sig}")
wr("depth_null_fdr.tsv", null_rows,
   ["pair", "species", "hap", "locus", "family", "ratio_mapq20", "log2", "z", "p", "q", "significant", "old_call"])
agree = Counter((r["significant"], r["old_call"]) for r in null_rows)
same_column = Counter(
    (r["significant"], 0.55 <= r["ratio_mapq20"] < 1.7) for r in null_rows
)
mixed_agree = sum(n for (sig, call), n in agree.items() if (sig == "no") == (call == "about_1x"))
matched_agree = sum(n for (sig, inside), n in same_column.items() if (sig == "no") == inside)
print(f"  fixed all-primary corridor vs MAPQ ≥20 BH labels: {mixed_agree}/{len(null_rows)}")
print(f"  fixed corridor and BH labels both on MAPQ ≥20: {matched_agree}/{len(null_rows)}")
print("  mixed-column calls (BH, all-primary label):", dict(agree))
print("  matched-column calls (BH, inside corridor):", dict(same_column))
for pair, _ in PAIRS:
    measured = [
        r for r in rd(ROOT / "verify" / pair / "hifi_depth.tsv")
        if r["kind"] == "locus" and r["depth_mapq20"] not in ("-", "")
    ]
    eligible = [r for r in null_rows if r["pair"] == pair]
    print(f"  {pair}: zero MAPQ ≥20 depth "
          f"{sum(float(r['depth_mapq20']) == 0 for r in measured)}/{len(measured)} measured, "
          f"{sum(r['ratio_mapq20'] == 0 for r in eligible)}/{len(eligible)} tested")

print()
print("=" * 70)
print("4. CONFOUNDERS: GC CONTENT AND LOCUS LENGTH vs DEPTH RATIO")
print("=" * 70)
conf_rows: list[dict] = []
for pair, sp in PAIRS:
    # hifi_depth.tsv stores the window as its region name, not the full "hap|region|contig:start-end" key
    wins = {w["region"]: w for w in rd(ROOT / "verify" / pair / "windows.tsv")}
    seqs: dict[str, str] = {}
    for h in ("hap1", "hap2"):
        seqs.update(read_fa(ROOT / "verify" / pair / f"windows_{h}.fa"))
    xs_gc, xs_len, ys = [], [], []
    for r in rd(ROOT / "verify" / pair / "hifi_depth.tsv"):
        if r["kind"] != "locus" or r["ratio_mapq20"] in ("-", "") or r["window_median"] in ("-", ""):
            continue
        if float(r["window_median"]) < 4 or r["window"] not in wins:
            continue
        w = wins[r["window"]]
        seq = seqs.get(w["name"], "")
        a, b = int(r["start"]) - int(w["start"]), int(r["end"]) - int(w["start"])
        sub = seq[max(0, a):max(0, b)].upper()
        if len(sub) < 200:
            continue
        gc = (sub.count("G") + sub.count("C")) / len(sub)
        v = float(r["ratio_mapq20"])
        if v <= 0:
            continue
        xs_gc.append(gc)
        xs_len.append(math.log10(len(sub)))
        ys.append(math.log2(v))
        conf_rows.append(dict(pair=pair, species=sp, locus=r["id"], family=r["family"], gc=round(gc, 4),
                              length=len(sub), log2_ratio=round(math.log2(v), 3)))
    if len(ys) > 10:
        rg, pg = pearsonr(xs_gc, ys)
        rl, pl = pearsonr(xs_len, ys)
        print(f"{sp:24s} n={len(ys):3d}  GC vs log2(ratio): r={rg:+.3f} p={pg:.3g}   "
              f"log10(length) vs log2(ratio): r={rl:+.3f} p={pl:.3g}")
wr("confounders_gc_length.tsv", conf_rows, ["pair", "species", "locus", "family", "gc", "length", "log2_ratio"])

print()
print("=" * 70)
print("5. NEAR-IDENTICAL PARALOGS vs THE REST")
print("=" * 70)
for pair, sp in PAIRS:
    hi, lo = [], []
    for r in rd(ROOT / "verify" / pair / "hifi_depth.tsv"):
        if r["kind"] != "locus" or r["ratio_mapq20"] in ("-", "") or r["window_median"] in ("-", ""):
            continue
        if float(r["window_median"]) < 4:
            continue
        v = float(r["ratio_mapq20"])
        pid = r["paralog_pident"]
        (hi if pid not in ("-", "") and float(pid) >= 95 else lo).append(v)
    if len(hi) >= 5 and len(lo) >= 5:
        u, p = mannwhitneyu(hi, lo, alternative="two-sided")
        print(f"{sp:24s} paralog>=95%: n={len(hi):3d} median ratio={st.median(hi):.2f} | "
              f"others: n={len(lo):3d} median={st.median(lo):.2f} | Mann-Whitney p={p:.3g}")
    else:
        print(f"{sp:24s} too few loci with a close paralog (n={len(hi)})")

print()
print("=" * 70)
print("6. WILSON 95% CONFIDENCE INTERVALS FOR THE PROPORTIONS IN THE TEXT")
print("=" * 70)
diff_counts = {p: Counter(r["refined_class"].split("(")[0].split(":")[0]
                          for r in rd(ROOT / "summary" / f"{p}.refined_differences.tsv")) for p, _ in PAIRS}
tot = sum(sum(c.values()) for c in diff_counts.values())
same_status = sum(c["artefact_sensitivity_or_threshold"] for c in diff_counts.values())
pairing_flags = Counter()
raw = Counter()
start_flags, rna_junctions, rna_coverage, atg = 0, 0, 0, 0
for pair, _ in PAIRS:
    for row in rd(ROOT / pair / "curation" / f"{pair}.difference_curation.tsv"):
        status = row["final_class"]
        if status.startswith("same_status"):
            pairing_flags["multiple"] += "multi_copy_landing" in status
            pairing_flags["outside"] += "counterpart_outside_paired_window" in status
            pairing_flags["either"] += ("multi_copy_landing" in status
                                        or "counterpart_outside_paired_window" in status)
    for row in rd(ROOT / pair / "out" / f"{pair}.haplotype_differences.tsv"):
        feature = row["feature"]
        if feature.startswith("locus_"):
            raw["locus"] += 1
        elif feature.startswith(("alignment_gap", "internal_indel")):
            raw["gap_or_indel"] += 1
        else:
            raw["other"] += 1
    flagged = {
        r["locus_id"] for r in rd(ROOT / pair / "curation" / f"{pair}.locus_curation.tsv")
        if r["pilot_status"] == "intact" and r["start_codon"] not in ("ATG", "-", "")
    }
    start_flags += len(flagged)
    atg += sum(r["locus_id"] in flagged and r["bat_verdict"].split("(")[0] == "ATG_from_residue_1"
               for r in rd(ROOT / "verify" / pair / "refseq_start.tsv"))
    for r in rd(ROOT / "verify" / pair / "rna_loci.tsv"):
        if r["locus_id"] in flagged:
            rna_junctions += r["rna_verdict"].startswith("upstream_exon_spliced_into_model_start")
            rna_coverage += r["rna_verdict"] == "transcript_extends_upstream_of_model_start"
print(f"  raw difference rows: {sum(raw.values())} total = {raw['locus']} locus + "
      f"{raw['gap_or_indel']} gap/indel + {raw['other']} other")
print(f"  same family/status after re-search: {same_status}/{tot}, "
      f"multiple landings {pairing_flags['multiple']}, outside paired window {pairing_flags['outside']}, "
      f"either flag {pairing_flags['either']}")
all_inside = sum(r["old_call"] == "about_1x" for r in null_rows)
for label, k, n in [("rows with same family/status on re-search", same_status, tot),
                    ("curated difference rows supported by reads", 14, 15),
                    ("all-primary depth ratios inside the fixed corridor", all_inside, len(null_rows)),
                    ("agreement with same MAPQ ≥20 ratios", matched_agree, len(null_rows)),
                    ("agreement with mixed depth columns", mixed_agree, len(null_rows)),
                    ("start-codon flags with upstream junctions", rna_junctions, start_flags),
                    ("start-codon flags with separate RNA coverage", rna_coverage, start_flags),
                    ("start-codon flags with RefSeq ATG at residue 1", atg, start_flags)]:
    lo, up = wilson(k, n)
    print(f"  {label:48s} {k:3d}/{n:3d} = {100 * k / n:5.1f}%  95% CI {100 * lo:.1f}-{100 * up:.1f}%")

print()
print("Tables written to", OUT)

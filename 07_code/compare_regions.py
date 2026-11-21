#!/usr/bin/env python3
"""Compare hap1 vs hap2 for each region window: extract (+-100 kb), minimap2 -x asm10 -c, gaps/inversions/liftover.
usage: compare_regions.py base_dir pair_tag hap1.fna hap2.fna acc1 acc2 seqrep1 seqrep2"""
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

base, pair = Path(sys.argv[1]), sys.argv[2]
fna = {"hap1": sys.argv[3], "hap2": sys.argv[4]}
seqrep = {"hap1": Path(sys.argv[7]), "hap2": Path(sys.argv[8])}
SAMTOOLS = str(base / "bin" / "samtools")
out = base / "out"
tmp = base / "mapping" / "cmp_tmp"
tmp.mkdir(exist_ok=True)
FLANK, GAPMIN, EDGE, MININV = 100_000, 5_000, 200_000, 5_000
PRESET = "asm10"

chrname = {}
asm_contig = {}  # (hap, chr_name) -> accession of the assembled chromosome
for h in seqrep:
    for r in json.loads(seqrep[h].read_text())["reports"]:
        chrname[(h, r["genbank_accession"])] = r.get("chr_name", "?")
        if r.get("role") == "assembled-molecule":
            asm_contig[(h, r.get("chr_name", "?"))] = r["genbank_accession"]


def read_tsv(p):
    L = Path(p).read_text().splitlines()
    k = L[0].split("\t")
    return [dict(zip(k, l.split("\t"))) for l in L[1:]]


regs = {h: read_tsv(out / f"{pair}.{h}.regions.tsv") for h in ("hap1", "hap2")}
loci = {h: read_tsv(out / f"{pair}.{h}.loci.tsv") for h in ("hap1", "hap2")}


def fetch(h, contig, s, e):
    """1-based inclusive; returns upper-case sequence"""
    r = subprocess.run([SAMTOOLS, "faidx", fna[h], f"{contig}:{s}-{e}"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


def extract(h, reg):
    clen = int(reg["contig_len"])
    s = max(1, int(reg["start"]) - FLANK)
    e = min(clen, int(reg["end"]) + FLANK)
    seq = fetch(h, reg["contig"], s, e)
    p = tmp / f'{reg["region_id"]}.fa'
    p.write_text(f'>{reg["region_id"]}\n' + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
    return p, s, e


def edge_flags(h, reg):
    clen = int(reg["contig_len"])
    s, e = int(reg["start"]), int(reg["end"])
    flags = []
    if s - EDGE < 1:
        flags.append(f"within200kb_of_contig_start({s - 1}bp)")
    if e + EDGE > clen:
        flags.append(f"within200kb_of_contig_end({clen - e}bp)")
    seq = fetch(h, reg["contig"], max(1, s - EDGE), min(clen, e + EDGE))
    runs = [(m.start() + max(1, s - EDGE), m.end() - 1 + max(1, s - EDGE)) for m in re.finditer(r"N+", seq)]
    if runs:
        flags.append("assembly_gap(N-run)_within200kb:" + ";".join(f"{a}-{b}" for a, b in runs[:5]))
    return flags


def run_mm(ref, qry):
    r = subprocess.run(["minimap2", "-x", PRESET, "-c", "-t", "2", str(ref), str(qry)], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"minimap2 failed: {r.stderr[-500:]}")
    recs = []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        tags = {t[:2]: t[5:] for t in f[12:]}
        recs.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": f[5], "tlen": int(f[6]),
                     "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "al": int(f[10]), "mapq": int(f[11]), "cg": tags.get("cg", "")})
    return recs, r.stderr


def union(iv):
    iv = sorted(iv)
    res = []
    for s, e in iv:
        if res and s <= res[-1][1]:
            res[-1][1] = max(res[-1][1], e)
        else:
            res.append([s, e])
    return res


def uncovered(cov, lo, hi):
    gaps, pos = [], lo
    for s, e in cov:
        if s > pos:
            gaps.append((pos, s))
        pos = max(pos, e)
    if hi > pos:
        gaps.append((pos, hi))
    return gaps


def blocks_of(rec):
    """list of (t0,t1,q0,q1,strand) aligned blocks (0-based half-open) + internal indels >= GAPMIN"""
    t, q = rec["ts"], (rec["qs"] if rec["strand"] == "+" else rec["qe"])
    bl, indels = [], []
    for n, op in re.findall(r"(\d+)([MIDNSH=X])", rec["cg"]):
        n = int(n)
        if op in "M=X":
            if rec["strand"] == "+":
                bl.append((t, t + n, q, q + n, "+"))
                q += n
            else:
                bl.append((t, t + n, q - n, q, "-"))
                q -= n
            t += n
        elif op == "D":
            if n >= GAPMIN:
                indels.append(("deletion_in_query", t, t + n, q))
            t += n
        elif op == "I":
            if n >= GAPMIN:
                indels.append(("insertion_in_query", t, t, q if rec["strand"] == "+" else q - n, n))
            q += n if rec["strand"] == "+" else -n
    return bl, indels


def lift(blocks, x, y, to_query=True):
    """map [x,y) from target->query (or query->target). returns (covered_fraction, lo, hi) in other coordinate"""
    cov, los, his = 0, [], []
    for t0, t1, q0, q1, st in blocks:
        a0, a1, b0, b1 = (t0, t1, q0, q1) if to_query else (q0, q1, t0, t1)
        s, e = max(x, a0), min(y, a1)
        if s >= e:
            continue
        cov += e - s
        if (st == "+"):
            los.append(b0 + (s - a0))
            his.append(b0 + (e - a0))
        else:
            los.append(b1 - (e - a0))
            his.append(b1 - (s - a0))
    if not los:
        return 0.0, None, None
    return cov / max(1, y - x), min(los), max(his)


# ---- 1. extract all windows ----
ext = {}
for h in ("hap1", "hap2"):
    for r in regs[h]:
        ext[r["region_id"]] = (h, r, *extract(h, r))

# ---- 2. pair windows of the same type, one-to-one, by aligned bp (min 20 kb and >=25% of the shorter extracted window) ----
MIN_FRAC = 0.25
bp_cache = {}


def aligned_bp(r1, r2):
    k = (r1["region_id"], r2["region_id"])
    if k not in bp_cache:
        recs, _ = run_mm(ext[k[0]][2], ext[k[1]][2])
        bp_cache[k] = sum(x["al"] for x in recs)
    return bp_cache[k]


def ext_len(r):
    return ext[r["region_id"]][4] - ext[r["region_id"]][3] + 1


assigned1, assigned2 = {}, {}
for stage in ("same_chr", "any_chr"):
    cands = []
    for r1 in regs["hap1"]:
        if r1["region_id"] in assigned1:
            continue
        for r2 in regs["hap2"]:
            if r2["region_id"] in assigned2 or r2["type"] != r1["type"]:
                continue
            same = chrname[("hap2", r2["contig"])] == chrname[("hap1", r1["contig"])]
            if (stage == "same_chr") != same:
                continue
            cands.append((aligned_bp(r1, r2), r1, r2))
    cands.sort(key=lambda c: -c[0])
    for bp, r1, r2 in cands:
        if r1["region_id"] in assigned1 or r2["region_id"] in assigned2:
            continue
        if bp >= max(20_000, MIN_FRAC * min(ext_len(r1), ext_len(r2))):
            assigned1[r1["region_id"]] = (r2, bp)
            assigned2[r2["region_id"]] = r1
pairs = [(r1, *assigned1.get(r1["region_id"], (None, 0))) for r1 in regs["hap1"]]
unpaired2 = [r for r in regs["hap2"] if r["region_id"] not in assigned2]

def landing(src_hap, src_reg, dst_hap):
    """for a window with no counterpart window: align its extraction to the same-named chromosome of the other haplotype"""
    dst_c = asm_contig.get((dst_hap, chrname[(src_hap, src_reg["contig"])]))
    if dst_c is None:
        return "no same-named chromosome in other haplotype"
    clen = int(next(l.split("\t")[1] for l in Path(fna[dst_hap] + ".fai").read_text().splitlines() if l.split("\t")[0] == dst_c))
    chrfa = tmp / f"{dst_hap}_{dst_c}.fa"
    with chrfa.open("wb") as fh:
        subprocess.run([SAMTOOLS, "faidx", fna[dst_hap], dst_c], stdout=fh, check=True)
    r = subprocess.run(["minimap2", "-x", PRESET, "-t", "2", str(chrfa), str(ext[src_reg["region_id"]][2])], capture_output=True, text=True, check=True)
    chrfa.unlink()
    recs = [l.split("\t") for l in r.stdout.splitlines()]
    if not recs:
        return f"no alignment to {dst_hap} chromosome {dst_c} ({clen} bp)"
    best = max(recs, key=lambda f: int(f[10]))
    ts, te = int(best[7]), int(best[8])
    near = sum(int(f[10]) for f in recs if int(f[7]) < te + 100_000 and int(f[8]) > ts - 100_000)
    dl = [l for l in loci[dst_hap] if l["contig"] == dst_c and int(l["start"]) <= te + 20_000 and int(l["end"]) >= ts - 20_000]
    return (f"best alignment {dst_hap} {dst_c}:{ts + 1}-{te} ({near} aligned bp within 100 kb of it); "
            f"{dst_hap} loci there: " + (",".join(f'{l["locus_id"]}:{l["family"]}/{l["status"]}' for l in dl) or "none"))



LOCFLANK = 3_000


def resolve(src_hap, loc, dst_fa, dst_start, dst_loci, skip):
    """second-chance check for a locus left unmatched by the window-level alignment: align the locus (+-3 kb) directly to the
    other haplotype's window. returns (covered_fraction, lifted_lo, lifted_hi, partner_or_None) in genomic coordinates of the other haplotype"""
    cs = int(loc["start"])
    ce = int(loc["end"])
    clen = int(next(l.split("\t")[1] for l in Path(fna[src_hap] + ".fai").read_text().splitlines() if l.split("\t")[0] == loc["contig"]))
    a, b = max(1, cs - LOCFLANK), min(clen, ce + LOCFLANK)
    qp = tmp / f'loc_{loc["locus_id"]}.fa'
    qp.write_text(f'>{loc["locus_id"]}\n' + "\n".join(fetch(src_hap, loc["contig"], a, b)[i:i + 80] for i in range(0, b - a + 1, 80)) + "\n")
    recs, _ = run_mm(dst_fa, qp)
    qp.unlink()
    ivs, tl, th = [], [], []
    for x in recs:
        bl_, _ = blocks_of(x)
        for t0, t1, q0, q1, st in bl_:
            o0, o1 = max(q0, cs - a), min(q1, ce - a + 1)
            if o1 > o0:
                ivs.append((o0, o1))
                if st == "+":
                    tl.append(t0 + (o0 - q0))
                    th.append(t0 + (o1 - q0))
                else:
                    tl.append(t1 - (o1 - q0))
                    th.append(t1 - (o0 - q0))
    cov = sum(e - s_ for s_, e in union(ivs)) / max(1, ce - cs + 1)
    if not tl:
        return cov, None, None, None
    lo, hi = min(tl) + dst_start, max(th) + dst_start
    partner = None
    for m in dst_loci:
        if m["locus_id"] in skip:
            continue
        if int(m["start"]) < hi + 5_000 and int(m["end"]) > lo - 5_000:
            partner = m
            break
    return cov, lo, hi, partner


# ---- 3. detailed comparison ----
diff_rows, pair_rows = [], []


def add(rid1, rid2, kind, h, contig, s, e, note):
    diff_rows.append([rid1, rid2, kind, h, contig, chrname.get((h if h in ("hap1", "hap2") else "hap1", contig), "?"), s, e, note])


for r1, r2, bp in pairs:
    id1 = r1["region_id"]
    h1, _, p1, s1, e1 = ext[id1]
    f1 = edge_flags("hap1", r1)
    if r2 is None:
        pair_rows.append([id1, "", r1["type"], r1["contig"], chrname[("hap1", r1["contig"])], s1, e1, e1 - s1 + 1] + [""] * 15 + ["|".join(f1) or "-", "", "no hap2 window aligns"])
        add(id1, "", "region_without_hap2_counterpart", "hap1", r1["contig"], r1["start"], r1["end"], "no hap2 window of same type aligned; " + landing("hap1", r1, "hap2"))
        continue
    id2 = r2["region_id"]
    _, _, p2, s2, e2 = ext[id2]
    f2 = edge_flags("hap2", r2)
    recs, err = run_mm(p1, p2)  # target = hap1 region, query = hap2 region
    L1, L2 = e1 - s1 + 1, e2 - s2 + 1
    cov1 = union([(x["ts"], x["te"]) for x in recs])
    cov2 = union([(x["qs"], x["qe"]) for x in recs])
    bl, indels = [], []
    for x in recs:
        b, i = blocks_of(x)
        bl += b
        indels += [(k, ) + tuple(v) for k, *v in i]
    # gaps > 5 kb uncovered between first and last aligned base
    g1 = [g for g in uncovered(cov1, cov1[0][0], cov1[-1][1]) if g[1] - g[0] >= GAPMIN] if cov1 else []
    g2 = [g for g in uncovered(cov2, cov2[0][0], cov2[-1][1]) if g[1] - g[0] >= GAPMIN] if cov2 else []
    # unaligned flanks / ends count as gaps too when region is mostly unaligned
    lead1, trail1 = (cov1[0][0], L1 - cov1[-1][1]) if cov1 else (L1, 0)
    for a, b in g1:
        add(id1, id2, "alignment_gap_hap1_only(>=5kb; hap1 sequence absent from hap2 alignment)", "hap1", r1["contig"], s1 + a, s1 + b, f"{b - a} bp")
    for a, b in g2:
        add(id1, id2, "alignment_gap_hap2_only(>=5kb; hap2 sequence absent from hap1 alignment)", "hap2", r2["contig"], s2 + a, s2 + b, f"{b - a} bp")
    for k, *v in indels:
        if k == "deletion_in_query":
            add(id1, id2, "internal_indel_in_alignment(>=5kb): hap1 sequence missing in hap2", "hap1", r1["contig"], s1 + v[0], s1 + v[1], f"{v[1] - v[0]} bp")
        else:
            add(id1, id2, "internal_indel_in_alignment(>=5kb): hap2 sequence missing in hap1", "hap2", r2["contig"], s2 + v[2], s2 + v[2] + v[3], f"{v[3]} bp")
    # strand structure: use alignments >= MININV aligned bp
    big = sorted([x for x in recs if x["al"] >= MININV], key=lambda x: x["ts"])
    plus = sum(x["al"] for x in big if x["strand"] == "+")
    minus = sum(x["al"] for x in big if x["strand"] == "-")
    dom = "+" if plus >= minus else "-"
    n_switch = sum(1 for a, b in zip(big, big[1:]) if a["strand"] != b["strand"])
    inv = [x for x in big if x["strand"] != dom]
    for x in inv:
        add(id1, id2, f"strand_switch_block(opposite to dominant orientation {dom})", "hap1", r1["contig"], s1 + x["ts"], s1 + x["te"], f'hap2 {s2 + x["qs"]}-{s2 + x["qe"]} strand {x["strand"]}; aligned {x["al"]} bp')
    # liftover of loci
    l1 = [l for l in loci["hap1"] if id1 in l["region_ids"].split(",")]
    l2 = [l for l in loci["hap2"] if id2 in l["region_ids"].split(",")]
    # one-to-one locus matching: lift each hap1 locus to hap2 coordinates, score overlap (Jaccard, +0.01 if same family), greedy by score
    cand, lifted1 = [], {}
    for l in l1:
        x, y = int(l["start"]) - s1, int(l["end"]) - s1 + 1
        fr, lo, hi = lift(bl, x, y, True)
        lifted1[l["locus_id"]] = (fr, lo, hi)
        if lo is None:
            continue
        for m in l2:
            ms_, me_ = int(m["start"]) - s2, int(m["end"]) - s2 + 1
            ov = min(hi, me_) - max(lo, ms_)
            if ov > 0:
                jac = ov / (max(hi, me_) - min(lo, ms_))
                cand.append((jac + (0.01 if m["family"] == l["family"] else 0), l["locus_id"], m["locus_id"]))
    cand.sort(reverse=True)
    pair_of, used = {}, set()
    for sc, a_, b_ in cand:
        if a_ in pair_of or b_ in used or sc < 0.05:
            continue
        pair_of[a_] = b_
        used.add(b_)
    by2 = {m["locus_id"]: m for m in l2}
    n_both_same = 0
    for l in l1:
        fr, lo, hi = lifted1[l["locus_id"]]
        if l["locus_id"] in pair_of:
            m = by2[pair_of[l["locus_id"]]]
            if m["family"] != l["family"] or m["status"] != l["status"]:
                add(id1, id2, "locus_in_both_but_call_differs", "both", l["contig"], l["start"], l["end"],
                    f'{l["locus_id"]} {l["family"]}/{l["status"]} vs {m["locus_id"]} {m["family"]}/{m["status"]}')
            else:
                n_both_same += 1
        else:
            cov_, lo_, hi_, partner = resolve("hap1", l, p2, s2, l2, used)
            if partner is not None:
                pair_of[l["locus_id"]] = partner["locus_id"]
                used.add(partner["locus_id"])
                if partner["family"] != l["family"] or partner["status"] != l["status"]:
                    add(id1, id2, "locus_in_both_but_call_differs", "both", l["contig"], l["start"], l["end"],
                        f'{l["locus_id"]} {l["family"]}/{l["status"]} vs {partner["locus_id"]} {partner["family"]}/{partner["status"]} (matched by direct locus alignment)')
                else:
                    n_both_same += 1
                continue
            why = (f"locus not alignable to hap2 window (direct locus alignment covers {cov_:.0%} of span)" if cov_ < 0.2 else
                   f"hap2 sequence aligned (direct alignment covers {cov_:.0%} of span, hap2 {lo_}-{hi_}) but no hap2 locus called there")
            add(id1, id2, "locus_hap1_only", "hap1", l["contig"], l["start"], l["end"], f'{l["locus_id"]} {l["family"]} {l["status"]}; {why}')
    for m in l2:
        if m["locus_id"] in used:
            continue
        x, y = int(m["start"]) - s2, int(m["end"]) - s2 + 1
        fr, lo, hi = lift(bl, x, y, False)
        cov_, lo_, hi_, partner = resolve("hap2", m, p1, s1, [l for l in l1 if l["locus_id"] not in pair_of], set())
        if partner is not None:
            pair_of[partner["locus_id"]] = m["locus_id"]
            used.add(m["locus_id"])
            if partner["family"] != m["family"] or partner["status"] != m["status"]:
                add(id1, id2, "locus_in_both_but_call_differs", "both", partner["contig"], partner["start"], partner["end"],
                    f'{partner["locus_id"]} {partner["family"]}/{partner["status"]} vs {m["locus_id"]} {m["family"]}/{m["status"]} (matched by direct locus alignment)')
            else:
                n_both_same += 1
            continue
        why = (f"locus not alignable to hap1 window (direct locus alignment covers {cov_:.0%} of span)" if cov_ < 0.2 else
               f"hap1 sequence aligned (direct alignment covers {cov_:.0%} of span, hap1 {lo_}-{hi_}) but no hap1 locus called there")
        add(id1, id2, "locus_hap2_only", "hap2", m["contig"], m["start"], m["end"], f'{m["locus_id"]} {m["family"]} {m["status"]}; {why}')
    pair_rows.append([id1, id2, r1["type"], r1["contig"], chrname[("hap1", r1["contig"])], s1, e1, L1, r2["contig"], chrname[("hap2", r2["contig"])], s2, e2, L2,
                      f"{sum(b - a for a, b in cov1) / L1:.3f}", f"{sum(b - a for a, b in cov2) / L2:.3f}", len(g1), len(g2), len(indels), len(inv), n_switch,
                      len(l1), len(l2), n_both_same, "|".join(f1) or "-", "|".join(f2) or "-", f"minimap2 -x {PRESET} -c; n_records={len(recs)}; dominant_strand={dom}; hap2_vs_hap1"])
for r2 in unpaired2:
    _, _, _, s2, e2 = ext[r2["region_id"]]
    pair_rows.append(["", r2["region_id"], r2["type"]] + [""] * 5 + [r2["contig"], chrname[("hap2", r2["contig"])], s2, e2, e2 - s2 + 1] + [""] * 10 + ["", "|".join(edge_flags("hap2", r2)) or "-", "no hap1 window aligns"])
    add("", r2["region_id"], "region_without_hap1_counterpart", "hap2", r2["contig"], r2["start"], r2["end"], "no hap1 window of same type aligned; " + landing("hap2", r2, "hap1"))

hdr = "region_hap1\tregion_hap2\ttype\thap1_contig\thap1_chr\thap1_extract_start\thap1_extract_end\thap1_len\thap2_contig\thap2_chr\thap2_extract_start\thap2_extract_end\thap2_len\thap1_frac_aligned\thap2_frac_aligned\tgaps_ge5kb_hap1\tgaps_ge5kb_hap2\tinternal_indels_ge5kb\tinverted_blocks\tstrand_switches\tloci_hap1\tloci_hap2\tloci_matched_identical_call\tflags_hap1\tflags_hap2\tmethod_note\n"
with (out / f"{pair}.region_pairs.tsv").open("w") as f:
    f.write(hdr)
    for r in pair_rows:
        f.write("\t".join(map(str, r + [""] * (26 - len(r)))) + "\n")
with (out / f"{pair}.haplotype_differences.tsv").open("w") as f:
    f.write("region_hap1\tregion_hap2\tfeature\thaplotype\tcontig\tchromosome\tstart\tend\tnote\n")
    for r in diff_rows:
        f.write("\t".join(map(str, r)) + "\n")
print(f"pairs={len(pair_rows)} diff_rows={len(diff_rows)} unpaired_hap2={len(unpaired2)}")

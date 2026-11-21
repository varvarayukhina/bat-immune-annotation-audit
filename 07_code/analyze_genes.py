#!/usr/bin/env python3
"""Parse miniprot GFF -> hits, loci (clustered by CDS overlap), status, regions (3 Mb windows).
usage: analyze_genes.py base_dir pair_tag hapLabel acc seqreport.json
Writes out/<pair_tag>.<hap>.{hits,loci,regions,families}.tsv"""
import json
import sys
from collections import defaultdict
from pathlib import Path

base, pair, hap, acc, seqrep = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], Path(sys.argv[5])
gff = base / "mapping" / f"{hap}.miniprot.gff"
out = base / "out"

# ---- thresholds (all stated in methods) ----
OUTS = 0.3          # keep hit if ms >= OUTS * best ms of that query over the whole genome (mimics --outs on a whole-genome run)
MIN_POS = 0.45      # min Positive fraction for a "supported" hit
MIN_COV = 0.30      # min fraction of query aligned for a "supported" hit
INTACT_COV = 0.90   # full-length: >=90% of query aligned
WINDOW = 3_000_000  # merge loci within 3 Mb into a region window
MAX_INTRON_OK = 15_000  # hits with an intron longer than this are excluded from locus building (miniprot default max is 200 kb; these gene families have short introns)
MAX_LOCUS_SPAN = 40_000 # a clustered locus longer than this is flagged uncertain (possible merge of adjacent genes)
MIN_IDENT_OK = 0.30     # evidence hit with identity below this -> uncertain

# ---- query table -> family ----
qinfo = {}
for line in (base / "queries" / "queries_final_table.tsv").read_text().splitlines()[1:]:
    f = line.split("\t")
    tag = f"{f[0]}|{f[1]}|{f[2]}|{f[3]}"
    qinfo[tag] = {"symbol": f[1], "queried": f[6].split(",")[0]}


def family_of(sym: str) -> str:
    s = sym.upper()
    for pre, fam in (("KLRC", "KLRC"), ("KLRD", "KLRD1"), ("KLRK", "KLRK1"), ("KLRB", "KLRB1"), ("KLRF", "KLRF"),
                     ("KLRG1", "KLRG1"), ("KLRG2", "KLRG2"), ("KLRA", "KLRA"), ("CLEC2D", "CLEC2D"), ("CLEC12A", "CLEC12A"),
                     ("CLEC4E", "CLEC4E"), ("KIR", "KIR"), ("LILR", "LILR"), ("NCR1", "NCR1"), ("NCR2", "NCR2"),
                     ("NCR3", "NCR3"), ("LAIR1", "LAIR1"), ("FCAR", "FCAR"), ("HLA-", "MHC-I"), ("H2-", "MHC-I"),
                     ("MR1", "MR1"), ("CD1", "CD1"), ("B2M", "B2M"), ("TAP1", "TAP1"), ("TAP2", "TAP2")):
        if s.startswith(pre):
            return fam
    raise SystemExit(f"no family for query symbol {sym}")


for t, d in qinfo.items():
    d["family"] = family_of(d["queried"])
REGION_OF = {}
for fam in "KLRC KLRD1 KLRK1 KLRB1 KLRF KLRG1 KLRA CLEC2D CLEC12A CLEC4E".split():
    REGION_OF[fam] = "NKC"
for fam in "LILR KIR NCR1 LAIR1".split():
    REGION_OF[fam] = "LRC"
for fam in "MHC-I MR1 CD1 B2M TAP1 TAP2".split():
    REGION_OF[fam] = "MHC-I"
# KLRG2 (not a CTLD NK receptor), FCAR, NCR2, NCR3 are queried but do not define a window (spec); they join a window only by position.

# ---- chromosome names ----
chrname = {}
for r in json.loads(seqrep.read_text())["reports"]:
    chrname[r["genbank_accession"]] = f'{r.get("chr_name", "?")}:{r.get("role", "?")}'

# ---- parse GFF ----
hits = []
cur = None
for line in gff.read_text().splitlines():
    if line.startswith("##PAF"):
        f = line.split("\t")
        tags = {t.split(":")[0]: t.split(":", 2)[2] for t in f[13:]}
        cur = {"query": f[1], "qlen": int(f[2]), "qs": int(f[3]), "qe": int(f[4]), "strand": f[5],
               "contig": f[6], "clen": int(f[7]), "AS": int(tags["AS"]), "np": int(tags["np"]),
               "fs": int(tags["fs"]), "st": int(tags["st"]), "cds": [], "stop": False}
        hits.append(cur)
    elif line.startswith("#") or not line.strip():
        continue
    else:
        f = line.split("\t")
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            cur["id"], cur["ms"] = at["ID"], int(f[5])
            cur["rank"], cur["ident"], cur["pos"] = int(at["Rank"]), float(at["Identity"]), float(at["Positive"])
            cur["start"], cur["end"] = int(f[3]), int(f[4])
        elif f[2] == "CDS":
            cur["cds"].append((int(f[3]), int(f[4])))
        elif f[2] == "stop_codon":
            cur["stop"] = True
for h in hits:
    c = sorted(h["cds"])
    h["max_intron"] = max([b[0] - a[1] - 1 for a, b in zip(c, c[1:])] or [0])
    h["cds_bp"] = sum(e - s_ + 1 for s_, e in c)
    h["cov"] = (h["qe"] - h["qs"]) / h["qlen"]
    h["family"] = qinfo[h["query"]]["family"]
    h["qsym"] = h["query"].split("|")[1]
    h["taxon"] = h["query"].split("|")[0]

best_ms = defaultdict(int)
for h in hits:
    best_ms[h["query"]] = max(best_ms[h["query"]], h["ms"])
for h in hits:
    h["rel"] = h["ms"] / best_ms[h["query"]]
    h["core"] = h["rel"] >= OUTS and h["pos"] >= MIN_POS and h["cov"] >= MIN_COV
    h["long_intron"] = h["max_intron"] > MAX_INTRON_OK
    h["supported"] = h["core"] and not h["long_intron"]
    dis = h["fs"] > 0 or h["st"] > 0
    h["disrupted"] = dis
    h["intact"] = h["supported"] and not dis and h["cov"] >= INTACT_COV and h["stop"] and h["max_intron"] <= MAX_INTRON_OK and h["ident"] >= MIN_IDENT_OK

# ---- hits table ----
with (out / f"{pair}.{hap}.hits.tsv").open("w") as f:
    f.write("hit_id\tcontig\tstart\tend\tstrand\tquery\tfamily\tqlen\tq_start\tq_end\tcov\tidentity\tpositive\tms\trel_to_best\tframeshifts\tinframe_stops\tmax_intron\tstop_codon_feature\tn_cds\tsupported\tlong_intron_excluded\tintact\n")
    for h in hits:
        f.write("\t".join(map(str, [h["id"], h["contig"], h["start"], h["end"], h["strand"], h["query"], h["family"], h["qlen"], h["qs"], h["qe"],
                                    f'{h["cov"]:.3f}', h["ident"], h["pos"], h["ms"], f'{h["rel"]:.3f}', h["fs"], h["st"], h["max_intron"], int(h["stop"]), len(h["cds"]),
                                    int(h["supported"]), int(h["core"] and h["long_intron"]), int(h["intact"])])) + "\n")

# ---- cluster supported hits into loci: single linkage on overlapping CDS exons, same contig and strand ----
# stage 1: hits with all introns <= MAX_INTRON_OK.  stage 2: long-intron hits that overlap no stage-1 CDS -> extra loci, always "uncertain".


def cluster(hl):
    par = list(range(len(hl)))

    def fnd(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    ivs = defaultdict(list)
    for i, h in enumerate(hl):
        for s_, e_ in h["cds"]:
            ivs[(h["contig"], h["strand"])].append((s_, e_, i))
    for lst in ivs.values():
        lst.sort()
        ce, ci = -1, None
        for s_, e_, i in lst:
            if ci is not None and s_ <= ce:
                par[fnd(i)] = fnd(ci)
                ce = max(ce, e_)
            else:
                ce, ci = e_, i
    gg = defaultdict(list)
    for i in range(len(hl)):
        gg[fnd(i)].append(hl[i])
    return list(gg.values())


sup = [h for h in hits if h["supported"]]
groups = [(g, "") for g in cluster(sup)]
cds1 = defaultdict(list)
for h in sup:
    for s_, e_ in h["cds"]:
        cds1[(h["contig"], h["strand"])].append((s_, e_))
for k in cds1:
    cds1[k].sort()


def overlaps_stage1(h):
    L = cds1.get((h["contig"], h["strand"]), [])
    return any(not (e1 < s2 or s1 > e2) for s1, e1 in L for s2, e2 in h["cds"] if abs(s1 - s2) < 400_000)


long_h = [h for h in hits if h["core"] and h["long_intron"] and not overlaps_stage1(h)]
groups += [(g, f"only_hits_with_intron>{MAX_INTRON_OK // 1000}kb") for g in cluster(long_h)]
loci = []
for g, flag in groups:
    best = max(g, key=lambda h: h["ms"])
    intact_hits = [h for h in g if h["intact"]]
    reason = ""
    if flag:
        status, ev, reason = "uncertain", best, flag
    elif intact_hits:
        status, ev = "intact", max(intact_hits, key=lambda h: h["ms"])
    elif any(h["disrupted"] and h["cov"] >= 0.5 for h in g) and not any((not h["disrupted"]) and h["cov"] >= 0.7 for h in g):
        status, ev = "pseudogene-like", max((h for h in g if h["disrupted"]), key=lambda h: h["ms"])
        reason = "frameshift/in-frame stop in alignment"
    else:
        status, ev = "partial", best
        rs = []
        if ev["cov"] < INTACT_COV:
            rs.append(f"query_cov<{INTACT_COV}")
        if not ev["stop"]:
            rs.append("no_terminal_stop_codon_in_alignment")
        if ev["ident"] < MIN_IDENT_OK:
            rs.append(f"identity<{MIN_IDENT_OK}")
            status = "uncertain"
        if ev["disrupted"]:
            rs.append("disrupted_but_low_cov")
        reason = ",".join(rs)
    span = max(e_ for h in g for _, e_ in h["cds"]) - min(s_ for h in g for s_, _ in h["cds"]) + 1
    if span > MAX_LOCUS_SPAN and status != "pseudogene-like":
        reason = (reason + "," if reason else "") + f"locus_span>{MAX_LOCUS_SPAN // 1000}kb(possible_merged_genes)"
        status = "uncertain"
    fam_best = defaultdict(int)
    for h in g:
        fam_best[h["family"]] = max(fam_best[h["family"]], h["ms"])
    ranked = sorted(fam_best.items(), key=lambda kv: -kv[1])
    second = ranked[1] if len(ranked) > 1 else ("", 0)
    ambiguous = int(second[1] >= 0.9 * ranked[0][1]) if second[0] else 0
    loci.append({
        "contig": best["contig"], "strand": best["strand"],
        "start": min(s for h in g for s, _ in h["cds"]), "end": max(e for h in g for _, e in h["cds"]),
        "family": best["family"], "best": best, "status": status, "ev": ev, "n_hits": len(g),
        "n_queries": len({h["query"] for h in g}), "second_family": second[0], "ambiguous": ambiguous,
        "queries": ",".join(sorted({h["qsym"] for h in g}))[:200], "reason": reason,
    })
loci.sort(key=lambda l: (l["contig"], l["start"]))
for i, l in enumerate(loci, 1):
    l["id"] = f"{hap}_L{i:03d}"
    l["region_type"] = REGION_OF.get(l["family"], "")

# ---- windows (3 Mb) per region type, built from defining families only ----
regions = []
for rtype in ("NKC", "LRC", "MHC-I"):
    mem = [l for l in loci if l["region_type"] == rtype]
    cur = []
    for l in mem:
        if cur and l["contig"] == cur[-1]["contig"] and l["start"] - max(x["end"] for x in cur) <= WINDOW:
            cur.append(l)
        else:
            if cur:
                regions.append((rtype, cur))
            cur = [l]
    if cur:
        regions.append((rtype, cur))
regions.sort(key=lambda r: (r[0], -len(r[1])))
reg_rows = []
counter = defaultdict(int)
for rtype, mem in regions:
    counter[rtype] += 1
    rid = f"{hap}_{rtype}_{counter[rtype]}"
    c = mem[0]["contig"]
    s, e = min(x["start"] for x in mem), max(x["end"] for x in mem)
    extra = [l for l in loci if l["contig"] == c and l["region_type"] == "" and l["start"] >= s and l["end"] <= e]
    allm = sorted(mem + extra, key=lambda l: l["start"])
    for l in allm:
        l.setdefault("regions", []).append(rid)
    reg_rows.append({"id": rid, "type": rtype, "contig": c, "chr": chrname.get(c, "?"), "start": s, "end": e, "n": len(allm),
                     "fams": ",".join(f"{k}x{v}" for k, v in sorted(defaultdict(int, {fam: sum(1 for l in allm if l['family'] == fam) for fam in {l['family'] for l in allm}}).items())),
                     "intact": sum(l["status"] == "intact" for l in allm), "partial": sum(l["status"] == "partial" for l in allm), "uncertain": sum(l["status"] == "uncertain" for l in allm),
                     "pseudo": sum(l["status"] == "pseudogene-like" for l in allm), "clen": mem[0]["best"]["clen"]})

with (out / f"{pair}.{hap}.loci.tsv").open("w") as f:
    f.write("haplotype\tassembly\tlocus_id\tregion_ids\tcontig\tchromosome\tstart\tend\tstrand\tfamily\tbest_query\tbest_ms\tbest_identity\tbest_positive\tbest_cov\tstatus\tevidence_query\tev_cov\tev_frameshifts\tev_inframe_stops\tev_stop_codon\tev_max_intron\tstatus_note\tn_hits\tn_queries\tsecond_family\tfamily_ambiguous\tqueries\n")
    for l in loci:
        b, ev = l["best"], l["ev"]
        f.write("\t".join(map(str, [hap, acc, l["id"], ",".join(l.get("regions", [])), l["contig"], chrname.get(l["contig"], "?"), l["start"], l["end"], l["strand"],
                                    l["family"], b["query"], b["ms"], b["ident"], b["pos"], f'{b["cov"]:.3f}', l["status"], ev["query"], f'{ev["cov"]:.3f}',
                                    ev["fs"], ev["st"], int(ev["stop"]), ev["max_intron"], l["reason"], l["n_hits"], l["n_queries"], l["second_family"], l["ambiguous"], l["queries"]])) + "\n")
with (out / f"{pair}.{hap}.regions.tsv").open("w") as f:
    f.write("region_id\ttype\tcontig\tchromosome\tstart\tend\tspan_bp\tcontig_len\tn_loci\tfamilies\tintact\tpartial\tpseudogene_like\tuncertain\n")
    for r in reg_rows:
        f.write("\t".join(map(str, [r["id"], r["type"], r["contig"], r["chr"], r["start"], r["end"], r["end"] - r["start"] + 1, r["clen"], r["n"], r["fams"], r["intact"], r["partial"], r["pseudo"], r["uncertain"]])) + "\n")
fam_tab = defaultdict(lambda: defaultdict(int))
for l in loci:
    fam_tab[l["family"]]["loci"] += 1
    fam_tab[l["family"]][l["status"]] += 1
with (out / f"{pair}.{hap}.families.tsv").open("w") as f:
    f.write("haplotype\tfamily\tregion_type\tn_loci\tintact\tpartial\tpseudogene_like\tuncertain\n")
    for fam, d in sorted(fam_tab.items()):
        f.write(f"{hap}\t{fam}\t{REGION_OF.get(fam, '')}\t{d['loci']}\t{d['intact']}\t{d['partial']}\t{d['pseudogene-like']}\t{d['uncertain']}\n")
print(f"{hap}: hits={len(hits)} supported={len(sup)} loci={len(loci)} regions={len(reg_rows)}")

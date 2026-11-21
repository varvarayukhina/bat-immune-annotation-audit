#!/usr/bin/env python3
"""Curation pass for one hap1/hap2 pair ("manual checks" of the automated pilot), run while the genome FASTA files still exist.
usage: curate.py work_dir pair_tag

Reads  work_dir/out/<tag>.hap{1,2}.loci.tsv, <tag>.haplotype_differences.tsv, <tag>.region_pairs.tsv,
       work_dir/genomes/hap{1,2}/*_genomic.fna (+ .fai), <root>/queries/queries_final.faa, <root>/ref/swissprot_hm.dmnd
Writes work_dir/curation/
  <tag>.locus_curation.tsv       one row per locus: re-alignment on the locus region (+-5 kb, all 99 queries, default settings),
                                 gene-model checks (start codon, splice-site dinucleotides, terminal stop, internal stops),
                                 back-search of the predicted protein against Swiss-Prot human+mouse (family verification),
                                 identity to the closest paralog of the same family in the same haplotype
  <tag>.locus_proteins.faa       predicted proteins (miniprot --trans) used for the back-search
  <tag>.backsearch.tsv           diamond blastp hits (top 5 per locus)
  <tag>.difference_curation.tsv  one row per locus-level haplotype difference, re-examined directly:
                                 locus +-10 kb mapped to the WHOLE other haplotype (minimap2 asm20), counterpart region re-searched with
                                 miniprot at default and at high sensitivity (-k5 -M0 -n2 -l4), variants inside the CDS of the better copy
  <tag>.gap_curation.tsv         one row per >=5 kb haplotype-specific gap/indel: where that sequence maps in the other haplotype
  <tag>.curation_summary.txt
All thresholds for "supported / intact / pseudogene-like / partial" are the ones of analyze_genes.py.
"""
import csv
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

work, tag = Path(sys.argv[1]).resolve(), sys.argv[2]
root = work.parent
out = work / "out"
cur = work / "curation"
tmp = cur / "tmp"
logs = cur / "logs"
for d in (cur, tmp, logs):
    d.mkdir(parents=True, exist_ok=True)
SAMTOOLS = str(root / "bin" / "samtools")
QUERIES = root / "queries" / "queries_final.faa"
DMND = root / "ref" / "swissprot_hm.dmnd"
HAPS = ("hap1", "hap2")
OTHER = {"hap1": "hap2", "hap2": "hap1"}
fna = {h: str(next((work / "genomes" / h).glob("*_genomic.fna"))) for h in HAPS}
fai = {h: {l.split("\t")[0]: int(l.split("\t")[1]) for l in Path(fna[h] + ".fai").read_text().splitlines()} for h in HAPS}
LOCFLANK, DIFFFLANK = 5_000, 10_000
# as the pilot, except -N/--outn raised: the pilot's cap of 100 secondary alignments applied per 250 Mb chunk (9 chunks), one small-region run
# over ~100 regions of tandem arrays would otherwise lose hits that the pilot kept
BASE = ["-j2", "-p", "0.3", "-N", "1000", "--outs=0.3", "--outc=0.3", "--outn=1000"]
SENS = ["-k5", "-M0", "-n2", "-l4"] + BASE
MIN_POS, MIN_COV, INTACT_COV, MAX_INTRON_OK, MIN_IDENT_OK = 0.45, 0.30, 0.90, 15_000, 0.30  # = analyze_genes.py
THREADS = "3"
NEAR_END, HOMOPOL = 50_000, 5


def log(*a):
    print(*a, flush=True)


waited = 0
while not DMND.exists() and waited < 3600:  # get_ref.sh may still be downloading Swiss-Prot when the first pair arrives here
    if waited == 0:
        log("waiting for", DMND)
    time.sleep(30)
    waited += 30
if not DMND.exists():
    sys.exit(f"reference database missing: {DMND}")


# ---------------- families ----------------
PREFIXES = (("KLRC", "KLRC"), ("KLRD", "KLRD1"), ("KLRK", "KLRK1"), ("KLRB", "KLRB1"), ("KLRF", "KLRF"), ("KLRG1", "KLRG1"),
            ("KLRG2", "KLRG2"), ("KLRA", "KLRA"), ("CLEC2D", "CLEC2D"), ("CLEC12A", "CLEC12A"), ("CLEC4E", "CLEC4E"), ("KIR", "KIR"),
            ("LILR", "LILR"), ("NCR1", "NCR1"), ("NCR2", "NCR2"), ("NCR3", "NCR3"), ("LAIR1", "LAIR1"), ("FCAR", "FCAR"), ("HLA-", "MHC-I"),
            ("H2-", "MHC-I"), ("MR1", "MR1"), ("CD1", "CD1"), ("B2M", "B2M"), ("TAP1", "TAP1"), ("TAP2", "TAP2"))


def family_of_query(sym: str) -> str:
    s = sym.upper()
    for pre, fam in PREFIXES:
        if s.startswith(pre):
            return fam
    raise SystemExit(f"no family for query symbol {sym}")


# stricter patterns for arbitrary Swiss-Prot gene names (CD14 is not CD1, HLA-DRB1 is class II, H2-Ab1 is class II)
BACK = (("KLRC", r"^KLRC\d"), ("KLRD1", r"^KLRD1$"), ("KLRK1", r"^KLRK1$"), ("KLRB1", r"^KLRB1"), ("KLRF", r"^KLRF\d"), ("KLRG1", r"^KLRG1$"),
        ("KLRG2", r"^KLRG2$"), ("KLRA", r"^KLRA\d"), ("CLEC2D", r"^CLEC2D$"), ("CLEC12A", r"^CLEC12A$"), ("CLEC4E", r"^CLEC4E$"), ("KIR", r"^KIR\d"),
        ("LILR", r"^LILR[AB]\d"), ("NCR1", r"^NCR1$"), ("NCR2", r"^NCR2$"), ("NCR3", r"^NCR3$"), ("LAIR1", r"^LAIR1$"), ("FCAR", r"^FCAR$"),
        ("MHC-I", r"^HLA-(?!D)"), ("MHC-II", r"^HLA-D"), ("MHC-I", r"^H2-(K\d|D\d|L$|Q\d|T\d|M\d|BL$)"), ("MHC-II", r"^H2-(A|E|DM|O)"),
        ("MR1", r"^MR1$"), ("CD1", r"^CD1[A-E]\d?$"), ("B2M", r"^B2M$"), ("TAP1", r"^TAP1$"), ("TAP2", r"^TAP2$"))


def back_family(sym: str) -> str:
    s = sym.upper()
    for fam, pat in BACK:
        if re.match(pat, s):
            return fam
    return "other"


qfam = {}
for line in (root / "queries" / "queries_final_table.tsv").read_text().splitlines()[1:]:
    f = line.split("\t")
    qfam[f"{f[0]}|{f[1]}|{f[2]}|{f[3]}"] = family_of_query(f[6].split(",")[0])
QSEQ: dict[str, str] = {}
_n = None
for line in QUERIES.read_text().splitlines():
    if line.startswith(">"):
        _n = line[1:].split()[0]
        QSEQ[_n] = ""
    elif _n:
        QSEQ[_n] += line.strip()

# ---------------- helpers ----------------
CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_aa):
    CODON[_b[_i // 16] + _b[_i // 4 % 4] + _b[_i % 4]] = _a
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(s: str) -> str:
    return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def read_tsv(p: Path) -> list[dict]:
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fetch(h: str, regions: list[str]) -> dict[str, str]:
    """regions 'contig:start-end' (1-based inclusive) -> upper-case sequence"""
    res: dict[str, list[str]] = {}
    for i in range(0, len(regions), 150):
        r = subprocess.run([SAMTOOLS, "faidx", fna[h], *regions[i:i + 150]], capture_output=True, text=True, check=True)
        name = None
        for line in r.stdout.splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                res[name] = []
            elif name is not None:
                res[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in res.items()}


def write_fa(path: Path, recs) -> None:
    with open(path, "w") as f:
        for n, s in recs:
            f.write(f">{n}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")


def miniprot(fa: Path, opts: list[str], prefix: str, queries: Path = QUERIES) -> list[dict]:
    """all 99 queries (or the given query file) against the records of fa; coordinates are 1-based on each record"""
    cmd = ["miniprot", "-t", THREADS, *opts, "--gff", "--trans", "-P", prefix, str(fa), str(queries)]
    r = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if r.returncode:
        sys.exit(f"miniprot failed: {' '.join(cmd)}\n{r.stderr[-800:]}")
    (logs / f"{prefix}{fa.stem}.gff").write_text(r.stdout)  # kept so that individual calls can be inspected later
    hits, h = [], None
    for line in r.stdout.splitlines():
        if line.startswith("##PAF"):
            f = line.split("\t")
            tags = {t.split(":")[0]: t.split(":", 2)[2] for t in f[13:]}
            h = {"query": f[1], "qlen": int(f[2]), "qs": int(f[3]), "qe": int(f[4]), "strand": f[5], "target": f[6], "tlen": int(f[7]),
                 "fs": int(tags.get("fs", 0)), "st": int(tags.get("st", 0)), "cs": tags.get("cs", ""), "cds": [], "stop": False, "prot": "",
                 "ms": 0, "ident": 0.0, "pos": 0.0, "start": 0, "end": 0}
            hits.append(h)
        elif line.startswith("##STA"):
            h["prot"] = line.split("\t")[1].strip()
        elif line.startswith("#") or not line.strip():
            continue
        else:
            f = line.split("\t")
            if f[2] == "mRNA":
                at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                h["ms"], h["ident"], h["pos"] = int(f[5]), float(at["Identity"]), float(at["Positive"])
                h["start"], h["end"] = int(f[3]), int(f[4])
            elif f[2] == "CDS":
                h["cds"].append((int(f[3]), int(f[4])))
            elif f[2] == "stop_codon":
                h["stop"] = True
    for h in hits:
        c = sorted(h["cds"])
        h["cds"] = c
        h["max_intron"] = max([b[0] - a[1] - 1 for a, b in zip(c, c[1:])] or [0])
        h["cov"] = (h["qe"] - h["qs"]) / h["qlen"]
        h["family"] = qfam[h["query"]]
        h["supported"] = h["pos"] >= MIN_POS and h["cov"] >= MIN_COV and h["max_intron"] <= MAX_INTRON_OK
        h["disrupted"] = h["fs"] > 0 or h["st"] > 0
        h["intact"] = h["supported"] and not h["disrupted"] and h["cov"] >= INTACT_COV and h["stop"] and h["ident"] >= MIN_IDENT_OK
    return hits


def classify(hs: list[dict]):
    """status of the hits of one family in one region; mirrors the locus rules of analyze_genes.py"""
    sup = [h for h in hs if h["supported"]]
    if not sup:
        return "none", None
    intact = [h for h in sup if h["intact"]]
    if intact:
        return "intact", max(intact, key=lambda h: h["ms"])
    if any(h["disrupted"] and h["cov"] >= 0.5 for h in sup) and not any((not h["disrupted"]) and h["cov"] >= 0.7 for h in sup):
        return "pseudogene-like", max((h for h in sup if h["disrupted"]), key=lambda h: h["ms"])
    best = max(sup, key=lambda h: h["ms"])
    return ("uncertain" if best["ident"] < MIN_IDENT_OK else "partial"), best


RANK = {"intact": 4, "partial": 3, "uncertain": 2, "pseudogene-like": 1, "none": 0}


def gene_model(h: dict, seq: str) -> dict:
    """start codon, splice-site dinucleotides (cs tag), terminal stop and internal stops of the hit's CDS on seq (1-based coords)"""
    spliced = "".join(seq[s - 1:e] for s, e in h["cds"])
    if h["strand"] == "-":
        spliced = rc(spliced)
    start = ("ATG" if spliced[:3] == "ATG" else f"non-ATG({spliced[:3]})") if h["qs"] == 0 else f"not_aligned(q_start={h['qs']})"
    introns = re.findall(r"~([a-z]{2})(\d+)([a-z]{2})", h["cs"])
    noncan = [f"{d}-{a}" for d, _, a in introns if (d, a) not in (("gt", "ag"), ("gc", "ag"), ("at", "ac"))]
    res = {"start_codon": start, "n_introns": len(introns), "noncanonical_splice": ";".join(noncan) or "-", "frame": "ok" if len(spliced) % 3 == 0 else "len%3!=0"}
    if h["fs"]:
        res.update(terminal_stop="na(frameshift)", internal_stops="na(frameshift)")
        return res
    prot = translate(spliced)
    if h["stop"]:
        res["terminal_stop"] = spliced[-3:] if spliced[-3:] in STOPS else f"last_codon_not_stop({spliced[-3:]})"
        res["internal_stops"] = str(prot[:-1].count("*"))
    else:
        res["terminal_stop"] = "no_stop_codon_feature"
        res["internal_stops"] = str(prot.count("*"))
    return res


def blocks_of(rec: dict):
    """aligned blocks (t0,t1,q0,q1,strand), 0-based half-open, from the cigar of a PAF record"""
    t, q = rec["ts"], (rec["qs"] if rec["strand"] == "+" else rec["qe"])
    bl = []
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
        elif op in "DN":
            t += n
        elif op == "I":
            q += n if rec["strand"] == "+" else -n
    return bl


def lift_q2t(blocks, x: int, y: int):
    """query interval [x,y) -> (covered fraction, target lo, target hi)"""
    cov, los, his = 0, [], []
    for t0, t1, q0, q1, st in blocks:
        s, e = max(x, q0), min(y, q1)
        if s >= e:
            continue
        cov += e - s
        if st == "+":
            los.append(t0 + (s - q0))
            his.append(t0 + (e - q0))
        else:
            los.append(t1 - (e - q0))
            his.append(t1 - (s - q0))
    if not los:
        return 0.0, None, None
    return cov / max(1, y - x), min(los), max(his)


def read_paf(p: Path) -> dict[str, list[dict]]:
    recs = defaultdict(list)
    for line in p.read_text().splitlines():
        f = line.split("\t")
        tags = {t[:2]: t[5:] for t in f[12:]}
        recs[f[0]].append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": f[5], "tlen": int(f[6]),
                           "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "al": int(f[10]), "cg": tags.get("cg", ""), "cs": tags.get("cs", "")})
    return recs


def minimap(target: str, query: Path, paf: Path, preset: str, whole_genome: bool) -> None:
    cmd = ["minimap2", "-x", preset, "-c", "--cs", "-N", "20", "-p", "0.5", "-t", THREADS]
    if whole_genome:
        cmd += ["-I", "300M"]  # multi-part index: ~1.5 GB RAM instead of ~6 GB for a 2 Gb genome
    cmd += [target, str(query)]
    with open(paf, "w") as fo, open(logs / (paf.stem + ".minimap2.log"), "w") as lg:
        subprocess.run(cmd, stdout=fo, stderr=lg, check=True)


def homopolymer(seq: str, pos0: int) -> int:
    """length of the run of identical bases that contains 0-based position pos0"""
    if pos0 < 0 or pos0 >= len(seq):
        return 0
    b, i, j = seq[pos0], pos0, pos0
    while i > 0 and seq[i - 1] == b:
        i -= 1
    while j + 1 < len(seq) and seq[j + 1] == b:
        j += 1
    return j - i + 1


def variants(ref_seq: str, ref_hit: dict, alt_seq: str, name: str) -> dict:
    """align alt_seq to ref_seq (minimap2 asm5, cs tag) and list the variants that fall inside the CDS exons of ref_hit
    (frameshift indels, stop-gain substitutions, indels in frame, splice-site changes)."""
    rfa, afa, paf = tmp / f"{name}_ref.fa", tmp / f"{name}_alt.fa", tmp / f"{name}.paf"
    write_fa(rfa, [("ref", ref_seq)])
    write_fa(afa, [("alt", alt_seq)])
    minimap(str(rfa), afa, paf, "asm5", False)
    recs = read_paf(paf).get("alt", [])
    exons = ref_hit["cds"]  # sorted ascending, 1-based inclusive
    spliced = "".join(ref_seq[s - 1:e] for s, e in exons)
    minus = ref_hit["strand"] == "-"
    if minus:
        spliced = rc(spliced)
    lens = [e - s + 1 for s, e in exons]

    def tpos(k: int, pos1: int) -> int:  # transcript position (0-based) of ref position pos1 inside exon k
        return (sum(lens[:k]) + (pos1 - exons[k][0])) if not minus else (sum(lens[k + 1:]) + (exons[k][1] - pos1))

    def exon_of(pos1: int):
        for k, (s, e) in enumerate(exons):
            if s <= pos1 <= e:
                return k
        return None

    items, n_cds, n_fs, n_stop, n_splice, aligned = [], 0, 0, 0, 0, 0
    for r in recs:
        aligned += r["al"]
        t = r["ts"]
        for tok in re.findall(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})", r["cs"]):
            if tok[0] == ":":
                t += int(tok[1:])
                continue
            if tok[0] == "~":
                t += int(re.search(r"\d+", tok).group())
                continue
            pos1 = t + 1
            k = exon_of(pos1)
            near = k is None and any(0 < s - pos1 <= 2 or 0 < pos1 - e <= 2 for s, e in exons)
            if tok[0] == "*":
                if k is not None:
                    n_cds += 1
                    p = tpos(k, pos1)
                    c0 = p // 3 * 3
                    codon = spliced[c0:c0 + 3]
                    alt_b = tok[2].upper() if not minus else tok[2].upper().translate(COMP)
                    new = codon[:p % 3] + alt_b + codon[p % 3 + 1:]
                    if len(codon) == 3 and new in STOPS and codon not in STOPS:
                        n_stop += 1
                        items.append(f"CDS{k + 1}:{tok[1].upper()}>{tok[2].upper()}@{pos1}->stop({new})")
                elif near:
                    n_splice += 1
                    items.append(f"splice_site:{tok[1].upper()}>{tok[2].upper()}@{pos1}")
                t += 1
            elif tok[0] == "+":  # insertion in alt relative to ref, between t-1 and t
                n = len(tok) - 1
                if k is not None or exon_of(pos1 - 1) is not None:
                    n_cds += 1
                    kk = k if k is not None else exon_of(pos1 - 1)
                    hp = homopolymer(ref_seq, t - 1) if len(set(tok[1:])) == 1 and ref_seq[t - 1:t].lower() == tok[1] else 0
                    if n % 3:
                        n_fs += 1
                        items.append(f"CDS{kk + 1}:ins{n}@{pos1}{'(homopolymer' + str(hp) + ')' if hp >= HOMOPOL else ''}")
                    else:
                        items.append(f"CDS{kk + 1}:inframe_ins{n}@{pos1}")
                elif near:
                    n_splice += 1
                    items.append(f"splice_site:ins{n}@{pos1}")
            elif tok[0] == "-":  # deletion of ref bases t .. t+n-1 in alt
                n = len(tok) - 1
                ks = {exon_of(p1) for p1 in range(pos1, pos1 + n)} - {None}
                if ks:
                    n_cds += 1
                    kk = min(ks)
                    hp = homopolymer(ref_seq, t) if len(set(tok[1:])) == 1 else 0
                    if n % 3:
                        n_fs += 1
                        items.append(f"CDS{kk + 1}:del{n}@{pos1}{'(homopolymer' + str(hp) + ')' if hp >= HOMOPOL else ''}")
                    else:
                        items.append(f"CDS{kk + 1}:inframe_del{n}@{pos1}")
                elif near:
                    n_splice += 1
                    items.append(f"splice_site:del{n}@{pos1}")
                t += n
    for p in (rfa, afa, paf):
        p.unlink(missing_ok=True)
    return {"items": ";".join(items[:8]) + (f";+{len(items) - 8}more" if len(items) > 8 else "") if items else "-",
            "n_cds_variants": n_cds, "n_frameshift_indels": n_fs, "n_stop_gains": n_stop, "n_splice": n_splice, "aligned_bp": aligned}


# ================= A. loci: re-alignment, gene model, proteins =================
log("A. loci")
loci = {h: read_tsv(out / f"{tag}.{h}.loci.tsv") for h in HAPS}
L = {(h, l["locus_id"]): l for h in HAPS for l in loci[h]}
loc_rows: dict[tuple, dict] = {}
proteins = []
for h in HAPS:
    regions = []
    for l in loci[h]:
        s, e, clen = int(l["start"]), int(l["end"]), fai[h][l["contig"]]
        a, b = max(1, s - LOCFLANK), min(clen, e + LOCFLANK)
        regions.append((l, f'{l["contig"]}:{a}-{b}', a))
    seqs = fetch(h, [r for _, r, _ in regions])
    fa = tmp / f"{h}_loci.fa"
    write_fa(fa, [(l["locus_id"], seqs[r]) for l, r, _ in regions])
    hits = miniprot(fa, BASE, f"{h}r_")
    byt = defaultdict(list)
    for x in hits:
        byt[x["target"]].append(x)
    for l, r, a in regions:
        lid = l["locus_id"]
        ls, le = int(l["start"]) - a + 1, int(l["end"]) - a + 1
        hs = [x for x in byt[lid] if x["start"] <= le and x["end"] >= ls]
        fam_hits = [x for x in hs if x["family"] == l["family"]]
        st_re, ev = classify(fam_hits)
        evq = [x for x in fam_hits if x["query"] == l["evidence_query"]]
        chosen = max(evq, key=lambda x: x["ms"]) if evq else ev
        row = {"haplotype": h, "locus_id": lid, "family": l["family"], "pilot_status": l["status"], "contig": l["contig"], "start": l["start"], "end": l["end"],
               "strand": l["strand"], "evidence_query": l["evidence_query"], "realign_status": st_re,
               "realign_consistent": "yes" if st_re == l["status"] else ("pilot_uncertain_by_span/identity_rule" if l["status"] == "uncertain" else "no"),
               "realign_query": chosen["query"] if chosen else "-", "realign_cov": f'{chosen["cov"]:.3f}' if chosen else "-",
               "realign_identity": chosen["ident"] if chosen else "-", "realign_frameshifts": chosen["fs"] if chosen else "-",
               "realign_inframe_stops": chosen["st"] if chosen else "-", "realign_stop_feature": int(chosen["stop"]) if chosen else "-",
               "start_codon": "-", "n_introns": "-", "noncanonical_splice": "-", "frame": "-", "terminal_stop": "-", "internal_stops": "-", "protein_len": 0,
               "curated_intact": "no", "curation_note": ""}
        if chosen:
            gm = gene_model(chosen, seqs[r])
            row.update(gm)
            row["protein_len"] = len(chosen["prot"])
            if chosen["prot"]:
                proteins.append((f'{h}|{lid}|{l["family"]}|{l["status"]}', chosen["prot"].replace("*", "")))
            if l["status"] == "intact":
                ok = st_re == "intact" and gm["start_codon"] == "ATG" and gm["noncanonical_splice"] == "-" and gm["terminal_stop"] in STOPS and gm["internal_stops"] == "0" and gm["frame"] == "ok"
                row["curated_intact"] = "yes" if ok else "no"
                if not ok:
                    why = []
                    if st_re != "intact":
                        why.append(f"realigned_as_{st_re}")
                    if gm["start_codon"] != "ATG":
                        why.append("start:" + gm["start_codon"])
                    if gm["noncanonical_splice"] != "-":
                        why.append("splice:" + gm["noncanonical_splice"])
                    if gm["terminal_stop"] not in STOPS:
                        why.append("stop:" + gm["terminal_stop"])
                    if gm["internal_stops"] != "0":
                        why.append("internal_stops:" + gm["internal_stops"])
                    if gm["frame"] != "ok":
                        why.append(gm["frame"])
                    row["curation_note"] = ",".join(why)
        else:
            # diagnostic: the same region searched alone with the evidence query alone (no competition, no relative-score filter)
            qseq = QSEQ.get(l["evidence_query"])
            if qseq:
                sfa, qfa_ = tmp / f"single_{h}_{lid}.fa", tmp / f"single_{h}_{lid}.faa"
                write_fa(sfa, [(lid, seqs[r])])
                write_fa(qfa_, [(l["evidence_query"], qseq)])
                single = miniprot(sfa, BASE, f"{h}q_", qfa_)
                sfa.unlink()
                qfa_.unlink()
                if single:
                    b = max(single, key=lambda x: x["ms"])
                    row["curation_note"] = (f"no_hit_of_this_family_on_realignment;alone_with_evidence_query:found(cov={b['cov']:.2f},fs={b['fs']},st={b['st']},"
                                            f"ident={b['ident']},ms={b['ms']},max_intron={b['max_intron']})")
                else:
                    row["curation_note"] = "no_hit_of_this_family_on_realignment;alone_with_evidence_query:none"
            else:
                row["curation_note"] = "no_hit_of_this_family_on_realignment"
        loc_rows[(h, lid)] = row
    log(f"  {h}: {len(regions)} loci, {len(hits)} region hits")
pfa = cur / f"{tag}.locus_proteins.faa"
write_fa(pfa, proteins)

# ================= B. back-search (Swiss-Prot human+mouse) and paralog identity =================
log("B. back-search")


def diamond(q: Path, db: str, k: int) -> dict[str, list[list[str]]]:
    o = tmp / f"{q.stem}.{Path(db).stem}.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", db, "-o", str(o), "-k", str(k), "--quiet", "--threads", THREADS, "--sensitive",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "stitle"], check=True)
    rows = defaultdict(list)
    for line in o.read_text().splitlines():
        f = line.split("\t")
        rows[f[0]].append(f)
    return rows


for r in loc_rows.values():
    r.update(back_symbol="-", back_entry="-", back_pident="-", back_evalue="-", back_family="-", back_second="-", family_check="no_protein",
             closest_paralog="-", paralog_pident="-")
if proteins:
    bs = diamond(pfa, str(DMND), 5)
    with open(cur / f"{tag}.backsearch.tsv", "w") as f:
        f.write("locus\tsubject\tpident\talign_len\tqlen\tslen\tevalue\tbitscore\tgene_name\tfamily_of_subject\n")
        for q, rows in bs.items():
            h, lid, fam, _ = q.split("|")
            parsed = []
            for x in rows:
                m = re.search(r"\bGN=(\S+)", x[8])
                sym = m.group(1) if m else x[1]
                parsed.append((x, sym, back_family(sym)))
                f.write("\t".join([q, x[1], x[2], x[3], x[4], x[5], x[6], x[7], sym, back_family(sym)]) + "\n")
            top, tsym, tfam = parsed[0]
            r = loc_rows[(h, lid)]
            r.update(back_symbol=tsym, back_entry=top[1], back_pident=top[2], back_evalue=top[6], back_family=tfam)
            second = next(((x, s, fm) for x, s, fm in parsed[1:] if fm != tfam), None)
            if second:
                r["back_second"] = f"{second[1]}({second[2]},{second[0][7]}/{top[7]})"
            if tfam == fam:
                r["family_check"] = "confirmed"
                if second and float(second[0][7]) >= 0.9 * float(top[7]):
                    r["family_check"] = f"confirmed_but_close_to_{second[2]}"
            elif second and second[2] == fam and float(second[0][7]) >= 0.9 * float(top[7]):
                r["family_check"] = f"ambiguous_{tfam}_vs_{fam}"
            else:
                r["family_check"] = f"relabel:{tfam}({tsym})"
    # paralogs: all-vs-all within haplotype and family
    subprocess.run(["diamond", "makedb", "--in", str(pfa), "-d", str(tmp / "loci_prot"), "--quiet"], check=True)
    pa = diamond(pfa, str(tmp / "loci_prot"), 200)
    for q, rows in pa.items():
        h, lid, fam, _ = q.split("|")
        best = None
        for x in rows:
            h2, lid2, fam2, _ = x[1].split("|")
            if h2 != h or fam2 != fam or lid2 == lid:
                continue
            if int(x[3]) < 0.5 * min(int(x[4]), int(x[5])):
                continue
            if best is None or float(x[2]) > float(best[2]):
                best = x
        if best:
            loc_rows[(h, lid)]["closest_paralog"] = best[1].split("|")[1]
            loc_rows[(h, lid)]["paralog_pident"] = best[2]
    for r in loc_rows.values():
        if r["family_check"] == "no_protein" and r["protein_len"]:
            r["family_check"] = "no_swissprot_hit"

LOC_COLS = ["haplotype", "locus_id", "family", "pilot_status", "contig", "start", "end", "strand", "evidence_query", "realign_status", "realign_consistent",
            "realign_query", "realign_cov", "realign_identity", "realign_frameshifts", "realign_inframe_stops", "realign_stop_feature", "start_codon", "n_introns",
            "noncanonical_splice", "frame", "terminal_stop", "internal_stops", "protein_len", "curated_intact", "curation_note", "back_symbol", "back_entry",
            "back_pident", "back_evalue", "back_family", "back_second", "family_check", "closest_paralog", "paralog_pident"]
with open(cur / f"{tag}.locus_curation.tsv", "w") as f:
    f.write("\t".join(LOC_COLS) + "\n")
    for r in loc_rows.values():
        f.write("\t".join(str(r[c]) for c in LOC_COLS) + "\n")

# ================= C. haplotype differences: direct re-examination =================
log("C. differences")
diffs = read_tsv(out / f"{tag}.haplotype_differences.tsv")
win = {}
for r in read_tsv(out / f"{tag}.region_pairs.tsv"):
    if r["region_hap1"] and r["hap1_extract_start"]:
        win[r["region_hap1"]] = ("hap1", r["hap1_contig"], int(r["hap1_extract_start"]), int(r["hap1_extract_end"]))
    if r["region_hap2"] and r["hap2_extract_start"]:
        win[r["region_hap2"]] = ("hap2", r["hap2_contig"], int(r["hap2_extract_start"]), int(r["hap2_extract_end"]))
tasks = []
seen: dict[tuple, int] = {}  # the pilot reports one event once per window pair that contains it; count duplicates instead of re-examining them
dup_count: Counter = Counter()
for i, d in enumerate(diffs):
    ft = d["feature"]
    if ft.startswith(("locus_hap1_only", "locus_hap2_only")):
        m = re.match(r"(hap[12]_L\d+) (\S+) (\S+);", d["note"])
        key = ("locus", m.group(1))
        if key in seen:
            dup_count[seen[key]] += 1
            continue
        seen[key] = i
        tasks.append({"i": i, "kind": "hap_only", "d": d, "src": (d["haplotype"], L[(d["haplotype"], m.group(1))])})
    elif ft == "locus_in_both_but_call_differs":
        m = re.match(r"(hap1_L\d+) (\S+)/(\S+) vs (hap2_L\d+) (\S+)/(\S+)", d["note"])
        key = ("pair", m.group(1), m.group(4))
        if key in seen:
            dup_count[seen[key]] += 1
            continue
        seen[key] = i
        tasks.append({"i": i, "kind": "call_differs", "d": d, "src": ("hap1", L[("hap1", m.group(1))]), "dst": ("hap2", L[("hap2", m.group(4))])})
    elif ft.startswith(("alignment_gap", "internal_indel")):
        key = ("gap", d["haplotype"], d["contig"], d["start"], d["end"])
        if key in seen:
            dup_count[seen[key]] += 1
            continue
        seen[key] = i
        tasks.append({"i": i, "kind": "gap", "d": d})
n_dup = sum(dup_count.values())


def region_of(h: str, contig: str, s: int, e: int, flank: int):
    a, b = max(1, s - flank), min(fai[h][contig], e + flank)
    return f"{contig}:{a}-{b}", a, b


# C1. whole-genome mapping of hap-only loci (+-10 kb) and gap sequences to the other haplotype
wg = {h: [] for h in HAPS}  # source hap -> (name, region string)
for t in tasks:
    d = t["d"]
    if t["kind"] == "hap_only":
        h, l = t["src"]
        reg, a, b = region_of(h, l["contig"], int(l["start"]), int(l["end"]), DIFFFLANK)
        t.update(name=f"diff{t['i']}", reg=reg, a=a, b=b, span=(int(l["start"]) - a, int(l["end"]) - a + 1))
        wg[h].append((t["name"], reg))
    elif t["kind"] == "gap":
        h = d["haplotype"]
        t.update(name=f"gap{t['i']}", reg=f'{d["contig"]}:{int(d["start"]) + 1}-{d["end"]}', hap=h)
        wg[h].append((t["name"], t["reg"]))
seqcache = {h: {} for h in HAPS}
paf = {}
for h in HAPS:
    if not wg[h]:
        continue
    seqs = fetch(h, [r for _, r in wg[h]])
    seqcache[h].update(seqs)
    qfa = tmp / f"{h}_to_{OTHER[h]}.fa"
    write_fa(qfa, [(n, seqs[r]) for n, r in wg[h]])
    pf = tmp / f"{h}_to_{OTHER[h]}.paf"
    log(f"  minimap2: {len(wg[h])} {h} sequences -> whole {OTHER[h]} genome")
    minimap(fna[OTHER[h]], qfa, pf, "asm20", True)
    paf.update(read_paf(pf))


def location_class(other: str, contig: str, lo: int, hi: int, rid_other: str) -> str:
    w = win.get(rid_other)
    if not w:
        return "no_paired_window"
    if w[1] == contig and lo <= w[3] and hi >= w[2]:
        return "in_paired_window"
    if w[1] == contig:
        return f"same_contig_outside_window({(min(abs(lo - w[3]), abs(w[2] - hi))) // 1000}kb_away)"
    return "other_contig"


# C2. counterpart regions -> miniprot at default and at high sensitivity; source regions -> high sensitivity
cp_regions = {h: [] for h in HAPS}   # hap of the counterpart -> (name, reg)
src_regions = {h: [] for h in HAPS}
for t in tasks:
    if t["kind"] == "hap_only":
        h, l = t["src"]
        o = OTHER[h]
        t["src_name"], t["cp_name"] = f"src{t['i']}", f"cp{t['i']}"  # distinct names: source and counterpart hits must not merge
        src_regions[h].append((t["src_name"], t["reg"]))
        recs = paf.get(t["name"], [])
        x, y = t["span"]
        qlen = t["b"] - t["a"] + 1
        scored = []
        for r in recs:
            fr, lo, hi = lift_q2t(blocks_of(r), x, y)
            scored.append((fr, r["al"], r, lo, hi))
        scored.sort(key=lambda s: (-s[0], -s[1]))
        t["n_alt"] = sum(1 for s in scored if s[0] >= 0.5)

        def flank_cov(contig):  # union over all records on that contig: covered fraction of the left and of the right 10 kb flank
            iv = []
            for r in recs:
                if r["t"] == contig:
                    iv += [(q0, q1) for _, _, q0, q1, _ in blocks_of(r)]
            iv.sort()
            merged = []
            for s, e in iv:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])

            def cov(a, b):
                return sum(max(0, min(b, e) - max(a, s)) for s, e in merged) / max(1, b - a)
            return (cov(0, x) if x > 0 else 0.0, cov(y, qlen) if qlen > y else 0.0)

        if scored and scored[0][3] is not None:
            fr, al, r, lo, hi = scored[0]
            t.update(cov=fr, rec=r, lo=lo, hi=hi, ident=r["nm"] / max(1, r["al"]), flanks=flank_cov(r["t"]))
            creg, ca, cb = region_of(o, r["t"], lo + 1, hi, DIFFFLANK)
            t.update(creg=creg, ca=ca, cb=cb)
            cp_regions[o].append((t["cp_name"], creg))
        else:
            best_any = max(recs, key=lambda r: r["al"], default=None)
            t.update(cov=0.0, rec=best_any, lo=None, hi=None, ident=(best_any["nm"] / max(1, best_any["al"])) if best_any else 0.0,
                     flanks=flank_cov(best_any["t"]) if best_any else (0.0, 0.0))
            if best_any:  # the flanks (or part of them) align although the locus itself does not -> look at that place for N-runs
                creg, ca, cb = region_of(o, best_any["t"], best_any["ts"] + 1, best_any["te"], DIFFFLANK)
                t.update(creg=creg, ca=ca, cb=cb)
                cp_regions[o].append((t["cp_name"], creg))
    elif t["kind"] == "call_differs":
        for key in ("src", "dst"):
            h, l = t[key]
            reg, a, b = region_of(h, l["contig"], int(l["start"]), int(l["end"]), DIFFFLANK)
            t[key + "_reg"], t[key + "_a"] = reg, a
            src_regions[h].append((f"{t['i']}_{key}", reg))
hits_sens, hits_base = {}, {}
for h in HAPS:
    allr = src_regions[h] + cp_regions[h]
    if not allr:
        continue
    need = [r for _, r in allr if r not in seqcache[h]]
    if need:
        seqcache[h].update(fetch(h, need))
    fa = tmp / f"{h}_diffregions.fa"
    write_fa(fa, [(n, seqcache[h][r]) for n, r in allr])
    for x in miniprot(fa, SENS, f"{h}S_"):
        hits_sens.setdefault(x["target"], []).append(x)
    if cp_regions[h]:
        fa2 = tmp / f"{h}_cpregions.fa"
        write_fa(fa2, [(n, seqcache[h][r]) for n, r in cp_regions[h]])
        for x in miniprot(fa2, BASE, f"{h}B_"):
            hits_base.setdefault(x["target"], []).append(x)


def status_in(hits_by_target: dict, name: str, fam: str, span=None):
    hs = [x for x in hits_by_target.get(name, []) if x["family"] == fam]
    if span:
        hs = [x for x in hs if x["start"] <= span[1] and x["end"] >= span[0]]
    return classify(hs)


def fmt_hit(h):
    return (f'{h["cov"]:.2f}', h["fs"], h["st"]) if h else ("-", "-", "-")


D_COLS = ["diff_index", "feature", "region_hap1", "region_hap2", "haplotype", "locus_id", "family", "pilot_status", "contig", "start", "end",
          "src_status_sensitive", "src_cov", "src_fs", "src_stops", "counterpart_hap", "counterpart_contig", "counterpart_start", "counterpart_end",
          "locus_span_aligned_frac", "alignment_identity", "flank_aligned_frac_left_right", "n_alt_landings", "counterpart_location",
          "counterpart_status_default", "counterpart_status_sensitive", "cp_cov", "cp_fs", "cp_stops", "flags", "variants_in_cds_of_better_copy",
          "n_frameshift_indels", "n_stop_gains", "final_class", "note"]
drows, grows = [], []
for t in tasks:
    d = t["d"]
    if t["kind"] == "gap":
        continue
    row = {c: "" for c in D_COLS}
    row.update(diff_index=t["i"], feature=d["feature"], region_hap1=d["region_hap1"], region_hap2=d["region_hap2"])
    if t["kind"] == "hap_only":
        h, l = t["src"]
        o = OTHER[h]
        fam = l["family"]
        row.update(haplotype=h, locus_id=l["locus_id"], family=fam, pilot_status=l["status"], contig=l["contig"], start=l["start"], end=l["end"], counterpart_hap=o)
        sst, sev = status_in(hits_sens, t["src_name"], fam, (t["span"][0] + 1, t["span"][1]))
        row.update(src_status_sensitive=sst)
        row["src_cov"], row["src_fs"], row["src_stops"] = fmt_hit(sev)
        row.update(locus_span_aligned_frac=f'{t["cov"]:.2f}', alignment_identity=f'{t["ident"]:.3f}', n_alt_landings=t["n_alt"],
                   flank_aligned_frac_left_right=f'{t["flanks"][0]:.2f}/{t["flanks"][1]:.2f}')
        flags = []
        srcseq = seqcache[h][t["reg"]]
        if int(l["start"]) <= NEAR_END or fai[h][l["contig"]] - int(l["end"]) <= NEAR_END:
            flags.append("src_within50kb_of_contig_end")
        if "N" in srcseq:
            flags.append(f"src_region_has_{srcseq.count('N')}N")
        cst, cev, bst = "none", None, "none"
        if t.get("creg"):
            cseq = seqcache[o][t["creg"]]
            rid_other = d["region_hap2"] if h == "hap1" else d["region_hap1"]
            if t["lo"] is not None:
                row.update(counterpart_contig=t["rec"]["t"], counterpart_start=t["lo"] + 1, counterpart_end=t["hi"],
                           counterpart_location=location_class(o, t["rec"]["t"], t["lo"], t["hi"], rid_other))
                cspan = (t["lo"] + 1 - t["ca"] + 1, t["hi"] - t["ca"] + 1)
            else:
                row.update(counterpart_contig=t["rec"]["t"], counterpart_start=t["rec"]["ts"] + 1, counterpart_end=t["rec"]["te"],
                           counterpart_location="flank_only:" + location_class(o, t["rec"]["t"], t["rec"]["ts"], t["rec"]["te"], rid_other))
                cspan = None
            cst, cev = status_in(hits_sens, t["cp_name"], fam, cspan)
            bst, _ = status_in(hits_base, t["cp_name"], fam, cspan)
            if "N" in cseq:
                flags.append(f"counterpart_region_has_{cseq.count('N')}N")
            if t["rec"]["ts"] <= NEAR_END or t["rec"]["tlen"] - t["rec"]["te"] <= NEAR_END:
                flags.append("counterpart_within50kb_of_contig_end")
        row.update(counterpart_status_default=bst, counterpart_status_sensitive=cst)
        row["cp_cov"], row["cp_fs"], row["cp_stops"] = fmt_hit(cev)
        row["flags"] = ";".join(flags) or "-"
        var = None
        if t["cov"] < 0.2:
            fl = t["flanks"]
            if fl[0] >= 0.5 and fl[1] >= 0.5:
                cls = f"locus_sequence_absent_from_{o}_between_aligned_flanks"
            elif t["rec"] is None:
                cls = f"no_alignment_of_locus_or_flanks_to_{o}"
            else:
                cls = f"locus_sequence_absent_from_{o}(flanks_partly_aligned)"
            if any(f.startswith("counterpart_region_has") or f.startswith("counterpart_within") for f in flags):
                cls += ";possible_assembly_gap_or_contig_end"
        elif cst == sst and cst != "none":
            cls = "same_status_at_sensitive_realignment(sensitivity_artefact_of_the_pilot_call)"
        elif cst == "none":
            cls = "sequence_aligned_but_no_gene_hit_in_counterpart(diverged_or_fragmented)"
        else:
            cls = f"status_differs_confirmed:{h}={sst}_vs_{o}={cst}"
            if sev and cev:
                ref_side = (srcseq, sev) if RANK[sst] >= RANK[cst] else (seqcache[o][t["creg"]], cev)
                alt_seq = seqcache[o][t["creg"]] if RANK[sst] >= RANK[cst] else srcseq
                var = variants(ref_side[0], ref_side[1], alt_seq, t["src_name"])
        if t["n_alt"] >= 2:
            cls += f";multi_copy_landing(n={t['n_alt']})"
        if row["counterpart_location"].startswith(("other_contig", "same_contig_outside")):
            cls += ";counterpart_outside_paired_window"
        if var:
            row.update(variants_in_cds_of_better_copy=var["items"], n_frameshift_indels=var["n_frameshift_indels"], n_stop_gains=var["n_stop_gains"])
        row["final_class"] = cls
        row["note"] = d["note"]
    else:  # call_differs
        (h1, l1), (h2, l2) = t["src"], t["dst"]
        fam = l1["family"]
        row.update(haplotype="both", locus_id=f'{l1["locus_id"]}|{l2["locus_id"]}', family=f'{l1["family"]}|{l2["family"]}',
                   pilot_status=f'{l1["status"]}|{l2["status"]}', contig=l1["contig"], start=l1["start"], end=l1["end"], counterpart_hap="hap2",
                   counterpart_contig=l2["contig"], counterpart_start=l2["start"], counterpart_end=l2["end"], counterpart_location="matched_locus_pair")
        s1 = (int(l1["start"]) - t["src_a"] + 1, int(l1["end"]) - t["src_a"] + 1)
        s2 = (int(l2["start"]) - t["dst_a"] + 1, int(l2["end"]) - t["dst_a"] + 1)
        st1, ev1 = status_in(hits_sens, f"{t['i']}_src", fam, s1)
        st2, ev2 = status_in(hits_sens, f"{t['i']}_dst", fam, s2)
        if st2 == "none" and l2["family"] != fam:  # the hap2 locus carries another family label: judge it by its own family too
            st2b, ev2b = status_in(hits_sens, f"{t['i']}_dst", l2["family"], s2)
            if st2b != "none":
                st2, ev2 = st2b, ev2b
        row.update(src_status_sensitive=st1, counterpart_status_sensitive=st2, counterpart_status_default="-")
        row["src_cov"], row["src_fs"], row["src_stops"] = fmt_hit(ev1)
        row["cp_cov"], row["cp_fs"], row["cp_stops"] = fmt_hit(ev2)
        row["flags"] = "-"
        row["locus_span_aligned_frac"], row["alignment_identity"], row["n_alt_landings"], row["flank_aligned_frac_left_right"] = "-", "-", "-", "-"
        seq1, seq2 = seqcache[h1][t["src_reg"]], seqcache[h2][t["dst_reg"]]
        var = None
        if st1 == st2 and st1 != "none":
            cls = "same_status_at_sensitive_realignment(call_difference_is_a_threshold/sensitivity_artefact)"
        elif st1 == "none" or st2 == "none":
            cls = f"gene_hit_missing_on_one_side_at_sensitive_realignment:hap1={st1},hap2={st2}"
        else:
            cls = f"status_differs_confirmed:hap1={st1}_vs_hap2={st2}"
        if ev1 and ev2 and st1 != st2:
            ref_seq, ref_hit, alt_seq = (seq1, ev1, seq2) if RANK[st1] >= RANK[st2] else (seq2, ev2, seq1)
            var = variants(ref_seq, ref_hit, alt_seq, f"cd{t['i']}")
            row.update(variants_in_cds_of_better_copy=var["items"], n_frameshift_indels=var["n_frameshift_indels"], n_stop_gains=var["n_stop_gains"])
        row["final_class"] = cls
        row["note"] = d["note"]
    drows.append(row)
with open(cur / f"{tag}.difference_curation.tsv", "w") as f:
    f.write("\t".join(D_COLS) + "\n")
    for r in drows:
        f.write("\t".join(str(r[c]) for c in D_COLS) + "\n")

# C3. gaps
G_COLS = ["diff_index", "feature", "haplotype", "contig", "start", "end", "length", "N_frac", "other_hap", "best_contig", "best_start", "best_end",
          "covered_frac", "identity", "n_landings_ge50pct", "location_vs_paired_window", "gap_class", "n_duplicate_rows_in_pilot_table"]
for t in tasks:
    if t["kind"] != "gap":
        continue
    d = t["d"]
    h = t["hap"]
    o = OTHER[h]
    seq = seqcache[h][t["reg"]]
    n = len(seq)
    recs = paf.get(t["name"], [])
    scored = []
    for r in recs:
        fr, lo, hi = lift_q2t(blocks_of(r), 0, n)
        scored.append((fr, r["al"], r, lo, hi))
    scored.sort(key=lambda s: (-s[0], -s[1]))
    nl = sum(1 for s in scored if s[0] >= 0.5)
    nfrac = seq.count("N") / max(1, n)
    rid_other = d["region_hap2"] if h == "hap1" else d["region_hap1"]
    if scored and scored[0][3] is not None:
        fr, al, r, lo, hi = scored[0]
        loc = location_class(o, r["t"], lo, hi, rid_other)
        if fr < 0.2:
            cls = f"absent_from_{o}"
        elif loc == "in_paired_window":
            cls = "present_in_paired_window(window_alignment_missed_it)"
        else:
            cls = f"present_elsewhere_in_{o}"
        if nl >= 5:
            cls += f";repetitive({nl}_landings)"
        best = (r["t"], lo + 1, hi, f"{fr:.2f}", f'{r["nm"] / max(1, r["al"]):.3f}')
    else:
        loc, cls, best = "-", f"absent_from_{o}", ("-", "-", "-", "0.00", "-")
    if nfrac >= 0.5:
        cls = f"N_run({nfrac:.2f})"
    grows.append([t["i"], d["feature"], h, d["contig"], d["start"], d["end"], n, f"{nfrac:.2f}", o, *best, nl, loc, cls, dup_count.get(t["i"], 0)])
with open(cur / f"{tag}.gap_curation.tsv", "w") as f:
    f.write("\t".join(G_COLS) + "\n")
    for r in grows:
        f.write("\t".join(map(str, r)) + "\n")

# ================= D. summary =================
lines = [f"CURATION SUMMARY {tag}", ""]
for h in HAPS:
    rows = [r for r in loc_rows.values() if r["haplotype"] == h]
    intact = [r for r in rows if r["pilot_status"] == "intact"]
    lines.append(f"{h}: {len(rows)} loci; realignment consistent with pilot status: {sum(r['realign_consistent'] == 'yes' for r in rows)}; "
                 f"pilot intact {len(intact)}, pass all gene-model checks (ATG, GT-AG/GC-AG, terminal stop, no internal stop): {sum(r['curated_intact'] == 'yes' for r in intact)}")
    fc = Counter(r["family_check"].split(":")[0].split("_")[0] for r in rows)
    lines.append("  back-search: " + ", ".join(f"{k}={v}" for k, v in sorted(fc.items())))
    rel = Counter(f'{r["family"]}->{r["back_family"]}({r["back_symbol"]})' for r in rows if r["family_check"].startswith("relabel"))
    if rel:
        lines.append("  relabel suggestions: " + "; ".join(f"{k} x{v}" for k, v in rel.most_common(12)))
    hi = Counter(r["family"] for r in rows if r["paralog_pident"] != "-" and float(r["paralog_pident"]) >= 98)
    if hi:
        lines.append("  loci with a >=98% identical paralog (same haplotype, same family): " + ", ".join(f"{k}={v}" for k, v in sorted(hi.items())))
lines.append("")
lines.append(f"pilot difference rows: {len(diffs)}; of these {n_dup} are the same event repeated in another window pair (not re-examined twice)")
lines.append(f"locus-level differences re-examined: {len(drows)}")
for k, v in Counter(r["final_class"].split(";")[0] for r in drows).most_common():
    lines.append(f"  {v:3d}  {k}")
lines.append(f"gaps/indels >=5 kb re-mapped (unique): {len(grows)}")
for k, v in Counter(r[-2].split(";")[0] for r in grows).most_common():
    lines.append(f"  {v:3d}  {k}")
(cur / f"{tag}.curation_summary.txt").write_text("\n".join(lines) + "\n")
log("\n".join(lines))
for p in tmp.glob("*.fa"):
    p.unlink()
log("done")

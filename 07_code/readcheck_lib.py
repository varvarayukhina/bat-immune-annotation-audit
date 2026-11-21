"""Shared helpers for the read-level checks (hifi_check.py, rna_check.py): windows, PAF loading, cs parsing, coordinate lifting."""
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

CS = re.compile(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})")
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
VAR = re.compile(r"^(CDS\d+|splice_site):(?:([ACGT])>([ACGT])@(\d+)(?:->stop\((\w+)\))?|(ins|del|inframe_ins|inframe_del)(\d+)@(\d+)(?:\(homopolymer(\d+)\))?)$")
RANK = {"intact": 4, "partial": 3, "uncertain": 2, "pseudogene-like": 1, "none": 0}
DIFFFLANK = 10_000


def rd(p: Path) -> list[dict]:
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fa(p: Path) -> dict[str, str]:
    seqs, name = {}, None
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def load_windows(vdir: Path):
    W = rd(vdir / "windows.tsv")
    for w in W:
        w["start"], w["end"] = int(w["start"]), int(w["end"])
    WSEQ = {}
    for h in ("hap1", "hap2"):
        WSEQ.update(read_fa(vdir / f"windows_{h}.fa"))
    return W, WSEQ


def window_of(W, hap, contig, s, e):
    """window containing [s,e] of contig in hap -> (window row, offset of s, offset of e), 1-based within the window; None if outside"""
    for w in W:
        if w["hap"] == hap and w["contig"] == contig and w["start"] <= s and e <= w["end"]:
            return w, s - w["start"] + 1, e - w["start"] + 1
    return None


def load_paf(paths, min_block: int, max_de: float, keep=None):
    """primary records with block >= min_block and divergence <= max_de, grouped by target and sorted by start;
    keep(target, ts, te) -> bool, if given, drops records at load time (memory: a 24 M-record PAF does not fit in 8 GB)"""
    recs, n_in, n_keep = defaultdict(list), 0, 0
    for paf in paths:
        with open(paf) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                n_in += 1
                tags = {t[:2]: t[5:] for t in f[12:]}
                if tags.get("tp") != "P" or int(f[10]) < min_block or float(tags.get("de", 0)) > max_de:
                    continue
                if keep is not None and not keep(f[5], int(f[7]), int(f[8])):
                    continue
                n_keep += 1
                recs[f[5]].append({"q": f[0], "ql": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "ts": int(f[7]), "te": int(f[8]),
                                   "mq": int(f[11]), "cs": tags.get("cs", "")})
    for v in recs.values():
        v.sort(key=lambda r: r["ts"])
    return recs, n_in, n_keep


def lift_table(ref_seq: str, alt_seq: str, tmp: Path, name: str):
    """align alt to ref (minimap2 asm5, cs) -> (lift: ref_pos1 -> alt_pos1 | None, strand, aligned_bp)"""
    rfa, afa = tmp / f"{name}_r.fa", tmp / f"{name}_a.fa"
    rfa.write_text(f">ref\n{ref_seq}\n")
    afa.write_text(f">alt\n{alt_seq}\n")

    def run():
        out = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", "-t", "2", str(rfa), str(afa)], capture_output=True, text=True, check=False).stdout
        return [ln.split("\t") for ln in out.splitlines()]

    recs_ = run()
    if not recs_:
        def none(_p):
            return None
        return none, "+", 0
    strand = max(recs_, key=lambda f: int(f[10]))[4]
    if strand == "-":  # re-align the reverse complement so that query coordinates run with the target
        afa.write_text(f">alt\n{rc(alt_seq)}\n")
        recs_ = run()
    blocks = []
    for f in recs_:
        if f[4] != "+":
            continue
        cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), "")
        blocks.append((int(f[7]), int(f[2]), cs, int(f[10])))
    blocks.sort(key=lambda b: -b[3])
    alen = sum(b[3] for b in blocks)
    n_alt = len(alt_seq)

    def conv(qq):
        return (n_alt - qq + 1) if strand == "-" else qq

    def lift(pos1):
        for ts, qs, cs, _ in blocks:
            t, q = ts, qs
            for tok in CS.findall(cs):
                if tok[0] == ":":
                    n = int(tok[1:])
                    if t < pos1 <= t + n:
                        return conv(q + (pos1 - t))
                    t += n
                    q += n
                elif tok[0] == "*":
                    if t + 1 == pos1:
                        return conv(q + 1)
                    t += 1
                    q += 1
                elif tok[0] == "+":
                    q += len(tok) - 1
                elif tok[0] == "-":
                    n = len(tok) - 1
                    if t < pos1 <= t + n:  # deleted in alt: junction = last alt base before the deletion
                        return conv(q)
                    t += n
                if t > pos1:
                    break
        return None

    return lift, strand, alen


def events(r: dict, off1: int, slop: int):
    """events of one read within target offsets [off1-slop, off1+slop] (1-based); None if the read does not span that interval;
    introns (~) are reported as 'intron' when they cover the site"""
    lo, hi = off1 - slop, off1 + slop
    if r["ts"] + 1 > lo or r["te"] < hi:
        return None
    t, ev = r["ts"], []
    for tok in CS.findall(r["cs"]):
        c = tok[0]
        if c == ":":
            t += int(tok[1:])
        elif c == "*":
            if lo <= t + 1 <= hi:
                ev.append("sub:" + tok[2].upper())
            t += 1
        elif c == "+":
            if lo <= t + 1 <= hi + 1:
                ev.append(f"ins{len(tok) - 1}")
        elif c == "-":
            n = len(tok) - 1
            if t + 1 <= hi and t + n >= lo:
                ev.append(f"del{n}")
            t += n
        elif c == "~":
            n = int(re.search(r"\d+", tok).group())
            if t + 1 <= hi and t + n >= lo:
                ev.append("intron")
            t += n
        if t > hi + 1:
            break
    return ev


def site_counts(recs_of_target: list, off1: int, slop: int, mapq: int = 20):
    c = {"n": 0, "n_mq": 0, "match": 0, "sub": 0, "ins": 0, "del": 0, "intron": 0, "detail": defaultdict(int)}
    for r in recs_of_target:
        if r["te"] < off1 - slop:
            continue
        if r["ts"] + 1 > off1 + slop:
            break
        ev = events(r, off1, slop)
        if ev is None:
            continue
        if "intron" in ev:  # spliced over the site: not informative for the allele
            c["intron"] += 1
            continue
        c["n"] += 1
        c["n_mq"] += r["mq"] >= mapq
        if not ev:
            c["match"] += 1
        else:
            for e in ev:
                c["detail"][e] += 1
            kinds = {e.split(":")[0][:3] for e in ev}
            for k in ("sub", "ins", "del"):
                if k in kinds:
                    c[k] += 1
    return c


def fmt_counts(c) -> str:
    return f"n={c['n']}(mq20 {c['n_mq']}) match={c['match']} sub={c['sub']} ins={c['ins']} del={c['del']} [{','.join(f'{k}x{v}' for k, v in sorted(c['detail'].items()))}]"


def variant_sides(r: dict):
    """for a status_differs row: (ref_hap, ref_locus(contig,s,e), ref_status, alt_hap, alt_locus, alt_status, classes)"""
    m = re.match(r"status_differs_confirmed:(hap[12])=([\w-]+)_vs_(hap[12])=([\w-]+)", r["final_class"].split(";")[0])
    hA, stA, hB, stB = m.groups()
    locA = (r["contig"], int(r["start"]), int(r["end"]))
    locB = (r["counterpart_contig"], int(r["counterpart_start"]), int(r["counterpart_end"]))
    if RANK[stA] >= RANK[stB]:
        return hA, locA, stA, hB, locB, stB, f"{hA}={stA},{hB}={stB}"
    return hB, locB, stB, hA, locA, stA, f"{hA}={stA},{hB}={stB}"

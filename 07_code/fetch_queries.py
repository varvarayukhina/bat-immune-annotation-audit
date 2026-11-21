#!/usr/bin/env python3
"""Fetch query proteins from NCBI Datasets v2. Throttled ~4 req/s, 5 retries."""
import io
import json
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)

NKC = "KLRC1 KLRC2 KLRC3 KLRC4 KLRD1 KLRK1 KLRB1 KLRF1 KLRF2 KLRG1 KLRG2 CLEC2D CLEC12A CLEC4E".split()
LRC = "KIR2DL1 KIR2DL4 KIR3DL1 KIR3DL3 KIR2DS4 LILRB1 LILRB2 LILRB4 LILRA2 LILRA5 NCR1 NCR2 NCR3 LAIR1 FCAR".split()
MHC = "HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G MR1 CD1A CD1D B2M TAP1 TAP2".split()
MOUSE_EXTRA = "Klra1 Klra2 Klra3 Klra4 Klra5 Klra7 Klra8 Klra9 Klrb1c H2-K1 H2-D1 H2-Q7 H2-T23 Cd1d1".split()
HUMAN = NKC + LRC + MHC

TAXA = {
    "Hsapiens": (9606, HUMAN),
    "Mmusculus": (10090, [s.capitalize() for s in HUMAN] + MOUSE_EXTRA),
    "Raegyptiacus": (9407, HUMAN),
}
# capitalize() turns "HLA-A" into "Hla-a"; H2-* mouse symbols are given verbatim in MOUSE_EXTRA.

_last = [0.0]


def get(url: str) -> bytes:
    for attempt in range(1, 6):
        wait = 0.25 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json, */*"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):  # not transient
                raise
            err = f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = repr(e)
        print(f"  retry {attempt}/5 {url[:100]}: {err}", file=sys.stderr)
        time.sleep(2 * attempt)
    raise RuntimeError(f"giving up after 5 attempts: {url}")


def parse_fasta(txt: str):
    name, seq = None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


rows = []  # per-gene table
missing = []  # symbols with nothing
fasta_all = []

for label, (taxid, symbols) in TAXA.items():
    genes = {}  # gene_id -> {"symbol", "type", "queried": set}
    answered = set()  # queried symbols (lower) that got any report
    for i in range(0, len(symbols), 20):
        batch = symbols[i : i + 20]
        url = f"{BASE}/gene/symbol/{','.join(batch)}/taxon/{taxid}?page_size=1000"
        try:
            d = json.loads(get(url))
        except urllib.error.HTTPError as e:
            print(f"{label} batch {batch}: HTTP {e.code} {e.read()[:200]!r}", file=sys.stderr)
            continue
        for rep in d.get("reports", []):
            g = rep["gene"]
            e = genes.setdefault(g["gene_id"], {"symbol": g["symbol"], "type": g.get("type"), "queried": set(), "desc": g.get("description", "")})
            for q in rep.get("query", []):
                e["queried"].add(q)
                answered.add(q.lower())
    for s in symbols:
        if s.lower() not in answered:
            missing.append((label, s, "symbol endpoint returned no gene"))
    for gid, e in genes.items():
        sym = e["symbol"]
        q = ",".join(sorted(e["queried"]))
        try:
            z = get(f"{BASE}/gene/id/{gid}/download?include_annotation_type=FASTA_PROTEIN&filename=g.zip")
            with zipfile.ZipFile(io.BytesIO(z)) as zf:
                names = [n for n in zf.namelist() if n.endswith("protein.faa")]
                if not names:
                    missing.append((label, q, f"gene {gid} {sym} ({e['type']}): zip has no protein.faa"))
                    continue
                txt = zf.read(names[0]).decode()
        except urllib.error.HTTPError as ex:
            missing.append((label, q, f"gene {gid} {sym}: HTTP {ex.code} on protein download"))
            continue
        recs = list(parse_fasta(txt))
        if not recs:
            missing.append((label, q, f"gene {gid} {sym}: empty protein.faa"))
            continue
        hdr, seq = max(recs, key=lambda r: len(r[1]))
        acc = hdr.split()[0]
        tag = f"{label}|{sym}|{gid}|{acc}"
        fasta_all.append((tag, seq))
        alias = "" if sym.lower() in {x.lower() for x in e["queried"]} else "ALIAS"
        rows.append((label, sym, gid, acc, len(recs), len(seq), q, alias, e["desc"]))
        print(f"{tag}\tlen={len(seq)}\tisoforms={len(recs)}\tqueried={q}\t{alias}", flush=True)

with (OUT / "queries.faa").open("w") as f:
    for tag, seq in fasta_all:
        f.write(f">{tag}\n")
        for j in range(0, len(seq), 80):
            f.write(seq[j : j + 80] + "\n")
with (OUT / "queries_table.tsv").open("w") as f:
    f.write("taxon\tsymbol\tgene_id\tprotein_acc\tn_isoforms\tlongest_len_aa\tqueried_as\talias_match\tdescription\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
with (OUT / "queries_missing.tsv").open("w") as f:
    f.write("taxon\tsymbol\treason\n")
    for m in missing:
        f.write("\t".join(m) + "\n")
print(f"\nqueries: {len(rows)}  missing: {len(missing)}")

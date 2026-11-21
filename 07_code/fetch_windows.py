#!/usr/bin/env python3
"""Fetch the pilot windows (+-100 kb, as compare_regions.py extracted them) from NCBI nuccore by efetch,
so the read-level checks do not need the whole genomes again.
usage: fetch_windows.py pair_dir out_dir  -> out_dir/windows_hap1.fa, windows_hap2.fa, windows.tsv"""
import csv
import sys
import time
import urllib.request
from pathlib import Path

pair, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
tag = pair.name
FLANK = 100_000
EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def efetch(acc: str, s: int, e: int) -> str:
    url = f"{EUT}?db=nuccore&id={acc}&seq_start={s}&seq_stop={e}&rettype=fasta&retmode=text"
    for i in range(4):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                txt = r.read().decode()
            if txt.startswith(">"):
                return "".join(txt.split("\n", 1)[1].split())
            raise ValueError(txt[:200])
        except Exception as ex:  # noqa: BLE001 - retry then fail loudly
            print(f"  retry {i} {acc}:{s}-{e}: {ex}", file=sys.stderr)
            time.sleep(3 * (i + 1))
    raise SystemExit(f"efetch failed for {acc}:{s}-{e}")


rows = list(csv.DictReader(open(pair / "out" / f"{tag}.region_summary.tsv"), delimiter="\t"))
meta = []
for h in ("hap1", "hap2"):
    recs = []
    for r in rows:
        if not r[f"{h}_window_start"]:  # window exists in the other haplotype only
            continue
        contig, ws, we = r[f"{h}_contig"], int(r[f"{h}_window_start"]), int(r[f"{h}_window_end"])
        a = max(1, ws - FLANK)
        b = we + FLANK  # efetch clips seq_stop at the sequence end
        seq = efetch(contig, a, b).upper()
        b = a + len(seq) - 1
        name = f"{h}|{r[f'region_{h}']}|{contig}:{a}-{b}"
        recs.append((name, seq))
        meta.append({"hap": h, "region": r[f"region_{h}"], "contig": contig, "start": a, "end": b, "len": len(seq), "name": name})
        print(f"{name} {len(seq)} bp", file=sys.stderr)
        time.sleep(0.4)
    with open(out / f"windows_{h}.fa", "w") as f:
        for n, s in recs:
            f.write(f">{n}\n" + "\n".join(s[i:i + 100] for i in range(0, len(s), 100)) + "\n")
with open(out / "windows.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(meta[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(meta)

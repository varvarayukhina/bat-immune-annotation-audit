#!/usr/bin/env python3
"""Run miniprot on a genome in chunks (bin-packed sequences, <= MAXBP each) with identical
parameters, because whole-genome indexing thrashed on this 8 GB machine. Concatenates GFFs.
usage: run_miniprot_chunked.py genome.fna label base_dir [MAXBP]"""
import subprocess
import sys
from pathlib import Path

genome, label, base = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
maxbp = int(sys.argv[4]) if len(sys.argv) > 4 else 250_000_000
SAMTOOLS = str(base / "bin" / "samtools")
queries = base / "queries" / "queries_final.faa"
out_gff = base / "mapping" / f"{label}.miniprot.gff"
tmp = base / "mapping" / "tmp_chunks"
tmp.mkdir(exist_ok=True)

seqs = []
for line in Path(f"{genome}.fai").read_text().splitlines():
    f = line.split("\t")
    seqs.append((f[0], int(f[1])))
seqs.sort(key=lambda x: -x[1])
chunks = []  # first-fit decreasing
for name, ln in seqs:
    for c in chunks:
        if c["bp"] + ln <= maxbp:
            c["names"].append(name)
            c["bp"] += ln
            break
    else:
        chunks.append({"names": [name], "bp": ln})
print(f"{label}: {len(seqs)} sequences -> {len(chunks)} chunks", flush=True)

with out_gff.open("w") as fout:
    for i, c in enumerate(chunks, 1):
        fa = tmp / f"{label}_c{i:02d}.fa"
        with fa.open("wb") as fh:
            subprocess.run([SAMTOOLS, "faidx", str(genome), *c["names"]], stdout=fh, check=True)
        cmd = ["miniprot", "-t4", "-j2", "-p", "0.3", "-N", "100", "--outs=0.3", "--outc=0.3",
               "--outn=200", "--gff", "-P", f"{label}c{i:02d}_", str(fa), str(queries)]
        with (base / "logs" / f"miniprot_{label}_c{i:02d}.log").open("w") as lg:
            r = subprocess.run(cmd, stdout=fout, stderr=lg)
        fa.unlink()
        print(f"chunk {i}/{len(chunks)} {c['bp']/1e6:.0f} Mb {len(c['names'])} seqs rc={r.returncode}", flush=True)
        if r.returncode != 0:
            sys.exit(f"miniprot failed on chunk {i}; see logs/miniprot_{label}_c{i:02d}.log")
print("done", flush=True)

#!/usr/bin/env python3
"""Higher-sensitivity miniprot (-k5 -M0 -n2 -l4) restricted to the default-run region windows +-FLANK; coordinates are lifted back to genome coordinates.
usage: run_sens_windows.py src_base sens_base hap genome.fna"""
import subprocess, sys
from collections import defaultdict
from pathlib import Path
src, dst, hap, fna = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
pair_src = sys.argv[5]
FLANK = 500_000
for d in ("mapping", "out", "logs"):
    (dst / d).mkdir(parents=True, exist_ok=True)
for n in ("queries", "bin"):
    if not (dst / n).exists():
        (dst / n).symlink_to(src / n)
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in Path(f"{fna}.fai").read_text().splitlines()}
byc = defaultdict(list)
for l in (src / "out" / f"{pair_src}.{hap}.regions.tsv").read_text().splitlines()[1:]:
    f = l.split("\t")
    byc[f[2]].append((max(1, int(f[4]) - FLANK), min(fai[f[2]], int(f[5]) + FLANK)))
wins = []
for c, iv in byc.items():
    for s, e in sorted(iv):
        if wins and wins[-1][0] == c and s <= wins[-1][2]:
            wins[-1][2] = max(wins[-1][2], e)
        else:
            wins.append([c, s, e])
print(hap, len(wins), "windows", sum(w[2] - w[1] + 1 for w in wins) / 1e6, "Mb", flush=True)
wfa = dst / "mapping" / f"{hap}_windows.fa"
with wfa.open("wb") as fh:
    subprocess.run([str(src / "bin" / "samtools"), "faidx", str(fna), *[f"{c}:{s}-{e}" for c, s, e in wins]], stdout=fh, check=True)
raw = dst / "mapping" / f"{hap}_windows.raw.gff"
with (dst / "logs" / f"miniprot_{hap}.log").open("w") as lg, raw.open("w") as fo:
    r = subprocess.run(["miniprot", "-t4", "-k5", "-M0", "-n2", "-l4", "-j2", "-p", "0.3", "-N", "100", "--outs=0.3", "--outc=0.3", "--outn=200",
                        "--gff", "-P", f"{hap}s_", str(wfa), str(src / "queries" / "queries_final.faa")], stdout=fo, stderr=lg)
if r.returncode:
    sys.exit(f"miniprot failed rc={r.returncode}")
# lift back
def lift(name):
    c, rng = name.rsplit(":", 1)
    return c, int(rng.split("-")[0]) - 1
out = (dst / "mapping" / f"{hap}.miniprot.gff").open("w")
for line in raw.read_text().splitlines():
    if line.startswith("##PAF"):
        f = line.split("\t")
        c, off = lift(f[6])
        f[6], f[7] = c, str(fai[c])
        f[8], f[9] = str(int(f[8]) + off), str(int(f[9]) + off)
        out.write("\t".join(f) + "\n")
    elif line.startswith("#"):
        out.write(line + "\n")
    else:
        f = line.split("\t")
        c, off = lift(f[0])
        f[0], f[3], f[4] = c, str(int(f[3]) + off), str(int(f[4]) + off)
        out.write("\t".join(f) + "\n")
out.close()
wfa.unlink(); raw.unlink()
print("done", flush=True)

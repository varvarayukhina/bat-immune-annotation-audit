import csv, json, sys
from collections import defaultdict
from pathlib import Path
base, pair, hap, seqrep = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
chr_ = {r["genbank_accession"]: r.get("chr_name","?") for r in json.loads(seqrep.read_text())["reports"]}
loci = [r for r in csv.DictReader(open(base/"out"/f"{pair}.{hap}.loci.tsv"), delimiter="\t")]
def fam(sym):
    s = sym.upper()
    for p in ("KLRB1","KLRG1","KLRA","LAIR1","NCR2","FCAR"):
        if s.startswith(p): return p
rows = []
for line in (base/"mapping"/f"{hap}.sens.miniprot.gff").read_text().splitlines():
    if not line.startswith("##PAF"): continue
    f = line.split("\t"); tags = {t.split(":")[0]: t.split(":")[2] for t in f[13:]}
    q = f[1].split("|")
    cov = (int(f[4])-int(f[3]))/int(f[2])
    rows.append(dict(fam=fam(q[1]), q=q[0][:2]+"|"+q[1], contig=f[6], s=int(f[8])+1, e=int(f[9]), cov=cov, AS=int(tags["AS"]), ms=int(tags["ms"]), pos=int(tags["np"])/max(1,int(f[2])), fs=int(tags["fs"]), st=int(tags["st"])))
out = base/"out"/f"{pair}.{hap}.sensitivity_check.tsv"
with out.open("w") as o:
    o.write("family\tquery\tchr\tcontig\tstart\tend\tquery_cov\tAS\tpositives_per_qlen\tframeshifts\tinframe_stops\toverlaps_called_locus\n")
    by = defaultdict(list)
    for r in rows:
        ov = [l["locus_id"]+":"+l["family"] for l in loci if l["contig"]==r["contig"] and int(l["start"])<=r["e"] and int(l["end"])>=r["s"]]
        by[r["fam"]].append((r, ov))
        o.write("\t".join(map(str,[r["fam"],r["q"],chr_[r["contig"]],r["contig"],r["s"],r["e"],f'{r["cov"]:.2f}',r["AS"],f'{r["pos"]:.2f}',r["fs"],r["st"],",".join(ov) or "-"]))+"\n")
for fm in ("KLRB1","KLRG1","KLRA","LAIR1","NCR2","FCAR"):
    L = by.get(fm, [])
    print(fm, "hits:", len(L))
    for r, ov in sorted(L, key=lambda x: -x[0]["AS"])[:4]:
        print("   ", r["q"], f'chr{chr_[r["contig"]]}:{r["s"]}-{r["e"]}', f'cov={r["cov"]:.2f} AS={r["AS"]} fs={r["fs"]} st={r["st"]}', "overlaps:", ",".join(ov) or "none")

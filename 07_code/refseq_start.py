#!/usr/bin/env python3
"""Start-codon check with bat RefSeq proteins (evidence-based models of an annotated relative) instead of human/mouse queries.
usage: refseq_start.py pair_dir verify_dir relative_protein.faa
For every pilot-intact locus (+-5 kb from the pilot windows) the proteins of the same family from the relative species are aligned
with miniprot; the best hit's query start and first codon are reported. A locus whose human/mouse query did not align from residue 1
but whose bat query does, with ATG, has a bat-type N-terminus, not a defect.
outputs: verify_dir/refseq_start.tsv, refseq_start_summary.txt"""
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readcheck_lib import load_windows, rc, rd, window_of  # noqa: E402

pair, vdir, faa = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
tag = pair.name
tmp = vdir / "tmp"
tmp.mkdir(exist_ok=True)
FLANK = 5000
# RefSeq product names -> pilot families (class II, ligands and 'minor histocompatibility' excluded)
FAM = [("KLRD1", r"natural killer cells antigen CD94"), ("KLRK1", r"NKG2-D"), ("KLRC", r"NKG2-[ABCEF]"), ("KLRB1", r"lectin-like receptor subfamily B member 1"),
       ("KLRF", r"lectin-like receptor subfamily F member"), ("KLRG1", r"lectin-like receptor subfamily G member 1"), ("KLRG2", r"lectin-like receptor subfamily G member 2"),
       ("KLRA", r"killer cell lectin-like receptor 2|lectin-like receptor subfamily A"), ("CLEC2D", r"C-type lectin domain family 2 member [AD]|early activation antigen CD69"),
       ("CLEC12A", r"C-type lectin domain family 12 member A"), ("CLEC4E", r"C-type lectin domain family 4 member E"), ("KIR", r"killer cell immunoglobulin-like receptor"),
       ("LILR", r"leukocyte immunoglobulin-like receptor subfamily [AB] member"), ("NCR1", r"natural cytotoxicity triggering receptor 1(?! ligand)"),
       ("NCR2", r"natural cytotoxicity triggering receptor 2(?! ligand)"), ("NCR3", r"natural cytotoxicity triggering receptor 3(?! ligand)"),
       ("LAIR1", r"leukocyte-associated immunoglobulin-like receptor 1"), ("FCAR", r"immunoglobulin alpha Fc receptor"),
       ("MHC-I", r"class I histocompatibility antigen(?!.*class II)"), ("MR1", r"class I-related gene protein"), ("CD1", r"glycoprotein CD1[a-e]"),
       ("B2M", r"beta-2-microglobulin"), ("TAP1", r"antigen peptide transporter 1"), ("TAP2", r"antigen peptide transporter 2")]

# 1. proteins of the pilot families from the relative
prot, name, keep = {}, None, False
with open(faa) as fh:
    for line in fh:
        if line.startswith(">"):
            desc = line[1:].strip()
            fam = next((f for f, pat in FAM if re.search(pat, desc, re.I)), None)
            keep = fam is not None and "LOW QUALITY" not in desc and "minor histocompatibility" not in desc
            if keep:
                name = f"{fam}|{desc.split()[0]}"
                prot[name] = []
        elif keep:
            prot[name].append(line.strip())
prot = {k: "".join(v) for k, v in prot.items()}
qfa = tmp / "relative_queries.faa"
with open(qfa, "w") as f:
    for k, s in prot.items():
        f.write(f">{k}\n{s}\n")
fam_n = Counter(k.split("|")[0] for k in prot)

# 2. regions of pilot-intact loci
W, WSEQ = load_windows(vdir)
L = [lc for lc in rd(pair / "curation" / f"{tag}.locus_curation.tsv") if lc["pilot_status"] == "intact"]
regs, rfa = {}, tmp / "intact_regions.fa"
with open(rfa, "w") as f:
    for lc in L:
        loc = window_of(W, lc["haplotype"], lc["contig"], int(lc["start"]), int(lc["end"]))
        if not loc:
            continue
        w, a, b = loc
        a, b = max(1, a - FLANK), min(w["end"] - w["start"] + 1, b + FLANK)
        seq = WSEQ[w["name"]][a - 1:b]
        rid = f'{lc["haplotype"]}_{lc["locus_id"]}'
        regs[rid] = (lc, seq)
        f.write(f">{rid}\n{seq}\n")

# 3. miniprot, GFF parse: best hit of the locus family per region
gff = subprocess.run(["miniprot", "-t", "3", "-j2", "-p", "0.3", "-N", "1000", "--outs=0.3", "--outc=0.3", "--outn=1000", "--gff", str(rfa), str(qfa)],
                     capture_output=True, text=True, check=True).stdout
hits = {}  # (region, mRNA id) -> dict
for line in gff.splitlines():
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if len(f) < 9:
        continue
    attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    if f[2] == "mRNA":
        tq = attr["Target"].split()
        hits[(f[0], attr["ID"])] = {"q": tq[0], "qs": int(tq[1]), "qe": int(tq[2]), "score": int(f[5]), "strand": f[6], "ident": float(attr.get("Identity", 0)),
                                   "s": int(f[3]), "e": int(f[4]), "cds": [], "fs": int(attr.get("Frameshift", 0)), "st": int(attr.get("StopCodon", 0))}
    elif f[2] == "CDS":
        hits[(f[0], attr["Parent"])]["cds"].append((int(f[3]), int(f[4])))
    elif f[2] == "stop_codon":
        hits[(f[0], attr["Parent"])]["stop_feature"] = True
rows = []
for rid, (lc, seq) in regs.items():
    fam = lc["family"]
    cands = [h for (r, _), h in hits.items() if r == rid and h["q"].split("|")[0] == fam and h["cds"]]
    row = {"haplotype": lc["haplotype"], "locus_id": lc["locus_id"], "family": lc["family"], "pilot_start_codon": lc["start_codon"], "curated_intact": lc["curated_intact"],
           "curation_note": lc["curation_note"], "relative_queries_of_family": fam_n.get(fam, 0)}
    if not cands:
        row.update(bat_query="-", bat_q_start="-", bat_identity="-", bat_start_codon="-", bat_verdict="no_hit_of_family" if fam_n.get(fam) else "no_relative_protein_for_family")
        rows.append(row)
        continue
    # prefer hits that start at residue 1; among them the best score
    best = max(cands, key=lambda h: (h["qs"] == 1, h["score"]))
    cds = sorted(best["cds"])
    # sequence of the first exon read in the gene's direction, plus 60 bp upstream, to look for the nearest in-frame ATG without a stop in between
    if best["strand"] == "+":
        s0 = cds[0][0] - 1
        first = seq[s0:s0 + 3]
        up = seq[max(0, s0 - 60):s0]
        down = seq[s0:cds[0][1]]
    else:
        e0 = cds[-1][1]
        first = rc(seq[e0 - 3:e0])
        up = rc(seq[e0:e0 + 60])
        down = rc(seq[cds[-1][0] - 1:e0])
    near = "-"
    if first != "ATG" and best["qs"] == 1:
        for k in range(1, len(up) // 3 + 1):  # upstream, in frame
            cod = up[len(up) - 3 * k:len(up) - 3 * k + 3]
            if cod in ("TAA", "TAG", "TGA"):
                break
            if cod == "ATG":
                near = f"in-frame_ATG_{k}_codons_upstream"
                break
        if near == "-":
            for k in range(1, min(15, len(down) // 3)):
                cod = down[3 * k:3 * k + 3]
                if cod in ("TAA", "TAG", "TGA"):
                    break
                if cod == "ATG":
                    near = f"in-frame_ATG_{k}_codons_downstream"
                    break
        if near == "-":
            near = "no_in-frame_ATG_within_20_codons_up/15_down"
    row.update(bat_query=best["q"].split("|")[1], bat_q_start=best["qs"], bat_identity=f'{best["ident"]:.2f}', bat_frameshifts=best["fs"], bat_stops=best["st"],
               bat_start_codon=first if best["qs"] == 1 else f"not_aligned(q_start={best['qs']})", nearest_inframe_ATG=near)
    if best["qs"] == 1 and first == "ATG":
        row["bat_verdict"] = "ATG_from_residue_1"
    elif best["qs"] == 1:
        row["bat_verdict"] = f"non_ATG({first};{near})"
    elif best["qs"] <= 10:
        row["bat_verdict"] = f"starts_within_10_residues(q_start={best['qs']})"
    else:
        row["bat_verdict"] = f"N_terminus_missing(q_start={best['qs']})"
    rows.append(row)
cols = ["haplotype", "locus_id", "family", "pilot_start_codon", "curated_intact", "curation_note", "relative_queries_of_family", "bat_query", "bat_q_start", "bat_identity",
        "bat_frameshifts", "bat_stops", "bat_start_codon", "nearest_inframe_ATG", "bat_verdict"]
with open(vdir / "refseq_start.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "-")) for c in cols) + "\n")
fail = [r for r in rows if r["curated_intact"] == "no" and r["curation_note"].startswith("start:")]
out = [f"{tag}: relative {faa.parts[-4] if len(faa.parts) > 4 else faa.name}: {len(prot)} proteins of pilot families ({', '.join(f'{k}={v}' for k, v in sorted(fam_n.items()))})",
       f"intact loci {len(L)}, in windows {len(rows)}; bat-query verdicts: " + ", ".join(f"{k}={v}" for k, v in Counter(r["bat_verdict"].split("(")[0] for r in rows).most_common()),
       f"loci failing the human/mouse start-codon check ({len(fail)}): " + ", ".join(f"{k}={v}" for k, v in Counter(r["bat_verdict"].split("(")[0] for r in fail).most_common())]
for r in fail:
    out.append(f'  {r["haplotype"]} {r["locus_id"]} {r["family"]} pilot:{r["pilot_start_codon"]} -> bat {r["bat_query"]} q_start={r["bat_q_start"]} ident={r["bat_identity"]} codon={r["bat_start_codon"]} -> {r["bat_verdict"]}')
(vdir / "refseq_start_summary.txt").write_text("\n".join(out) + "\n")
print("\n".join(out))

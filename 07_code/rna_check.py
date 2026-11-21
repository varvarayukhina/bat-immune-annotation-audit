#!/usr/bin/env python3
"""RNA-seq evidence for the pilot loci (expression, splice junctions, upstream exons at problematic start codons) and for the
disruptive variants (which allele the transcripts of another individual / a related species carry).
usage: rna_check.py pair_dir verify_dir [--map]
  --map : stream verify_dir/rna/*.fastq.gz (mates as single-end; truncated prefixes allowed) through
          minimap2 -x splice:sr against windows_hap1.fa and windows_hap2.fa separately -> verify_dir/rna/hap{1,2}.paf
outputs: verify_dir/rna_loci.tsv, rna_sites.tsv, rna_summary.txt
Filters: primary, aligned block >= 50 bp, divergence <= 0.12 (cross-species reads), within 5 kb of a pilot locus."""
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readcheck_lib import CS, DIFFFLANK, VAR, fmt_counts, lift_table, load_paf, load_windows, rd, site_counts, variant_sides, window_of  # noqa: E402

pair, vdir = Path(sys.argv[1]), Path(sys.argv[2])
tag = pair.name
rna = vdir / "rna"
tmp = vdir / "tmp"
tmp.mkdir(exist_ok=True)
MIN_BLOCK, MAX_DE = 50, 0.12

if "--map" in sys.argv:
    fqs = sorted(str(p) for p in rna.glob("*.fastq.gz"))
    for h in ("hap1", "hap2"):
        paf = rna / f"{h}.paf"
        if paf.exists():
            continue
        loop = " ".join(f"gzip -dc '{f}' 2>/dev/null | paste - - - - | awk -F'\\t' 'NF==4 && length($2)==length($4) && length($2)>0';" for f in fqs)
        # no pipefail: gzip exits 1 on a truncated prefix although every complete record was emitted; minimap2's own status decides
        cmd = (f"( {loop} ) | tr '\\t' '\\n' | nice -n 10 minimap2 -x splice:sr -t 2 -c --cs --secondary=no '{vdir / f'windows_{h}.fa'}' - > '{paf}.part' "
               f"2> '{rna / f'{h}.minimap2.log'}'; st=${{PIPESTATUS[2]}}; [[ $st -eq 0 ]] && mv '{paf}.part' '{paf}'; exit $st")
        r = subprocess.run(["bash", "-c", cmd], check=False)
        if r.returncode != 0 or not paf.exists():
            sys.exit(f"minimap2 failed for {h} (exit {r.returncode}), see {rna / f'{h}.minimap2.log'}")
        print(f"mapped {h}: {sum(1 for _ in open(paf))} PAF records", flush=True)

W, WSEQ = load_windows(vdir)
# only reads within 5 kb of a pilot locus are used below; dropping the rest at load keeps memory flat
# (pair1: a blood library put 24 M reads on the first 50 kb of a chromosome end, 50 kb away from the nearest locus)
BIN, IV = 1000, defaultdict(set)
for lc in rd(pair / "curation" / f"{tag}.locus_curation.tsv"):
    loc = window_of(W, lc["haplotype"], lc["contig"], int(lc["start"]), int(lc["end"]))
    if loc:
        w, a, b = loc
        IV[w["name"]].update(range(max(0, a - 5000) // BIN, (b + 5000) // BIN + 1))


def near_locus(t, ts, te):
    bins = IV.get(t)
    return bool(bins) and any(x in bins for x in range(ts // BIN, te // BIN + 1))


R = {h: load_paf([rna / f"{h}.paf"], MIN_BLOCK, MAX_DE, near_locus) for h in ("hap1", "hap2")}
recs = {h: R[h][0] for h in R}


def introns(r):
    """(start1, end1) of every intron of a read, 1-based window coordinates"""
    t, out = r["ts"], []
    for tok in CS.findall(r["cs"]):
        if tok[0] == ":":
            t += int(tok[1:])
        elif tok[0] == "*":
            t += 1
        elif tok[0] == "-":
            t += len(tok) - 1
        elif tok[0] == "~":
            n = int(re.search(r"\d+", tok).group())
            out.append((t + 1, t + n))
            t += n
    return out


# ---------- loci ----------
lrows = []
for lc in rd(pair / "curation" / f"{tag}.locus_curation.tsv"):
    h = lc["haplotype"]
    loc = window_of(W, h, lc["contig"], int(lc["start"]), int(lc["end"]))
    row = {"haplotype": h, "locus_id": lc["locus_id"], "family": lc["family"], "pilot_status": lc["pilot_status"], "strand": lc["strand"], "start_codon": lc["start_codon"],
           "curated_intact": lc["curated_intact"], "n_introns_model": lc["n_introns"]}
    if not loc:
        row["rna_verdict"] = "outside_windows"
        lrows.append(row)
        continue
    w, a, b = loc
    rr = [r for r in recs[h].get(w["name"], []) if r["te"] >= a and r["ts"] + 1 <= b]
    jc = Counter(j for r in rr for j in introns(r) if j[1] >= a - 5000 and j[0] <= b + 5000)
    plus = lc["strand"] == "+"
    first_s, first_e = (a, a + 60) if plus else (b - 60, b)
    # an intron that ends (acceptor) just before/inside the model's first 60 bp, coming from upstream -> an exon exists upstream of the model start
    up_j = sum(n for (s, e), n in jc.items() if (plus and first_s - 10 <= e <= first_e and s < first_s - 10) or (not plus and first_s <= s <= first_e + 10 and e > first_e + 10))
    cov_first = sum(1 for r in rr if r["ts"] + 1 <= first_e and r["te"] >= first_s)
    cov_up = sum(1 for r in recs[h].get(w["name"], []) if (plus and r["ts"] + 1 <= a - 1 and r["te"] >= a - 100) or (not plus and r["ts"] + 1 <= b + 100 and r["te"] >= b + 1))
    row.update(rna_reads=len(rr), rna_spliced_reads=sum(1 for r in rr if "~" in r["cs"]), rna_junctions_ge2=sum(1 for n in jc.values() if n >= 2),
               reads_over_model_start=cov_first, reads_100bp_upstream_of_model=cov_up, junction_reads_into_model_start=up_j)
    if len(rr) < 3:
        v = "not_expressed_or_not_covered"
    elif up_j >= 2:
        v = "upstream_exon_spliced_into_model_start(real_start_is_upstream)"
    elif cov_up >= 3 and cov_first >= 3:
        v = "transcript_extends_upstream_of_model_start"
    elif cov_first >= 3:
        v = "transcript_starts_near_model_start"
    else:
        v = "expressed_but_5prime_end_not_covered"
    row["rna_verdict"] = v
    lrows.append(row)
LCOLS = ["haplotype", "locus_id", "family", "pilot_status", "strand", "start_codon", "curated_intact", "n_introns_model", "rna_reads", "rna_spliced_reads", "rna_junctions_ge2",
         "reads_over_model_start", "reads_100bp_upstream_of_model", "junction_reads_into_model_start", "rna_verdict"]
with open(vdir / "rna_loci.tsv", "w") as f:
    f.write("\t".join(LCOLS) + "\n")
    for r in lrows:
        f.write("\t".join(str(r.get(c, "-")) for c in LCOLS) + "\n")

# ---------- variant sites ----------
srows = []
for r in rd(pair / "curation" / f"{tag}.difference_curation.tsv"):
    if not r["final_class"].startswith("status_differs_confirmed") or r["variants_in_cds_of_better_copy"] in ("-", ""):
        continue
    ref_h, ref_loc, _rs, alt_h, alt_loc, _as, classes = variant_sides(r)
    a_ref, a_alt = max(1, ref_loc[1] - DIFFFLANK), max(1, alt_loc[1] - DIFFFLANK)
    wr = window_of(W, ref_h, ref_loc[0], a_ref, ref_loc[2] + DIFFFLANK)
    wa = window_of(W, alt_h, alt_loc[0], a_alt, alt_loc[2] + DIFFFLANK)
    base = {"diff_index": r["diff_index"], "locus": r["locus_id"], "family": r["family"], "classes": classes, "better_hap": ref_h, "alt_hap": alt_h}
    if not wr or not wa:
        srows.append({**base, "variant": r["variants_in_cds_of_better_copy"], "rna_verdict": "outside_windows"})
        continue
    lift, strand, _ = lift_table(WSEQ[wr[0]["name"]][wr[1] - 1:wr[2]], WSEQ[wa[0]["name"]][wa[1] - 1:wa[2]], tmp, f"rlift{r['diff_index']}")
    for item in r["variants_in_cds_of_better_copy"].split(";"):
        mv = VAR.match(item)
        if not mv:
            continue
        _ex, _rb, _ab, spos, _stop, ikind, ilen, ipos, hp = mv.groups()
        pos1 = int(spos or ipos)
        n = int(ilen) if ilen else 1
        if ikind and n > 50:
            continue
        slop = max(3, (int(hp) if hp else 0) + 2) if ikind else 2
        cr = site_counts(recs[ref_h].get(wr[0]["name"], []), wr[1] + pos1 - 1, slop, mapq=1)
        alt_p = lift(pos1)
        row = {**base, "variant": item, "rna_at_better_copy": fmt_counts(cr), "rna_at_alt_copy": "lift_failed"}
        if alt_p is not None:
            ca = site_counts(recs[alt_h].get(wa[0]["name"], []), wa[1] + alt_p - 1, slop, mapq=1)
            row["rna_at_alt_copy"] = fmt_counts(ca)
            if ikind:
                alt_allele_at_ref, better_allele_at_alt = (cr["ins"], ca["del"]) if ikind.endswith("ins") else (cr["del"], ca["ins"])
            else:
                alt_allele_at_ref, better_allele_at_alt = cr["sub"], ca["sub"]
            tot = cr["n"] + ca["n"]
            if tot < 3:
                v = "no_RNA_coverage"
            else:
                # reads matching the alt copy or showing the alt allele on the better copy carry the disruptive allele
                carry_alt = ca["match"] + alt_allele_at_ref
                carry_better = cr["match"] + better_allele_at_alt
                v = f"transcripts:{carry_alt}_with_disruptive_allele,{carry_better}_with_intact_allele"
            row["rna_verdict"] = v
        else:
            row["rna_verdict"] = "lift_failed"
        srows.append(row)
SCOLS = ["diff_index", "locus", "family", "classes", "better_hap", "alt_hap", "variant", "rna_at_better_copy", "rna_at_alt_copy", "rna_verdict"]
with open(vdir / "rna_sites.tsv", "w") as f:
    f.write("\t".join(SCOLS) + "\n")
    for r in srows:
        f.write("\t".join(str(r.get(c, "-")) for c in SCOLS) + "\n")

# ---------- summary ----------
out = [f"{tag}: RNA PAF records hap1 {R['hap1'][1]} kept {R['hap1'][2]}; hap2 {R['hap2'][1]} kept {R['hap2'][2]}"]
for st in ("intact", "partial", "pseudogene-like", "uncertain"):
    c = Counter(r["rna_verdict"].split("(")[0] for r in lrows if r["pilot_status"] == st)
    out.append(f"{st}: " + ", ".join(f"{k}={v}" for k, v in c.most_common()))
sc = [r for r in lrows if r["pilot_status"] == "intact" and r["curated_intact"] == "no" and r["start_codon"] not in ("ATG", "-", "")]
out.append(f"intact loci with a start-codon problem ({len(sc)}): " + ", ".join(f"{k}={v}" for k, v in Counter(r["rna_verdict"].split("(")[0] for r in sc).most_common()))
for r in sc:
    out.append(f'  {r["haplotype"]} {r["locus_id"]} {r["family"]} {r["start_codon"]} reads={r.get("rna_reads")} spliced={r.get("rna_spliced_reads")} over_start={r.get("reads_over_model_start")} '
               f'upstream100={r.get("reads_100bp_upstream_of_model")} junction_into_start={r.get("junction_reads_into_model_start")} -> {r["rna_verdict"]}')
fam = defaultdict(lambda: [0, 0])
for r in lrows:
    fam[r["family"]][0] += 1
    fam[r["family"]][1] += r["rna_verdict"] not in ("not_expressed_or_not_covered", "outside_windows")
out.append("expressed loci by family (any RNA reads >= 3): " + ", ".join(f"{k}={v[1]}/{v[0]}" for k, v in sorted(fam.items())))
out.append("variant sites: " + ", ".join(f"{k}={v}" for k, v in Counter(r["rna_verdict"].split(":")[0] for r in srows).most_common()))
for r in srows:
    out.append(f'  {r["locus"]} {r["family"]} {r["classes"]} {r["variant"]} | better {r["better_hap"]}: {r.get("rna_at_better_copy")} | alt {r["alt_hap"]}: {r.get("rna_at_alt_copy")} -> {r["rna_verdict"]}')
(vdir / "rna_summary.txt").write_text("\n".join(out) + "\n")
print("\n".join(out))

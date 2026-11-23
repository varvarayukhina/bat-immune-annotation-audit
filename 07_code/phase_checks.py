#!/usr/bin/env python3
"""Two checks on whether reported haplotype differences could be phasing artefacts.

Check A, counterpart identity. A switch error puts the same parental sequence in both
haplotypes, so a locus and its counterpart would be near-identical. Summarising the
alignment_identity column of the difference curation by family therefore bounds how much
of each species' difference set a switch error could explain.

Check B, within-array paralogue identity for KLRD1. Redundant calls on one locus would be
near-identical to each other; a real array of paralogues would not. All pairwise protein
identities among the KLRD1 models are computed from the deposited model proteins, and the
cross-haplotype pairs are reported separately because those are the allelic pairs.

Scope: check B is run for KLRD1 only, the family where a published count conflicts with the
annotation. Running it on MHC-I would mean thousands of pairwise alignments in pure Python.
"""

from __future__ import annotations

import csv
import itertools
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAIRS = ("pair1", "pair2", "pair3")
SPECIES = {"pair1": "Myotis sodalis", "pair2": "Rhinolophus clivosus", "pair3": "Eidolon dupreanum"}
OUT = ROOT / "07_code" / "derived"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    buf: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def identity(a: str, b: str, gap: int = -2) -> float:
    """Global alignment identity: matches divided by alignment length."""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = gap * i
    for j in range(1, m + 1):
        score[0][j] = gap * j
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = score[i], score[i - 1]
        for j in range(1, m + 1):
            row[j] = max(
                prev[j - 1] + (1 if ai == b[j - 1] else -1),
                prev[j] + gap,
                row[j - 1] + gap,
            )
    i, j, same, length = n, m, 0, 0
    while i > 0 and j > 0:
        if score[i][j] == score[i - 1][j - 1] + (1 if a[i - 1] == b[j - 1] else -1):
            same += a[i - 1] == b[j - 1]
            length += 1
            i -= 1
            j -= 1
        elif score[i][j] == score[i - 1][j] + gap:
            length += 1
            i -= 1
        else:
            length += 1
            j -= 1
    length += i + j
    return same / length if length else 0.0


def check_a() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pair in PAIRS:
        by_family: dict[str, list[float]] = {}
        for record in read_tsv(ROOT / "05_tables" / "curation" / f"{pair}.difference_curation.tsv"):
            raw = record.get("alignment_identity", "")
            if raw in ("", "-", "NA"):
                continue
            by_family.setdefault(record["family"], []).append(float(raw))
        for family, values in sorted(by_family.items()):
            values.sort()
            rows.append(
                {
                    "pair": pair,
                    "species": SPECIES[pair],
                    "family": family,
                    "n_pairs": len(values),
                    "min_identity": f"{values[0]:.4f}",
                    "median_identity": f"{statistics.median(values):.4f}",
                    "max_identity": f"{values[-1]:.4f}",
                    "n_ge_0.99": sum(1 for v in values if v >= 0.99),
                }
            )
    return rows


def check_b(family: str = "KLRD1") -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    pairs_out: list[dict[str, object]] = []
    summary: dict[str, dict[str, object]] = {}
    for pair in PAIRS:
        proteins = read_fasta(ROOT / "05_tables" / "curation" / f"{pair}.locus_proteins.faa")
        models = {k: v for k, v in proteins.items() if f"|{family}|" in k}
        if len(models) < 2:
            summary[pair] = {"n_models": len(models), "n_pairs": 0}
            continue
        values: list[float] = []
        for a, b in itertools.combinations(sorted(models), 2):
            value = identity(models[a], models[b])
            values.append(value)
            hap_a, hap_b = a.split("|")[0], b.split("|")[0]
            pairs_out.append(
                {
                    "pair": pair,
                    "species": SPECIES[pair],
                    "family": family,
                    "model_a": a.replace("|", "_"),
                    "model_b": b.replace("|", "_"),
                    "cross_haplotype": "yes" if hap_a != hap_b else "no",
                    "identity": f"{value:.4f}",
                }
            )
        cross = [
            float(r["identity"])
            for r in pairs_out
            if r["pair"] == pair and r["cross_haplotype"] == "yes"
        ]
        summary[pair] = {
            "n_models": len(models),
            "n_pairs": len(values),
            "min": min(values),
            "median": statistics.median(values),
            "max": max(values),
            "n_ge_095": sum(1 for v in values if v >= 0.95),
            "n_cross_ge_095": sum(1 for v in cross if v >= 0.95),
        }
    return pairs_out, summary


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    rows_a = check_a()
    write_tsv(OUT / "counterpart_identity_by_family.tsv", rows_a)
    print("== check A: counterpart identity by family ==")
    for pair in PAIRS:
        subset = [r for r in rows_a if r["pair"] == pair]
        if not subset:
            continue
        worst = max(subset, key=lambda r: float(str(r["max_identity"])))
        total_hi = sum(int(str(r["n_ge_0.99"])) for r in subset)
        total_n = sum(int(str(r["n_pairs"])) for r in subset)
        print(
            f"  {SPECIES[pair]:22} n={total_n:3}  counterparts >=0.99: {total_hi:3}"
            f"  highest: {worst['family']} at {worst['max_identity']}"
        )

    rows_b, summary = check_b()
    write_tsv(OUT / "klrd1_pairwise_identity.tsv", rows_b)
    print("\n== check B: pairwise protein identity among KLRD1 models ==")
    for pair in PAIRS:
        s = summary.get(pair, {})
        if not s.get("n_pairs"):
            print(f"  {SPECIES[pair]:22} {s.get('n_models', 0)} model(s), no pairs")
            continue
        print(
            f"  {SPECIES[pair]:22} {s['n_models']} models, {s['n_pairs']} pairs  "
            f"min={s['min']:.3f} median={s['median']:.3f} max={s['max']:.3f}  "
            f">=0.95: {s['n_ge_095']} (cross-haplotype {s['n_cross_ge_095']})"
        )

    print(f"\nwrote {OUT / 'counterpart_identity_by_family.tsv'}")
    print(f"wrote {OUT / 'klrd1_pairwise_identity.tsv'}")


def demo() -> None:
    assert identity("ACDEF", "ACDEF") == 1.0
    assert identity("ACDEF", "ACDEG") == 0.8
    assert 0.0 <= identity("ACDEF", "WWWWW") < 0.3
    assert identity("ACDEF", "ACDE") == 0.8  # one gap in a 5-column alignment
    print("identity() self-check passed")


if __name__ == "__main__":
    import sys

    if "--demo" in sys.argv:
        demo()
    else:
        main()

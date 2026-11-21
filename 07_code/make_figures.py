#!/usr/bin/env python3
"""Figures for the short report, drawn from the verification tables of the pilot.

Run: uv run --with matplotlib --with scipy python manuscript/make_figures.py
Writes manuscript/figures/fig1..fig5 (png + pdf, 300 dpi) and prints the numbers used in the captions.

Design rules applied (dataviz skill): validated categorical palette, emphasis instead of
rainbow, hairline solid grid, thin marks with a 2 px surface gap, selective direct labels,
text in ink tokens and never in a series colour, no dual axes, no number on every point.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)
PAIRS = [("pair1", "Myotis sodalis"), ("pair2", "Rhinolophus clivosus"), ("pair3", "Eidolon dupreanum")]

# validated categorical slots (light surface #fcfcfb) + neutrals
BLUE, ORANGE, AQUA, YELLOW, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6"
BLUE_LT = "#9ec5f4"          # same hue, lighter step: the two ends of a dumbbell
INK, INK2, MUTED = "#15171a", "#52514e", "#8b8f96"
GRID, GREY, SURF = "#e3e3df", "#c7c9cc", "#ffffff"

_have = {f.name for f in font_manager.fontManager.ttflist}
SANS = next((f for f in ("Helvetica Neue", "Avenir Next", "Inter", "Source Sans 3", "Arial") if f in _have),
            "DejaVu Sans")
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [SANS, "DejaVu Sans"], "font.size": 8.5,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.facecolor": SURF,
    "axes.facecolor": SURF, "figure.facecolor": SURF,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
    "axes.edgecolor": GREY, "axes.linewidth": 0.8, "axes.labelcolor": INK2, "axes.labelsize": 8.5,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "xtick.major.size": 3, "ytick.major.size": 0, "xtick.major.width": 0.8,
    "legend.frameon": False, "legend.fontsize": 7.8, "legend.handlelength": 1.1,
    "lines.solid_capstyle": "round",
})
W1, W2 = 3.46, 7.2  # single- and double-column widths in inches


def rd(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def strip(v: str) -> str:
    return v.split("(")[0].split(":")[0]


def panel(ax, letter: str, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="bottom")


def save(fig, name: str) -> None:
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)


def grid_x(ax) -> None:
    ax.grid(axis="x", color=GRID, linewidth=0.7, linestyle="-")
    ax.set_axisbelow(True)


def grid_y(ax) -> None:
    ax.grid(axis="y", color=GRID, linewidth=0.7, linestyle="-")
    ax.set_axisbelow(True)


# ------------------------------------------------------------------ Fig 1: observed haplotype differences
# Per family, the number of loci called in hap1 against hap2 of the SAME animal.
# Dumbbell: one row per family that differs; families that agree are summarised in the subtitle.
fam_counts: dict[str, dict[str, dict[str, int]]] = {}
for pair, _ in PAIRS:
    d: dict[str, dict[str, int]] = defaultdict(dict)
    for r in rd(ROOT / pair / "out" / f"{pair}.family_summary.tsv"):
        d[r["family"]][r["hap"]] = int(r["loci"])
    fam_counts[pair] = d

fig, axes = plt.subplots(1, 3, figsize=(W2, 3.1), sharex=True)
fig1_note = []
xmax1 = max(max(v.get("hap1", 0), v.get("hap2", 0)) for d in fam_counts.values() for v in d.values())
for ax, (pair, sp) in zip(axes, PAIRS):
    fams = fam_counts[pair]
    diff = sorted(((f, v.get("hap1", 0), v.get("hap2", 0)) for f, v in fams.items() if v.get("hap1", 0) != v.get("hap2", 0)),
                  key=lambda t: abs(t[1] - t[2]))
    same = len(fams) - len(diff)
    fig1_note.append((sp, len(fams), same, [(f, a, b) for f, a, b in diff]))
    ys = range(len(diff))
    for y, (f, a, b) in zip(ys, diff):
        ax.plot([a, b], [y, y], color=GREY, lw=1.6, zorder=1, solid_capstyle="round")
        ax.scatter([a], [y], s=34, color=BLUE_LT, zorder=3, edgecolor=SURF, linewidths=1.2)
        ax.scatter([b], [y], s=34, color=BLUE, zorder=3, edgecolor=SURF, linewidths=1.2)
        ax.annotate(f"{abs(a - b)}", (max(a, b), y), textcoords="offset points", xytext=(7, -2.6),
                    fontsize=7.4, color=INK2)
    ax.set_yticks(list(ys), [f for f, _, _ in diff], fontsize=7.6)
    ax.set_xlim(0, xmax1 * 1.16)
    ax.set_ylim(-0.8, max(len(diff) - 0.2, 1.2))
    ax.set_title(f"{sp}\n{len(diff)} of {len(fams)} families differ", style="italic", fontsize=8.4,
                 color=INK, pad=7)
    grid_x(ax)
axes[1].set_xlabel("Loci called in the region, per haplotype of one individual")
axes[0].legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=BLUE_LT, mec=SURF, label="haplotype 1"),
                        Line2D([], [], marker="o", ls="", ms=6, mfc=BLUE, mec=SURF, label="haplotype 2")],
               loc="lower right", bbox_to_anchor=(1.02, -0.02))
save(fig, "fig1_haplotype_counts")

# ------------------------------------------------------------------ Fig 2: fate of the differences
CLASSES = [
    ("artefact_sensitivity_or_threshold", "Same family and status\non renewed search", BLUE),
    ("status_differs_no_disruptive_variant_found", "Status differs, no candidate\nsequence change found", YELLOW),
    ("status_differs_with_disruptive_variant", "Status differs, candidate\nsequence change found", ORANGE),
    ("status_differs_large_indel_in_gene", "Large deletion in a gene model", AQUA),
    ("locus_absent_from_other_haplotype", "Sequence absent from\npartner assembly", VIOLET),
]
fig2_counts = {p: Counter(strip(r["refined_class"]) for r in rd(ROOT / "summary" / f"{p}.refined_differences.tsv"))
               for p, _ in PAIRS}
fig, ax = plt.subplots(figsize=(W2, 2.25))
ypos = list(range(len(PAIRS)))
left = [0.0] * len(PAIRS)
for key, label, colour in CLASSES:
    vals = [fig2_counts[p][key] for p, _ in PAIRS]
    ax.barh(ypos, vals, left=left, height=0.46, color=colour, label=label, edgecolor=SURF, linewidth=1.6)
    for i, v in enumerate(vals):
        if v >= 6:  # only where the number fits inside the block with padding
            ax.text(left[i] + v / 2, i, str(v), ha="center", va="center", color="white",
                    fontsize=8, fontweight="semibold")
    left = [a + b for a, b in zip(left, vals)]
ax.set_yticks(ypos, [f"{sp}\nn = {int(sum(fig2_counts[p].values()))}" for p, sp in PAIRS], style="italic", fontsize=8.2)
ax.invert_yaxis()
ax.set_xlabel("Locus-level differences reported between the two haplotypes of one individual")
grid_x(ax)
ax.legend(loc="upper left", bbox_to_anchor=(1.015, 1.05), labelspacing=0.85)
save(fig, "fig2_differences")

# ------------------------------------------------------------------ Fig 3: depth
depth_tests = rd(ROOT / "manuscript" / "derived" / "depth_null_fdr.tsv")
zero_floor = 0.013
fig = plt.figure(figsize=(W2, 3.5))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.65], wspace=0.48)
axes = [fig.add_subplot(gs[0, 0])]
axes.extend(fig.add_subplot(gs[0, i], sharey=axes[0]) for i in (1, 2))
for k, (ax, (pair, sp)) in enumerate(zip(axes, PAIRS)):
    rows = [r for r in depth_tests if r["pair"] == pair]
    ax.axhspan(0.55, 1.7, color=BLUE_LT, alpha=0.27, lw=0)
    ax.axhline(1.0, color=GREY, lw=0.9)
    zeros, outliers = 0, 0
    for i, r in enumerate(rows):
        x = (i % 15) * 0.058 - 0.42
        v = float(r["ratio_mapq20"])
        out = r["significant"] == "yes"
        outliers += out
        zeros += v == 0
        ax.scatter([x], [v if v else zero_floor], s=13 if out else 9, zorder=3 if out else 2,
                   color=ORANGE if out else MUTED, edgecolor=SURF, linewidths=0.8 if out else 0.5)
        if v >= 4:
            suffix = "*" if pair == "pair2" and r["locus"] == "hap1_L008" else ""
            ax.annotate(f"{r['family']} {v:.0f}×{suffix}", (x, v), textcoords="offset points",
                        xytext=(0, 5), ha="center", fontsize=7, color=INK)
    print(f"FIG3a {pair}: {len(rows)} evaluated, {outliers} with nominal BH q < 0.05, {zeros} zeros among evaluated")
    ax.set_yscale("log")
    ax.set_xlim(-0.61, 0.61)
    ax.set_xticks([])
    ax.set_title(f"{sp}\n{len(rows)} loci", style="italic", fontsize=8.4, color=INK, pad=7)
    grid_y(ax)
    ax.minorticks_off()
    ax.set_ylim(0.01, 116)
    if k:
        ax.tick_params(labelleft=False)
axes[0].set_ylabel("Normalized locus depth (MAPQ ≥20)")
axes[0].set_yticks([zero_floor, 0.05, 0.2, 0.5, 1, 2, 5, 10, 30, 60],
                   ["0", "0.05", "0.2", "0.5", "1", "2", "5", "10", "30", "60"])

axb = fig.add_subplot(gs[0, 3])
cn = rd(ROOT / "manuscript" / "derived" / "copies_implied_vs_assembled.tsv")
mx = max(max(int(r["assembled_loci"]), float(r["implied_copies"])) for r in cn) * 1.12
axb.plot([0, mx], [0, mx], color=GREY, lw=0.9, zorder=1)
exceptions = []
for r in cn:
    a, b = int(r["assembled_loci"]), float(r["implied_copies"])
    far = abs(b - a) > 2
    if far:
        exceptions.append(r)
    withdrawn = (r["pair"], r["family"], r["haplotype"]) == ("pair2", "MHC-I", "hap1")
    color = ORANGE if far else MUTED
    axb.scatter([a], [b], s=27 if far else 13, marker="D" if withdrawn else "o",
                facecolor=SURF if withdrawn else color, edgecolor=color if withdrawn else SURF,
                linewidths=1 if withdrawn else 0.5, zorder=3 if far else 2)
    if b - a > 2:
        label = r["family"] + ("*" if withdrawn else "")
        axb.annotate(label, (a, b), textcoords="offset points", xytext=(4, 3), fontsize=7, color=INK)
axb.annotate("CD1", (10, 4.4), textcoords="offset points", xytext=(6, -2), fontsize=7, color=INK)
axb.set_xlabel("Annotated loci with usable depth")
axb.set_ylabel("Sum of normalized depths")
axb.set_xlim(-1.8, mx)
axb.set_ylim(-1.8, mx)
axb.spines["left"].set_visible(True)
axb.spines["left"].set_color(GREY)
axb.grid(color=GRID, linewidth=0.7)
axb.set_axisbelow(True)
axb.tick_params(axis="y", length=3)
panel(axes[0], "a", x=-0.34)
panel(axb, "b", x=0)
fig.legend(handles=[Patch(facecolor=BLUE_LT, alpha=0.27, label="0.55–1.7 corridor (a)"),
                    Line2D([], [], marker="o", ls="", color=ORANGE, label="BH q < 0.05 (a)"),
                    Line2D([], [], marker="o", ls="", color=ORANGE, label="Absolute deviation >2 (b)")],
           loc="lower center", bbox_to_anchor=(0.5, -0.09), ncol=3, columnspacing=1.25)
save(fig, "fig3_depth")

# ------------------------------------------------------------------ Fig 4: start codons
fig4 = {}
for pair, _ in PAIRS:
    cur = rd(ROOT / pair / "curation" / f"{pair}.locus_curation.tsv")
    flagged = [r["locus_id"] for r in cur if r["pilot_status"] == "intact" and r["start_codon"] not in ("ATG", "-", "")]
    rs = {r["locus_id"]: strip(r["bat_verdict"]) for r in rd(ROOT / "verify" / pair / "refseq_start.tsv")}
    rna = {r["locus_id"]: strip(r["rna_verdict"]) for r in rd(ROOT / "verify" / pair / "rna_loci.tsv")}
    protein = Counter(rs.get(i) for i in flagged)
    transcript = Counter(rna.get(i) for i in flagged)
    fig4[pair] = {
        "intact": sum(r["pilot_status"] == "intact" for r in cur), "flagged": len(flagged),
        "atg": protein["ATG_from_residue_1"], "near_start": protein["starts_within_10_residues"],
        "other_protein": len(flagged) - protein["ATG_from_residue_1"] - protein["starts_within_10_residues"],
        "junction": transcript["upstream_exon_spliced_into_model_start"],
        "coverage": transcript["transcript_extends_upstream_of_model_start"],
    }

fig, (axa, axb) = plt.subplots(1, 2, figsize=(W2, 2.5), sharex=True, sharey=True, layout=None)
fig.subplots_adjust(wspace=0.20)
ypos = list(range(len(PAIRS)))
for ax, title, categories in (
    (axa, "Relative’s RefSeq protein", (("atg", BLUE, "ATG at query residue 1"),
                                      ("near_start", YELLOW, "Starts at query residue 6"),
                                      ("other_protein", GRID, "Other alignment outcome"))),
    (axb, "RNA from relatives", (("junction", ORANGE, "Upstream splice junction"),
                                 ("coverage", BLUE_LT, "Separate coverage counts"))),
):
    left = [0] * len(PAIRS)
    for key, color, description in categories:
        counts = [fig4[pair][key] for pair, _ in PAIRS]
        ax.barh(ypos, counts, left=left, height=0.51, color=color, edgecolor=SURF, linewidth=1.0,
                label=description)
        for y, x, count in zip(ypos, left, counts):
            if count:
                ax.text(x + count / 2, y, str(count), ha="center", va="center", fontsize=7,
                        color=SURF if color == BLUE else INK)
        left = [x + count for x, count in zip(left, counts)]
    ax.set_title(title, fontsize=9, pad=9)
    ax.set_xlabel("Flagged models")
    ax.set_xlim(0, 26)
    ax.set_xticks(range(0, 26, 5))
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.22))
    grid_x(ax)
axa.set_yticks(ypos, [f"{sp}\n{fig4[pair]['flagged']} flagged" for pair, sp in PAIRS],
               style="italic", fontsize=8.2)
axa.invert_yaxis()
panel(axa, "a")
panel(axb, "b")
save(fig, "fig4_starts")

# ------------------------------------------------------------------ Fig 5: transcripts
fig = plt.figure(figsize=(W2, 3.9))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.42)
axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
expr: dict[str, dict[str, tuple[int, int]]] = {}
for pair, _ in PAIRS:
    t, e = Counter(), Counter()
    for r in rd(ROOT / "verify" / pair / "rna_loci.tsv"):
        t[r["family"]] += 1
        if r["rna_verdict"] not in ("not_expressed_or_not_covered", "outside_windows"):
            e[r["family"]] += 1
    expr[pair] = {f: (e[f], t[f]) for f in t}
fams = sorted({f for p in expr for f in expr[p]}, reverse=True)
MARK = {"pair1": "o", "pair2": "s", "pair3": "^"}
for k, (pair, sp) in enumerate(PAIRS):
    dy = (k - 1) * 0.22
    for f in fams:
        if f not in expr[pair]:
            continue
        got, tot = expr[pair][f]
        frac = got / tot
        full = frac >= 0.999
        axa.scatter([frac], [fams.index(f) + dy], marker=MARK[pair], s=20,
                    facecolor=SURF if full else ORANGE, edgecolor=GREY if full else ORANGE,
                    linewidths=1.0, zorder=3)
axa.set_yticks(range(len(fams)), fams, fontsize=7.2)
axa.set_xlim(-0.08, 1.1)
axa.set_xticks([0, 0.5, 1], ["0", "50", "100"])
axa.set_ylim(-0.9, len(fams) - 0.1)
axa.set_xlabel("Models with ≥3 RNA alignments (%)")
axa.grid(color=GRID, linewidth=0.7)
axa.set_axisbelow(True)
axa.legend(handles=[Line2D([], [], marker=MARK[p], ls="", ms=5.5, mfc=SURF, mec=GREY, label=sp)
                    for p, sp in PAIRS] +
                   [Line2D([], [], marker="o", ls="", ms=5.5, mfc=ORANGE, mec=ORANGE, label="not all loci covered")],
           loc="upper center", bbox_to_anchor=(0.5, -0.15), prop={"size": 7.2}, labelspacing=0.45)

kl = rd(ROOT / "manuscript" / "derived" / "klrd1_junction_mapq.tsv")
idx = list(range(len(kl)))
axb.barh(idx, [int(r["junction_alignments"]) for r in kl], height=0.58, color=BLUE_LT,
         label="Before MAPQ filter")
axb.barh(idx, [int(r["junction_alignments_mapq20"]) for r in kl], height=0.58, color=ORANGE,
         label="After MAPQ ≥20")
axb.axvline(2, color=GREY, linewidth=0.9, linestyle="--")
axb.set_yticks(idx, [f'{r["haplotype"]} {r["locus_id"].split("_")[-1]}' for r in kl], fontsize=7)
axb.invert_yaxis()
axb.set_xlabel("Upstream junction alignments")
axb.set_title("KLRD1 models, Eidolon dupreanum", style="italic", fontsize=8.4, color=INK, pad=7)
grid_x(axb)
axb.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), labelspacing=0.5, ncol=1)
panel(axa, "a", x=-0.30)
panel(axb, "b", x=-0.24)
save(fig, "fig5_rna")

# ------------------------------------------------------------------ caption numbers
print("FIG1 observed haplotype differences per family")
for sp, nfam, same, diff in fig1_note:
    print(f"  {sp}: {len(diff)} of {nfam} families differ between the haplotypes; {same} identical; "
          f"largest {', '.join(f'{f} {a}v{b}' for f, a, b in sorted(diff, key=lambda t: -abs(t[1]-t[2]))[:3])}")
print("FIG2", {p: dict(fig2_counts[p]) for p, _ in PAIRS})
print(f"FIG3b: {len(cn) - len(exceptions)} of {len(cn)} groups within two units")
for r in exceptions:
    print(f"   {r['species']:24s} {r['family']:8s} {r['haplotype']} "
          f"{r['assembled_loci']:>3s} -> {r['implied_copies']:>5s}")
print("FIG4", fig4)
print("FIG5 silent or partial families:",
      {p: {f: expr[p][f] for f in expr[p] if expr[p][f][0] < expr[p][f][1]} for p, _ in PAIRS})
print("FIG5 KLRD1 pair3: models", len(kl), "with >=2 junction alignments",
      sum(int(r["junction_alignments"]) >= 2 for r in kl),
      "after MAPQ >=20", sum(int(r["junction_alignments_mapq20"]) >= 2 for r in kl))
print("font:", SANS)

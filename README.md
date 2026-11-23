# Auditing immune gene-family annotations across phased assemblies of three bat species

Code, parameter files, derived result tables and figure code for the manuscript
[`01_manuscript.md`](01_manuscript.md).

The study asks how much the number of annotated immune-receptor loci in a bat genome depends on
which haplotype of a phased assembly is annotated and on which annotation pass is used. Three
haplotype-resolved assemblies are annotated twice each, every count difference between the two
haplotypes is curated by locus re-alignment, and the surviving differences are checked against
PacBio HiFi reads from the same individual, RNA-seq from relatives and RefSeq proteins of
annotated relatives.

## Reproducing the figures and statistics

Everything downstream of the alignments runs from the tables in this repository, with no raw data:

```sh
./reproduce.sh
```

It needs a Python with `matplotlib` and `scipy`. If they are not installed, let `uv` supply them
without installing anything:

```sh
PYTHON="uv run --with matplotlib --with scipy python3" ./reproduce.sh
```

The script writes the derived tables to `07_code/derived/` and figures 1–5 to `07_code/figures/`.
All five regenerated tables reproduce byte for byte against the deposited copies in
`05_tables/derived/`: `copies_implied_vs_assembled`, `depth_null_fdr`,
`haplotype_count_differences`, `counterpart_identity_by_family` and `klrd1_pairwise_identity`.

One step does not run from this repository: section 4 of `analysis_extra.py` correlates the depth
ratio with GC content and locus length, so it reads the window sequences
(`verify/<pair>/windows_<hap>.fa`, 53 MB of assembly sequence). Those are not deposited here;
rebuild them with `07_code/fetch_windows.py` from the public assemblies, or use the deposited
result in `05_tables/derived/confounders_gc_length.tsv`. The script reports this and continues.

Rebuilding the tables themselves from raw data needs the public assemblies and reads listed below,
roughly 16 GB on disk, and the driver scripts `07_code/run_pair.sh` and `07_code/verify_driver.sh`.
The paths at the top of those scripts point at a local working directory and have to be set first.

## What is here

| Path | Contents |
|---|---|
| `01_manuscript.md` | Manuscript text, figure legends, Table 1, references |
| `02_methods_full.txt` | Full methods, including every command line and parameter |
| `04_numbers/` | Raw stdout of the analysis steps, and `digest.txt`, the per-species accounting of every count quoted in the text |
| `05_tables/summary/` | Per-family locus counts per haplotype, for the default and the sensitive pass, and the curated difference classes |
| `05_tables/curation/` | Per-locus curation, difference curation, gap curation and back-search tables, one set per species, plus the model protein sequences (`*.locus_proteins.faa`) |
| `05_tables/readcheck/` | HiFi depth per locus, variant-site read support, transcript evidence, start-codon checks, window definitions |
| `05_tables/derived/` | Per-family count differences, summed normalized depth vs assembled loci, the depth screen with BH-adjusted values, the KLRD1 junction check, the GC/length confounder check, counterpart identity per family and all pairwise KLRD1 protein identities |
| `05_tables/queries/` | The query sets used for annotation |
| `06_figures/` | Figures 1–5 as submitted, PNG and PDF at 300 dpi |
| `07_code/` | All analysis and figure code |

Every number in the manuscript is produced by the code in `07_code/` from these tables.
`04_numbers/digest.txt` is the place to start when checking one.

## Software

macOS 13, 4 cores, 8 GB RAM. Python 3.14.8; `miniprot` 0.18-r281; `minimap2` 2.31-r1302;
`diamond blastp` 2.1.16; `seqkit` 2.12.0 (`samtools faidx` is served by a `seqkit faidx` wrapper,
see `02_methods_full.txt` §1). Figures need `matplotlib`, the statistics need `scipy`.

## Primary data

All primary data are public and are not duplicated here.

Assemblies: GCA_057945255.1, GCA_057945245.1 (*Myotis sodalis*, mMyoSod1); GCA_054645015.1,
GCA_054645025.1 (*Rhinolophus clivosus*, mRhiCli1); GCA_059749095.1, GCA_059749105.1
(*Eidolon dupreanum*, mEidDup1).

PacBio HiFi: SRR38296513, SRR36951662, SRR39556278.

RNA-seq: SRR1657060, SRR13417526, SRR36602444, SRR27278001.

## Known limitations

The manuscript states these in full. The two that matter most when reusing these tables:

1. HiFi reads were mapped to a reference made of the selected windows of both haplotypes, not to
   the complete haplotypes. A read whose true locus is absent from that reference is not discarded
   but forced onto the most similar window, which inflates local depth in places. The depth columns
   should be read with that in mind, and the largest excess in summed depth reported by the screen
   is withdrawn in the manuscript on those grounds.
2. Splice-site labels in the curation tables are positional: a variant is labelled `splice_site`
   when it falls within 2 bp of an aligned exon boundary. The label does not test whether a donor
   or acceptor motif is destroyed.

## Licence

MIT, see [LICENSE](LICENSE). If you use the derived tables, a citation to the manuscript is
appreciated.

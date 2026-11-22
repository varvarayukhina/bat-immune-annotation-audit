# Auditing immune gene-family annotations across phased assemblies of three bat species

**Running title:** An audit of immune gene-family counts in phased bat genomes

**Author:** Varvara Iukhina¹

¹ Independent researcher

**Corresponding author:** Varvara Iukhina, juchinawarja@gmail.com

**Keywords:** haplotype-resolved assembly; gene family expansion; natural killer cell receptors; MHC class I; Chiroptera; copy number variation; annotation artefact

**Word count:** main text ~3,600; 5 figures; 1 table; 19 references

---

## Abstract

Immune gene-family counts inform comparisons across bat species. We examined how these counts
depend on the haplotype and annotation pass in paired assemblies of *Myotis sodalis*,
*Rhinolophus clivosus* and *Eidolon dupreanum*, with one individual per species. We audited
locus-level calls in natural killer receptor, leukocyte receptor and MHC class I regions using
renewed annotation searches, same-individual PacBio HiFi reads and RNA-seq and RefSeq proteins
from relatives. In the default pass, the largest difference was seven annotated CD1 loci in
*M. sodalis* (16 versus 23). The sums of absolute per-family differences were 21, 8 and 6
across the three individuals. In 91 of 117 flagged rows, renewed searches found a model of
the same family and status in the investigated partner region. Multiple landings or
out-of-window matches were flagged in 61 of these 91 rows. Fifteen rows carried candidate
sequence differences, and 14 met the study's read-support rules. A localized flanking
pile-up drove the largest excess in summed normalized depth. Of 50 start-codon flags,
28 had a canonical ATG at residue 1 in a relative's RefSeq alignment and 14 had upstream
RNA junction evidence. Together, the cases identify where paralogue assignment, mapping
and gene-model completion can change the interpretation of annotated locus counts.

---

## 1. Introduction

Bats are studied as reservoirs of viruses that rarely make them ill, and gene-family size is a
common explanatory variable in that literature. Expanded or contracted natural killer (NK) receptor
genes, MHC class I loci, interferon loci and related families are reported as species or lineage
characteristics [1–5], and functional work on single pathways, such as dampened STING-dependent
interferon activation [6] or attenuated NLRP3 inflammasome responses [7], is often read together
with these counts. Chromosome-level and haplotype-resolved bat assemblies are now routine [2,3,8],
and the counts drawn from them increasingly carry comparative weight.

Two steps in that chain can create differences that are not in the genome. The first is gene
calling: alignment-based annotation of a tandem array depends on score thresholds, on which query
proteins are used and on how secondary alignments are filtered, so the same array can yield
different locus counts in two runs or two assemblies. The MHC class I region of bats is organised
differently from the human region and contains lineage-specific duplications [9], and the
killer-cell receptor haplotypes of humans, of comparable complexity, were resolved only when long
reads made full haplotypes available [10]; both are warnings that short-query annotation of such
arrays is fragile. The second step is assembly: near-identical copies collapse into one and
heterozygous segments are sometimes represented twice [11,12], so the assembled sequence need not
carry the true number of copies. Pangenome work in humans has shown how much structural variation a
single reference hides [13].

Paired haplotypes allow an annotation workflow to be inspected within one individual.
Those haplotypes can differ in their gene content through structural variation. Errors shared
by the two assemblies or gene annotations can also leave their counts unchanged.
The magnitude of a count difference therefore has no calibrated relationship to the error
of the method. We compared annotated counts across three pairs of bat haplotypes and two
annotation passes, then used regional realignment and read evidence to audit selected
differences and identify ambiguous assignments.

## 2. Materials and methods

### 2.1 Assemblies and regions

Haplotype pairs were taken from the public assemblies of *Myotis sodalis* (mMyoSod1,
GCA_057945255.1 / GCA_057945245.1), *Rhinolophus clivosus* (mRhiCli1, GCA_054645015.1 /
GCA_054645025.1) and *Eidolon dupreanum* (mEidDup1, GCA_059749095.1 / GCA_059749105.1). Analysis
windows covered the natural killer complex (NKC), the leukocyte receptor complex (LRC) and MHC
class I regions with ±100 kb flanks: 39 windows (24.4 Mb) for *M. sodalis*, 24 (15.6 Mb) for
*R. clivosus* and 22 (15.2 Mb) for *E. dupreanum*, across both haplotypes.

### 2.2 Gene calling and haplotype comparison

Ninety-nine curated proteins of the target families (human, mouse and *Rousettus aegyptiacus*) were
aligned to each haplotype with miniprot 0.18-r281 [14]. Loci were assembled from hits and classified
as intact, partial, pseudogene-like or uncertain by the presence of frameshifts, premature stop
codons and non-canonical splice sites in the aligned coding sequence. The same pipeline was run
twice per haplotype. The default pass searched each whole haplotype in chunks of at most 250 Mb with
`-N 100 --outn 200`; the sensitive pass searched the selected windows and added `-k5 -M0 -n2 -l4`.
The two passes therefore differ in the reference searched as well as in the settings, and the cap on
secondary alignments applies per chunk in the first and once per window set in the second, so their
difference is not attributable to sensitivity alone. The separate cap `-N 1000 --outn 1000` was used
only in the curation step (§2.3). Haplotype regions were compared with minimap2 2.31-r1302
(`-x asm10`) [15], giving the locus-level and gap-level difference lists used below (205, 136 and
148 loci, respectively).

### 2.3 Curation of the reported differences

Every reported difference was re-examined. The locus and ±5 kb were extracted and searched again
with all queries at sensitive settings, so that a locus missing only because of a relative-score
cutoff could be recovered. Provisional family labels were checked with DIAMOND blastp 2.1.16 [16]
against reviewed Swiss-Prot entries of human (20,429) and mouse (17,285), release 2026_04 [17]. Loci
reported in one haplotype only were mapped with 10 kb flanks onto the whole partner haplotype, as
were all gaps ≥5 kb. A variant was labelled as affecting a splice site when it fell within 2 bp of
an aligned exon boundary; the rule is positional and does not test whether the canonical GT or AG
dinucleotide is destroyed, so such labels identify position and not consequence. The curation code
was validated on a synthetic two-haplotype genome carrying
four real query genes with introns and three planted events (a 1 bp deletion in exon 2 of KLRD1, a
deleted LILRB2 locus and a 6 kb insertion): all three were recovered and placed correctly.

### 2.4 Read-level check

For each assembly, a prefix of a PacBio run from the same BioSample as the assembly was streamed
through minimap2 (`-x map-hifi`) against the windows of both haplotypes: SRR38296513 (Revio, 18 GB)
for *M. sodalis*, SRR36951662 (Sequel II, 20 GB) for *R. clivosus* and SRR39556278 (Revio, 11.3 GB)
for *E. dupreanum*, yielding 2,572,740 / 705,323 / 386,399 alignments, of which 46,477 / 10,492 /
8,346 passed the filters (primary alignment, aligned block ≥2 kb, divergence ≤0.02). Median depth
per window was 3–19×, 3–7× and 4–9×. Locus depth is reported relative to the median depth of its
window, counting only reads with MAPQ ≥20 unless stated otherwise, and windows with a median below
4× were excluded from comparative depth analyses. Family-level sums and counts use the same
per-locus window threshold.

For exploratory screening, the log2 of positive MAPQ-filtered depth ratios was summarized
by its median and a robust scale (1.4826 × median absolute deviation) within each species.
Zero ratios were assigned a standardized score of −8. We calculated two-sided normal tail
areas and applied the Benjamini–Hochberg procedure. The resulting q values label an exploratory
screen with an uncalibrated false discovery rate. Two confounders were examined by Pearson
correlation of log2(ratio) against GC content and log10(locus length). Loci with a paralogue
of ≥95% protein identity in the same haplotype were compared with the rest by the Mann–Whitney
test.

For candidate sequence differences, the position was lifted to the comparison locus
(minimap2 `-x asm5 --cs`). We counted events within ±2 bp for substitutions and
±(homopolymer length + 2) bp for indels. Under the pilot rule, a candidate was read-supported
if ≥3 reads covered the alternative site and at most max(1, 20% of the covering reads)
carried a nearby event compatible by type with the other assembly. The rule uses the total
covering-read count and records the MAPQ ≥20 count separately. A high rate of reverse
events flagged candidates for possible assembly error. These checks provide provisional
support conditional on the window alignments and the locus pairing.

### 2.5 Transcript check

No RNA-seq exists for the assembled individuals, so RNA-seq of relatives was mapped to the same
windows with minimap2 (`-x splice:sr`): SRR1657060 (*Myotis myotis*, blood), SRR13417526
(*R. clivosus*, heart of an adult male, BioSample SAMN17088689) together with SRR36602444
(*Rhinolophus ferrumequinum*, spleen), and SRR27278001 (*Eidolon helvum*, blood). The *R. clivosus*
evidence therefore combines heart of that species with spleen of a congener, and no single tissue
covers it. Alignments were kept if primary, aligned block ≥50 bp, divergence ≤0.12
and within 5 kb of a locus. Because transcripts come from other animals, allelic statements were not
drawn from them; they were used for transcription, splice structure and the position of the first
exon.

### 2.6 Start-codon check

For every locus called intact, RefSeq proteins of the same family from an annotated relative
(*M. myotis*, *R. ferrumequinum*, *R. aegyptiacus*) were aligned to the locus with miniprot, and the
query start and first codon of the best hit were recorded.

### 2.7 Statistics and reproducibility

Wilson 95% intervals are descriptive summaries of row counts. Rows can share windows or
paralogues, which limits inference beyond the observed data. All analyses ran on a laptop (macOS 13,
4 cores, 8 GB RAM). Scripts, parameter files, result tables and the code that draws every figure are
deposited at https://github.com/varvarayukhina/bat-immune-annotation-audit; each number reported
below is produced by that code from the deposited tables.

## 3. Results

### 3.1 Within one individual, families differ by up to seven loci (Fig. 1)

Comparative studies state gene-family size as a count per species. Counted the same way, the two
haplotypes of one animal already disagree. In *M. sodalis*, 7 of 19 families differed between hap1
and hap2, by 21 loci in total and by as much as 7 in one family (CD1, 16 against 23; MHC class I, 34
against 40; LILR, 20 against 23). In *R. clivosus*, 4 of 22 families differed, by 8 loci in total
and at most 3 (MHC class I, 18 against 15). In *E. dupreanum*, 3 of 22 families differed, by 6 loci
in total and at most 3 (MHC class I, 19 against 16).

Raising the sensitivity of the search did not remove the disagreement; it redistributed it. At
sensitive settings the number of families that differ rose to 10, 8 and 6, while the summed
difference changed to 20, 12 and 10 loci. Five families per species changed their apparent
difference between the two passes, although no family reversed the direction of the difference.
The two passes differ in the reference searched as well as in the settings (§2.2).
The observed changes describe their combined effect on the reported counts.

### 3.2 Locus-level results after renewed search (Fig. 2)

The automated comparison reported 117 locus-level rows (54, 13 and 50). Renewed searches in
the investigated partner regions returned a model of the same family and status in 91 rows
(77.8%). Of those 91, 43 allowed multiple landings and 38 had a match outside the paired
window (61 with at least one flag). Same-family matches in these regions can represent a
different paralogue. Eleven rows retained different statuses without an identified candidate
sequence change in the examined alignment. Coverage-dependent classification remains one
possible explanation for those rows.

Fifteen rows carried candidate sequence differences. Thirteen involved a candidate small
sequence change, one a 32,236 bp deletion overlapping an exon of a LILR model in *M. sodalis*
hap1, and one a fragment from *M. sodalis* hap1 absent between aligned flanks in the hap2
assembly. The last fragment, hap1_L006, carried a pilot CLEC2D label. Its top back-search
match was CLEC2A, leaving the family assignment unresolved. The 15 rows are distributed
as 10, 2 and 3 across the individuals.

### 3.3 Read support for candidate sequence differences

Fourteen of the 15 rows met the pilot read-support rules. In *M. sodalis*, the small changes
with read support include substitutions and short indels in candidate LILR, MHC class I
and CD1 models. Five reads spanned the assembly junction of the 32 kb deletion. Nine reads
spanned the fragment at hap1_L006, supporting its presence on that assembly. The comparison
to the flanked hap2 assembly and the provisional family label were assessed separately.

In *E. dupreanum*, all three rows met the pilot support rule. Two comparisons involve the same
KLRD1 model, hap1_L027 (CM185962.1:24,025,356–24,029,142). One describes a premature stop
within that model. The other records a substitution at position 24,025,354, two bases before
the model boundary, marked by the position-only splice rule. Its effect on splicing is
unresolved. Counting distinct candidate loci gives read support at 13 of 14 by the pilot
criteria. In *R. clivosus*, the MHC class I candidate met the read-support rule. The NCR2
candidate remained undetermined with only two reads covering the alternative site.

### 3.4 Depth diagnostics depend on the reference and threshold (Fig. 3)

Read depth per locus, expressed relative to the median depth of its window, placed 406 of the 472
loci with usable coverage (86.0%, 95% CI 82.6–88.9%) in the fixed 0.55–1.7× corridor using all
primary reads. Of the 17 additional loci, 13 lay outside the mapped windows and four fell in
windows with median depth below 4×. At the same all-read thresholds, 33 loci were at or above
1.7×, 31 were between 0.2× and 0.55×, and two were below 0.2×. High and low ratios identify regions
for which mapping ambiguity or assembly copy-number errors merit inspection [12]. Across the
476 loci for which depth could be computed, 23 had zero coverage after the MAPQ ≥20 filter
(2, 16 and 5 per individual). Of these, 20 passed the window coverage criterion and enter Fig. 3a
(2, 13 and 5).

The depth screen flagged 45 of 472 loci at nominal q <0.05 (13, 19 and 13). The robust
centre of log2(ratio) was +0.06, +0.05 and −0.17 with respective scales of 0.48, 0.58 and 0.55.
The empirical fit and the assigned score for zero ratios make the screen exploratory. The two criteria can only be compared on one depth column, and as
implemented they do not share one: the corridor is applied to depth from all primary reads, the
empirical null to depth from reads at MAPQ ≥20. Applied both to the MAPQ-filtered column, they agree
for 415 of 472 loci (87.9%, 95% CI 84.7–90.6%); no significant locus then lies inside the corridor,
and all 57 disagreements are loci outside the corridor that the test does not call. Comparing the
corridor call as published against the test, which mixes the two columns, gives 395 of 472 (83.7%,
95% CI 80.1–86.7%), with 28 significant loci inside the corridor and 49 outside it (28 above 1.7×,
21 below 0.55×) that are not significant. These comparisons describe how the classifications change with the depth column and
screening rule.

We also compared summed normalized depth with the annotated locus count per family and
haplotype (Fig. 3b) after excluding individual loci from windows with median depth below 4×.
All three CD1 loci of *R. clivosus* hap2 fail this threshold, so that group is absent from the
depth comparison. The absolute deviation was at most two units in 108 of 115 groups. Four
groups had excess depth and three had deficits. The count and depth sum were 18 and 61.5 for
MHC class I in *R. clivosus* hap1, 20 and 50.2 for LILR in *M. sodalis* hap1, and 40 and 51.1
for MHC class I and 23 and 29.9 for LILR in *M. sodalis* hap2. The three deficits involved CD1
(3 against 0.1 in *R. clivosus* hap1 and 10 against 5.2 and 4.4 in *E. dupreanum* hap1 and hap2).
The largest excess is shown as a diagnostic of the original method and excluded from
copy-number interpretation.

The largest excess came from the hap1_L008 depth value in *R. clivosus*. In the saved PAF, the
3 kb immediately left of the candidate locus in genomic coordinates had a mean aligned-base
depth of 835× at MAPQ ≥20. The six CDS segments of the best HLA-B query alignment averaged
3.37×. This localized peak drove the high normalized depth of the padded interval. Its origin
remains undetermined. Summed depth in this window therefore gives no supported estimate of the
number of MHC-I copies. The other three excesses and three CD1 deficits also require
locus-specific validation.

Two confounders were tested. GC content was not associated with the depth ratio in any individual
(r = −0.10, +0.12 and +0.03; p = 0.18, 0.19 and 0.70). Locus length was associated in *M. sodalis*
(r = +0.35, p = 3.6 × 10⁻⁷) but not in the other two (p = 0.48 and 0.38), so in that individual part
of the spread reflects locus size and the depth statements there are correspondingly weaker. Loci
with a paralogue of ≥95% protein identity in the same haplotype did not show higher depth ratios
than the rest (Mann-Whitney p = 0.72 and 0.96) and in *E. dupreanum* tended to show lower ones
(median 0.46 against 0.89, p = 0.06): with a MAPQ ≥20 filter, near-identical copies lose reads
rather than accumulate them. The high-depth outliers are therefore not simply the loci with an
obvious close paralogue, and copy number in families whose members are nearly identical (CD1 in
*E. dupreanum* and *R. clivosus*) cannot be resolved by this approach at all.

### 3.5 RefSeq and RNA checks on start-codon flags (Fig. 4)

Across the three individuals, 122 loci were called intact, and 50 of them carried a start-codon flag
because the human or mouse query did not align from residue 1 or the first codon was not ATG. Using
RefSeq proteins of an annotated relative instead placed a canonical ATG at residue 1 in 28 of the 50
(56.0%, 95% CI 42.3–68.8%). Two further loci, both KLRD1 models of *E. dupreanum* hap2, aligned
from residue 6 and form a separate category in Fig. 4a.

All 50 carried transcript evidence, but of two kinds that differ in strength. In 14 (two in
*R. clivosus*, twelve in *E. dupreanum*; 28.0%, 95% CI 17.5–41.7%) at least two reads were spliced
from an upstream exon into the first part of the model, which supports an upstream splice
connection on individual molecules. In the remaining 36 (72.0%, 95% CI 58.3–82.5%) the evidence is two
coverage counts — at least three reads over the model start and at least three in the 100 bp above
it — which are counted independently and need not come from the same molecule; they establish that
the region above the start is transcribed, not that it is contiguous with the model.

### 3.6 RNA alignments across the KLRD1 models (Fig. 5)

All 13 pilot-intact KLRD1 models in *E. dupreanum* carried a start flag. Ten began with TCG or
TCT in the pilot annotation and three aligned from a later query residue. All 19 KLRD1 models
overlapped RNA alignments from *E. helvum* (67–209 alignments per model). Eighteen had at least
two upstream junction alignments (2–32). When we applied MAPQ ≥20 to the same PAF inputs and
reran the original junction-counting code, 10 models met that criterion (Fig. 5b). Among read
identifier and length pairs whose alignments overlapped at least one KLRD1 model, 395 of 456
hap1 pairs also occurred in the separate hap2 mapping, which contained 454 pairs. The shared
alignments support family-level transcription in the donor species. Individual models may
receive signal from the same read.

Transcription is otherwise near-universal in these windows (Fig. 5a): of the families represented,
only KLRG2 (all three individuals) and NCR2 (*M. sodalis*, *E. dupreanum*) had no RNA evidence in
the tissues used, and FCAR in *R. clivosus* reached half of its loci. Because the reads come from
other animals and distribute over near-identical copies, this should be read as transcription of the
family in the tissues sampled — blood, spleen and heart — rather than expression of an individual
copy; single-cell work in
*R. aegyptiacus* shows how unevenly such families are expressed across immune cell types [18].

## 4. Discussion

In these three individuals, default annotated counts differed between haplotypes by up to seven
loci per family. Both annotation passes and both haplotypes contributed to variation in the
reported counts. Renewed searches frequently found a model of the original family and status,
with many matches requiring further positional or paralogue validation. The study therefore
illustrates the sensitivity of these family counts to annotation choices and the availability
of evidence for individual loci.

Phased assemblies make it possible to report the within-individual differences alongside
between-species counts. Same-individual reads can test candidate absences, sequence changes and
junctions. Locus-wide normalized depth highlights where those checks are needed. Its classifications
change with the depth column and threshold, and a localized peak can dominate a family-level sum.
For candidate start-codon truncations, a relative's RefSeq protein and RNA splice junctions
give different kinds of evidence about the gene model.

The RNA data from *E. helvum* support KLRD1 family transcription, with substantial sharing of
alignments across the two *E. dupreanum* reference haplotypes. One candidate disruptive stop
was observed in a KLRD1 model in the assembled individual. The nearby substitution received a position-only splice label and requires separate
functional validation. Excess depth in LILR and MHC-I and low depth in
CD1 identify cases where resolving paralogy and mappability matters for copy counting.
Comparative claims about NK receptor and MHC class I repertoires in bats [1,4,5,9] would benefit
from such audits. Functional assays, such as the study of MR1-restricted T cells in
*Pteropus alecto* [19], address the function of the resulting gene models.

### 4.1 Limitations

Reads were mapped to extracted windows rather than whole genomes, and this is the largest limitation
of the depth analysis. Reads from elsewhere in the genome can align to similar sequence in the
available windows. Among primary
alignments with an aligned block ≥5 kb against the animal's own windows, those at divergence ≥0.10
cover 37.5% of the read in *E. dupreanum* at a mean MAPQ of 3.7 and make up 72.5% of the reads,
whereas those below 0.01 divergence cover 99.2% of the read at a mean MAPQ of 41.2; in
*R. clivosus* the corresponding population sits at 0.05–0.10 divergence, covers 28.5% of the read at
a mean MAPQ of 9.8 and makes up 54.7% of the reads. The divergence ≤0.02 and MAPQ ≥20 filters remove
that population. Close matches to sequences outside the windows could account for the pile-up in
§3.4. The origin of that signal would be tested by competitive remapping of the same reads against
both complete haplotypes. The MAPQ ≥20 filter also removes reads from near-identical copies and
leaves copy number in such families unresolved. HiFi depth was 4–9× per
window for *E. dupreanum* and 3–7× for *R. clivosus*, so a locus covered by three to five reads is
at the limit of what can be decided, and one difference (NCR2 in *R. clivosus*) stayed undecided for
that reason. In *M. sodalis* the depth ratio correlated with locus length, so part of the spread
there is not copy number. RNA came from other individuals and, for two species, other species of the
same genus, so it supports transcription and splice structure but not alleles; MHC class I alleles in
particular cannot be checked this way, and within a near-identical array the same transcript read is
counted at several copies, so per-copy transcription is not established. Non-canonical splice sites
were read off the dinucleotides at the aligned intron boundaries, and in *R. clivosus* 17 of the 19
such calls are the single pattern GT–AC, 7 of them in MHC class I loci of one haplotype; a pattern
that uniform is more likely to be systematic in the extraction than biological, and the status calls
resting on it should be treated as provisional. Prefixes of the read runs were used rather than
whole runs. The synthetic validation covered three planted events. Sensitivity across
variant classes remains unmeasured. With one individual per species, variation in gene
content within species and error rates for the workflow remain outside the scope of this study.

## 5. Conclusions

Annotated immune-family counts differed by up to seven loci between the paired haplotypes of
one bat. Gene-model searches, genomic position, read coverage and upstream RNA junctions
provided different kinds of evidence when inspecting those counts. Whole-locus depth flagged
a strong signal whose peak lay next to the coding sequence. Renewed searches also identified
cases where locus correspondence remained ambiguous. Reporting both haplotypes and the
evidence for individual gene models can make these uncertainties visible in comparative
immunogenomics.

## Declarations

**Funding.** This work received no specific grant from any funding agency in the public, commercial
or not-for-profit sectors. All computations were run on the author's personal laptop.

**Conflicts of interest.** The author declares no competing interests.

**Ethics approval.** Not applicable: the study analyses publicly archived sequence data only; no new
animal material was collected.

**Data availability.** All primary data are public. Assemblies: GCA_057945255.1, GCA_057945245.1,
GCA_054645015.1, GCA_054645025.1, GCA_059749095.1, GCA_059749105.1. Sequencing runs: SRR38296513,
SRR36951662, SRR39556278 (PacBio HiFi); SRR1657060, SRR13417526, SRR36602444, SRR27278001
(RNA-seq).

**Code and derived-data availability.** Analysis scripts, parameter files, all result tables (locus
curation, difference curation, read depth, variant-site support, transcript evidence, start-codon
checks, per-family count differences, normalized depth sums, exploratory depth screen) and the figure code are
deposited at https://github.com/varvarayukhina/bat-immune-annotation-audit under the MIT licence.

**Author contributions.** Varvara Iukhina designed the study, performed the analyses, verified all
results against the primary tables, prepared the figures, wrote the manuscript and approves the
submitted version.

**Use of assistive technology.** An AI coding assistant was used to edit the analysis code and to
improve the visual design of the figures. The author designed the study, performed the analyses,
verified every number against the primary result tables and is responsible for all content.

## References

1. Pavlovich SS, Lovett SP, Koroleva G, et al. The Egyptian Rousette genome reveals unexpected
   features of bat antiviral immunity. *Cell*. 2018;173(5):1098–1110.e18.
   doi:10.1016/j.cell.2018.03.070
2. Jebb D, Huang Z, Pippel M, et al. Six reference-quality genomes reveal evolution of bat
   adaptations. *Nature*. 2020;583(7817):578–584. doi:10.1038/s41586-020-2486-3
3. Morales AE, Dong Y, Brown T, et al. Bat genomes illuminate adaptations to viral tolerance and
   disease resistance. *Nature*. 2025;638(8050):449–458. doi:10.1038/s41586-024-08471-0
4. Moreno Santillán DD, Lama TM, Gutierrez Guerrero YT, et al. Large-scale genome sampling reveals
   unique immunity and metabolic adaptations in bats. *Molecular Ecology*. 2021;30(23):6449–6467.
   doi:10.1111/mec.16027
5. Scheben A, Mendivil Ramos O, Kramer M, et al. Long-read sequencing reveals rapid evolution of
   immunity- and cancer-related genes in bats. *Genome Biology and Evolution*. 2023;15(9):evad148.
   doi:10.1093/gbe/evad148
6. Xie J, Li Y, Shen X, et al. Dampened STING-dependent interferon activation in bats. *Cell Host &
   Microbe*. 2018;23(3):297–301.e4. doi:10.1016/j.chom.2018.01.006
7. Ahn M, Anderson DE, Zhang Q, et al. Dampened NLRP3-mediated inflammation in bats and implications
   for a special viral reservoir host. *Nature Microbiology*. 2019;4(5):789–799.
   doi:10.1038/s41564-019-0371-3
8. Rhie A, McCarthy SA, Fedrigo O, et al. Towards complete and error-free genome assemblies of all
   vertebrate species. *Nature*. 2021;592(7856):737–746. doi:10.1038/s41586-021-03451-0
9. Ng JH, Tachedjian M, Wang LF, Baker ML. Evolution and comparative analysis of the bat MHC-I
   region. *Scientific Reports*. 2016;6:21256. doi:10.1038/srep21256
10. Roe D, Vierra-Green C, Pyo CW, et al. Revealing complete complex KIR haplotypes phased by
    long-read sequencing technology. *Genes and Immunity*. 2017;18(3):127–134.
    doi:10.1038/gene.2017.10
11. Cheng H, Jarvis ED, Fedrigo O, et al. Haplotype-resolved assembly of diploid genomes without
    parental data. *Nature Biotechnology*. 2022;40(9):1332–1335. doi:10.1038/s41587-022-01261-x
12. Guan D, McCarthy SA, Wood J, Howe K, Wang Y, Durbin R. Identifying and removing haplotypic
    duplication in primary genome assemblies. *Bioinformatics*. 2020;36(9):2896–2898.
    doi:10.1093/bioinformatics/btaa025
13. Liao WW, Asri M, Ebler J, et al. A draft human pangenome reference. *Nature*.
    2023;617(7960):312–324. doi:10.1038/s41586-023-05896-x
14. Li H. Protein-to-genome alignment with miniprot. *Bioinformatics*. 2023;39(1):btad014.
    doi:10.1093/bioinformatics/btad014
15. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*.
    2018;34(18):3094–3100. doi:10.1093/bioinformatics/bty191
16. Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using
    DIAMOND. *Nature Methods*. 2021;18(4):366–368. doi:10.1038/s41592-021-01101-x
17. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids
    Research*. 2025;53(D1):D609–D617. doi:10.1093/nar/gkae1010
18. Friedrichs V, Toussaint C, Schäfer A, et al. Landscape and age dynamics of immune cells in the
    Egyptian rousette bat. *Cell Reports*. 2022;40(10):111305. doi:10.1016/j.celrep.2022.111305
19. Leeansyah E, Hey YY, Sia WR, et al. MR1-restricted T cells with MAIT-like characteristics are
    functionally conserved in the pteropid bat *Pteropus alecto*. *iScience*. 2020;23(12):101876.
    doi:10.1016/j.isci.2020.101876

> **Note on references 14–17.** Software and database citations are given in the form the publishers
> list; verify volume, page and year against the publisher record before submission, and apply the
> reference style of the target journal. References 1–13, 18 and 19 were verified against Europe PMC
> on 12 November 2026.

## Figure legends

**Figure 1.** *The two haplotypes of one individual, counted family by family.* Each row is a gene
family whose locus count differs between hap1 (light) and hap2 (dark) of the same animal; the number
beside each pair is the difference. Families that agree are not drawn and are counted in the panel
title. Counts are from the default annotation run; at sensitive settings the number of families that
differ rises to 10, 8 and 6, and the summed difference becomes 20, 12 and 10 loci. Panels share the
x-axis.

**Figure 2.** *Results of renewed searches and sequence inspection.* The automated comparison
reported 54, 13 and 50 locus-level rows. In 35, 10 and 46 rows the partner-region search recovered
the same family and status. In 9, 1 and 1 rows the model statuses differed without an identified
candidate sequence change. In 8, 2 and 3 rows the differing calls had a candidate small
sequence change. The remaining two rows describe a large deletion overlapping a LILR model and a
fragment absent from the partner assembly with an unresolved CLEC-family label.

**Figure 3.** *Normalized depth and annotated locus counts.* (**a**) Each point represents one
of 472 loci from windows with median depth ≥4×. The ordinate is MAPQ ≥20 locus depth divided by
the all-primary-read window median. The shaded band shows the fixed 0.55–1.7× corridor on this
column. Orange marks the exploratory screen at nominal BH q < 0.05. The zero row contains 2, 13 and 5 loci, respectively.
The asterisk marks the *R. clivosus* hap1_L008 pile-up described in §3.4. (**b**) Each point
compares the number of annotated loci with usable depth and the sum of their normalized
depths for a family-haplotype group. The diagonal shows equality. Seven groups deviate by more than two units, four above
and three below, including the withdrawn MHC-I estimate marked by an open diamond.

**Figure 4.** *Follow-up of 50 start-codon flags among 122 pilot-intact models.* (**a**) Alignment
to a relative's RefSeq protein gave a canonical ATG at query residue 1 in 28 models, alignment
from query residue 6 in two, and another alignment outcome in 20. (**b**) RNA alignments from
relatives provided upstream junction evidence for 14 models and separate coverage over and
above the model start for 36. Each panel classifies the same flagged models by its own
criterion.

**Figure 5.** *RNA alignment evidence from relatives.* (**a**) Fraction of annotated models per
family with at least three RNA alignments. Reference species identify the assemblies. RNA
sources are in Table 1. (**b**) Upstream junction alignments per KLRD1 model in *E. dupreanum*
reference windows, using RNA from *E. helvum*. Light bars show the original counts and orange
bars show the subset retained at MAPQ ≥20. The dashed line marks two alignments. Eighteen
models meet this threshold before MAPQ filtering and 10 after it. Alignments can be shared
between models and haplotypes.

## Table 1

*Material and the audit, per individual.*

| | *Myotis sodalis* | *Rhinolophus clivosus* | *Eidolon dupreanum* |
|---|---|---|---|
| Assembly (hap1 / hap2) | GCA_057945255.1 / GCA_057945245.1 | GCA_054645015.1 / GCA_054645025.1 | GCA_059749095.1 / GCA_059749105.1 |
| Windows (both haplotypes) | 39 (24.4 Mb) | 24 (15.6 Mb) | 22 (15.2 Mb) |
| Loci called | 205 | 136 | 148 |
| Families differing between haplotypes (default / sensitive) | 7 of 19 / 10 of 20 | 4 of 22 / 8 of 22 | 3 of 22 / 6 of 21 |
| Summed per-family difference, loci (default / sensitive) | 21 / 20 | 8 / 12 | 6 / 10 |
| Largest difference in one family | 7 (CD1, 16 vs 23) | 3 (MHC-I, 18 vs 15) | 3 (MHC-I, 19 vs 16) |
| Automated locus differences | 54 | 13 | 50 |
| Same family and status on renewed search | 35 | 10 | 46 |
| Status differs, no candidate sequence change found | 9 | 1 | 1 |
| Candidate sequence differences (rows / distinct loci) | 10 / 10 | 2 / 2 | 3 / 2 |
| Meets pilot read-support criteria (rows / distinct loci) | 10 / 10 | 1 / 1 | 3 / 2 |
| HiFi run (prefix used) | SRR38296513 (18 GB) | SRR36951662 (20 GB) | SRR39556278 (11.3 GB) |
| Alignments kept / total | 46,477 / 2,572,740 | 10,492 / 705,323 | 8,346 / 386,399 |
| Median depth per window | 3–19× | 3–7× | 4–9× |
| Loci with all-primary depth ratios ≥1.7× / 0.2–0.55× | 17 / 13 | 10 / 9 | 6 / 9 |
| Loci with zero MAPQ ≥20 coverage / loci with measured depth | 2 / 200 | 16 / 132 | 5 / 144 |
| Exploratory depth flags at nominal BH q <0.05 | 13 of 199 | 19 of 129 | 13 of 144 |
| RNA-seq runs (source) | SRR1657060 (*M. myotis*, blood) | SRR13417526 (*R. clivosus*, heart); SRR36602444 (*R. ferrumequinum*, spleen) | SRR27278001 (*E. helvum*, blood) |
| Intact loci / flagged for start codon | 44 / 25 | 28 / 8 | 50 / 17 |
| Flags with a canonical ATG at residue 1 from a relative's protein | 20 | 6 | 2 |
| Flags with an upstream-exon junction into the model start | 0 | 2 | 12 |
| Flags with coverage over and above the model start only | 25 | 6 | 5 |

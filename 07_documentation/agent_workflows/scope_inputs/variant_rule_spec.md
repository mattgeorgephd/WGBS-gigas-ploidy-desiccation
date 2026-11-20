# C-to-T variant exclusion "as in Venkataraman et al. 2022": what was done, what BS-Snper does, and a specification for R

Scoping note, 2026-11-20. Read only with respect to the repository. Everything downloaded is in
`downloads/` of this scratch directory; my scripts are in `code/`.

Labels used below: **[verified]** = read in a primary source or computed by me, with the source or
command given; **[assumed]** = not checked.

---

## 1. What Venkataraman et al. 2022 did

### 1.1 The article (Europe PMC full text of PMC9351233, `downloads/article/pmc.xml`)

Methods, "SNP identification" [verified, quoted]:

> "SNP variants across all samples and in individuals were identified with BS-Snper [67]. This program
> takes into consideration nucleotide identity on the opposite strand: a thymine paired with an adenine
> indicates a C- > T SNP, while a thymine paired with a guanine is a bisulfite-converted unmethylated
> cytosine. Sorted deduplicated BAM files generated with Bismark v.0.22.3 and SAMtools v.1.10 [65, 66]
> were merged (samtools merge) to identify SNPs. Variants were also identified in individual sorted
> deduplicated BAM files. Default BS-Snper settings were used, except for a minimum coverage of 5x to be
> consistent with methylation analyses. The C- > T SNPs were used to exclude loci of interest."

Methods, "Global methylation characterization": "using loci with at least 5x coverage that did not
overlap with C- > T SNPs". Methods, DML: "CpG dinucleotides that did not overlap with C- > T SNPs with at
least a 50% methylation difference ... and a q-value < 0.01". Results: "13,234,183 unique SNPs were found
in individual or merged BAM files, including 300,278 C- > T SNPs"; "11,238,223 CpG loci ... After C- > T
SNP removal, 10,939,918 CpG loci"; "1599 CpG loci were identified as potential DML ... 315 (19.7%)
overlapped with C- > T SNPs, and thus were removed ... A total of 1284 DML".

### 1.2 The code (OSF project ygctb, "Oyster Gonad Methylation"; files via api.osf.io)

Files read: `code/07-BS-SNPer.ipynb` (OSF id f56eh), `code/09-General-Methylation-Landscape.ipynb`
(am4b5) and `.Rmd` (u725h), `code/10-Genomic-Location-of-DML.ipynb` (qdez6) and `.Rmd` (3xwjr),
`code/06-methylKit.Rmd` (u46xj), `code/05-bismark-roslin.sh` (8gskp), and
`output/06-methylKit/DML/DML-pH-50-Cov5-All.csv`. Text renderings with cell numbers:
`downloads/osf/07.txt`, `09-General-Methylation-Landscape.txt`, `10-Genomic-Location-of-DML.txt`.

Step by step [verified, notebook 07 unless stated]:

1. **Merge.** Cell 18: `samtools merge .../merged-sorted-deduplicated.bam` of the 8
   `zr3616_*_bismark_bt2_pe.deduplicated.sorted.bam`.
2. **BS-Snper on the merged file.** Cell 24 (comment line: `#Defaults: --minhetfreq 0.1 --minhomfreq
   0.85 --minquali 15 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20`):
   `perl BS-Snper.pl --fa .../cgigas_uk_roslin_v1_genomic-mito.fa --input .../merged-sorted-deduplicated.bam
   --output .../SNP-candidates.txt --methcg ... --methchg ... --methchh ... --mincover 5 > .../SNP-results.vcf`.
   The VCF header (cell 25 output) records the parameters actually used: `--minhetfreq 0.1 --minhomfreq
   0.85 --minquali 15 --mincover 5 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20`,
   `##bssnperVersion=1.1`, and the fixed line `##Bisulfite=directional>`.
3. **BS-Snper on each of the 8 BAM files.** Cell 33, same command per file, `--mincover 5`.
4. **Keep C/T at CpG motifs.** Cells 28 and 37: `intersectBed -u -a <vcf> -b
   cgigas_uk_roslin_v1_fuzznuc_CGmotif.gff | grep "C<TAB>T" > ..._CT-SNPs.tab`. There is **no filter on
   FILTER, QUAL, GT or DP**. The grep matches REF = C with ALT = T or ALT = "T,G" (genotype G/T at a
   reference C). It does not match REF = G with ALT = A (the same change on the other strand), nor ALT =
   "A,T".
5. **Union.** Cell 52: `cat *CT-SNPs.tab >> all-CT-SNPs.tab` (merged + 8 individual, 685,762 rows);
   cell 53: `awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort | uniq > unique-CT-SNPs.tab` (300,278 rows).
6. **BED.** Notebook 10, cell 18: `awk '{print $1"\t"$2"\t"$2}' unique-CT-SNPs.tab > unique-CT-SNPs.bed`
   (start = end = the 1-based VCF POS; a zero-length BED record).
7. **Exclusion from the CpG sets of the global characterization.** Notebook 09, cell 32: `subtractBed -a
   <each 5x bedgraph and the union-average bedgraph> -b unique-CT-SNPs.bed > *.NO-SNPs` (bedtools v2.30.0,
   cell 15). All later counts of notebook 09 and the "all 5x CpGs" background of the genomic-location test
   are built from the `NO-SNPs` files (09 Rmd line 49; notebook 09 cells near output lines 991 to 1137).
8. **Exclusion from the DML list, after the test.** Notebook 10, cell 21: `intersectBed -u -a
   DML-pH-50-Cov5-All.csv.bed -b unique-CT-SNPs.bed` (315 rows); cell 22: `subtractBed` (1,284 rows).
   The methylKit input, the united table, the model fits and the SLIM q-values were **not** changed: in
   `06-methylKit.Rmd` the variants appear only in heatmaps of DML (lines 593 to 698), and the DML were
   called with `calculateDiffMeth(..., overdispersion = "MN", test = "Chisq")` (line 211) and
   `getMethylDiff(difference = 50, qvalue = 0.01)` (line 222) on all united CpGs.

### 1.3 Checks of their outputs (my computation)

I downloaded the nine C/T tables and `unique-CT-SNPs.tab` from
`https://gannet.fish.washington.edu/spartina/project-gigas-oa-meth/output/BS-Snper/` (line counts equal
the notebook outputs: 122,306 merged; 69,845 to 73,664 per animal; sum 685,762), and ran
`code/analyse_ct.py` (output `ct_analysis.txt`):

| Quantity | Value |
|---|---|
| Rows of `unique-CT-SNPs.tab` (the article's "300,278 C->T SNPs") | 300,278 |
| Distinct positions | **300,032** (246 positions appear twice, once with ALT T and once with ALT "T,G") |
| Positions that are not the C of a CpG on the plus strand (checked against `cgigas_uk_roslin_v1_genomic-mito.fa` of the pipeline cache) | 0 |
| Positions with at least one PASS record | 251,627 |
| Positions supported only by FILTER = Low records | 48,405 (16.1%) |
| Positions that never had 5 or more opposite-strand (Crick) reads in any of the 9 calls | 48,466 |
| Positions that never had more than 1 opposite-strand T read | 21,119 |
| Positions from the merged call / only from individual calls / only from the merged call | 122,306 / 177,726 / 2,583 |
| DML (of 1,599) whose C is a listed position | **315** (reproduces the article) |
| ... counting only positions with at least one PASS record | 266 |
| ... counting only the merged call | 192 |
| ... resting only on Low records | 49 |

Also: the article's "11,238,223 CpG loci" is the line count of `union_5x.bedgraph` including its header
line (notebook 09 cell 20); the file holds 11,238,222 CpGs, and 11,238,222 - 10,939,918 = 298,304 rows
were removed.

### 1.4 The rule Venkataraman et al. actually applied (derived from 1.2 and 1.3)

A CpG was excluded if, in **any** of the nine BS-Snper runs (the merged BAM of all 8 animals, or one
animal's BAM), BS-Snper wrote a VCF record at the **C of that CpG on the plus strand** with REF = C and
ALT = T or "T,G", **whatever its FILTER (PASS or Low), QUAL or depth**. The "minimum coverage of 5x" only
changed the PASS/Low label (section 2.6) and therefore had no effect on which loci were excluded. A
variant at the G of the CpG (a C-to-T change on the minus strand, written by BS-Snper as G>A) did **not**
exclude the CpG. Exclusion was applied to the CpG sets of the global characterization and to the DML list
after the test, not to the data the test was run on.

Coordinates [verified]: the coverage rows (theirs and this project's) are `chr, start, end` with start =
0-based position of the C and end = start + 2 (Bismark `coverage2cytosine --merge_CpG --zero_based`; their
`05-bismark-roslin.sh`, this project's `07_documentation/upstream_processing/README.md` and
`02_R_functions/coverage_io.R` lines 4-5; checked against the genome with `code/check_coords.py`: rows
`649 651`, `1547 1549`, `1571 1573` of NC_047559.1 and their DML rows `7867 7869` etc. cover "CG"). The SNP
BED record `POS POS` is widened by bedtools 2.30.0 to `[POS-1, POS+1)` (bedtools2 v2.30.0
`src/utils/FileRecordTools/Records/Record.cpp` lines 241-248, `adjustZeroLength()`), so it overlaps a CpG
row `[C0, C0+2)` when POS is C0, C0+1 or C0+2 (1-based C-1, C, G). Since every listed POS is a C followed
by a G, only POS = C0 + 1 (the C) can match. Computed both ways, the DML overlap is 315 (section 1.3).

**In R, for a coverage row (chr, start0): excluded if (chr, start0 + 1) is in the variant set.**

---

## 2. What BS-Snper does

Source: `https://raw.githubusercontent.com/hellbelly/BS-Snper/master/` (BS-Snper.pl, main.c,
sam_funcs.c, README.txt; md5 in facts). The Perl header says "Version: 1.1, Date: 2018-07-27"; the study's
VCF headers say `bssnperVersion=1.1` [verified]. That the master branch equals the copy downloaded in April
2021 is [assumed]: the GitHub API and commit feed were not reachable from this session. Paper: Gao et al.
2015, Bioinformatics 31:4006 (PMC4673977, read in full; its Supplementary Text and Table S1 were not read;
the code is the authority used here).

### 2.1 Parameters and defaults [verified, BS-Snper.pl lines 55-62]

| Option | Default | Where used | Effect |
|---|---|---|---|
| `--minhetfreq` | 0.1 | C: candidate filter (`vSnpRate`, sam_funcs.c line 129, strict `>`); Perl: heterozygous PASS (`>=`, on a value rounded to 3 decimals) | |
| `--minhomfreq` | 0.85 | Perl: homozygous PASS | |
| `--minquali` | 15 | C: bases with quality < 15 become N (sam_funcs.c line 1003); Perl: mean quality of the variant base >= 15 for PASS (always true after the masking) | |
| `--mincover` | **10** (README usage line and code); the 2022 study used 5 | C: only the methylation output files (`nLayerMin`, printMeth); Perl: PASS/Low label | does not remove records |
| `--maxcover` | 1000 | passed to the C program as `nLayerMax`, which is **never used** (only in the three function signatures, sam_funcs.c lines 13, 16, 801) | none |
| `--minread2` | 2 | Perl: variant reads >= 2 for PASS | |
| `--errorate` | 0.02 | C: per-read mismatch filter (`vSnpPerBase`, sam_funcs.c lines 1014-1016, 1063) | |
| `--mapvalue` | 20 | C: reads with MAPQ < 20 dropped (sam_funcs.c lines 183-185) | |
| `--regions-file`, `--variants-only` | unset | Perl only; parsed but their flag function is not used in output | none seen |

Call: `rrbsSnp $fasta $bam $output $methcg $methchg $methchh $minquali $mincover $maxcover $minhetfreq
$errorate $mapvalue` (BS-Snper.pl line 109; main.c maps them to vQualMin, nLayerMin, nLayerMax, vSnpRate,
vSnpPerBase, mapqThr).

### 2.2 Read handling (C program, `snpProcess`, sam_funcs.c 801-1332) [verified]

- Reads: skipped if unmapped (FLAG 0x4) or MAPQ < mapvalue. No other FLAG is checked (secondary, QC-fail
  and duplicate flags are not looked at). Overlapping mates are **both counted**.
- Strand and mate: `strand = '+'` if FLAG 0x10 unset, else `'-'`; `r12 = 2` if FLAG 0x80 set, else 1
  (lines 165-176).
- **Watson class** = (`+`, r12 1) or (`-`, r12 2), i.e. FLAG 99 and 147 in a Bismark paired-end BAM.
  **Crick class** = (`-`, r12 1) or (`+`, r12 2), i.e. FLAG 83 and 163. This is the directional-library
  assumption (read 1 forward = top strand).
- CIGAR: **only the last CIGAR operation is kept** (line 188: `sprintf(record->cigar, ...)` inside the
  loop overwrites the buffer each time; confirmed with `code/cigar_loop_test.c`, which prints `7M` for
  `66M7I7M`). The read is then laid out from POS with that one operation: for `66M7I7M` only the first 7
  bases are used; for `15M1D34M` the first 34 bases are placed at POS..POS+33, so bases 16-34 are shifted
  by one position. Reads without indels are unaffected. M copies bases, D writes N, I skips read bases,
  S writes N at reference positions (bowtie2 end-to-end alignments of Bismark have no S).
- Base quality: bases with Phred < minquali (15) become N; bases at a reference N become N (lines
  999-1009).
- Per-read mismatch filter: `iread = (int)(len * errorate + 0.5)`, at least 1 (`len` = laid-out length);
  `cnt` = number of non-N bases that differ from the reference, **not counting** T at a reference C in a
  Watson read and A at a reference G in a Crick read; the read is dropped entirely if `cnt > iread`
  (lines 1012-1063). For 100 to 149 laid-out bases iread is 2 or 3.
- Counting (lines 1063-1290), at each laid-out position, in plus-strand letters as stored in the BAM:
  - Watson read: T at a reference C is counted as **w_C** (bisulfite conversion folded into C); any
    other base X is counted as w_X.
  - Crick read: A at a reference G is counted as **c_G**; any other base X is counted as c_X.
  - Sums of base qualities per (class, base) are kept; MAPQ sums per class too.
  - Consequence at a reference C: w_T is always 0; the Watson reads cannot separate C from T; the Crick
    reads show C (genotype C) or T (genotype T, because the bottom strand carries A, which bisulfite does
    not change).

### 2.3 Candidate filter (`printSnp`, sam_funcs.c 58-139) [verified]

With m = sum of the 8 counts, for a reference C (lines 99-111):
C>A = (w_A + c_A)/m; **C>T = c_T/(c_C + c_T)**; C>G = (w_G + c_G + c_A)/m. For a reference G:
**G>A = w_A/(w_A + w_G)** (Watson reads only); G>T = (w_T + c_T)/m; G>C = (w_T + w_C + c_C)/m.
A position is written to the candidates file if any of its three rates is **> minhetfreq** (0.1). A 0/0
rate is NaN and never passes. No depth limit at this stage.

Candidates file columns (lines 131-137): chrom, 1-based pos, ref base, Watson counts A,T,C,G, Crick counts
A,T,C,G, mean base quality Watson A,T,C,G, mean base quality Crick A,T,C,G (rounded, `(unsigned)(sum/n +
0.5)`; written as 0 when n = 0 in the study's files), mean MAPQ Watson, mean MAPQ Crick. The header line
the C code writes (line 856) names other columns and is skipped by the Perl script.

### 2.4 Genotype (Perl `sub Bayes`, BS-Snper.pl 1329-1621) [verified]

For each candidate line, diploid genotypes are scored. For a **reference C** the read numbers are:
A = w_A; **T = c_T**; **C = w_C + c_C** (so every Watson read, converted or not, supports C); G = w_G + c_G
(c_A is not used). Error rates e_X = min over the two classes of 10^(-meanQ/10) (an absent base has meanQ
0, e = 1). A genotype is scored only if its bases were seen (for example `ct` needs a T and a C). Log
likelihoods (nn = log multinomial coefficient, common to all):

- cc = nn + C ln(1-e_C) + (A+T+G) ln(e_C/3)
- ct = nn + (C+T) ln((1-e_T)/2) + (A+G) ln(e_CT/3), with e_CT the Crick T error
- tt = nn + T ln(1-e_T) + (A+C+G) ln(e_T/3)
- gt = nn + (T+G) ln((1-e_T)/2) + (A+C) ln(e_T/3)
- cg, gg, aa, at, ac, ag analogously (lines 1385-1458).

Genotypes whose score is exactly 0 are dropped (line 1465). Priors for a reference C (lines 1487-1498):
cc 0.985, ct 0.000667, tt 0.00033, ac 0.00017, cg 0.00017, aa 0.000083, gg 0.000083, at 1.1e-7, gt 1.1e-7,
ag 2.78e-8. Genotype = largest log likelihood + log prior. QUAL uses 2.7 instead of e and does not affect
the genotype. Post-rules (lines 1580-1606): **CT with c_T = 0 becomes CC; CT with c_T of 1 or 2 becomes NN
(nothing written)**; TT with c_C > 1 and c_C/(c_C+c_T) > 0.02 becomes NN; AG with w_A = 0 becomes GG, with
w_A of 1 or 2 NN; AA with w_G > 1 and w_G/(w_A+w_G) > 0.02 becomes NN.

Effect, computed with a port of this code (`code/bayes_refC.py`, base quality 37): CT is chosen over CC at
a reference C when roughly 8.9 x T > 7.3 + 0.69 x C; TT only when there is essentially no C support,
which requires no Watson coverage. Smallest number of opposite-strand T reads that gives a call:

| Crick C reads \ Watson reads | 0 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|
| 0 | 1 (TT) | 3 | 3 | 3 | 4 |
| 5 | 3 | 3 | 3 | 3 | 5 |
| 20 | 3 | 3 | 4 | 4 | 6 |
| 50 | 6 | 6 | 6 | 7 | 8 |

### 2.5 Records written (Perl `sub genotype`, lines 193-1320) [verified]

For a reference C:

| Genotype | ALT | GT | depth (DP) and PASS condition |
|---|---|---|---|
| CT (lines 999-1023, comment "only use crick") | T | 0/1 | DP = c_C + c_T; PASS if DP >= mincover, c_T >= minread2, mean Crick T quality >= minquali, c_T/DP (3 decimals) >= minhetfreq; else Low |
| TT (lines 695-727) | T | 1/1 | DP = c_T + c_C + w_T + w_C; PASS if DP >= mincover, c_T + w_T >= minread2, quality, (c_T + w_T)/DP >= minhomfreq; else Low |
| GT (lines 1119-1150) | T,G | 1/2 | PASS if depth, both alleles >= minread2 and >= minhetfreq; else Low |
| AT | A,T | 1/2 | (not matched by the study's grep) |
| CC, NN | nothing | | |

Every non-reference genotype is written, as PASS or Low. FORMAT is `GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR`, where
BSD = the 8 counts (Watson A,T,C,G, Crick A,T,C,G), BSQ = their mean qualities, ALFR = reference and
alternative fractions. In the study's tables a heterozygous-labelled record often has ALFR 0.000,1.000
(for example `0/1 ... BSD 0,0,19,0,0,5,0,0`), because Watson reads count as C in the genotype step.

### 2.6 Directional assumption and the non-directional libraries [verified]

BS-Snper assigns strands from FLAG 0x10 and 0x80 as for a directional library, and prints the fixed
string `##Bisulfite=directional>` (BS-Snper.pl line 127). Both studies are non-directional (`--non_directional`).
Bismark 0.22.3, the version of both studies, swaps the R1/R2 FLAG values of the two complementary strands
"so that we don't end up with discordant pairs" (bismark 0.22.3 lines 8305-8350): OT 99/147, CTOT 147/99,
OB 83/163, CTOB 163/83 (read 1/read 2). CTOT reads carry top-strand information and CTOB bottom-strand
information, so with these FLAG values BS-Snper's Watson class is OT + CTOT and its Crick class is OB +
CTOB, which is the correct grouping. The rule therefore holds for these non-directional Bismark BAMs,
provided the default FLAG values (not `--old_flag`) were used, which the job scripts show [verified:
neither job script has `--old_flag`]. Bismark's `XG:Z:CT` (top) / `XG:Z:GA` (bottom) tag gives the same
grouping without relying on FLAG [assumed equivalent; to be checked on a sample of reads].

---

## 3. Specification of an R implementation

Goal: for a given set of CpGs and a set of BAM files, produce the set of CpGs that Venkataraman et al.'s
procedure would exclude, with the counts behind each call. BS-Snper itself is not run.

### 3.1 Inputs

- `bam_files`: the 22 deduplicated, coordinate-sorted, indexed Bismark BAMs of the June 2021 run
  (`/panopea/061021-big/`, produced by `deduplicate_bismark` then `samtools sort` and `samtools index`,
  `061021-big_job.sh` lines 70-113). Their presence on gannet is [assumed] and has to be checked.
- `genome`: `cgigas_uk_roslin_v1_genomic-mito.fa` (pipeline cache).
- `cpgs`: data frame `chr, start0` as in the coverage files (C at 1-based start0 + 1, G at start0 + 2).
  For the DML stage the tested CpGs suffice; for the global characterization all CpGs of the coverage
  union at the configured threshold.
- Parameters, with BS-Snper defaults: `minhetfreq = 0.1`, `minhomfreq = 0.85`, `minquali = 15`,
  `mincover = 5` (as in the study; only used for the FILTER label), `minread2 = 2`, `errorate = 0.02`,
  `mapvalue = 20`. `maxcover` is not needed (unused by BS-Snper).
- Switches (each documented; defaults reproduce the 2022 study): `use_filter = "any"` (as in the study;
  alternative `"PASS"`), `positions = "C_only"` (as in the study; alternative `"C_and_G"`, which adds G>A
  at the G), `calls = c("pooled", "per_library")` (both, as in the study), `count_both_mates = TRUE`,
  `cigar = "correct"` (alternative `"bs_snper"` to reproduce the last-operation behaviour).

### 3.2 Per-read processing (per BAM, per chromosome, in chunks)

Read with `Rsamtools::scanBam` / `GenomicAlignments::readGAlignments` (fields: flag, pos, mapq, cigar,
seq, qual, tags XG and MD).

1. Drop FLAG 0x4 and MAPQ < `mapvalue`. Keep secondary, QC-fail and duplicate flags (BS-Snper does not
   check them; Bismark sets none of them, to be confirmed on the files).
2. Class: Watson if (0x10 unset and 0x80 unset) or (0x10 set and 0x80 set); otherwise Crick. Cross-check
   against XG (CT = Watson, GA = Crick) on a sample; stop if they disagree.
3. Lay out SEQ and QUAL on the reference with the full CIGAR (`GenomicAlignments::sequenceLayer`,
   `pileLettersAt`), or with the last CIGAR operation only when `cigar = "bs_snper"`.
4. Set to N every base with Phred < `minquali` and every base at a reference N.
5. Mismatch count per read: non-N bases that differ from the reference, except Watson T at reference C and
   Crick A at reference G. Drop the read if count > max(1, floor(L x errorate + 0.5)), L = laid-out length
   (reference span with deletions). The MD tag of Bismark lists the reference base at every mismatch and
   avoids fetching the reference for every read.
6. At each target position (the C, and the G when `positions = "C_and_G"`): Watson T at reference C is
   counted as C; Crick A at reference G as G; any other base as itself. Keep per position and library:
   `w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G` and the sums of base qualities for the same 8 cells
   (sums, not means, so that libraries can be added).

Pooled ("merged") call: the 8 counts and 8 quality sums **summed over the libraries**. Because every
BS-Snper filter acts on single reads (MAPQ, base quality, per-read mismatches) and every later step uses
only these sums, this equals the result of running on a merged BAM, without creating one.

### 3.3 Per-position decision (port of printSnp, Bayes and genotype)

For each position and each call (pooled, and each library):
1. Candidate: rates of section 2.3 for the reference base; keep if any > `minhetfreq` (strict).
2. Mean qualities: floor(sum/n + 0.5), 0 when n = 0.
3. Bayes genotype and post-rules exactly as in section 2.4 (port all 10 genotypes and the priors for the
   reference base; ties are resolved arbitrarily by Perl's hash order, so record ties).
4. Record ALT, GT, DP and FILTER as in section 2.5.
5. Variant at a C of a CpG if REF = C and ALT is "T" or "T,G" (CT, TT or GT). With `positions =
   "C_and_G"`, also REF = G and ALT "A" or starting with "A," at the G (the mirror: AG, AA, AT genotypes
   at a reference G).
6. With `use_filter = "PASS"`, keep only PASS records.

### 3.4 Outputs

- `variant_calls.tsv.gz`: one row per (position, call) written: `chr, start0, pos1, ref, call`
  (`pooled` or library id), `w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, q_w_A ... q_c_G`, `genotype, alt, gt,
  dp, alt_freq, filter, qual`.
- `excluded_cpgs.tsv.gz`: one row per CpG of `cpgs`: `chr, start0, excluded` (TRUE if any written variant
  under the switches), `n_libraries_called, called_pooled, any_pass, max_crick_depth, max_crick_T`.
- A log with the numbers of section 1.3 for this project (positions, Low-only share, per-library
  counts, DML removed by stage).

### 3.5 Where to apply it (as in the 2022 study)

- Remove excluded CpGs from the CpG sets of the global methylation characterization and from the
  background of the genomic-location tests.
- Remove excluded CpGs from the DML lists **after** the test, leaving the united table, the model fits
  and the q-values unchanged. Report how many DML were removed per contrast.
- Toggle off by default (the author's decision): the default run must not change.

---

## 4. Where an R implementation differs from BS-Snper, and what it costs

| Point | BS-Snper | R (proposed) | Likely effect on calls |
|---|---|---|---|
| CIGAR | last operation only (bug) | full CIGAR (or emulate) | R uses all bases of reads with indels and places them correctly: slightly more depth and fewer spurious mismatches near indels; BS-Snper's misplaced bases can create false variants or push reads over the mismatch limit. Small share of reads; emulation option gives exactness |
| Overlapping mates | both mates counted | both counted (default) | match; counting once would lower DP and T counts and remove some calls that rest on one fragment |
| Base quality | Phred < 15 to N before everything | same | match if implemented as specified; note `Rsamtools::pileup` defaults to `min_base_quality = 13` and `max_depth = 250`, which must not be used as defaults |
| Mean base quality | rounded integer means per class and base | same rounding | match |
| MAPQ | < 20 dropped | same | match |
| Per-read mismatch filter | as 2.2 | same, on full CIGAR | small differences for indel reads |
| Pooled call | merged BAM | sum of per-library counts | identical by construction (2.2, 3.2) |
| FLAG checks | only 0x4 | only 0x4 | Rsamtools pileup silently drops secondary/QC-fail/duplicate; use scanBam to match |
| Floating point | C float for rates and the read limit; Perl doubles | R doubles | differences only at exact ties with 0.1 or with x.5 read limits |
| Genotype ties | Perl hash order (random) | deterministic | negligible |
| Positions | whole genome, then CpG motifs | CpG positions only | same result for the C of CpGs |
| Validation | none here | see section 6 | without a check against BS-Snper output an R port can carry an undetected difference |

Not reproducible exactly without running BS-Snper: nothing in the decision logic; the differences are
the CIGAR behaviour (selectable) and arithmetic edge cases. The cost of choosing correct CIGAR handling is
a small, unmeasured difference in calls near indels.

---

## 5. Diploids and triploids

BS-Snper's model is diploid (genotypes of two alleles). For a triploid, one variant copy of three gives an
expected T share of 1/3 among the opposite-strand reads and two copies give 2/3. Both are above
`minhetfreq` (0.1), so the frequency thresholds do not stop such a site, and both are below `minhomfreq`
(0.85), so a PASS call is labelled heterozygous (0/1) whatever the dosage. The constraints that bind are
the Bayes step and the rule that a CT call needs at least 3 opposite-strand T reads. Probability of a call
at opposite-strand depth n with the same number of same-strand reads (base quality 37, no sequencing
error; `code/bayes_refC.py`; first value = the 2022 rule, any FILTER; second = PASS with mincover 5):

| Allele share | n = 3 | n = 5 | n = 8 | n = 10 | n = 15 | n = 20 | n = 30 |
|---|---|---|---|---|---|---|---|
| 1 of 3 (triploid) | 0.04 / 0.00 | 0.21 / 0.21 | 0.53 / 0.53 | 0.70 / 0.70 | 0.92 / 0.92 | 0.94 / 0.94 | 0.96 / 0.96 |
| 1 of 2 (diploid heterozygote) | 0.12 / 0.00 | 0.50 / 0.50 | 0.86 / 0.86 | 0.95 / 0.95 | 1.00 | 1.00 | 1.00 |
| 2 of 3 (triploid) | 0.30 / 0.00 | 0.79 / 0.79 | 0.98 / 0.98 | 1.00 | 1.00 | 1.00 | 1.00 |
| homozygous | 1.00 / 0.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

So, per library, a single-copy allele in a triploid is found less often than a heterozygous allele in a
diploid at the same depth (0.21 against 0.50 at 5 opposite-strand reads, 0.70 against 0.95 at 10). The
libraries here have 6 to 12 reads on average over the genome (issues.md, I-49, item "What would settle
it"); about half of the reads at a position come from each strand class, so a typical C has about 3 to 6
opposite-strand reads [assumed equal split], where per-library detection is low for every genotype except
homozygotes. Masking
is therefore incomplete and its sensitivity depends on ploidy, which matters for a ploidy contrast.

Pooled call over 22 libraries: the candidate rule needs a pooled T share > 0.1, and at pooled depths of 40
to 200 opposite-strand reads the Bayes step needs about 15% to 18% (7 of 40, 14 of 90, 30 of 200). With
equal depth per library that is an allele carried by about 7 of 22 animals as diploid heterozygotes or
about 11 of 22 as single-copy triploids. Rare alleles are found only by the per-library calls.

---

## 6. Validation without running BS-Snper

The 2022 study's BAMs are on gannet (`/spartina/project-gigas-oa-meth/output/bismark-roslin/zr3616_*_R1_val_1_val_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam`,
9.6 GB each, `.bai` 1.8 MB, byte ranges accepted [verified by HEAD request]) and their C/T tables are in
`downloads/gannet_bssnper/`. Running the R port on one small sequence (for example NC_001276.1 or one
scaffold) of the 8 BAMs, fetched by range requests, and comparing position by position with the merged and
per-animal tables (counts in BSD, genotype, FILTER) checks the port against real BS-Snper output, once with
`cigar = "bs_snper"` (should match exactly) and once with `cigar = "correct"`.

---

## 7. Open points for the author

1. Plus-strand C only (exactly as in 2022) or both C and G of each CpG. The coverage files merge both
   strands, so the 2022 rule leaves CpGs whose minus-strand C is a variant. Proposed: the toggle reproduces
   2022 (C only); C and G as a second, separately labelled setting.
2. Any FILTER (as in 2022; 16% of their excluded positions rest only on Low records) or PASS only.
3. Pooled call over all 22 libraries (as in 2022, where the 8 animals were one experiment) or per
   experiment; this project has two experiments sequenced separately (I-01).
4. Whether the BAMs of the 22 libraries are still available.

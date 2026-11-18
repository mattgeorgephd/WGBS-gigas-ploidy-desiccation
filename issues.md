# Issue log

Project: *Epigenetic divergence in the thermal stress response of diploid and triploid Pacific oysters*.
Opened 2026-11-18 on branch `reproducibility-audit`, from a full read of the manuscript draft (Google Doc as of 2026-11-18), every script and notebook in the repository at commit `523d734`, and the project directory on gannet (843 files, 217 GB).

**How to read this file.** Each issue has an ID, a severity, the evidence, why it matters, and a proposed resolution. Issues that need a decision or information from the authors say so. The plan that these issues feed into is [`07_documentation/reproducibility_plan.md`](07_documentation/reproducibility_plan.md).

**Paths.** File paths in the evidence are the paths at commit `523d734`, before the reorganization. After the reorganization the same file is at `08_archive/<same path>` unless the plan lists it as moved to `03_input_files/`.

**Severity.**
- **S1**: a main conclusion of the manuscript depends on it.
- **S2**: it changes a reported number or a figure, or materially changes how a result should be read.
- **S3**: wording, labels, small inconsistencies.
- **R**: repository and reproducibility engineering; no effect on results by itself.

**Status values.** `open`, `needs decision` (authors), `needs information` (authors), `fix planned` (phase in the plan), `resolved`.

**How things were verified.** Statements marked *verified* were checked directly: by reading the code, by recomputing from files in the repository or on gannet, or by rerunning the step in R 4.3.3 with methylKit 1.28.0. Anything not verified is labelled as an inference or as undetermined.

---

## Summary

| ID | Sev. | Area | Issue | Status |
|---|---|---|---|---|
| I-01 | S1 | Design | Heat treatment is completely confounded with experiment and sequencing batch in the WGBS data | needs decision |
| I-02 | S1 | Statistics | "Interaction-DML" are the overlap of two main-effect lists, not a tested interaction | needs decision |
| I-03 | S1 | Statistics | GO enrichment treats every transcript isoform as an independent gene; the enrichment results are artifacts of isoform number | needs decision |
| I-04 | S2 | Design | Four of the ten heat libraries also received a 45 C acute heat shock; not disclosed, and mislabelled in one spreadsheet | needs decision |
| I-05 | S2 | Code | The 99.9th percentile coverage filter described in the Methods was never applied (argument name typo) | needs decision |
| I-06 | S2 | Statistics | Differential methylation model has no overdispersion correction although groups are biological replicates | needs decision |
| I-07 | S2 | Text vs code | Ploidy-DML are described as control-only ("baseline") differences but were computed across all 22 libraries | open |
| I-08 | S2 | Design | Sequence variants were not masked; ploidy groups and cohorts differ genetically | open |
| I-09 | S2 | Interpretation | DML are spatially clustered; one 2.2 kb locus supplies 39 of the 200 "interaction" loci and all of the reported lncRNA signal | open |
| I-10 | S2 | Statistics | ELISA conclusions for ploidy and interaction depend on the transformation, which was chosen by a test on the wrong quantity | needs decision |
| I-11 | S2 | Figure | Figure 1A drops three of the eight diploid heat samples through its axis limits | fix planned (3) |
| I-12 | S2 | Data | ELISA calculation rests on unverified inputs (DNA mass per well, standards, extrapolation) | needs information |
| I-13 | S2 | Data | Figure 1B values come from a hand-built sheet that mixes two alignment runs and carries four wrong sample names | fix planned (3) |
| I-14 | S2 | Number | The "all CpG" genomic background in Figure 3A and Results 3.3 was not computed from the analysed CpGs | fix planned (4) |
| I-15 | S2 | Number | The 1x genomic distribution in Results 3.1 double counts exons and comes from an unidentified CpG set | fix planned (4) |
| I-16 | S3 | Number | CpG counts: the 1x value is defined differently from the 5x, 10x, 20x values; 480,479 is a typo | fix planned (3) |
| I-17 | S2 | Statistics | The chi-squared test of genomic location has spreadsheet errors and is not reported in Results | fix planned (4) |
| I-18 | S3 | Method | Genomic feature percentages are shares of overlap events, not of loci | open |
| I-19 | S3 | Figure | Figure 3A: hypo and hyper counts swapped; Figure 2D legend says "desiccation" | fix planned (4) |
| I-20 | S2 | Figure | Figure 5 plots row counts of join tables, not the quantities its panel titles name | needs decision |
| I-21 | S2 | Text vs code | "Enriched" terms behind Figure 9, Table 1, and Section 3.5 are the top ten by p-value, not p < 0.01 | needs decision |
| I-22 | S2 | Figure | Figures 6 to 8 depend on the REVIGO web server; inputs and labels cannot be regenerated | fix planned (4) |
| I-23 | S3 | Figure | Figure 2C to 2G: PCA on selected loci is circular; loci are matched by start coordinate only | fix planned (4) |
| I-24 | S2 | Text | The "CAD" gene is misidentified in the Discussion and Conclusions | open |
| I-25 | S3 | Number | Chromosome percentages are hand-computed and do not match the DML lists exactly | fix planned (4) |
| I-26 | S3 | Figure | Figure S1 B and C were built from hand-transcribed histogram values | fix planned (3) |
| I-27 | S3 | Text | Supplementary lists are cited but not assembled | fix planned (6) |
| I-28 | S2 | Text | Methods describe several steps differently from what was run | open |
| I-29 | S2 | Text | qPCR and "mitochondrial abundance" claims appear only in the Discussion and are not supported by the data in the repository | needs decision |
| I-30 | S3 | Text | Manuscript structure: empty abstract, numbering gaps, orphan section, about 60 citations missing from the reference list | open |
| I-31 | R | Code | No analysis script runs as committed | fix planned (3, 4) |
| I-32 | R | Code | Non-R and manual steps (Python, bash, bedtools, spreadsheets, PowerPoint, REVIGO) | fix planned (3, 4) |
| I-33 | R | Data | Inputs on a retired server; one annotation source has no surviving copy | fix planned (2) |
| I-34 | R | Repo | Repository hygiene: duplicates, lock files, misplaced files, 613 MB | fix planned (1) |
| I-35 | R | Data | gannet directory: about 200 GB of redundant intermediates, one incomplete upload, one mislabelled sample set | needs decision |
| I-36 | R | Data | Sample metadata are incomplete for the control libraries | needs information |
| I-37 | R | Environment | Software versions are not recorded anywhere | needs information |

---

## Design and statistics

### I-01 (S1) Heat treatment is completely confounded with experiment and sequencing batch in the WGBS data

**Where.** Methods 2 and 2.2; every result that uses the heat contrast: heat-DML, "interaction-DML", Figure 1B, Figure 2, Figure 3, Figures 5 to 9, Table 1, Discussion 4.2 to 4.9.

**Evidence (verified).**
- `metadata.xlsx`, sheet `list`: the 10 heat libraries are all Zymo project zr3534 (`dataset = ronit`); the 12 control libraries are all Zymo project zr3644 (`dataset = yaamini`). `code/2_WGBS_Methylkit-all-gigas.R:149-150` builds the heat covariate as `rep("yes", 10), rep("no", 12)` in that file order. The author has confirmed that the controls are animals from a separate experiment.
- The two sets were sequenced as different projects and differ systematically as libraries (Bismark summary of the June 2021 alignment, `07_documentation/upstream_processing/061021-big_bismark_summary_report.txt`):

  | | Heat (zr3534, n = 10) | Control (zr3644, n = 12) |
  |---|---|---|
  | Read pairs, mean (range) | 77.2 M (66.1 to 90.6 M) | 49.6 M (36.2 to 58.1 M) |
  | Duplicate alignments, mean (range) | 18.3% (16.5 to 20.6%) | 11.6% (10.2 to 12.6%) |
  | CpGs with at least 10x per library | 4.41 to 5.51 M | 1.98 to 4.00 M |
  | Non-CpG methylation (CHH), mean | 0.51% | 0.50% |

  Read depth, duplication rate, and per-library 10x coverage do not overlap between the two sets. Non-CpG methylation is equal, so bisulfite conversion was comparable.
- The control libraries are not the heat experiment's own controls. The sample guide in `README.md` lists control animals D01 to D08 and T01 to T08 for the heat experiment; none was sequenced (`gene_expression/data/zr3534_wgbs_info.csv` lists only heat-treated animals).

**Why it matters.** When a batch variable is perfectly aligned with the variable of interest, no statistical adjustment can separate them (Leek et al. 2010, Nat Rev Genet 11:733, doi:10.1038/nrg2825). Every difference between the two experiments is an equally valid explanation for a "heat-DML": source population and family, age, season, holding conditions, tissue handling, DNA extraction, library batch, flow cell, and depth. The 2,351 heat-DML, the overlap set built from them, and the percent CpG methylation interaction in Figure 1B (F(1,18) = 14.49) cannot be attributed to heat on the present evidence. The interaction in Figure 1B is, numerically, opposite-sign ploidy differences in the two experiments (heat set: diploid minus triploid = +0.86 points; control set: -0.65 points).

**What is not affected.** The ELISA (N = 32) used the heat experiment's own control animals on one plate, so its heat effect is free of this confound (see I-10 and I-12 for its own caveats). The ploidy contrast is balanced across the two experiments (5 + 6 diploid against 5 + 6 triploid), so batch does not align with ploidy.

**Options.**
1. Sequence the heat experiment's own controls (D01 to D08, T01 to T08 ctenidia), if tissue or DNA remains. This is the only route to an unconfounded heat contrast.
2. Keep the data and reframe: report the contrast as "heat-treated cohort versus control cohort", state the confound in Methods and in the limitations, and remove mechanistic interpretation of heat-DML.
3. Shift the weight of the paper to what the design supports: the ploidy contrast, tested within each cohort and reported where it replicates in both, and the ELISA for the heat response.

The pipeline will produce the within-cohort ploidy analyses and the replication table under any of these options.

**Needs from authors.** Whether option 1 is possible; otherwise which framing. Also the provenance of the control animals (I-36).

---

### I-02 (S1) "Interaction-DML" are the overlap of two main-effect lists, not a tested interaction

**Where.** Methods 2.5; Results 3.2, 3.3, 3.5; Figures 2G, 3A, 4, 8, 9; Table 1; Discussion 4.3 and Conclusions.

**Evidence (verified).** `code/2_WGBS_Methylkit-all-gigas.R:400-413` builds the set by merging the ploidy-DML and heat-DML tables on their row names, then takes the rows and the `meth.diff` values from the heat table. Recomputed from `tables/DML-getMethylDiff-*.csv`: the 200 "interaction" loci are exactly the intersection of the 1,038 ploidy-DML and the 2,351 heat-DML, and their `meth.diff` column equals the heat contrast. The three classes in the text are sign combinations of the two main effects: 28 with both positive, 49 with both negative, 123 with opposite signs.

**Why it matters.** A locus that differs by ploidy and also differs by heat shows two additive effects. It says nothing about whether the response to heat differs between ploidies, which is what the manuscript claims ("methylation responses to heat differed between diploids and triploids"; "thermal exposure can actively reverse baseline ploidy-associated methylation differences"; "contrasting methylation patterns between diploids and triploids following heat exposure"). "Opposite sign" here means only that triploids are higher than diploids while the heat cohort is lower than the control cohort, or the reverse. Inferring an interaction from two separate significance calls is a known error (Gelman and Stern 2006, Am Stat 60:328, doi:10.1198/000313006X152649; Nieuwenhuis et al. 2011, Nat Neurosci 14:1105, doi:10.1038/nn.2886).

**Options.**
1. Rename the set to what it is (for example "loci shared by the ploidy and heat contrasts") and rewrite the interpretation.
2. Fit a model with a ploidy by heat term. In methylKit this can be done by using the product of the two indicators as the tested variable with both main effects as covariates; DSS fits the general design directly (Park and Wu 2016, Bioinformatics 32:1446, doi:10.1093/bioinformatics/btw026).
3. Both: option 1 as description, option 2 as the test.

Under I-01, a formally tested interaction is still ploidy by cohort, and the text must say so.

**Needs from authors.** Choice among the options. The pipeline will implement the legacy overlap set (to reproduce the manuscript numbers) and the formal test.

---

### I-03 (S1) GO enrichment treats every transcript isoform as an independent gene

**Where.** Methods 2.6; Results 3.4 and 3.5; Figures 5 to 9; Table 1; most of Discussion 4.1 to 4.3 and 4.6 to 4.9.

**Evidence (verified).**
- `code/5_WGBS_Functional-enrichement-analysis-BOTH.R:186-206`: the topGO universe is the list of transcript IDs in `geneid2go-union1x-BOTH.tab`, and a transcript is flagged as "of interest" if its gene contains a DML. Every isoform of a gene with one DML is therefore counted as a separate significant gene carrying the same GO terms.
- The RefSeq annotation has 63,341 mRNAs for 31,371 protein-coding genes; 613 genes have 10 or more isoforms.
- The top terms are single genes repeated. Recomputed from `GOterms/BOTH-*/all-*-EnrichedGO-DML-withTranscript.csv` (transcripts annotated to each top term) and the DML annotation tables on gannet (`analyses/GOterm-annotation/BOTH/*.GOAnnot`):

  | DML set, ontology | Example top term | Annotated / significant transcripts | Genes behind the significant transcripts |
  |---|---|---|---|
  | Interaction, BP | pyrimidine nucleobase catabolic process, p = 1.8e-21 | 12 / 12 | 1: dihydropyrimidinase (12 isoforms, one DML) |
  | Interaction, BP | nuclear matrix organization, p = 1.8e-21 | 12 / 12 | 1: nesprin-1 (12 isoforms) |
  | Heat, BP | iron ion import across plasma membrane, p < 1e-30 | 73 / 51 | 2: LAP (50 isoforms) and LOC117680672 |
  | Ploidy, BP | gephyrin clustering involved in postsynaptic density assembly, p = 4.3e-19 | 13 / 13 | 1: neurexin-1 (13 isoforms) |

  Across the three DML sets, the ten top BP terms rest on DML in 7 genes (ploidy), 4 genes (heat), and 3 genes (interaction). Of all 90 top-ten terms (three sets, three ontologies), 70 rest on two or fewer genes and 43 on exactly one.
- No multiple testing correction: 12,720 BP terms were tested per DML set with a raw p < 0.01 cutoff (`:1009-1063`), using the classic Fisher test that ignores the GO hierarchy (`:255`).
- The background does not match the Methods. Methods: genes with CpGs at 10x in all samples. Code: transcripts overlapping a CpG with at least 1x in any sample (`code/4_WGBS_GOterm-Annotation-BOTH.ipynb`, cells 55 to 74), which is nearly every annotated transcript.
- Annotation bug: cell 10 of the same notebook rebuilds transcript IDs with `sed 's/^/XM_/'`, which relabels the 157 `NM_` (curated RefSeq) and 8,059 `XR_` queries in the blastx file as `XM_`; the `NM_` transcripts then fail to join and lose their GO terms.

**Why it matters.** The reported p-values (1e-17 to below 1e-30) measure how many isoforms a few genes have, not functional enrichment. The gene lists in Table 1 and Figure 9 are therefore lists of isoform-rich genes that contain at least one DML, and the functional narrative built on them (DNA replication and meiosis for ploidy; MAPK, ER overload, RNA processing for heat; tRNA modification and "alarm system" thresholds for interaction; AK9, CAD, nesprin-1, LAP, dihydropyrimidinase as mechanisms) has no statistical support from this analysis. Two further biases compound it: gene-set tests on methylation data are biased toward genes with many assayed CpGs unless that is modelled (Geeleher et al. 2013, Bioinformatics 29:1851, doi:10.1093/bioinformatics/btt311; Maksimovic et al. 2021, Genome Biol 22:173, doi:10.1186/s13059-021-02388-x), and a background that includes genes that could never be detected inflates enrichment (Timmons et al. 2015, Genome Biol 16:186, doi:10.1186/s13059-015-0761-7).

**Proposed resolution.** Test at gene level (one entry per gene, GO terms pooled over its transcripts); use as background only genes containing at least one of the CpGs that were tested; account for the number of tested CpGs per gene; control the false discovery rate; use a hierarchy-aware algorithm (Alexa et al. 2006, Bioinformatics 22:1600, doi:10.1093/bioinformatics/btl140). The legacy analysis will be reproduced unchanged alongside it so the two can be compared. It is likely that few or no terms survive; the manuscript should be prepared for that outcome.

**Needs from authors.** Agreement to make the gene-level analysis primary.

---

### I-04 (S2) Four of the ten heat libraries also received a 45 C acute heat shock

**Where.** Methods 2 and 2.2 (not mentioned); all heat contrasts.

**Evidence (verified).** `gene_expression/data/zr3534_wgbs_info.csv`: zr3534_4, 5, 9, 10 are animals D19, D20, T19, T20 with `Heat_Stress = yes`. The sample guide in `README.md` defines these as 24 h at 27 C followed by 1 h at 45 C. `metadata.xlsx` sheet `list3` agrees. In `data/bismark-percent-mCPG-summary.xlsx` the same four libraries are typed as D14-C, D15-C, T14-C, T15-C, which are different animals.

**Why it matters.** The "heat" group in the WGBS data is six emersion-only animals plus four that also had an acute heat shock. The manuscript describes one treatment. The ELISA used emersion-only animals (D11 to D18, T11 to T18), so the two assays do not cover the same treatment either.

**Options.** Disclose and keep pooled; drop the four (leaving three per ploidy); or include the heat shock as a covariate. The pipeline will report the sensitivity of the DML lists to this choice.

**Needs from authors.** Choice, and confirmation of the animal identities.

---

### I-05 (S2) The 99.9th percentile coverage filter described in the Methods was never applied

**Where.** Methods 2.5 ("removing CpG in the 99.9th percentile of coverage"); the CpG set behind every WGBS result.

**Evidence (verified).** `code/2_WGBS_Methylkit-all-gigas.R:115-121` and `code/1_coverage_analysis.R:145-224` call `filterByCoverage(..., lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9)`. The methylKit arguments are named `hi.count` and `hi.perc`. The misspelled names fall into the function's `...` argument and are ignored, so no upper filter runs. Recomputed from the 22 coverage files on gannet: CpGs with at least 10x in all libraries number **480,579 with no upper filter**, exactly the manuscript value, and **454,555 with the filter applied as described**.

**Why it matters.** About 26,000 analysed CpGs (5.4%) are the high-coverage sites the filter was meant to remove, typically repeats, collapsed duplications, and the mitochondrial genome. DML counts will change when the filter is applied; the size of the change will be reported by the pipeline.

**Options.** Apply the filter and update all downstream numbers, or keep the analysis as run and correct the Methods. Applying it is the more defensible choice.

**Needs from authors.** Choice. The pipeline implements both.

---

### I-06 (S2) No overdispersion correction in the differential methylation model

**Where.** Methods 2.5; all DML.

**Evidence (verified).** `code/2_WGBS_Methylkit-all-gigas.R:154-162` calls `calculateDiffMeth()` with a covariate and all other arguments at their defaults, which means `overdispersion = "none"`: a binomial logistic regression in which reads from different animals are treated as exchangeable draws.

**Why it matters.** With biological replicates, variation among animals exceeds binomial sampling variation, and a model without overdispersion is anticonservative. The methylKit authors' own benchmark shows a markedly higher false positive rate without the correction and recommends `overdispersion = "MN"` with a chi-squared test (Wreczycka et al. 2017, J Biotechnol 261:105, doi:10.1016/j.jbiotec.2017.08.007). The 20% difference threshold limits but does not remove the problem.

**Proposed resolution.** Report a sensitivity analysis with the correction; if the authors agree, make it primary.

**Needs from authors.** Whether to change the primary model.

---

### I-07 (S2) Ploidy-DML are described as control-only differences but were computed across all 22 libraries

**Where.** Results 3.2 ("under control conditions"), 3.4 ("baseline epigenetic differences"), 3.5; Figure 6 caption ("before stress treatment"); Figure 9 caption ("controls only"); Discussion 4.1 ("prior to environmental stress exposure").

**Evidence (verified).** `code/2_WGBS_Methylkit-all-gigas.R:143-156`: ploidy is tested on all 22 libraries with heat as a covariate. A control-only ploidy analysis exists in the legacy outputs (`yaamini`, 12 libraries) but is a different DML list and is not the one reported.

**Resolution.** Change the text to describe a ploidy main effect adjusted for cohort, or report the control-only analysis. The pipeline produces both and the overlap between them.

---

### I-08 (S2) Sequence variants were not masked

**Where.** Methods 2.4 and 2.5.

**Evidence.** No step in the repository identifies or removes C/T polymorphisms. The Methods cite Venkataraman et al. (2022) for methodology; that project included a C/T SNP identification step with BS-Snper (verified: its OSF archive, doi:10.17605/OSF.IO/YGCTB, contains `code/07-BS-SNPer.ipynb`), which has no counterpart here. Diploid and triploid oysters in each cohort come from different crosses, and the two cohorts are different stocks (inference from the design; to be confirmed by the authors).

**Why it matters.** In bisulfite data a C-to-T variant is indistinguishable from an unmethylated cytosine (Liu et al. 2012, Genome Biol 13:R61, doi:10.1186/gb-2012-13-7-r61). Genetic differences between stocks can therefore appear as DML with large effect sizes, which is the class the 20% threshold selects for. In triploids allele dosage adds intermediate values.

**Proposed resolution.** Flag DML whose per-animal methylation values sit at genotype-like levels; if BAM files are still available (gannet `/panopea/061021-big/`), call variants and mask them. To be scoped after the legacy reproduction.

---

### I-09 (S2) DML are spatially clustered; one locus supplies the reported lncRNA signal

**Where.** Results 3.3 ("pronounced enrichment in lncRNA loci (15.5%)"), 3.4, Figure 3A; Discussion.

**Evidence (verified, recomputed from the DML lists and `cgigas_uk_roslin_v1_gene.gff`).**
- 29% of ploidy-DML, 33% of heat-DML, and 30% of overlap loci have another DML of the same list within 1 kb.
- A single 2.2 kb region, NC_047563.1:39,633,323-39,635,520, contains 53 ploidy-DML, 39 heat-DML, and 39 of the 200 overlap loci. It overlaps two lncRNA genes (LOC105318567, LOC105341636) and one protein-coding gene (LOC105317102). These 39 loci are 39 of the 42 overlap loci counted in lncRNA features.

**Why it matters.** The lncRNA "enrichment" among interaction loci, and the link drawn to "gene silencing by RNA", describe one locus, not a genome-wide pattern. Counting neighbouring CpGs as independent loci also overstates the number of independent differences throughout. A block of many adjacent CpGs that differs between groups of animals from different stocks is also what a structural or copy-number variant looks like (see I-08).

**Proposed resolution.** Report regions as well as loci (number of distinct regions and genes), and describe this locus explicitly.

---

### I-10 (S2) ELISA conclusions for ploidy and interaction depend on the transformation

**Where.** Results 3.1 (5-mC statistics); Figure 1A letters; Discussion 4.5.

**Evidence (verified in R).** `code/5-mC & 5mCpG plots.Rmd:84-101` applies an arcsine square-root transformation when a Shapiro-Wilk test on the pooled raw values is significant. Normality is an assumption about model residuals, not about the pooled response; the residuals of the untransformed model pass (p = 0.59).

| Scale | Ploidy | Heat | Ploidy x heat | Residual normality | Equal variance |
|---|---|---|---|---|---|
| Untransformed | p = 0.033 | p = 0.0014 | p = 0.036 | p = 0.59 | p = 0.043 |
| Arcsine square root (manuscript) | p = 0.057 | p = 0.0009 | p = 0.058 | p = 0.83 | p = 0.33 |
| Log | p = 0.12 | p = 0.0008 | p = 0.11 | p = 0.78 | p = 0.65 |

Residual normality is Shapiro-Wilk on model residuals; equal variance is Brown-Forsythe across the four groups. Three versions of this analysis exist in the repository (untransformed, Johnson, arcsine), and an older figure prints the untransformed p-values.

**Why it matters.** The heat effect on global 5-mC is robust. The ploidy effect and the interaction are not: they are significant only on the untransformed scale, which is the one that fails the equal-variance check, and are not significant on either scale that meets both model assumptions (arcsine 0.057 and 0.058; log 0.12 and 0.11). The arcsine analysis in the manuscript is therefore defensible, although it was selected by a test on the wrong quantity. The statement that diploids increased while triploids "showed no significant change" compares a significant with a non-significant result, which does not establish a difference (Gelman and Stern 2006); it needs the interaction, and the interaction is the fragile term. The Tukey letters are computed on the transformed scale and the reported means are not back-transformed.

**Proposed resolution.** Choose the scale from diagnostics of the fitted model, state the choice and the reason in the Methods, report the sensitivity to scale, and report effect sizes with confidence intervals on the original scale. Rewrite Discussion 4.5 so that it does not rest on the interaction.

**Needs from authors.** Agreement on the scale; see also I-12, which could change the values themselves.

---

### I-11 (S2) Figure 1A drops three of the eight diploid heat samples

**Evidence (verified).** `code/5-mC & 5mCpG plots.Rmd:60` sets `ylim(0.5, 5)`. In ggplot2, scale limits remove out-of-range observations before statistics are computed (see the ggplot2 reference for `lims()`; `coord_cartesian()` is the zoom that keeps the data). D14 (5.21), D15 (5.22), and D16 (6.96) are removed, so the diploid heat box is drawn from five animals: plotted median 2.99 against a true median of 3.98.

**Resolution.** Plot all points with an axis that contains them. Phase 3.

---

### I-12 (S2) ELISA calculation rests on unverified inputs

**Evidence.**
- *Verified.* Percent 5-mC is computed with a fixed 25 ng of DNA per well (`ELISA/GlobalDNAMeth_Polyploids.R:117-120`). The Methods say samples were diluted to 5 ng/uL and "4 mL" was loaded, which is 20 ng if uL is meant. A sheet named `Brconc.` in `ELISA/docs/MethylFlashELISA_StdCrv_and_Data_Analysis_Feb2019.xlsx` lists per-sample "DNA added" values of 10 to 22 ng. If those are the true inputs, recomputation gives heat p = 0.04 to 0.05 and no ploidy or interaction effect.
- *Verified.* The standard curve is fitted to the 0, 0.1, 0.2, 0.5, and 1% standards; the 5% and 10% standards are dropped as saturating. The Methods list 0.1, 0.2, 0.5, 1, 2, and 5%. Four diploid heat samples (D13 to D16) have absorbance above the highest standard used in the fit, so their values are extrapolated.
- *Verified.* Three wells are flagged `contam` with no explanation; two samples have a single well.
- *Verified.* The plate is one plate with the four groups in contiguous column strips, so plate position is aligned with group.

**Needs from authors.** The DNA mass loaded per sample well and per standard well; what the `Brconc.` sheet represents; why wells were flagged.

---

### I-13 (S2) Figure 1B values come from a hand-built sheet that mixes two alignment runs

**Evidence (verified).** `data/bismark-percent-mCPG-summary.xlsx` holds typed read counts copied from per-library Bismark alignment reports in `data/PE_reports/`. Nine of the ten heat reports are from the March 2021 alignment (Bismark 0.21.0); the tenth and all twelve control reports are from the June 2021 alignment (Bismark 0.22.3). The counts are from the alignment step, before deduplication. Four sample names are wrong (I-04) and the control libraries have no names. The differential methylation analysis itself used June 2021 coverage files for all 22 libraries (the project's coverage files on gannet have the same sizes as the June 2021 outputs, and the two that were checksummed, zr3534_7 and zr3644_14, are identical), so the mix affects only this figure.

**Resolution.** Recompute from the June 2021 summary for all 22 libraries (now in `07_documentation/upstream_processing/`). The values change in the second decimal. The confound in I-01 applies to this panel regardless.

---

## Numbers that do not trace to the data

### I-14 (S2) The "all CpG" genomic background was not computed from the analysed CpGs

**Where.** Results 3.3, first sentence; Figure 3A, first bar; the chi-squared test.

**Evidence (verified).** `analyses/DML_locations-all-datasets.xlsx`, sheet `data`: the column used as the 10x background is the cell-by-cell sum of two other columns, the control-only CpG set and the heat-only CpG set (for example exon UTR: 47,587 + 3,877 = 51,464). The heat-only counts come from a bedtools run that aborted on a malformed record (`code/3_WGBS_DML_feature_location-new.ipynb`, cell 22: "Error: Invalid record in file DML/ALL-CpG-Cov10-Ronit.bed"), so they cover only the start of the file. Neither column is the 480,579 CpGs that were tested.

**Resolution.** Compute the background from the analysed CpG set. The percentages in Results 3.3 and the enrichment and depletion statements will be regenerated.

### I-15 (S2) The 1x genomic distribution double counts exons

**Where.** Results 3.1, last paragraph ("exon UTR (27.3%)").

**Evidence (verified).** Same workbook: the 1x "exon UTR" count, 832,055, is the sum of the exon count (660,217) and the UTR count (171,838) from `analyses/DML_locations_genome_feature.xlsx`, so coding sequence is counted twice (once as CDS, once inside exon). The underlying CpG set totals about 1.9 million, which is neither the 12.8 million nor the 5.9 million 1x set; it matches the size of the heat-only 10x set.

**Resolution.** Recompute with one category definition at every coverage level.

### I-16 (S3) CpG counts

**Evidence (verified).** The 5x, 10x, and 20x counts (2,297,236; 480,579; 29,884) are CpGs covered in all 22 libraries. The corresponding 1x count in the code and in `analyses/coverage_summary.xlsx` is 5,852,939. The manuscript's 12,794,074 is on a different definition (CpGs covered in at least one library) and its source is not in the repository. Section 3.2 and Figure 3A say 480,479; the value is 480,579. The percentages of genomic CpGs (96.6%, 3.6%) need a stated denominator.

**Resolution.** Report both definitions, labelled; compute the genome total from the reference sequence.

### I-17 (S2) Chi-squared test of genomic location

**Evidence (verified).** `analyses/DML_locations-all-datasets.xlsx`, sheet `chi-square`: the expected counts and n = 1,299 computed for ploidy are reused for the heat and interaction tests; the interaction test reuses the heat test's contributions; df is 6 for eight categories; the categories overlap (I-18). The Methods describe the test; the Results report no statistic.

**Resolution.** Test association with a method that respects overlapping features and clustered loci, or drop the test from the Methods.

### I-18 (S3) Feature percentages are shares of overlap events

**Evidence (verified).** The eight feature tracks overlap (a transposable element inside an intron is counted in both). For ploidy-DML the eight counts sum to 1,299 for 1,038 loci. The percentages in Results 3.3 use the sum as the denominator.

**Resolution.** State the definition, or assign each locus to one category by a fixed priority.

### I-19 (S3) Figure labels

Figure 3A labels the bars "hypo (28)" and "hyper (49)"; the data have 28 loci positive in both contrasts and 49 negative in both, as the text says. Figure 2D's legend says "desiccation". Figure 3A says 480,479.

### I-20 (S2) Figure 5 does not plot what its panels name

**Evidence (verified).** The values are typed into `code/gene-plots.Rmd:61-66` and match `GOterms/summary_table.xlsx`, whose entries in turn match row counts noted in comments of the enrichment script. "GO terms" (for example 39,688 for heat BP) is the number of rows in a table of DML x transcript x term combinations, not a number of terms. "Unique transcripts" and "gene products" count transcripts annotated to the ten top terms. "Unique genes" (1, 3, 1 for BP) disagrees with Table 1, which lists six genes for ploidy BP.

**Options.** Replace with defined quantities computed by the pipeline (genes containing DML, by set), or drop the figure.

### I-21 (S2) "Enriched" means top ten

**Evidence (verified).** `code/5_WGBS_Functional-enrichement-analysis-BOTH.R:321-340` takes rows 1 to 10 of each topGO table. The Methods say p < 0.01. Between 36 and 856 terms per set and ontology have p < 0.01, so Figure 9, Table 1, and the percentages in Section 3.5 (for example "CC 49%, MF 31%, BP 20%") describe ten terms each. The "Category" column of Table 1 is hand-assigned.

**Resolution.** Tied to I-03.

### I-22 (S2) Figures 6 to 8 depend on the REVIGO web server

**Evidence (verified).** Only one REVIGO export exists (`figures/revigo_figure/Revigo-PLOIDY-BP.R`, without p-values); the heat and interaction figures exist as images only. The labelled terms are a hand-picked handful out of 341 to 856 plotted terms, and Results 3.4 narrates those labels.

**Resolution.** Regenerate with `rrvgo` from the pipeline's term tables, labelling cluster representatives by rule. Tied to I-03.

### I-23 (S3) Figure 2C to 2G

**Evidence (verified).** The PCAs in panels C and D use only the loci selected for differing between the groups shown, so separation is guaranteed and cannot support the statement that the DML "captured biologically meaningful variation" (Kriegeskorte et al. 2009, Nat Neurosci 12:535, doi:10.1038/nn.2303). `code/2_WGBS_Methylkit-all-gigas.R:823, 896, 945` look up each DML by start coordinate alone, without the chromosome, so a locus can be matched to a different chromosome's CpG with the same coordinate. The heat maps carry no sample annotation.

**Resolution.** Match on chromosome and position; annotate samples; reword the text.

### I-24 (S2) The CAD gene is misidentified

**Evidence (verified).** Table 1 lists gene 105319299, "CAD protein". Its best blastx hit in the project's own annotation is UniProt P27708 (human CAD: carbamoyl-phosphate synthetase 2, aspartate transcarbamylase, dihydroorotase), the pyrimidine synthesis enzyme, and it sits under pyrimidine metabolism GO terms. Discussion 4.3 and the Conclusions treat it as caspase-activated DNase, "a terminal effector of programmed cell death" (citing Nagata et al. 2003), which is a different protein (DFFB; in this genome LOC105341521, "DNA fragmentation factor subunit beta"). Table 1 assigns it the category "Apoptosis".

**Resolution.** Correct the text and the table. The apoptosis argument in 4.3, 4.6, and 4.7 does not follow from this locus.

### I-25 (S3) Chromosome distribution

The percentages in Results 3.3 (7: 16%, 3: 15.5%, 10: 12.9%, 6: 12.9%) have no code behind them; the chromosome figure is drawn from a hand-built workbook on gannet (`DML/ALL_gigas/DML-gene-distribution/DML-gene-distribution.xlsx`). Recomputed from the DML lists: 15.7%, 15.1%, 12.5%, 12.5% of all DML, with 203 of 3,589 on unplaced scaffolds. The legacy code keeps "the first ten" levels of a table, which is correct only while no DML falls on the mitochondrial genome. The four chromosomes hold 57% of DML, which the text calls "the majority" while also saying there is no clustering. "Tracked gene density" is asserted without a test. The "chromosomes" are RefSeq linkage groups NC_047559.1 to NC_047568.1.

### I-26 (S3) Figure S1

Panels B and C were built from histogram bar heights read off 44 per-sample plots and typed into `analyses/coverage_summary.xlsx`. The correlation range (0.87 to 0.93) and the PCA have no saved output. All will be recomputed.

### I-27 (S3) Supplementary lists

Section 7 cites lists of DML, GO terms, and genes that are not assembled. The pipeline will write them.

---

## Manuscript text

### I-28 (S2) Methods describe several steps differently from what was run

| Methods text | What the files show | Source |
|---|---|---|
| Reads trimmed with Cutadapt | File names carry `fastp-trim.20201202`; no trimming command or parameters are recorded in the repository | PE reports; `old_WGBS/WGBS/sbatch_scripts/` (a Bismark job script despite its name) |
| Bismark v0.22.3 | Correct for the June 2021 run used for DML; the March 2021 run (0.21.0) fed part of Figure 1B | `07_documentation/upstream_processing/` |
| Bowtie 2-2.3.4; SAMtools v1.10 | bowtie2 2.3.4.1; samtools 1.9 | job scripts |
| "alignment scores must exceed 90%" | `--score_min L,0,-0.6`, a minimum alignment score function, not a percentage | job scripts |
| Methylation calls extracted, then BAMs deduplicated | Deduplication first, then extraction | job scripts |
| CpGs in the 99.9th percentile of coverage removed | Not applied (I-05) | methylKit script |
| GO background: genes with CpGs at 10x in all samples | Transcripts with a CpG at 1x in any sample (I-03) | annotation notebook |
| Enriched terms: Fisher p < 0.01 | Top ten for Table 1 and Figure 9; p < 0.01 for Figures 6 to 8 (I-21) | enrichment script |
| Ploidy and treatment "as covariates" | Each factor tested with the other as covariate; no interaction term (I-02) | methylKit script |
| Chi-squared contingency test | Spreadsheet with errors; not reported (I-17) | workbook |
| One heat treatment | Four libraries also had a 45 C acute heat shock (I-04) | sample sheet |
| Controls from a mesocosm experiment "investigating ... gonad methylation" | Control libraries are ctenidia from Zymo project zr3644, a separate diploid/triploid experiment with ambient and low pH arms; details to be confirmed (I-36) | `metadata.xlsx` |
| ELISA: 5 ng/uL, "4 mL"; standards 0.1 to 5% | 25 ng assumed; standards 0 to 1% used (I-12) | ELISA script |
| BioProject PRJNAXXXX | zr3534 is PRJNA678408 in the sample sheet; zr3644 not recorded | `zr3534_wgbs_info.csv` |
| methylKit v1.17.4, topGO v2.34.0 | Not recorded anywhere (I-37) | none |

### I-29 (S2) qPCR and "mitochondrial abundance" claims

**Where.** Discussion 4.2, 4.5, 4.7, Conclusions. No Methods or Results.

**Evidence (verified by recomputation; details in `07_documentation/audit/elisa_qpcr_audit.md`).**
- *Mitochondrial abundance.* The repository has COX1 Cq values for 32 animals and nothing that expresses them against a nuclear reference. Group means differ by less than 0.3 cycles; the ploidy by treatment term has p = 0.070 and no pairwise comparison is significant. The data do not show an increase in diploids and a decrease in triploids.
- *HSP90.* About twofold higher in triploids, with p = 0.035 or p = 0.068 depending on how one failed sample is handled; the clearer effect is induction by heat in both ploidies with no interaction.
- Actin is the only reference gene, and its Cq ranges from 21.9 to 35.3. Non-detects are imputed as Cq 45 in some scripts. One script has unresolved merge-conflict markers.

**Options.** Add the qPCR to Methods and Results with an analysis that the data support, or remove the statements.

### I-30 (S3) Structure

- The abstract is a placeholder URL; keywords are empty.
- Section numbering jumps from 4.3 to 4.5. A section "3.5 Ploidy-Dependent Methylation Responses to Heat Stress" sits after the Conclusions. Section 3.5 contains two near-duplicate paragraphs on heat-DML genes.
- Sections 4.5 to 4.9 contain Markdown syntax (`##`, `**`) from pasted text.
- About 60 author-year citations in the text have no entry in the reference list: 24 are Paperpile citations whose bibliography was not refreshed, the rest are plain-text citations in Sections 4.3 to 4.9. Each of the latter should be checked against its source before it is kept; the descriptions of some (for example parental conditioning in geoduck) should be compared with what those papers report.
- Figure 9 is cited, but the figure files number two different figures as 6.
- Adenylate kinase 9 is discussed as an interaction locus but is not in Table 1. The gene (LOC105318231) does contain two of the overlap loci; Table 1 covers only genes under the ten top BP terms.
- The Discussion refers to "stress-resilient diploids and stress-susceptible triploids"; this study measured no survival or performance trait.
- Limitations (4.8) do not mention I-01, I-02, or I-04.

---

## Repository and reproducibility

### I-31 (R) No analysis script runs as committed

Absolute paths (`E:/...`, `/home/mngeorge/...`, `~/Documents/GitHub/...`); `setwd()` and `rstudioapi` throughout; scripts that depend on saved workspaces of up to 8.9 GB; objects used but never created (for example `differentialMethylationStatsINTERACTIONCov10`); notebook cells that were interrupted or failed; `gene_expression/scripts/dct_analysis.R` contains merge-conflict markers; `ELISA/GlobalDNAMeth_Polyploids.R` stops on undefined objects.

### I-32 (R) Non-R and manual steps

Feature overlap and GO annotation are in Python and bash notebooks that call bedtools. Seven workbooks sit between steps as hand-edited intermediates (`coverage_summary.xlsx`, `DML_difference_threshold_comparison.xlsx`, `DML_locations-all-datasets.xlsx`, `DML_locations_genome_feature.xlsx`, `summary_table.xlsx`, `DML-gene-distribution.xlsx`, `DML_summary.xlsx`), and two more hold the Figure 1 inputs. Figures are assembled, and statistics letters typed, in PowerPoint.

### I-33 (R) Inputs on a retired server

`owl.fish.washington.edu` does not respond. It held the exon and CDS tracks (also on OSF, doi:10.17605/OSF.IO/YGCTB, and now being brought into GitHub) and the UniProt-to-GO table `uniprot-SP-GO.sorted` used to build the gene-to-GO map, for which no copy has been found. The derived map (`geneid2go*.tab`) survives on gannet and will be used as the frozen annotation; the alternative is to rebuild from current UniProt, which would change the annotation.

### I-34 (R) Repository hygiene

613 MB; `old_WGBS/` duplicates code and data; `tables/` duplicates parts of `GOterms/`; qPCR analysis sits under `ELISA/`; the WGBS sample sheet and sequencing reports sit under `gene_expression/`; Office lock files and `.DS_Store` are committed; the README links to paths that no longer exist; no license; no lockfile.

**Needs from authors.** Whether regenerable legacy outputs above 10 MB may be dropped from the working tree (they remain in history and on gannet). The largest: `old_WGBS/WGBS/DML/ALL-CpG-Cov10.csv` (81 MB), the heat-only GO product tables under `GOterms/RONIT-*` (about 150 MB), two 26 MB TIFFs.

### I-35 (R) gannet directory

Of 217 GB, the pipeline needs the 22 coverage files (9.7 GB). About 117 GB are GO annotation joins of every covered CpG with every transcript and term (files of 19 to 33 GB), 43 GB are bedgraphs, 18 GB are saved workspaces. `bed_files/BOTH/` ends mid-listing (an incomplete upload). `bed_files/YAAMINI/` holds zr3644_1 to 12, which per `metadata.xlsx` are diploid control and diploid low-pH animals, not the 12 controls (1 to 6 and 13 to 18); if those bedgraphs fed the control-only background, it used the wrong animals. A stray 16 KB file named `name` sits at the top level. This session has read-only access; any cleanup is for the owner.

### I-36 (R) Sample metadata for the control libraries

No animal IDs, source experiment description, holding conditions, or accessions for zr3644 exist in the repository. `metadata.xlsx` labels them ctenidia and labels the other half of that project "OA".

**Needs from authors.** Experiment, location, dates, holding temperature, whether the animals were submerged, age and stock, relationship to the heat-treated animals, tissue handling, and SRA accessions.

### I-37 (R) Software versions

No script records `sessionInfo()`. The versions in the Methods (methylKit 1.17.4, topGO 2.34.0) cannot be confirmed; topGO 2.34.0 dates from 2018 and may have been carried over from an earlier paper. The audit environment is R 4.3.3 with Bioconductor 3.18 (methylKit 1.28.0, topGO 2.54.0); the CpG count reproduces exactly under it.

**Needs from authors.** The R version and operating system to target for the lockfile.

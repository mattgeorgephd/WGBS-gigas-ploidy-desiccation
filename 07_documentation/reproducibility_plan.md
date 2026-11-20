# Reproducibility plan

Project: *Epigenetic divergence in the thermal stress response of diploid and triploid Pacific oysters* (George, Jain, Roberts).
Branch: `reproducibility-audit`. Plan written 2026-11-18 after a full read of the manuscript draft, every script and notebook in this repository, and the project directory on gannet.

The companion document is [`issues.md`](../issues.md) at the repository root. This plan says what will be built and in what order; `issues.md` records every problem found, the evidence for it, and what has to be decided or fixed.

---

## 1. Goal

A pipeline, written only in R, that lives in this repository and regenerates every number, table, and figure cited in the manuscript, starting from the per-sample Bismark coverage files on gannet and the small raw data files kept in the repository. Each script reads the outputs of the script before it, so any step can be rerun on its own and any manuscript value can be traced to one script and one input.

A second goal follows from the audit: several manuscript values cannot be regenerated because they came from hand-edited spreadsheets or from steps that did not do what the Methods describe. For those, the pipeline computes the value correctly and reports the difference next to the manuscript value, so the manuscript can be revised against a known baseline.

## 2. Rules

These hold for everything added on this branch.

1. **R only.** No Python, shell, Jupyter, bedtools, or web tools in the pipeline. The legacy notebooks are translated: bedtools interval operations become `GenomicRanges` overlaps, shell joins become `data.table` joins, the REVIGO web step becomes `rrvgo`.
2. **GitHub first for data.** A file that the pipeline needs and that fits in GitHub lives in GitHub. The working ceiling is 50 MB per file, after gzip where the format allows it (GFF, BED, tab-delimited text). Only files above that ceiling stay on gannet. Every gannet-hosted input is listed in `03_input_files/gannet_manifest.csv` with its URL, size, and md5, and scripts verify the checksum after download. Nothing small exists only on gannet.
3. **One source of truth per value.** No hand-typed numbers in scripts and no hand-edited spreadsheets between steps. A value shown in a figure is computed by the script that draws the figure.
4. **Paths through `here::here()`.** No `setwd()`, no `rstudioapi`, no absolute paths. Scripts run from RStudio or from `Rscript` on any operating system.
5. **Each step declares its inputs and outputs** in a header block and writes into its own folder under `04_output/`. A step never writes outside its folder (exceptions: the download cache, and `05_manuscript/`, which step 13 assembles), and never reads a file it did not declare.
6. **Legacy behaviour and corrected behaviour are both available where they differ.** The switch is in `run_config.R`, the default is stated there, and both sets of outputs carry the mode in their file names. This is what lets the manuscript values be reproduced and corrected values be compared with them.
7. **Every reproduced value is checked.** `06_diagnostics/` compares regenerated outputs with the legacy outputs and with the values quoted in the manuscript, and writes a pass/differs table. A difference is never silently accepted.
8. **Nothing is deleted without being listed.** Legacy files move to `08_archive/` with their original relative paths. Files proposed for removal are listed in `issues.md` first.
9. **Naming.** The treatment is "heat" (or "heat-treated") everywhere in new code, outputs, and text. The legacy label "desiccation" survives only inside `08_archive/`. The legacy dataset labels map as follows: `ronit` = heat-treated animals only (Zymo project zr3534), `yaamini` = control animals only (Zymo project zr3644), `both` = the 22-library comparison reported in the manuscript.

## 3. Target layout

Modelled on `Coastal-Rec-Crab-BSS`: numbered stage folders in pipeline order, a README in each, a single configuration file and a single orchestrator at the root.

```text
WGBS-gigas-ploidy-desiccation/
├── 01_pipeline/          Numbered R scripts, one per analysis step, run in order
├── 02_R_functions/       Helper functions sourced by every step (paths, gannet fetch, plotting, overlaps)
├── 03_input_files/       Inputs kept in GitHub, plus the manifest of inputs fetched from gannet
│   ├── metadata/           Sample sheets (one row per library, one row per animal)
│   ├── elisa/              Raw plate reads, plate map, standard curve inputs
│   ├── qpcr/               Raw Cq tables
│   ├── bismark_reports/    Alignment and summary reports for all 22 libraries (one consistent run)
│   ├── genome_features/    Roslin v1 feature tracks, gzipped (gene, mRNA, exon, CDS, UTR, intron, ...)
│   ├── annotation/         blastx hits, transcript to gene map, GO mapping
│   ├── gannet_manifest.csv URL, size, md5 for every input that stays on gannet
│   └── gannet_cache/       Download cache (git-ignored)
├── 04_output/            One subfolder per pipeline step; regenerated by the pipeline
├── 05_manuscript/        Figures and tables as numbered in the manuscript, plus the value-by-value traceability table
├── 06_diagnostics/       Checks: regenerated vs legacy outputs vs manuscript values
├── 07_documentation/     This plan, audit notes, upstream (HPC) processing scripts, methods notes
├── 08_archive/           Everything from before the reorganization, original relative paths preserved
├── run_config.R          Single control surface: thresholds, sample sets, legacy/corrected switches, cache location
├── run_pipeline.R        Orchestrator: runs the steps in order, or a chosen subset
├── issues.md             Issue log
├── README.md
├── renv.lock, .Rprofile, renv/   Package versions (added once the package set is final)
└── WGBS-gigas-ploidy-desiccation.Rproj
```

### Where the existing files go

| Now | After | Note |
|---|---|---|
| `code/`, `analyses/`, `GOterms/`, `figures/`, `plots/`, `tables/`, `stats/`, `old_WGBS/`, `presentations/`, `ELISA/`, `gene_expression/`, `data/` | `08_archive/<same path>` | Moved with `git mv`, so history follows each file |
| `ELISA/docs/20190206_*ReadingSTRIGG.xls`, `ELISA/docs/20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx` | `03_input_files/elisa/raw/` (the workbook as `20190131_RonitsSamplesforGlobalDNAMeth.xlsx`, without its apostrophe) | Raw absorbance and plate map for the N = 32 run; converted to CSV by a script so the pipeline has no `.xls` dependency |
| `gene_expression/data/qpcr_ct_values/*.csv` | `03_input_files/qpcr/raw/` | Raw Cq values |
| `gene_expression/data/zr3534_wgbs_info.csv`, `metadata.xlsx` | `03_input_files/metadata/raw/zr3534_wgbs_info.csv`, `03_input_files/metadata/raw/legacy_metadata.xlsx` | Kept as raw inputs; `build_wgbs_libraries.R` merges them into `03_input_files/metadata/wgbs_libraries.csv` (treatment, ploidy, acute heat shock flag, source experiment, Zymo project, accession) |
| `old_WGBS/WGBS/GFF/*` (11 tracks) | `03_input_files/genome_features/*.gz` | Gzipped; the uncompressed copies leave the tree |
| gannet `GFF/`, OSF `genome-feature-files/` (CDS, exon) | `03_input_files/genome_features/*.gz` | CDS and exon tracks were only on the retired `owl` server and on OSF; brought into GitHub |
| gannet `data/20210605-cgigas-roslin-mito-blastx.outfmt6.txt`, `analyses/GOterm-annotation/` small files | `03_input_files/annotation/` | Brought into GitHub |
| gannet `/panopea/061021-big/` Bismark reports and job scripts | `03_input_files/bismark_reports/`, `07_documentation/upstream_processing/` | Brought into GitHub |
| gannet `cov_files/` (22 files, about 440 MB each) | stay on gannet | Listed in the manifest; fetched on demand |

Office lock files (`~$*.xlsx`, `~$*.pptx`) and `.DS_Store` are removed rather than archived.

## 4. Pipeline steps

Each row is one script in `01_pipeline/`. "Reproduces" names the manuscript items the step is responsible for.

| Step | Script | Reads | Writes (under `04_output/<step>/`) | Reproduces |
|---|---|---|---|---|
| 01 | `01_fetch_inputs.R` | `gannet_manifest.csv` | verified cache, `input_inventory.csv` | none (inputs) |
| 02 | `02_alignment_qc.R` | Bismark reports | per-library read, mapping, duplication, and methylation table; batch comparison | Methods 2.4; Fig. 1B inputs |
| 03 | `03_global_methylation.R` | ELISA plate reads, plate map, step 02 table | per-sample %5-mC from raw absorbance, model tables, post hoc tests | Results 3.1 (ELISA, %mCpG); Fig. 1 |
| 04 | `04_cpg_coverage.R` | 22 coverage files, reference genome, feature tracks | CpG counts by coverage threshold and sample set, the same CpG sets by genome feature track, CpG sites in the reference, methylation distributions, sample correlations, global PCA | Results 3.1; Fig. S1 |
| 05 | `05_differential_methylation.R` | 22 coverage files | united 10x table, per-CpG statistics, DML lists by threshold, BED files | Results 3.2; Fig. 2A, 2B |
| 06 | `06_dml_figures.R` | step 05 | PCA, heat maps, paired-difference plot | Fig. 2C to 2G; Fig. 4 |
| 07 | `07_genome_features.R` | step 05, feature tracks | feature overlap counts for the analysed CpG set and each DML set, association tests, per-chromosome counts | Results 3.3; Fig. 3 |
| 08 | `08_gene_annotation.R` | step 05, blastx, GO maps, feature tracks | transcript and gene annotation tables, DML to gene tables, per-gene counts of tested and eligible CpGs, the legacy annotation tables | Supplementary lists |
| 09 | `09_go_enrichment_legacy.R` | step 08, legacy universe, GO.db 3.13.0 | the legacy enrichment tables, Figure 5 and Figure 9 values, Table 1, and a record of what the legacy procedure did | Results 3.4, 3.5; Fig. 5, Fig. 9; Table 1 |
| 10 | `10_go_enrichment_gene_level.R` | steps 08 and 09 | gene-level enrichment per DML set, ontology, universe, and mode; each term of the draft retested | Replacement candidates for the same items |
| 11 | `11_semantic_similarity.R` | steps 09 and 10 | semantic similarity plots and cluster tables | Fig. 6, 7, 8 |
| 12 | `12_qpcr.R` | qPCR Cq tables | relative expression tables and tests | Discussion statements on HSP90 and COX1 (see issue I-29) |
| 13 | `13_manuscript_values.R` | all steps | `05_manuscript/` figures, tables, and `manuscript_values.csv` (manuscript value, regenerated value, status) | Every number in the text |

Design notes:

- **Memory.** The legacy analysis saved an 8.9 GB workspace. The new steps read one coverage file at a time and keep only filtered tables, so steps 04 and 05 run within 16 GB.
- **Step 05 is the expensive one** (about 10.5 GB of downloads the first time, then 18 minutes of model fitting on four cores, about an hour on one). Its outputs, including the united 10x methylation table (gzipped), are committed, so steps 06 to 13 can be rerun from the repository; step 09 also needs the GO.db 3.13.0 source package (31 MB), which step 01 or step 09 fetches from the Bioconductor archive into the cache.
- **Legacy and corrected modes** (rule 6). As built: the modes of `run_config.R` cover the upper coverage filter (issue I-05) and the lost first line of each coverage file (I-38), and steps 04 to 08 and 10 write both. The GO enrichment was split into two steps instead of two modes, because the legacy procedure differs from a sound one in too many ways to share a script: step 09 repeats it exactly (I-03, I-21, I-39, I-40, I-43, and the 2021 ontology release), and step 10 is the gene-level analysis. The formally tested interaction (I-02) was run as a diagnostic and is not a pipeline step: analysis A2 of `06_diagnostics/dml_sensitivity.R` lists 16,971 CpGs without the overdispersion correction and 2 with it (`06_diagnostics/dml_sensitivity/dml_counts.csv`; I-02, I-06).
- **The plan had twelve steps; the pipeline has thirteen**, because of that split. Steps 10 to 12 of the first version of this plan are now steps 11 to 13.

## 5. Phases and commits

Each phase ends in one or more pushed commits on `reproducibility-audit`.

| Phase | Content | State |
|---|---|---|
| 0 | Audit; `issues.md`; this plan; gannet manifest; upstream processing scripts | done (`d6a8112`, `503bf22`) |
| 1 | Reorganization into the numbered layout; README per folder; lock files removed | done (`117be05`) |
| 2 | `run_config.R`, `run_pipeline.R`, `02_R_functions/`, input manifest, inputs brought into GitHub under the 50 MB rule | done (`879969b`, `6a6a3c8`) |
| 3 | Steps 01 to 05, each validated against legacy outputs (CpG counts, DML lists row by row) | done; all 25 CpG counts and the three DML lists reproduced exactly |
| 4 | Steps 06 to 10, validated against legacy outputs | done; the feature overlap counts of the DML sets, the nine GO tables, Figure 5, Figure 9, and Table 1 reproduced exactly. The 18 genome-feature percentages of Results 3.3 follow from the counts; the chromosome percentages of the draft are not reproduced (I-25). Semantic similarity (Figures 6 to 8) moved to phase 6 |
| 5 | Corrected analyses, with side-by-side comparison tables | coverage filter and gene-level GO enrichment done; the interaction test, the overdispersion correction, the analyses without the acute heat-shock libraries and the within-cohort contrasts have been run as sensitivity analyses (`06_diagnostics/dml_sensitivity.R`, full run of 2026-11-19; results in `issues.md`, I-06, and `06_diagnostics/README.md`), and whether any of them becomes the reported analysis needs the decisions listed in section 7 (I-06 supports reporting no list of loci), as does variant masking; for step 10 (I-41) still to add: an exact test beside the Wallenius value, the logarithm of the CpG count as the bias variable, remapping of secondary GO IDs, the expectation adjusted for CpG count beside the unadjusted one, and in `06_diagnostics/` a label-permutation control, a spike-in positive control and the matched and permuted location comparisons |
| 6 | Steps 11 to 13; `05_manuscript/`; traceability table; list of manuscript edits implied by the results | steps 11 to 13 and `05_manuscript/` done; the traceability table holds 545 statements (352 match, 146 differ, 47 not regenerable); the list of edits is the set of rows that differ, with their notes and issues |
| 7 | `renv.lock`, final READMEs, clean checkout test from a fresh container | |

## 6. Validation targets

Legacy mode is accepted for a step when these are matched exactly:

| Quantity | Legacy value | Source |
|---|---|---|
| CpGs with at least 1x, 3x, 5x, 10x, 20x in all 22 libraries | 5,852,939; 3,714,898; 2,297,236; 480,579; 29,884 | `analyses/coverage_summary.xlsx`, sheet `total_cpg` |
| Ploidy DML at 20% difference, q < 0.01 | 1,038 (442 hyper, 596 hypo) | `tables/DML-getMethylDiff-PLOIDY-Cov10-20.csv` |
| Heat DML | 2,351 (1,150 hyper, 1,201 hypo) | `tables/DML-getMethylDiff-HEAT-Cov10-20.csv` |
| Loci in both lists | 200 (28 both hyper, 49 both hypo, 123 opposite sign) | `tables/DML-getMethylDiff-INTERACTION-Cov10-20.csv` |
| DML counts at 5% to 50% thresholds | table of 60 values (10 thresholds, 3 DML sets, hyper and hypo) | `analyses/DML_difference_threshold_comparison.xlsx` |
| Feature overlaps of each DML set (13 tracks) | for example ploidy: gene 981, CDS 371, intron 515, lncRNA 65, TE 194, intergenic 34 | `analyses/DML_locations-all-datasets.xlsx`, sheet `data` |
| topGO tables | 9 files, term by term | `GOterms/BOTH-*/*FisherTestResults.csv` |
| ELISA and %mCpG model tables | F and p as in the manuscript | `stats/*.csv` |

**State on 2026-11-19.** All of these are matched. `06_diagnostics/check_legacy_outputs.R` records 331 comparisons in `06_diagnostics/legacy_check_results.csv`; 329 match exactly. The two that differ are explained by the legacy files, not by the pipeline: the %mCpG interaction F (14.43 regenerated, 14.49 legacy; the legacy sheet mixed two alignment runs, issue I-13) and one threshold count where the legacy workbook has two typed values swapped (I-19). Beyond the targets above, the comparison also covers the values of Figure 5 and Figure 9 and the gene lists behind Table 1. Matching the topGO tables required four things, all part of step 09: the ontology release of 2021 (GO.db 3.13.0), the legacy universe file with its gap (I-39), a replica of the file-reading fault of the legacy script (I-40), and the universe names taken with their repeats (I-43).

## 7. Decisions needed from the authors

These do not block phases 1 to 4. They decide what phase 5 builds and what the manuscript can claim. Details and evidence are in `issues.md` under the issue numbers given.

1. **Control and heat libraries come from different experiments and sequencing batches (I-01).** Do tissue or DNA from the heat experiment's own control animals (D01 to D08, T01 to T08) still exist, and is sequencing them possible? If not, how should the heat contrast be framed? The sensitivity analyses make this the decisive question: inside the heat experiment diploids and triploids differ in overall methylation level and at more CpGs than any regrouping of its libraries gives; inside the control experiment the number of such CpGs is what regroupings give and the level differs in the opposite direction; only the heat experiment's own controls can show whether heat has anything to do with it (I-06, I-49). Because libraries with a raised level occur in 2 of the 12 untreated libraries here, that needs more than a few animals per group (I-49).
2. **Four of the ten heat libraries also received a 45 C acute heat shock (I-04).** Keep them pooled and disclose it, drop them, or model the heat shock as a covariate? Measured since: the four are not distinguishable from the other six heat libraries (55 of 100 regroupings give at least as many CpGs above 20 points), and dropping them changes 43% of the ploidy list and 24% of the heat list without making either more reliable. The evidence favours disclosing and keeping them. **Answered by the author on 2026-11-19:** all heat-treated animals count as one treatment group, so the four libraries stay pooled, as in the draft; the 24 h emersion that Methods 2 describe is not confirmed (I-04).
3. **"Interaction" loci (I-02).** Keep the overlap set under a different name, fit a formal ploidy by heat model, or both? The formal model has been run: without an overdispersion correction it calls 16,971 CpGs, of which 99 are among the draft's 200; with the correction it calls 2. The interaction it detects is a difference between cohorts in the ploidy contrast, most of which sits in five libraries (I-49); a smaller part remains without them and is the difference between two diploid libraries (D11, D12) and the triploids of the heat cohort (I-02). A renamed overlap set would still be built from the ploidy list, which I-06 finds cannot be reported as loci that differ by ploidy. Is the set to be dropped?
4. **GO enrichment (I-03, I-39, I-40, I-41).** The gene-level analysis is built (step 10). Compared with the genes that hold an eligible CpG (mean methylation 10% to 90%), and adjusted for the number of such CpGs per gene, no specific GO category is established in any DML set: at the configured settings no term passes; under exact tests the near-universal term "organelle" passes for the heat set, and one group of related terms (clathrin and vesicle coat, heat set) sits at the threshold and is unresolved. The draft's enrichment results also rest on an incomplete gene universe and on an input table that was mostly skipped when read. Two ways forward: report that no GO category is established, with the power limits and the one unresolved group named (I-41, option 2), and replace Figures 5 to 9 and Table 1 with a descriptive table of genes that hold DML; or design a different functional analysis (new work, for example a region-level statistic within one cohort; the sensitivity analyses limit what that can show, because among the controls the ploidy contrast gives no more CpGs above 20 points than regroupings do, I-07, and the difference among the heat-treated libraries cannot be attributed to ploidy, I-49). Which?
5. **Upper coverage filter (I-05).** Apply the 99.9th percentile filter the Methods describe, or keep the analysis as run (480,579 CpGs) and change the Methods? If the filter is applied, the Methods must also say over which CpGs the percentile is taken: 454,555 CpGs remain when it is taken per library over CpGs with at least 2 reads (the import threshold of the legacy DML script, used by the pipeline), and 450,065 to 469,453 under the three other readings listed in I-05. The DML counts differ by at most 7 (ploidy) and 19 (heat) between these readings.
6. **ELISA (I-10, I-12).** What DNA mass was loaded per well, and which scale of analysis should be reported?
7. **qPCR (I-29, I-46, I-47).** The pipeline step exists (step 12). The mitochondrial abundance statements are not supported by the data; the HSP90 statements are partly supported and depend on actin, the only reference gene, which differs between groups. Remove the qPCR statements from the Discussion, or keep HSP90 as an exploratory observation with Methods, fold changes with intervals, and its limits stated?
8. **Provenance of the control animals (I-01, I-36).** Source experiment, location, holding conditions, age, and relationship to the heat-treated cohort; animal IDs and SRA accessions for zr3644.
9. **Computing environment (I-37).** Install comments, workspace headers and the ontology release of the archived tables point to R 4.1.0 on Windows with Bioconductor 3.13 for the legacy enrichment; can the authors confirm this? Is there a record of the methylKit version (candidates 1.18.0 and 1.20.0, less likely 1.16.1)? What were the three topGO tables whose names end in `-new` made for (GO.db 3.16.0, not on gannet)? Which REVIGO version, database and settings were used for Figures 6 to 8? Which R version and operating system should the lockfile target?
10. **Bulky legacy files (I-34).** May regenerable legacy outputs above 10 MB be dropped from the working tree (they stay in git history and on gannet)?
11. **Location of DML (I-41).** Report where DML lie descriptively and drop the enrichment statements? Against all analysed CpGs the gene-body odds ratio is 6 to 7; against CpGs with a mean methylation between 10% and 90% it is 2.2 to 2.5; against CpGs matched for methylation level and between-animal variation it is 1.0 to 1.1, and CpG lists made with shuffled ploidy labels show the same odds ratio as the real DML (median 2.65 over all 232,848 labellings of the ploidy labels within cohort; `06_diagnostics/feature_odds_label_permutation.csv`). No location preference attributable to ploidy or heat is shown. A small excess in coding sequence remains after matching (odds ratio 1.2), but CpG lists made with shuffled ploidy labels show the same excess (median 1.58 against eligible CpGs, 1.17 after the same matching), so it comes from the selection rule and not from ploidy; for heat no such control is possible. (An earlier version of this item proposed the 10% to 90% background as the fair comparison; verification showed it is not.)
12. **Functional annotation (I-42).** The GO terms are a July 2017 UniProt snapshot. Keep it frozen, or refresh it from current UniProt (one download, then frozen in the repository)? Only matters if a functional analysis stays in the paper.
13. **Pipeline outputs in git.** Step 05 commits six tables of 14 to 22 MB each (about 100 MB in total) so that later steps run without the coverage files. This follows the GitHub-first rule but makes the repository heavier with every rerun that changes them. Keep them in git, or move them to gannet and fetch them like the coverage files?
14. **Semantic similarity threshold (I-22).** Step 11 uses rrvgo's default threshold of 0.7; REVIGO's "allowed similarity" of 0.7 corresponds to 0.3 in rrvgo. Keep 0.7 or use 0.3? One value in `run_config.R`. Only matters if a semantic similarity figure stays in the paper, which I-48 argues against.
15. **Differential methylation model (I-06).** The draft's model has no overdispersion term, and its test statistic is inflated nearly fourfold (about threefold by variation among libraries, the rest by the scaling of the counts; `06_diagnostics/dml_sensitivity/dispersion_summary.csv`). Under it the ploidy list (1,038 CpGs) cannot be told apart in length from lists of relabelled libraries. methylKit's corrected test leaves 2 ploidy CpGs and 16 heat CpGs (6 with the coverage filter), but it applies no correction where its dispersion estimate is not above 1, and relabelled libraries get such calls too. Report no list of loci and say that these data cannot name ploidy or heat loci; keep the uncorrected lists as description only, with the relabelling reference beside each ploidy count; or commission a different method (beta-binomial or region-level), checked against relabelled libraries, knowing that it removes neither the confound of I-01 nor the structure of I-49? The benchmark of Wreczycka et al. 2017 has since been read in full (I-06, correction paragraph): it rates the default test lowest in specificity and best in overall accuracy, in comparisons where half of the CpGs truly differ. The paper on which the draft's Methods 2.5 are modelled, Venkataraman et al. 2022, applied the correction with a chi-squared test, a 50% threshold and the removal of CpGs at C-to-T variants, and the draft's sentence on the model omits the correction: was it left out on purpose? **Answered by the author on 2026-11-19:** the correction was left out on purpose. The pipeline keeps the draft's model; the options named here stay as notes, and none is implemented without the author's approval (I-06).
16. **Animals and library batches (I-49).** What were the source and the crosses of the diploid and the triploid animals of each experiment? Were the triploids of the heat experiment one family, and did its diploids come from the same stock as its triploids and as the animals of the control experiment? In what order and batches were the libraries prepared, and were the libraries numbered 1 and 2 of each project (zr3534_1, zr3534_2, zr3644_1, zr3644_2) handled differently from the rest? Is anything recorded about animals D11, D12, D13, D19 and D20 and about the animals behind zr3644_14 and zr3644_15 (size, sex, gonad stage, condition, position in the sampling order, the time between the acute heat shock and sampling for D19 and D20, amount and quality of the DNA)? Five of these libraries (D13, D19, D20, zr3644_14 and zr3644_15) carry most of the ploidy by heat interaction of Figure 1B. What was the bisulfite conversion efficiency of each library (the E. coli amplicons the Methods mention)? Were zr3534_9 and zr3534_10 converted or prepared separately from zr3534_1 to 8? How were the triploids of each experiment produced (mating of tetraploid males with diploid females, or induction), and how many sires and dams contributed to each group? Is the sex of the 22 animals recorded? For both projects: extraction method, library kit, conversion protocol and number of PCR cycles; was the ctenidium tissue fixed or embedded (Methods 2.2 use the wording of a paper that extracted DNA from histology blocks; I-28)? Was the 24 h emersion applied to the heat-treated animals in one container and one run (I-01)? Does tissue or DNA of the control animals of the heat experiment remain (decision 1)? Methods 2.2 says one lane of a single flowcell, and the read names show two runs on two flowcells, one per experiment (`06_diagnostics/fastq_read_names.csv`; I-28): is that as the authors remember it? **Answered by the author on 2026-11-19:** how the triploids were produced and from how many parents is not known at present (to be added); sex was not recorded; the ctenidia were flash frozen and stored at -80 C before extraction, not fixed or embedded; the 24 h heat emersion is not confirmed. The other questions of this item are open.
17. **Two values of Results 3.1 that could not be traced (I-16).** Which file or command gave the 12,794,074 CpG loci at 1x, and for which libraries? Which libraries, coverage threshold and script gave the Pearson range 0.87 to 0.93?

### Decided

| Date | Decision | By |
|---|---|---|
| 2026-11-19 | Keep `07_documentation/agent_workflows/` (workflow scripts, the agents' raw reports, the verifying agents' own scripts) in the repository for now. To be revisited before the branch is merged. | M. George |
| 2026-11-19 | **Scope.** The goal is a reproducible pipeline that produces the results of the paper and notes where improvements would help. No improvement is implemented without the author's approval. Consequences drawn by the session: no further sensitivity analysis, diagnostic or alternative analysis is started; the remaining code work is limited to defects that affect running the pipeline or reproducing the paper's results, and to documentation; what was already built beyond the paper's analysis (corrected mode, the gene-level enrichment of step 10, the scripts of `06_diagnostics/`) stays as it is, as evidence behind the notes of `issues.md`, until the author says what is to become of it. | M. George |
| 2026-11-20 | **Coverage: the pipeline is to run with 5x or with 10x as the standard, chosen by a toggle.** 10x is the draft's requirement; 5x is that of Venkataraman et al. 2022. To be built in the toggles pass. | M. George |
| 2026-11-20 | **Answers on the variant exclusion and the old records.** (1) The exclusion of loci at C-to-T variants is a toggle, off by default, so that the default run still reproduces the draft. (2) Variants are called by R code only, applying the rule of BS-Snper to the same sorted, deduplicated alignments; BS-Snper itself is not run. (3) The records of earlier agent work under `07_documentation/agent_workflows/` that mention the F test are left as they are (historical records). | M. George |
| 2026-11-20 | **Loci at C-to-T variants are to be excluded, as in Venkataraman et al. 2022.** That study called variants from the sorted, deduplicated alignments with BS-Snper (minimum coverage 5x, other settings default) and removed C-to-T variant loci from the CpGs it characterised and from its candidate DML. This project has no variant calls; the alignments of the 2021 run are the source (`07_documentation/upstream_processing/`). Open, to be put to the author: whether the exclusion is the standard of every run, which changes the DML counts of the draft, or an option; and whether variants are called with BS-Snper outside R or by R code that applies the same principle. | M. George |
| 2026-11-20 | **Remove every mention of the F test** of the DML model (the author's reply to the question whether the counts of that test should stay in `06_diagnostics/dml_sensitivity/`). The sensitivity analysis keeps one model with the correction, the chi-squared one. Read by the working session as covering the F test of methylKit's DML model, not the F statistics of the draft's own analyses of variance (ELISA, qPCR), which are results of the paper. | M. George |
| 2026-11-20 | **The repeat track keeps its content and its numbers and is called "repetitive elements (RepeatMasker, all classes)"** in Methods 2.3, Results 3.1 and 3.3 and Figure 3A, with one sentence on its composition (option 1 of the answer of 2026-11-19). Nothing is recomputed. | M. George |
| 2026-11-20 | **Wording of the difference threshold.** The author asked why the session wrote "above 20 points" where the draft says "at least 20%". Answer given: the same threshold; methylKit's `getMethylDiff` keeps loci with an absolute difference strictly above the value (`abs(meth.diff) > difference`), and the difference is in percentage points of methylation. With "at least", the heat list of the draft would hold 2,352 loci instead of 2,351. The documentation is to say "a difference of more than 20 percentage points (the draft: at least 20%)" where precision matters. | M. George (question) |
| 2026-11-20 | **The option for the overdispersion correction offers the chi-squared test only**, as in Venkataraman et al. 2022 (`overdispersion = "MN"`, `test = "Chisq"`; BMC Genomics 23:556, Methods, "an overdispersion correction with a chi-squared test"). The F test is not coded as an option of the pipeline. Read by the working session as a decision about the test alone: the thresholds of the draft (q < 0.01, difference above 20 points), its coverage requirement and its covariates stay as they are when the option is on; the other settings in which the 2022 paper differs from the draft (5x coverage, a difference of 50%, exclusion of loci at C-to-T variants) are not part of the option. To be confirmed by the author if meant otherwise. | M. George |
| 2026-11-19 | **What was built beyond the paper's analysis stays in place and is off by default.** The default run produces only the paper's results; the rest runs on request and is labelled as not adopted, so that the evidence behind the notes stays reproducible. The author asks for toggles to turn each part on, and is particularly interested in the overdispersion correction of the DML model. To be built after the hardening pass: toggles in `run_config.R` and on the command line of `run_pipeline.R` for the coverage filter of the Methods (now "corrected mode"), for the DML model with methylKit's overdispersion correction, for the gene-level GO enrichment (step 10) and for the diagnostics. | M. George |
| 2026-11-19 | The overdispersion correction of methylKit was left out of the DML model on purpose. The pipeline keeps the draft's model (I-06). | M. George |
| 2026-11-19 | All heat-treated animals count as one treatment group (heated). The four libraries that also received the acute heat shock stay pooled with the other six, as in the draft (I-04; the session's reading of the answer). The 24 h heat emersion that Methods 2 describe is not confirmed. | M. George |
| 2026-11-19 | Information: sex of the 22 animals was not recorded. The ctenidia were flash frozen and stored at -80 C before DNA extraction, not fixed or embedded. How the triploids were produced, and from how many parents, is not known at present and will be added. | M. George |
| 2026-11-19 | The text of the manuscript draft is not committed. The repository is public and the draft is unpublished. `03_input_files/manuscript/README.md` records which copy of the Google Doc the line numbers of the values map refer to (date, size, checksum, line ranges), and each row carries the verbatim clause it concerns. | M. George |

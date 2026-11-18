# 01_pipeline

The analysis, as numbered R scripts that run in order. Each script reads the inputs in `03_input_files/` and the outputs of earlier steps in `04_output/`, and writes only into its own folder `04_output/<step>/`.

For the project overview see the [root README](../README.md); for the design and the order of work see [`07_documentation/reproducibility_plan.md`](../07_documentation/reproducibility_plan.md).

## Steps

| Step | Script | What it does | Manuscript items | Status |
|---|---|---|---|---|
| 01 | `01_fetch_inputs.R` | Downloads the files listed in `03_input_files/gannet_manifest.csv` into the cache and verifies their checksums | none | done |
| 02 | `02_alignment_qc.R` | Per-library read, mapping, duplication, and methylation table from the Bismark reports; comparison of the two sequencing batches | Methods 2.4; Figure 1B inputs | done |
| 03 | `03_global_methylation.R` | %5-mC from raw ELISA absorbance; models and post hoc tests for %5-mC and %mCpG | Results 3.1; Figure 1 | done |
| 04 | `04_cpg_coverage.R` | CpG counts by coverage threshold and sample set; methylation distributions; sample correlations; whole-data PCA | Results 3.1; Figure S1 | planned |
| 05 | `05_differential_methylation.R` | Coverage filtering, normalization, union across libraries, per-CpG tests for ploidy and heat, DML lists by difference threshold | Results 3.2; Figure 2A, 2B | planned |
| 06 | `06_dml_figures.R` | PCA and heat maps of DML; paired-difference plot of loci in both lists | Figure 2C to 2G; Figure 4 | planned |
| 07 | `07_genome_features.R` | Overlap of the analysed CpGs and each DML set with the genome feature tracks; association tests; per-chromosome counts | Results 3.3; Figure 3 | planned |
| 08 | `08_gene_annotation.R` | Transcript, gene, and GO annotation tables; DML-to-gene tables; the tested-gene background | Supplementary lists | planned |
| 09 | `09_go_enrichment.R` | GO enrichment per DML set and ontology; DML within enriched genes | Results 3.4, 3.5; Figures 5 and 9; Table 1 | planned |
| 10 | `10_semantic_similarity.R` | Semantic similarity clustering and plots of enriched terms | Figures 6 to 8 | planned |
| 11 | `11_qpcr.R` | Relative expression from raw Cq values | Discussion (see issue I-29) | planned |
| 12 | `12_manuscript_values.R` | Collects every number quoted in the manuscript, compares it with the regenerated value, and assembles `05_manuscript/` | all | planned |

The status column is updated as each script is added and validated.

## How a step is written

1. A header block states the purpose, the inputs with their locations, the outputs, and the manuscript items the step is responsible for.
2. The step sources `run_config.R` and every file in `02_R_functions/`, and loads its packages through one helper, so a missing package is reported by name before any work starts.
3. All paths are built with `here::here()`. There is no `setwd()`, no `rstudioapi`, and no absolute path. A script runs the same from RStudio (Source) and from a terminal (`Rscript 01_pipeline/05_differential_methylation.R`).
4. Sample groups come from `03_input_files/metadata/wgbs_libraries.csv`, never from positions in a list.
5. No number is typed into a script and no spreadsheet is edited between steps.
6. Where the legacy analysis and a corrected analysis differ, the script produces the one selected in `run_config.R`, and the mode is part of the output file names.
7. The step ends by writing its session information to `04_output/<step>/session_info.txt`.

## Running

From the repository root:

```sh
Rscript run_pipeline.R              # all steps in order
Rscript run_pipeline.R --steps 5:7  # a subset
```

or source an individual script in RStudio with the project (`WGBS-gigas-ploidy-desiccation.Rproj`) open.

Steps 04 and 05 need the 22 coverage files (9.1 GB, downloaded once by step 01) and about 12 GB of memory. Every later step runs from files that are in the repository.

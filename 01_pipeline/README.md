# 01_pipeline

The analysis, as numbered R scripts that run in order. Each script reads the inputs in `03_input_files/` and the outputs of earlier steps in `04_output/`, and writes only into its own folder `04_output/<step>/`, with two exceptions: step 01 (and step 09, if step 01 was skipped) fills the download cache, and step 13 writes `05_manuscript/`.

For the project overview see the [root README](../README.md); for the design and the order of work see [`07_documentation/reproducibility_plan.md`](../07_documentation/reproducibility_plan.md).

## Steps

| Step | Script | What it does | Manuscript items | Status |
|---|---|---|---|---|
| 01 | `01_fetch_inputs.R` | Downloads the files listed in `03_input_files/gannet_manifest.csv` into the cache and verifies their checksums; installs the ontology release of the legacy enrichment into the cache; writes the size and md5 of every input to `input_inventory.csv` | none | done |
| 02 | `02_alignment_qc.R` | Per-library read, mapping, duplication, and methylation table from the Bismark reports; comparison of the two sequencing batches | Methods 2.4; Figure 1B inputs | done |
| 03 | `03_global_methylation.R` | %5-mC from raw ELISA absorbance; models and post hoc tests for %5-mC and %mCpG | Results 3.1; Figure 1 | done |
| 04 | `04_cpg_coverage.R` | CpG counts by coverage threshold and sample set; the overlap of each of those CpG sets with the genome feature tracks (`cpg_counts_by_feature_track.csv`), with the percentages computed the two ways the draft computed them (issue I-15); CpG sites in the reference; methylation distributions; library correlations; whole-data PCA | Results 3.1; Figure S1 | done |
| 05 | `05_differential_methylation.R` | Coverage filtering, normalization, restriction to the CpGs covered in all libraries, per-CpG tests for ploidy and heat, DML lists by difference threshold, loci in both lists. DML are counted as methylKit counts them, with `abs(meth.diff) > threshold` on the unrounded difference; the tables carry 15 significant digits, so a count retaken from a written table differs at three loci of the heat contrast whose difference is exactly 15 or 20 points in the counts and is held one rounding error away from it: NC_047559.1:17861107 and NC_047562.1:49679999 (written -15, held as -15.000000000000002; counted at 15 points, which gives the 5,910 loci of the legacy script and not 5,908) and NW_022994958.1:20934 (written -20, held as -19.999999999999996; not a DML at 20 points) | Results 3.2; Figure 2A, 2B | done |
| 06 | `06_dml_figures.R` | Draws Figure 2 (DML counts and percent hypermethylated by threshold, from the table of step 05; PCA and heat maps of the DML) and Figure 4 (both differences of each locus in both lists). A panel whose DML set is too small to draw says so, and the step still finishes | Figure 2; Figure 4 | done |
| 07 | `07_genome_features.R` | Overlap of the analysed CpGs and each DML set with the 13 genome feature tracks; one exclusive class per CpG; Fisher tests per track against all analysed CpGs and against eligible CpGs, and a chi-squared test on the classes against all analysed CpGs; counts per linkage group; DML grouped into regions. The track called TE is every RepeatMasker repeat, read as the legacy analysis read it (`03_input_files/genome_features/README.md`). A set without loci gets NA in place of a test | Results 3.3; Figure 3 | done |
| 08 | `08_gene_annotation.R` | Transcript, gene, and GO tables; for each configured mode every (DML, gene) overlap and per-gene counts of tested CpGs, eligible CpGs, and DML, with a summary row per mode and DML set; in legacy mode the legacy annotation tables, regenerated line for line | Supplementary lists; inputs of steps 09 and 10 | done |
| 09 | `09_go_enrichment_legacy.R` | The GO enrichment exactly as run for the draft (transcripts as genes, legacy universe, legacy reader, GO.db 3.13.0), with a record of what that procedure did. Legacy mode only: when `legacy` is not a configured mode the step says so and ends without reading or writing anything | Results 3.4, 3.5; Figures 5 and 9; Table 1 | done |
| 10 | `10_go_enrichment_gene_level.R` | Gene-level GO enrichment against tested genes and against eligible genes; Fisher and Wallenius tests; false discovery rate; each term of the draft retested. Known defect, left as it is by the scope decision of 2026-11-19: GO IDs that are secondary IDs of current terms are dropped, and with them replaced by their primary IDs the step calls one term (organelle; heat set, CC) where it now calls none (`issues.md`, I-41) | Replacement candidates for the same items | done |
| 11 | `11_semantic_similarity.R` | Semantic similarity of GO terms in R (rrvgo, Rel measure) in place of the REVIGO web server: similarity matrix, clusters and a deterministic layout for the nine legacy term lists (p < 0.01), for the three lists the draft's figures were in fact drawn from (every BP term annotated to a gene with a DML), and for the gene-level terms (none at present); a check of how far the clusters depend on the order of the terms; each surviving REVIGO export matched to its term list and redrawn; every GO term named in Results 3.4 traced. Legacy mode only: when `legacy` is not a configured mode the step says so and ends without reading or writing anything | Figures 6 to 8; Results 3.4; Methods 2.6, last sentence (issues I-22, I-48) | done |
| 12 | `12_qpcr.R` | Per-animal Cq from the replicate-level qPCR files, checked against the hand-built consolidated file; the legacy delta Cq ANOVA and Tukey tests under both treatments of non-detects, checked against the knitted legacy output; delta Cq models by complete cases and by censored regression, with fold changes, confidence intervals, q-values and leave-one-out ranges; a simulation of how often the censored test rejects a true null; checks of the reference gene and of the normalization; raw Cq of COX1; every qPCR statement of the draft with its numbers and a verdict | Discussion 4.2, 4.5, 4.7; Conclusions (issues I-29, I-46, I-47) | done |
| 13 | `13_manuscript_values.R` | Evaluates the map of the manuscript's quantitative statements (545 rows), compares each with the value that legacy mode regenerates and, when `corrected` is configured, with the value of corrected mode, and assembles `05_manuscript/` (figures, Table 1, supplementary lists). When `legacy` is not a configured mode the step says so and ends without reading or writing anything | all | done |

The status column is updated as each script is added and validated. "Done" means the script runs from a clean R process and, in legacy mode, its outputs agree with the legacy outputs wherever a legacy output exists (`06_diagnostics/legacy_check_results.csv`).

## Modes

`modes` in `run_config.R` is `"legacy"`, `"corrected"` or both (the setting in the repository); any other value stops every script with a message. Steps 04 to 08 and step 10 write one set of outputs for each mode listed, with the mode in the file name. Every step prints the configured modes under its banner, and `run_pipeline.R` prints them at the start and at the end of a run.

| Mode | Coverage filter | Coverage files |
|---|---|---|
| `legacy` | No upper filter: the legacy call misnamed the argument, so the filter never ran (issue I-05) | First line of each file lost, as in the legacy import (issue I-38) |
| `corrected` | CpGs at or above the 99.9th percentile of coverage in a library are removed from that library, as the Methods describe | Every line read |

Legacy mode reproduces the numbers of the draft. Step 09 exists only for the legacy mode; it repeats the legacy GO procedure, defects included, so that the GO numbers of the draft can be traced. When `legacy` is not among the modes of `run_config.R`, step 08 writes no legacy annotation table and no legacy row in `annotation_summary.csv` (tables of an earlier legacy run stay in the folder as they were), and `run_pipeline.R` does not start step 09; started alone, step 09 says so, reads and writes nothing, and ends with exit status 0. Step 09 stops with an error when `legacy` is configured and the summary of step 08 has no legacy rows, that is when step 08 last ran without the legacy mode. Step 10 applies one procedure to the DML of both modes.

What a run does with each setting:

| `modes` | Steps 04 to 08 | Steps 09 and 11 | Step 10 | Step 13 |
|---|---|---|---|---|
| `c("legacy", "corrected")` | both sets of outputs | run | both modes | both columns of `manuscript_values.csv`; all items copied |
| `"legacy"` | the legacy outputs | run | legacy mode alone; the tables that span the modes hold legacy rows only | the legacy values; the corrected columns are empty, and the items of corrected mode are not copied |
| `"corrected"` | the corrected outputs | not started by `run_pipeline.R`; started alone, they say so and end with status 0 | not started by `run_pipeline.R` | not started by `run_pipeline.R`; started alone, it says so and ends with status 0 |

Steps 01 to 03 and 12 do not depend on the mode. The single-mode runs were tested on steps 02, 03 and 06 to 13, with the outputs of steps 04 and 05 taken from a run with both modes; steps 04 and 05 themselves were not run with one mode. Reading the code, step 05 loops over the configured modes without naming one, and step 04 names legacy mode in one statement only, a console print of the legacy shares after its two count tables, which runs only when legacy mode is configured. Steps 09, 11 and 13, started alone, end in the same way; they call one helper for it (`require_mode()` in `02_R_functions/logging.R`). Step 10 has no check of its own: started alone without the legacy mode it finishes and gives every term of the draft the status "no eligible gene carries it", which is wrong, so `run_pipeline.R` leaves it out; a check inside the step is one of the changes that wait for the author's decision. Outputs of a mode that is not configured are neither written nor removed: files of an earlier run stay in `04_output/` and `05_manuscript/` as they were.

## Run times

Measured on 4 cores and 15 GB of memory, with the downloads already in the cache.

| Step | Minutes | Notes |
|---|---|---|
| 01 | 1 to 2 | First run downloads 10.5 GB (plus the 31 MB GO.db archive) |
| 02, 03 | under 1 | |
| 04 | 12 | About 9.5 GB of memory (peak resident size) |
| 05 | 18 | `WGBS_CORES=4`; about 10 GB of memory. `WGBS_CORES=2`: 37 to 45 minutes and about 7.5 GB (the R process and its two workers, shared pages counted once). Single core: about an hour. About half a minute when both fits come from the cache |
| 06, 07, 08 | under 1 each | |
| 09 | 3 | |
| 10 | 3 | |
| 11 | 5 | About 1.6 GB of memory |
| 12 | under 1 | About 35 seconds, 24 of them the simulation |
| 13 | under 1 | About 10 seconds for 545 statements and 67 copied files |

## How a step is written

1. A header block states the purpose, the inputs with their locations, the outputs, and the manuscript items the step is responsible for.
2. The step sources `run_config.R` and every file in `02_R_functions/`, and loads its packages through one helper, so a missing package is reported by name before any work starts. The helper also asks for R.utils whenever a step names data.table, because `fread()` reads the `.gz` tables through it. `run_config.R` checks its settings as it is sourced, so a setting that cannot work stops the step at this point too.
3. All paths are built with `here::here()`. There is no `setwd()`, no `rstudioapi`, and no absolute path. A script runs the same from RStudio (Source) and from a terminal (`Rscript 01_pipeline/05_differential_methylation.R`).
4. Sample groups come from `03_input_files/metadata/wgbs_libraries.csv`, never from positions in a list.
5. No number is typed into a script and no spreadsheet is edited between steps. Settings are in `run_config.R`, with one list for a step that has many (`semantic_similarity` for step 11, `qpcr` for step 12).
6. Where the legacy analysis and a corrected analysis differ, the script produces both (the modes listed in `run_config.R`), and the mode is part of the output file names.
7. The step ends by writing its session information to `04_output/<step>/session_info.txt`.

## Running

```sh
Rscript run_pipeline.R                   # all steps in order
Rscript run_pipeline.R --steps 5:7       # steps 05, 06 and 07
Rscript run_pipeline.R --steps 2:3,6:10  # steps 02, 03 and 06 to 10
Rscript run_pipeline.R --list            # show the steps and exit
```

or source an individual script in RStudio with the project (`WGBS-gigas-ploidy-desiccation.Rproj`) open.

`run_pipeline.R` can be started from any working directory: it takes the repository to be the folder it lies in and starts each step there. It reads and checks `run_config.R` before the first step, stops at the first step that fails, and ends with a table of the steps, their minutes and their result. The selection after `--steps` is one word without spaces; a selection it cannot read (`--steps 3, 12`, an empty selection, a second `--steps`, a range that runs downward, a step that does not exist) stops the run before any step starts. A single step is started with the working directory inside the repository (`Rscript 01_pipeline/12_qpcr.R`); see How paths work in the root README.

Steps 04 and 05 need the 22 coverage files and the reference genome (10.5 GB, downloaded once by step 01), and about 9.5 GB (step 04) and about 10 GB (step 05 with `WGBS_CORES=4`) of memory. Step 05 keeps its two model fits in `04_output/05_differential_methylation/cache/` and reuses a fit only while the settings, sample sheet, coverage files, methylKit version and fitting code it was made from are unchanged; otherwise it fits again. Every later step runs from files that are in the repository, except that step 09 also needs the GO.db 3.13.0 package (31 MB, fetched from the Bioconductor archive into the cache by step 01, or by step 09 itself if step 01 was skipped).

Two environment variables are read by `run_config.R`: `WGBS_CACHE_DIR` (where downloads are kept; default `03_input_files/gannet_cache/`) and `WGBS_CORES` (cores for the model fits of step 05; default 1). `WGBS_CORES` must be a whole number of 1 or more; set to anything else, an empty value included, it stops every script at once. On Windows leave it unset: the fits fork their workers, which Windows cannot do, and a larger value is replaced by 1 with a message.

Steps 09 and 10 load different releases of GO.db and must not share an R session; `run_pipeline.R` starts each step in a fresh process, and both steps stop with a message if the wrong release is already loaded. Step 11 uses the GO.db of the R installation, as step 10 does, and stops if the legacy release is loaded; it checks that the outputs of steps 08, 09 and 10 exist before it starts.

No step reads `08_archive/`. The legacy files that steps 11 and 12 redraw or compare with are copies under `03_input_files/`: the REVIGO exports and three Fisher tables of another analysis in `annotation/legacy/revigo_exports/` and `annotation/legacy/fisher_other_analysis/` (step 11), and the knitted legacy qPCR output in `qpcr/legacy/` (step 12). The README of each folder says where a file came from. The checks in `06_diagnostics/` do read the archive.

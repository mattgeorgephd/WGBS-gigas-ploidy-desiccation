# WGBS-gigas-ploidy-desiccation

DNA methylation in diploid and triploid Pacific oysters (*Crassostrea gigas*) after an aerial heat exposure: global 5-mC by ELISA and locus-level methylation by whole-genome bisulfite sequencing (WGBS).

- **Manuscript (draft):** [*Epigenetic divergence in the thermal stress response of diploid and triploid Pacific oysters*](https://docs.google.com/document/d/17mcGDI-TWmU4vgBXmiXmeofe4qEuFH5inBKBHhG9tzg/edit), M. N. George, R. Jain, S. B. Roberts
- **Lead:** Matt George
- **Status:** The repository is being rebuilt as a reproducible, R-only pipeline on the branch `reproducibility-audit`. The audit, the folder reorganization, and pipeline steps 01 to 13 are done: the WGBS, ELISA and qPCR figures and tables of the draft are regenerated in R from the raw inputs, and the outputs of the legacy analysis are reproduced in 329 of 331 comparisons (`06_diagnostics/legacy_check_results.csv`). Step 13 ties 545 quantitative statements of the draft to pipeline outputs: 352 are reproduced, 146 are not, and 47 concern laboratory or design details that no pipeline can regenerate (`05_manuscript/manuscript_values.csv`). Reproducing a number does not make it sound: see the next paragraph. Still open: the lockfile and a final run from a clean clone; the confirmed findings of the six code reviews were applied in the code pass of 2026-11-20 (changes that would alter an analysis are left as notes for the author). The legacy code is kept in `08_archive/`.

> **Scope (decision of the author, 2026-11-19).** The goal is a reproducible pipeline that produces the results of the paper, with notes on where improvements would help; no improvement is implemented without the author's approval. The default run of the pipeline (`Rscript run_pipeline.R`) is the paper's analysis and writes nothing else. The optional analyses (the coverage filter that the Methods describe, the 5x standard, the overdispersion correction, the exclusion of CpGs at C-to-T variants, the gene-level enrichment of step 10) and the scripts of `06_diagnostics/` were built to test statements of the audit or on the author's request; they are the evidence behind the notes in `issues.md`, are not adopted analyses, are off by default, and write only under folders named `not_adopted` (section Optional analyses, below).
>
> **Resuming the work.** [`07_documentation/session_handoff.md`](07_documentation/session_handoff.md) says what is done, what is in progress, and what comes next, and how to set up a fresh machine.
>
> **Where the work stands.** Read [`issues.md`](issues.md) first: it lists 50 problems, eight of which bear on the main conclusions. (1) The heat treatment is confounded with experiment and sequencing batch in the WGBS data. (2) The "interaction" loci are the overlap of two lists, not a tested interaction. (3) The GO enrichment counts transcript isoforms as genes. (4) Its gene universe is missing 36% of the annotated transcripts. (5) A file-reading fault skipped 76% of its ploidy input table and 24% of its heat input table. (6) Where DML lie, and what the genes holding them do, is explained by where methylation and between-animal variation sit in the genome: against a matched or label-permuted background DML show no preference for gene bodies, and no specific GO category is established. (7) Figures 6 to 8 were not drawn from enriched terms at all, but from the GO terms attached to genes with a DML, as a faulty file read returned them. (8) The differential methylation model has no term for variation among libraries, and its test statistic is inflated nearly fourfold (about threefold by that variation, the rest by the scaling of the counts): the ploidy list cannot be told apart in length from lists made from relabelled libraries, and methylKit's corrected test (chi-squared, the setting of the paper the Methods follow) leaves 14 ploidy CpGs and 185 heat CpGs, from a test that does not hold its level on these libraries: with relabellings that keep little of the real grouping it gives 1.6% and 0.31% of p-values below 0.01 and 0.001 (1% and 0.1% expected), and relabelled libraries average 1.6 calls with p below 2.1e-8, where 0.01 are expected. Separately, one relabelling of the libraries in eight gives a corrected ploidy list at least as long as the real one. One further finding (I-49, rated S2) changes how Results 3.1 is read. Most of the ploidy by heat interaction in overall methylation (Figure 1B) rests on five libraries with a raised level, two of them from untreated animals. Inside the heat experiment every diploid library departs from a profile that its triploids share with most control libraries, and the ploidy difference seen in one experiment is not found in the other, where the overall level differs in the opposite direction; whether heat, the source of the animals or the preparation of the libraries is the cause cannot be decided without sequencing that experiment's own controls. The statements of the log written up to the third pass of 2026-11-19 were re-derived by independent agents instructed to refute them; what that changed is listed at the end of `issues.md`. Point (8) and the finding on Figure 1B come from the fourth pass. Their text has been through three checks, which reproduced the numbers of the tables, corrected seven values quoted wrongly, withdrew one reading, narrowed twenty-four and lowered the rating of I-49 from S1 to S2; the independent review of the second script reproduced its tables and found no wrong calculation; the check of the citations read the benchmark that I-06 cites and corrected three cited statements of the log (I-02, I-06, I-08). The plan and the decisions needed from the authors are in [`07_documentation/reproducibility_plan.md`](07_documentation/reproducibility_plan.md).

---

## The study in brief

Adult diploid and triploid oysters were either held as controls or exposed to air at 27 C for 24 h ("heat"; older files call this treatment "desiccation"). Ctenidia were sampled.

| Assay | Samples | Design |
|---|---|---|
| Global 5-mC ELISA | 32 animals, 8 per ploidy by treatment group | Controls and heat-treated animals from the same experiment, one plate |
| WGBS | 22 libraries: 5 diploid and 5 triploid heat-treated, 6 diploid and 6 triploid controls | Heat-treated libraries are Zymo project zr3534 (the heat experiment); control libraries are Zymo project zr3644 (animals from a separate experiment). Four of the ten heat-treated animals also received 1 h at 45 C. |
| qPCR | 40 animals, 9 to 12 genes | By their numbers the 32 animals of the ELISA and the 8 animals given the acute heat shock; to be confirmed by the authors (I-29) |

Sample names in the heat experiment (D = diploid, T = triploid; suffix `-C` = ctenidia):

| Animals | Treatment |
|---|---|
| D01 to D08, T01 to T08 | Control |
| D09, D10, T09, T10 | Control, then 1 h acute heat shock at 45 C |
| D11 to D18, T11 to T18 | Heat (24 h aerial emersion at 27 C) |
| D19, D20, T19, T20 | Heat, then 1 h acute heat shock at 45 C |

## Repository layout

Numbered folders follow the order of the work: scripts, helpers, inputs, outputs, manuscript products, checks, documentation, archive. Each folder has a `README.md` with a file inventory.

```text
WGBS-gigas-ploidy-desiccation/
├── 01_pipeline/          Numbered R scripts, one per analysis step
├── 02_R_functions/       Helper functions sourced by every step
├── 03_input_files/       Inputs kept in GitHub, and the manifest of inputs fetched from gannet
├── 04_output/            Results, one subfolder per pipeline step
├── 05_manuscript/        Figures and tables as numbered in the manuscript; value-by-value traceability table
├── 06_diagnostics/       Checks of regenerated results against legacy results and manuscript values
├── 07_documentation/     Plan, audit notes, upstream processing records, protocols
├── 08_archive/           Everything from before the reorganization, original paths preserved
├── issues.md             Issue log: what is wrong, the evidence, what needs deciding
├── README.md             This file
├── run_config.R          Settings of every step, checked each time a script reads them
├── run_pipeline.R        Runs the steps in order
├── .gitattributes        Keeps git from converting line endings in the data files
├── .gitignore
└── WGBS-gigas-ploidy-desiccation.Rproj
```

| Folder | Contents | README |
|---|---|---|
| `01_pipeline/` | The thirteen analysis steps and how a step is written | [01_pipeline/README.md](01_pipeline/README.md) |
| `02_R_functions/` | Shared helpers | [02_R_functions/README.md](02_R_functions/README.md) |
| `03_input_files/` | Sample sheets, ELISA plate reads, qPCR Cq values, Bismark reports, genome feature tracks, annotation tables, and the gannet manifest | [03_input_files/README.md](03_input_files/README.md) |
| `04_output/` | Tables and figures written by the pipeline | [04_output/README.md](04_output/README.md) |
| `05_manuscript/` | Final figures, tables, supplementary lists, and the table tying each manuscript number to its source | [05_manuscript/README.md](05_manuscript/README.md) |
| `04_output/not_adopted/`, `05_manuscript/not_adopted/` | The outputs of the optional analyses, which the paper did not adopt | [04_output/not_adopted/README.md](04_output/not_adopted/README.md) |
| `06_diagnostics/` | Reproduction checks | [06_diagnostics/README.md](06_diagnostics/README.md) |
| `07_documentation/` | Plan, audit, provenance | [07_documentation/README.md](07_documentation/README.md) |
| `08_archive/` | Legacy code, outputs, figures, slides, with a guide to which legacy script produced which manuscript item | [08_archive/README.md](08_archive/README.md) |
| `tests/` | Tests of the settings layout (settings, paper bytes, optional runs, step 13); not run by the pipeline | headers of the scripts |

At the root: `run_config.R` is the single control surface (thresholds, the optional analyses, cache location, and the settings of each step) and `run_pipeline.R` runs the numbered steps in order, each in a fresh R process. Every script reads `run_config.R` first, and the file checks the settings of the analysis: a value that cannot work (an unknown optional analysis, a change of one of the five frozen settings of the paper, a `WGBS_CORES` that is empty or not a whole number, a threshold, tolerance or seed out of range) stops the script at once with a message that names the setting. The sizes, colours and label placement of the figures, and the facts it records about the draft and the legacy files, are not checked; the header of `run_config.R` lists what is not. `.gitattributes` tells git to store and check out the tables and every input exactly as they are, so that a checkout on Windows holds the same bytes. `renv.lock` records the package versions (step 1 of Running the pipeline).

### How paths work

Scripts build every path with `here::here()`, which looks for the repository root (the folder with the `.Rproj` file) upward from the working directory. A step therefore runs unchanged from RStudio, with the project open, or from a terminal whose working directory is anywhere inside the repository, on any operating system, and no step calls `setwd()`. Two folders of the archive, `08_archive/gene_expression/` and `08_archive/ELISA/`, hold project files of their own: a step started from inside one of them does not find the repository and stops with an error when it looks for `run_config.R`. `run_pipeline.R` can be started from any directory: it takes the repository to be the folder it lies in, starts every step with that folder as working directory, and stops with an error when it finds no steps there. If a numbered folder is renamed, update the matching string in `02_R_functions/paths.R` and the places that cannot use it: the line of each script that sources `02_R_functions/`, `run_pipeline.R` (`01_pipeline`), the two build scripts under `03_input_files/`, step 13 and the diagnostics scripts (which name `04_output` and `06_diagnostics`), the `source` column of `03_input_files/manuscript/manuscript_items_map.csv`, two expressions of `manuscript_values_map.csv` in the same folder (rows D4.8-002 and R3.4-007 name `03_input_files`), and the patterns of `.gitignore` and `.gitattributes`.

## Where the data live

**GitHub first.** Any input that fits in GitHub is in `03_input_files/` (ceiling: 50 MB per file, after gzip where possible). That covers everything except the per-library coverage files and the reference genome.

**gannet for the large files.** The 22 Bismark coverage files (about 445 MB each, 9.8 GB in total) stay at <https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/cov_files/>, and the reference genome (656 MB) at `/panopea/Cg-roslin/`. They are listed with md5 checksums in [`03_input_files/gannet_manifest.csv`](03_input_files/gannet_manifest.csv), downloaded once by the first pipeline step into a git-ignored cache, and verified.

**One file from Bioconductor.** Step 09 reproduces the legacy GO enrichment, which matches the archived tables only under the ontology release it was run on. That release (the R package GO.db 3.13.0, 31 MB) is fetched from the Bioconductor archive into the same cache; its address and md5 are in `run_config.R`.

Other locations, for reference:

| What | Where |
|---|---|
| Alignment outputs (BAM, reports) for all libraries | gannet `/panopea/061021-big/` (June 2021 run, used for the analysis) and `/panopea/030521-ronrosM/` (March 2021 run, superseded) |
| Genome feature tracks, original copies | OSF, doi:[10.17605/OSF.IO/YGCTB](https://doi.org/10.17605/OSF.IO/YGCTB) |
| Reference genome | Roslin v1, RefSeq GCF_902806645.1, with mitochondrial genome NC_001276.1; gannet `/panopea/Cg-roslin/` |
| Raw reads | NCBI SRA BioProject PRJNA678408 for the heat-treated libraries; the accession for the control libraries is not recorded here (issue I-36) |
| Legacy intermediates and saved workspaces | gannet project directory; inventory in `07_documentation/audit/gannet_manifest_2026-11-18.tsv` |

## Running the pipeline

Steps 01 to 13 run end to end (see the status column in [`01_pipeline/README.md`](01_pipeline/README.md)).

1. Install R 4.3.3 and the packages. `renv.lock` at the root records the 201 packages the pipeline ran with (R 4.3.3, Bioconductor 3.18, methylKit 1.28.0, among them the optional analyses' Rsamtools and GenomicAlignments). The project does not activate `renv` (there is no `.Rprofile`), so a library is restored explicitly, from the repository root:
   ```r
   install.packages("renv")
   renv::restore(lockfile = "renv.lock", library = .libPaths()[1], prompt = FALSE)
   ```
   On Linux several packages need system libraries first (libcurl, libxml2, libssl, zlib, libbz2, liblzma, libglpk; the `apt-get` line of `07_documentation/cloud_environment_setup.sh` lists them). That script is how the machine of the recorded runs was set up (Ubuntu 24.04: Ubuntu's R packages and Bioconductor 3.18, not pinned), and `renv.lock` is a snapshot of the library it produced. A restore from the lockfile on a clean machine has not yet been tested.
2. Run `Rscript run_pipeline.R` for everything or `Rscript run_pipeline.R --steps 6:10` for a subset (the selection is one word without spaces, such as `5:7`, `3,12` or `2:3,6:10`; anything else on the command line stops the run before a step starts). Or open `WGBS-gigas-ploidy-desiccation.Rproj` and source the numbered scripts in `01_pipeline/` one at a time, each in a fresh R session.

The first run downloads 10.5 GB (plus the 31 MB GO.db archive). Step 04 needs about 9.5 GB of memory and step 05 about 10 GB with `WGBS_CORES=4`; together they take about half an hour on four cores; every other step takes a few minutes or less and runs from files in the repository. To use the outputs already in the repository and skip the heavy steps, start at step 06.

### Optional analyses (not adopted)

`Rscript run_pipeline.R` with no switch is the paper: it writes the paper's outputs and nothing else. Each optional analysis is off by default and is turned on by a switch (or the same setting in `run_config.R`); the switches combine, and one run computes one analysis, whose outputs go to `04_output/not_adopted/<id>/` and `05_manuscript/not_adopted/`. Every log names the run kind and the analysis; `Rscript run_pipeline.R --plan <switches>` prints the plan of a run without running it. The five settings of the paper (`min_coverage`, `import_min_coverage`, `upper_coverage_percentile`, `dml_qvalue`, `dml_difference`) are frozen in `run_config.R`. Details: [`04_output/not_adopted/README.md`](04_output/not_adopted/README.md).

| Analysis | Adopted in the paper? | Switch of `run_pipeline.R` | Setting in `run_config.R` (list `optional`) | Id or token | Outputs | What it changes | Issue |
|---|---|---|---|---|---|---|---|
| The paper | yes | none: `Rscript run_pipeline.R` | every setting at its default | `legacy` | `04_output/<step>/`, `05_manuscript/` | at least 10 reads per CpG in each of the 22 libraries, no upper coverage filter, first line of each coverage file dropped, logistic regression without overdispersion correction, DML at q < 0.01 and a difference of more than 20 points; steps 01 to 09 and 11 to 13 | |
| Coverage filter of the Methods | no | `--coverage-filter` | `coverage_filter_of_methods = TRUE` | `corrected` | `04_output/not_adopted/corrected/` | CpGs at or above the 99.9th percentile of coverage in a library removed; every line of each coverage file read | I-05, I-38 |
| 5x standard | no | `--coverage 5` | `coverage_standard = 5L` | `cov5` | `04_output/not_adopted/cov5/` | 5 reads per CpG in each library (2,297,236 CpGs tested instead of 480,579) | Venkataraman et al. 2022 |
| Overdispersion correction | no | `--overdispersion` | `overdispersion_chisq = TRUE` | `mnchisq` | `04_output/not_adopted/mnchisq/` | methylKit `overdispersion = "MN"`, `test = "Chisq"`; no other test is offered | I-06 |
| C-to-T variant exclusion | no | `--exclude-ct-variants` | `exclude_ct_variants = TRUE` | `ctexcl` | `04_output/not_adopted/ctexcl/` | CpGs listed by `01_pipeline/optional/ct_variant_loci.R` removed from the CpG sets of step 04 and, after the test, from the tested CpGs and the DML (fits, p-values and q-values unchanged); the run refuses to start without the list | I-08 |
| Difference of 50 points | no | `--difference 50` | `dml_difference_standard = 50L` | `d50` | `04_output/not_adopted/d50/` | DML at q < 0.01 and an absolute difference of more than 50 points instead of 20 (methylKit `getMethylDiff`, `difference = 50`); applied after the test, so the fit of the analysis without `d50` is reused (author's decision of 2026-11-21) | Venkataraman et al. 2022 |
| Any combination of the five | no | the switches together | the settings together | the tokens joined by `-` in the order above | `04_output/not_adopted/<id>/`, for example `corrected-cov5-ctexcl/`, `mnchisq-d50/` | each of its parts | |
| Gene-level GO enrichment | no | `--gene-level-go`, alone (on the paper's DML) or with any of the above | `gene_level_go = TRUE` | | `04_output/not_adopted/<id>/10_go_enrichment_gene_level/` and `.../11_semantic_similarity/` (the paper's DML: `<id>` = `legacy`) | step 10 and the gene-level part of step 11 | I-41 |
| Diagnostics | no | `--diagnostics` (all) or `--diagnostics name1,name2` | `diagnostics = "all"` or names | | `06_diagnostics/` | the checks and controls of the audit, listed in `06_diagnostics/diagnostics.csv` | |
| Side by side | no | made by every optional run; the table of all analyses by `Rscript 01_pipeline/optional/compare_analyses.R` | | | `04_output/not_adopted/comparisons/` (`analysis_comparison.md`), `05_manuscript/not_adopted/` | the tables that list analyses, from every analysis folder on disk, the paper first; the manuscript values of every analysis, computed as the paper computed them and with the corrected recipe of the audit ("not run" where absent) | |

## Conventions

- The treatment is "heat" or "heat-treated" in all new code, outputs, and text. "Desiccation" survives only in `08_archive/` and in this repository's name.
- Legacy dataset labels: `ronit` = the heat-treated animals alone; `yaamini` = the control animals alone; `both` or `all-gigas` = the 22-library comparison reported in the manuscript.
- Diploid is drawn in blue and triploid in orange-red throughout.

## Related material

- Research poster by Ronit Jain: `08_archive/presentations/Triploid-Diploid-Oyster-Poster_RonitJain.pptx`
- Earlier manuscript draft by Ronit Jain: [Google Doc](https://docs.google.com/document/d/1nwY9I3pVzF5Xlfdzb7SdO1QKaNp5qyj0DuPCCa_YXi8/edit?usp=sharing)
- Global DNA methylation analysis by Shelly Trigg (2019), the origin of `08_archive/ELISA/GlobalDNAMeth_Polyploids.*`: [script](https://github.com/shellytrigg/C_gigas/blob/master/Polyploids/GlobalDNAMeth_Polyploids.R), [rendered analysis](https://github.com/shellytrigg/C_gigas/blob/master/Polyploids/GlobalDNAMeth_Polyploids.md)

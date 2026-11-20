# 04_output/not_adopted

**Nothing in this folder is a result of the paper.** It holds the outputs of the optional analyses that the paper did not adopt. The paper's outputs are in `04_output/<step>/` and `05_manuscript/` outside the folders named `not_adopted`: a path that contains `not_adopted` is not the paper's.

Each optional analysis is off by default. It is turned on by a switch of `run_pipeline.R` or by a setting of `run_config.R` (list `optional`); the switches combine. One run computes one analysis, named by its id: the tokens of the switches that are on, in the order `corrected`, `cov5`, `mnchisq`, `ctexcl`, joined by `-` (16 ids, `legacy` being the paper's). `Rscript run_pipeline.R --plan <switches>` prints what a run would do.

| Analysis | Adopted in the paper? | Switch of `run_pipeline.R` | Setting in `run_config.R` (list `optional`) | Id or token | Outputs | What it changes | Issue |
|---|---|---|---|---|---|---|---|
| The paper | yes | none: `Rscript run_pipeline.R` | every setting at its default | `legacy` | `04_output/<step>/`, `05_manuscript/` | at least 10 reads per CpG in each of the 22 libraries, no upper coverage filter, first line of each coverage file dropped, logistic regression without overdispersion correction, DML at q < 0.01 and a difference of more than 20 points; steps 01 to 09 and 11 to 13 | |
| Coverage filter of the Methods | no | `--coverage-filter` | `coverage_filter_of_methods = TRUE` | `corrected` | `04_output/not_adopted/corrected/` | CpGs at or above the 99.9th percentile of coverage in a library removed; every line of each coverage file read | I-05, I-38 |
| 5x standard | no | `--coverage 5` | `coverage_standard = 5L` | `cov5` | `04_output/not_adopted/cov5/` | 5 reads per CpG in each library (2,297,236 CpGs tested instead of 480,579) | Venkataraman et al. 2022 |
| Overdispersion correction | no | `--overdispersion` | `overdispersion_chisq = TRUE` | `mnchisq` | `04_output/not_adopted/mnchisq/` | methylKit `overdispersion = "MN"`, `test = "Chisq"`; no other test is offered | I-06 |
| C-to-T variant exclusion | no | `--exclude-ct-variants` | `exclude_ct_variants = TRUE` | `ctexcl` | `04_output/not_adopted/ctexcl/` | CpGs listed by `01_pipeline/optional/ct_variant_loci.R` removed from the CpG sets of step 04 and, after the test, from the tested CpGs and the DML (fits, p-values and q-values unchanged); the run refuses to start without the list | I-08 |
| Any combination of the four | no | the switches together | the settings together | the tokens joined by `-` in the order above | `04_output/not_adopted/<id>/`, for example `corrected-cov5-ctexcl/` | each of its parts | |
| Gene-level GO enrichment | no | `--gene-level-go`, alone (on the paper's DML) or with any of the above | `gene_level_go = TRUE` | | `04_output/not_adopted/<id>/10_go_enrichment_gene_level/` and `.../11_semantic_similarity/` (the paper's DML: `<id>` = `legacy`) | step 10 and the gene-level part of step 11 | I-41 |
| Diagnostics | no | `--diagnostics` (all) or `--diagnostics name1,name2` | `diagnostics = "all"` or names | | `06_diagnostics/` | the checks and controls of the audit, listed in `06_diagnostics/diagnostics.csv` | |
| Side by side | no | made by every optional run | | | `04_output/not_adopted/comparisons/`, `05_manuscript/not_adopted/` | the tables that list analyses, from every analysis folder on disk, the paper first; the manuscript values of every analysis ("not run" where absent) | |

Examples:

```sh
Rscript run_pipeline.R                                                       # the paper, nothing else
Rscript run_pipeline.R --coverage-filter                                     # -> not_adopted/corrected/
Rscript run_pipeline.R --overdispersion                                      # -> not_adopted/mnchisq/
Rscript run_pipeline.R --coverage-filter --overdispersion --gene-level-go    # -> not_adopted/corrected-mnchisq/, with steps 10 and 11
Rscript run_pipeline.R --gene-level-go                                       # step 10 on the paper's DML -> not_adopted/legacy/
Rscript run_pipeline.R --diagnostics cohort_separation                       # one script of 06_diagnostics/
Rscript run_pipeline.R --plan --coverage 5 --exclude-ct-variants             # print the plan only
WGBS_ANALYSIS=mnchisq Rscript 01_pipeline/05_differential_methylation.R      # one step alone, for one analysis
```

## Layout

```text
04_output/not_adopted/
  README.md                           this file
  <id>/                               one optional analysis, for example corrected/, mnchisq/, corrected-mnchisq/
    analysis.csv                      stamp: settings of the analysis, md5 of its fit and of the variant list, commit, time, steps done
    04_cpg_coverage/ ... 09_go_enrichment_legacy/    the files of the paper's steps, <id> in place of "legacy" in the names
                                      (dml_ploidy_mnchisq.csv, figure_2_mnchisq.png); the paper's names without a mode keep their names
    10_go_enrichment_gene_level/      only with --gene-level-go
    11_semantic_similarity/           only with --gene-level-go: the gene-level part
  legacy/                             --gene-level-go on the paper's DML: steps 10 and 11 (gene level) only
  ct_variants/                        the C-to-T variant list of 01_pipeline/optional/ct_variant_loci.R (tables/; cache/ is not tracked)
  comparisons/                        the tables that list analyses side by side, rebuilt by every optional run from every
                                      analysis folder on disk, the paper first; analyses_included.csv lists the stamps used
  13_manuscript_values/tables/        expression_errors_by_analysis.csv, manuscript_items_copied.csv of the optional part of step 13
```

Rules that keep the paper's outputs safe: an optional run never opens a file of a paper folder for writing (`output_dir()` of `02_R_functions/paths.R` stops such a run, except for the git-ignored `cache/` folders), and the paper run never writes here. In an optional run the steps that do not depend on the analysis (01 to 03 and 12) are not started, and three tables of step 04 and step 08 that do not depend on it (`library_coverage.csv`, `genome_cpg_total.csv`; `gene_go_map.csv.gz`, `transcript_table.csv.gz`, `go_map_id_repair.csv`) are read from the paper folder and not written again.

Column `mode` of every table holds the analysis id (the paper's tables keep the value `legacy`).

## What is here now

| Folder | Made by | Steps |
|---|---|---|
| `corrected/` | the corrected mode of the pipeline before the settings layout (committed outputs moved here on 2026-11-20 with `git mv`); rebuilt on 2026-11-20 by `--coverage-filter` and `--coverage-filter --gene-level-go`, which gave every moved table, figure and step 10 file back byte for byte (454,555 CpGs, 952 ploidy and 2,273 heat DML, 159 in both lists), except that the eight `legacy_` columns of `annotation_summary.csv`, empty before, are now filled | 04 to 11 |
| `legacy/` | step 10 and the gene-level part of step 11 of the paper's DML, before the settings layout (moved); rebuilt by `--gene-level-go` | 10, 11 (gene level) |
| `mnchisq/` | `--overdispersion` (and `--gene-level-go`), 2026-11-20: 480,579 CpGs tested, 14 ploidy and 185 heat DML, 4 in both lists; the two lists are exactly those of `06_diagnostics/dml_sensitivity/dml_A1_*_mn_chisq_legacy.csv`, and every locus is in the paper's lists | 04 to 11 |
| `corrected-mnchisq/` | `--coverage-filter --overdispersion`, 2026-11-20: 454,555 CpGs tested, 14 ploidy and 157 heat DML, 4 in both lists | 04 to 09 |
| `comparisons/` | the tables that listed the legacy and the corrected mode side by side before the settings layout, now rebuilt by every optional run from every analysis folder | |
| `ct_variants/` | `01_pipeline/optional/ct_variant_loci.R` (see `01_pipeline/optional/README.md`) | |

Other analysis folders appear when their run has been made. The first optional analysis whose fit is not cached refits in step 05 (about 18 minutes at `WGBS_CORES=4`, 37 to 45 at 2); an analysis that adds `ctexcl` to another reuses that one's fit (`methylkit_<fit id>.rds` in `04_output/05_differential_methylation/cache/`). The per-CpG tables of an analysis at 5x (`*cov5*`) are larger than the 50 MB limit of the repository and are not tracked (`.gitignore`).

The values of the manuscript for every analysis are in `05_manuscript/not_adopted/manuscript_values_by_analysis.csv`, written by step 13 in an optional run, with "not run" for the analyses that have no folder.

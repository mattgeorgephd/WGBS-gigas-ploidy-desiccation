# 04_output/not_adopted

**Nothing in this folder is a result of the paper.** It holds the outputs of the optional analyses that the paper did not adopt. The paper's outputs are in `04_output/<step>/` and `05_manuscript/` outside the folders named `not_adopted`: a path that contains `not_adopted` is not the paper's.

Each optional analysis is off by default. It is turned on by a switch of `run_pipeline.R` or by a setting of `run_config.R` (list `optional`); the switches combine. One run computes one analysis, named by its id: the tokens of the switches that are on, in the order `corrected`, `cov5`, `mnchisq`, `ctexcl`, `d50`, joined by `-` (32 ids, `legacy` being the paper's). `Rscript run_pipeline.R --plan <switches>` prints what a run would do.

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

Examples:

```sh
Rscript run_pipeline.R                                                       # the paper, nothing else
Rscript run_pipeline.R --coverage-filter                                     # -> not_adopted/corrected/
Rscript run_pipeline.R --overdispersion                                      # -> not_adopted/mnchisq/
Rscript run_pipeline.R --overdispersion --difference 50                      # -> not_adopted/mnchisq-d50/ (reuses the fit of mnchisq)
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
                                      (step 04 removes it and the folders of steps 05 and later before it writes; a table
                                      counts, in step 13 and in comparisons/, only while the stamp lists its step)
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

Rules that keep the paper's outputs safe: an optional run never opens a file of a paper folder for writing (`output_dir()` of `02_R_functions/paths.R` stops such a run, except for the git-ignored `cache/` folders), and the paper run never writes here. In an optional run the steps that do not depend on the analysis (01 to 03 and 12) are not started, and the tables of step 04 and step 08 that only the paper run writes (`library_coverage.csv`, `genome_cpg_total.csv`; `gene_go_map.csv.gz`, `transcript_table.csv.gz`, `go_map_id_repair.csv`) are read from the paper folder and not written again. They do not depend on the analysis, with one exception: under `ctexcl` the libraries lose the listed CpGs, so the analysis's own `library_coverage.csv` would differ (CpGs per library, first line, medians, cutoffs); it is not written, and step 13 reads the paper's for every analysis.

Column `mode` of every table holds the analysis id (the paper's tables keep the value `legacy`).

## What is here now

| Folder | Made by | Steps |
|---|---|---|
| `corrected/` | the corrected mode of the pipeline before the settings layout (committed outputs moved here on 2026-11-20 with `git mv`); rebuilt on 2026-11-20 by `--coverage-filter` and `--coverage-filter --gene-level-go`, which gave every moved table, figure and step 10 file back byte for byte (454,555 CpGs, 952 ploidy and 2,273 heat DML, 159 in both lists), except that the eight `legacy_` columns of `annotation_summary.csv`, empty before, are now filled | 04 to 11 |
| `legacy/` | step 10 and the gene-level part of step 11 of the paper's DML, before the settings layout (moved); rebuilt by `--gene-level-go` | 10, 11 (gene level) |
| `mnchisq/` | `--overdispersion` (and `--gene-level-go`), 2026-11-20: 480,579 CpGs tested, 14 ploidy and 185 heat DML, 4 in both lists; the two lists are exactly those of `06_diagnostics/dml_sensitivity/dml_A1_*_mn_chisq_legacy.csv`, and every locus is in the paper's lists | 04 to 11 |
| `corrected-mnchisq/` | `--coverage-filter --overdispersion`, 2026-11-20: 454,555 CpGs tested, 14 ploidy and 157 heat DML, 4 in both lists | 04 to 09 |
| `ctexcl/` | `--exclude-ct-variants`, 2026-11-20, after the variant screen (no refit; the paper's fit): of the 480,579 tested CpGs 28,420 are listed and removed, leaving 452,159; 717 ploidy and 1,667 heat DML (1,038 and 2,351 in the paper), 139 in both lists (200); `ct_exclusion_summary_ctexcl.csv` (issues.md, I-08) | 04 to 09 |
| `cov5/` | `--coverage 5`, 2026-11-20 (step 04) and 2026-11-21 (steps 05 to 09 and 13): 2,297,236 CpGs tested, 9,942 ploidy DML (4,408 hyper, 5,534 hypo) and 21,361 heat DML, 1,784 in both lists; step 05 fitted by sequence in 107 minutes at `WGBS_CORES=2`, peak memory 9.3 GB. Step 13 cannot evaluate R3.5-046 and D4.9-012, because the 5x Table 1 lists two interaction CpGs of gene 105340590 (`13_manuscript_values/tables/expression_errors_by_analysis.csv`) (issues.md, I-51) | 04 to 09 |
| `cov5-ctexcl/` | `--coverage 5 --exclude-ct-variants`, 2026-11-21 (reuses the 5x fit): of the 2,297,236 tested CpGs 144,219 are listed and removed, leaving 2,153,017; 7,046 ploidy and 15,530 heat DML, 1,219 in both lists (issues.md, I-08 and I-51) | 04 to 09 |
| `corrected-cov5-ctexcl/` | `--coverage-filter --coverage 5 --exclude-ct-variants`, 2026-11-21: 5x fit `corrected-cov5` by sequence (88 minutes); of the 2,265,759 tested CpGs 142,790 are listed and removed, leaving 2,122,969; 6,902 ploidy and 15,376 heat DML (9,791 and 21,200 before the exclusion), 1,148 in both lists | 04 to 09 |
| `corrected-cov5-mnchisq-ctexcl/` | `--coverage-filter --coverage 5 --overdispersion --exclude-ct-variants`, 2026-11-21, the settings of Venkataraman et al. 2022 for coverage, the test and variant removal, except their 50-point threshold and the full CIGAR of the variant rule (I-08): 5x fit `corrected-cov5-mnchisq` by sequence (89 minutes, peak memory 9.04 GB); 142,790 CpGs removed, leaving 2,122,969; 96 ploidy and 875 heat DML (126 and 1,133 before the exclusion), 19 in both lists | 04 to 09 |
| `mnchisq-ctexcl/` | `--overdispersion --exclude-ct-variants --gene-level-go` (combination C2), 2026-11-21, the fit of `mnchisq` from the cache: 452,159 CpGs tested, 12 ploidy and 142 heat DML (14 and 185 before the exclusion), 4 in both lists; step 10 calls no GO term against the eligible universe | 04 to 11 |
| `corrected-ctexcl/`, `corrected-mnchisq-ctexcl/`, `corrected-cov5/`, `corrected-cov5-mnchisq/`, `cov5-mnchisq/`, `cov5-mnchisq-ctexcl/` and their `-d50` versions | the six combinations not run before, at 20 and 50 points, 2026-11-21 (author's request); all with cached fits except `cov5-mnchisq`, a new 5x fit by sequence (resumed from its saved parts after a container restart); counts in the comparison below | 04 to 09 |
| `d50/`, `ctexcl-d50/`, `mnchisq-d50/`, `mnchisq-ctexcl-d50/`, `corrected-d50/`, `corrected-mnchisq-d50/`, `cov5-d50/`, `cov5-ctexcl-d50/`, `corrected-cov5-ctexcl-d50/`, `corrected-cov5-mnchisq-ctexcl-d50/` | `--difference 50` added to the switches of each analysis above, 2026-11-21, each with the fit of its analysis without `d50` from the cache (10 to 14 minutes per run); the last is the analysis closest to the settings of Venkataraman et al. 2022 (their stated settings, except the variant rule, which counts reads with their full CIGAR, I-08). DML counts in the comparison below | 04 to 09 |
| `comparisons/` | the tables that listed the legacy and the corrected mode side by side before the settings layout, now rebuilt by every optional run from every analysis folder | |
| `ct_variants/` | `01_pipeline/optional/ct_variant_loci.R` (see `01_pipeline/optional/README.md`), finished 2026-11-20: of the 13,245,836 CpG sites of the reference 407,217 listed (3.07%) (`tables/ct_variant_log.txt`) | |

Other analysis folders appear when their run has been made. The first optional analysis whose fit is not cached refits in step 05 (22 to 24 minutes per fit at 10x and `WGBS_CORES=2`, as measured on 2026-11-20; 107 minutes at 5x and `WGBS_CORES=2`, peak memory 9.3 GB, fitted by sequence, measured on 2026-11-21); an analysis that adds `ctexcl` or `d50` (or both) to another reuses that one's fit (`methylkit_<fit id>.rds` in `04_output/05_differential_methylation/cache/`; the fit id is the analysis id without `ctexcl` and `d50`).

### The option `d50`

Author's decision of 2026-11-21: "The 2022 paper also used a 50% cut off. All variations of the runs should also be run with this cutoff as well as the 20% I use in the paper." Venkataraman et al. 2022 (BMC Genomics 23:556) called DML with methylKit `getMethylDiff` at `difference = 50` and q < 0.01. With `--difference 50` step 05 keeps the CpGs with q < 0.01 and an absolute difference of more than 50 points (strict, as `getMethylDiff` applies it); the paper's setting `dml_difference = 20` of `run_config.R` stays frozen, and the paper's run never changes. `--difference 20` names the paper's value and does not turn off `dml_difference_standard = 50L` set in `run_config.R` (the same rule as `--coverage 10`).

- **Fit.** The threshold is applied after the test, so a `d50` analysis reuses the fit of its analysis without `d50` (for `d50` itself the paper's fit, `methylkit_legacy.rds`); nothing of the threshold enters the fingerprint of a fit.
- **Step 04** of a `d50` analysis is that of its analysis without `d50` (same CpG sets), run again under its own id: 7 to 9.5 minutes and up to 9.2 GB at 10x. A `d50` run must therefore not run beside a 5x fit.
- **What changes besides the DML lists.** The "eligible" CpGs of steps 07, 08 and 10 (mean methylation from half the threshold to 100 minus that: 25% to 75% at 50 points, 10% to 90% at 20; with the paper's fit 30,315 of 480,579 CpGs, 6.3%, at 50 points against 65,480, 13.6%, at 20), and in step 13 the value of `config$dml_difference` in the expressions of both recipes (50), so the statements that name the threshold of 20 points, or apply it to group means, are evaluated at 50. `dml_counts_by_threshold_<id>.csv` (5 to 50 points) is the same as that of the analysis without `d50`.
- **Few or no DML.** At 50 points the lists are short (with the paper's fit: 1 ploidy DML, 4 heat DML, no locus in both lists). Steps 06 to 10 then write their tables with their columns and the rows there are, and replace a panel that cannot be drawn with a note. Step 09 does not run the Fisher test of a set without a transcript of interest in the legacy universe, and step 10 does not test a set none of whose DML genes carries a term of the ontology; both log it and write the set's tables without rows. Step 13 evaluates the statements as for any other analysis: a value about a locus or gene not in the lists is "not comparable: no value in this analysis", and an expression that fails because a DML set is too small is "not evaluable: too few DML (<message>)" (column `not_evaluable` of the summary), not an error. A set is too small when it has fewer than 2 DML in `dml_counts_by_mode.csv` of the analysis, or when step 09 could not test it (its Fisher table has no row); the expression must have read a table of that set (or, when every set has fewer than 2 DML, any table of steps 05 to 10 of the analysis). An expression that read `go_gene_level_<set>_<id>.csv` of a set step 10 did not test (no row) is not evaluable even when it gives a value. No error message is matched, so a failure with another cause in such an expression is also counted as not evaluable: its message is kept in the status and in `04_output/not_adopted/13_manuscript_values/tables/expression_errors_by_analysis.csv` (column `state`: `error` or `not evaluable: too few DML`), and step 13 logs the count of such statements on a line of its own. In `legacy_top10_at_gene_level.csv` (step 10) and `named_terms_results_3_4_gene_level.csv` (step 11) a term of an untested set has the status "set not tested: no DML gene carries a term of the ontology". Statements that read only `go_gene_level_summary.csv` (for example R3.4-007, "no term of the ploidy set is enriched") cannot be told apart and keep their value. Tested on 2026-11-21 in a scratch copy: `d50` and `mnchisq-d50` have 8 such statements under the paper recipe (R3.2-029 and R3.2-031, clusterings that need 2 loci; F6-002, F8-002 and D4.1-003, -005, -007, -010, which look up GO terms in Fisher tables without rows) and 2 under the corrected recipe (R3.2-029, R3.2-031); with `--gene-level-go` the corrected recipe of `d50` has 33, the 31 more reading the step 10 tables of the ploidy and shared sets, which step 10 did not test. None is an error; real errors (as the two of `cov5`, I-51) stay errors. With an unrelated fault injected into two statements that read tables of the ploidy set (an unknown object), the analyses at 20 points report it as an error and the `d50` analyses as not evaluable, both with the message `object '...' not found`.

#### Commands for the 50-point version of every analysis run so far

Of the 16 combinations, 10 were run first at 20 points (on 2026-11-20 and 2026-11-21); the other six, and the 50-point version of each, follow the same pattern (add `--coverage-filter`, `--coverage 5`, `--overdispersion`, `--exclude-ct-variants` and `--difference 50` as needed). The first 10 are `legacy` (the paper), `corrected`, `mnchisq`, `corrected-mnchisq`, `ctexcl`, `mnchisq-ctexcl`, `cov5`, `cov5-ctexcl`, `corrected-cov5-ctexcl` and `corrected-cov5-mnchisq-ctexcl`; their 50-point versions are `d50`, `corrected-d50`, `mnchisq-d50`, `corrected-mnchisq-d50`, `ctexcl-d50`, `mnchisq-ctexcl-d50`, `cov5-d50`, `cov5-ctexcl-d50`, `corrected-cov5-ctexcl-d50` and `corrected-cov5-mnchisq-ctexcl-d50`. Run them from the repository root, one at a time, never beside a 5x fit (step 04 needs about 9.2 GB). Before each, `--plan` with the same switches shows in its line "Fit of step 05" whether the fit is in the cache; if it is not, step 05 fits it (22 to 24 minutes at 10x, 88 to 107 minutes at 5x) or stops while another run holds the lock of that fit.

```sh
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --difference 50                                     # d50 (fit legacy)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage-filter --difference 50                   # corrected-d50 (fit corrected)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --overdispersion --difference 50                    # mnchisq-d50 (fit mnchisq)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage-filter --overdispersion --difference 50  # corrected-mnchisq-d50 (fit corrected-mnchisq)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --exclude-ct-variants --difference 50               # ctexcl-d50 (fit legacy)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --overdispersion --exclude-ct-variants --difference 50  # mnchisq-ctexcl-d50 (fit mnchisq)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage 5 --difference 50                        # cov5-d50 (fit cov5)
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage 5 --exclude-ct-variants --difference 50  # cov5-ctexcl-d50 (fit cov5)
# once methylkit_corrected-cov5.rds is in the cache (the fit of corrected-cov5-ctexcl):
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage-filter --coverage 5 --exclude-ct-variants --difference 50   # corrected-cov5-ctexcl-d50
# the analysis closest to Venkataraman et al. 2022 at 20 points (a new 5x fit by sequence, corrected-cov5-mnchisq, 89 minutes), then its d50 version:
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage-filter --coverage 5 --overdispersion --exclude-ct-variants                   # corrected-cov5-mnchisq-ctexcl
WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R --coverage-filter --coverage 5 --overdispersion --exclude-ct-variants --difference 50  # corrected-cov5-mnchisq-ctexcl-d50
```

Add `--gene-level-go` to a command for steps 10 and 11 (gene level) of that analysis.

#### For the author (left as they are until the author decides):
- Step 04 runs again for every `d50` id, although its outputs are those of the analysis without `d50` under another id (about 9 minutes and 9.2 GB per run). Reusing the step 04 folder of that analysis would need a change of step 04.
- Eleven step 13 rows apply `config$dml_difference` to differences of group means rather than to the DML lists (R3.5-029, R3.5-031, R3.5-033, R3.5-047, D4.2-004, D4.3-005, D4.3-016, D4.9-008, D4.9-010, D4.9-012 and D4.9-014) and so take 50 under `d50`, as every row that reads the threshold does (19 rows); whether these eleven should keep 20 (the threshold of the draft's text) is open (issues.md, I-52).

Step 05 fits every optional analysis by sequence (decision of the author, 2026-11-20), except `ctexcl` alone, which reuses the paper's fit in one piece: the import, the coverage normalisation and the united table are made genome-wide and saved, the tests are then run one part of the genome at a time (each of the ten linkage groups; the mitochondrial genome and the unplaced scaffolds in groups of at most 50,000 CpGs), each part is saved as soon as it is done in `04_output/05_differential_methylation/cache/parts_<fit id>/` (git-ignored), and the SLIM q-values are computed once over all CpGs when every part is in. A stopped run loses only the part it was fitting; the same command started again reuses the saved import and parts and fits the rest. The results are identical to those of a fit in one piece (`tests/check_fit_by_sequence.R`; the paper's own fit is still made in one piece). The fits of `mnchisq`, `corrected` and `corrected-mnchisq`, made in one piece before, are still used. One run at a time per fit: `parts_<fit id>/lock` names the running process, a second run stops while it lives, and takes the lock over after it ended. Before rerunning a fit that was killed, list leftover processes with `pgrep -af 05_differential` and stop them (a forked import or worker can outlive its parent and keep its memory). Saved files of another fingerprint are moved to `parts_<fit id>/stale/` (only the latest such set is kept).

## Comparison of all analyses

`Rscript 01_pipeline/optional/compare_analyses.R` (from the repository root) writes `comparisons/analysis_comparison.csv` (every count, rates per 1,000 tested CpGs, ratios to the paper and a column `venkataraman_2022_settings`) and `comparisons/analysis_comparison.md` from the tables every run writes; it computes nothing new from the data. The table below is a copy of `comparisons/analysis_comparison.md` as of 2026-11-22 (32 analyses: the 16 combinations of the four optional settings at 20 and at 50 points); after a new run, run the script again, and the file is authoritative. Only `legacy` is the paper's analysis. A DML has SLIM q < 0.01 and an absolute difference of more than the threshold (percentage points). Upper filter: the 99.9th percentile coverage filter that the draft describes but did not apply (I-05). Shared: loci in both lists; opposite: different directions in the two lists. V2022: the settings that the Methods of Venkataraman et al. 2022 (BMC Genomics 23:556) state for coverage, the test, variant removal and the threshold; the variant rule counts reads with their full CIGAR (I-08).

| Analysis | Coverage | Upper filter | Overdispersion | C-to-T loci | Threshold | CpGs tested | Ploidy DML (hyper/hypo) | Heat DML (hyper/hypo) | Shared (both hyper/both hypo/opposite) | Opposite % |
|---|---|---|---|---|---|---|---|---|---|---|
| legacy (paper) | 10x | no | none | kept | 20 | 480,579 | 1,038 (442/596) | 2,351 (1,150/1,201) | 200 (28/49/123) | 61.5 |
| ctexcl | 10x | no | none | excluded | 20 | 452,159 | 717 (302/415) | 1,667 (826/841) | 139 (16/26/97) | 69.8 |
| mnchisq | 10x | no | MN, chi-squared | kept | 20 | 480,579 | 14 (9/5) | 185 (95/90) | 4 (1/2/1) | 25.0 |
| mnchisq-ctexcl | 10x | no | MN, chi-squared | excluded | 20 | 452,159 | 12 (7/5) | 142 (75/67) | 4 (1/2/1) | 25.0 |
| corrected | 10x | yes | none | kept | 20 | 454,555 | 952 (430/522) | 2,273 (1,089/1,184) | 159 (27/49/83) | 52.2 |
| corrected-ctexcl | 10x | yes | none | excluded | 20 | 427,349 | 635 (292/343) | 1,596 (772/824) | 100 (15/26/59) | 59.0 |
| corrected-mnchisq | 10x | yes | MN, chi-squared | kept | 20 | 454,555 | 14 (9/5) | 157 (84/73) | 4 (1/2/1) | 25.0 |
| corrected-mnchisq-ctexcl | 10x | yes | MN, chi-squared | excluded | 20 | 427,349 | 12 (7/5) | 118 (65/53) | 4 (1/2/1) | 25.0 |
| cov5 | 5x | no | none | kept | 20 | 2,297,236 | 9,942 (4,408/5,534) | 21,361 (10,090/11,271) | 1,784 (318/485/981) | 55.0 |
| cov5-ctexcl | 5x | no | none | excluded | 20 | 2,153,017 | 7,046 (3,157/3,889) | 15,530 (7,265/8,265) | 1,219 (228/311/680) | 55.8 |
| cov5-mnchisq | 5x | no | MN, chi-squared | kept | 20 | 2,297,236 | 126 (58/68) | 1,191 (511/680) | 27 (10/12/5) | 18.5 |
| cov5-mnchisq-ctexcl | 5x | no | MN, chi-squared | excluded | 20 | 2,153,017 | 96 (48/48) | 923 (382/541) | 19 (8/8/3) | 15.8 |
| corrected-cov5 | 5x | yes | none | kept | 20 | 2,265,759 | 9,791 (4,380/5,411) | 21,200 (9,987/11,213) | 1,711 (312/482/917) | 53.6 |
| corrected-cov5-ctexcl | 5x | yes | none | excluded | 20 | 2,122,969 | 6,902 (3,133/3,769) | 15,376 (7,168/8,208) | 1,148 (222/308/618) | 53.8 |
| corrected-cov5-mnchisq | 5x | yes | MN, chi-squared | kept | 20 | 2,265,759 | 126 (58/68) | 1,133 (494/639) | 27 (10/12/5) | 18.5 |
| corrected-cov5-mnchisq-ctexcl | 5x | yes | MN, chi-squared | excluded | 20 | 2,122,969 | 96 (48/48) | 875 (368/507) | 19 (8/8/3) | 15.8 |
| d50 | 10x | no | none | kept | 50 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| ctexcl-d50 | 10x | no | none | excluded | 50 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| mnchisq-d50 | 10x | no | MN, chi-squared | kept | 50 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| mnchisq-ctexcl-d50 | 10x | no | MN, chi-squared | excluded | 50 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| corrected-d50 | 10x | yes | none | kept | 50 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| corrected-ctexcl-d50 | 10x | yes | none | excluded | 50 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| corrected-mnchisq-d50 | 10x | yes | MN, chi-squared | kept | 50 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| corrected-mnchisq-ctexcl-d50 | 10x | yes | MN, chi-squared | excluded | 50 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| cov5-d50 | 5x | no | none | kept | 50 | 2,297,236 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |
| cov5-ctexcl-d50 | 5x | no | none | excluded | 50 | 2,153,017 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |
| cov5-mnchisq-d50 | 5x | no | MN, chi-squared | kept | 50 | 2,297,236 | 1 (0/1) | 63 (30/33) | 0 (0/0/0) |  |
| cov5-mnchisq-ctexcl-d50 | 5x | no | MN, chi-squared | excluded | 50 | 2,153,017 | 1 (0/1) | 48 (21/27) | 0 (0/0/0) |  |
| corrected-cov5-d50 | 5x | yes | none | kept | 50 | 2,265,759 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |
| corrected-cov5-ctexcl-d50 | 5x | yes | none | excluded | 50 | 2,122,969 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |
| corrected-cov5-mnchisq-d50 | 5x | yes | MN, chi-squared | kept | 50 | 2,265,759 | 1 (0/1) | 62 (30/32) | 0 (0/0/0) |  |
| corrected-cov5-mnchisq-ctexcl-d50 (V2022) | 5x | yes | MN, chi-squared | excluded | 50 | 2,122,969 | 1 (0/1) | 47 (21/26) | 0 (0/0/0) |  |

## What is tracked

An optional analysis commits its results, not its tables of every tested CpG. For every analysis folder except `corrected/`, `05_differential_methylation/tables/united_cpg_<id>.csv.gz` and `diffmeth_<contrast>_<id>.csv.gz` (13 to 22 MB each at 10x, up to 98 MB at 5x) are git-ignored (decision of 2026-11-20, pending the author). Rerun step 05 for the analysis to make them again: `WGBS_ANALYSIS=<id> Rscript 01_pipeline/05_differential_methylation.R` (half a minute with the fit in the cache, else 22 to 24 minutes); steps 06 to 08 stop with that command when the united table is absent. Step 05 writes `united_at_dml_<id>.csv`, the rows of the united table at the DML of either list (tracked), and step 13 reads that, so the optional values in `05_manuscript/not_adopted/` are reproduced from a clone. `gene_table_<id>.csv.gz` (step 08) and `ct_excluded_cpgs_<id>.csv.gz` (step 05, `ctexcl`) stay tracked. The per-CpG tables of `corrected/` predate the settings layout and are read by seven scripts of `06_diagnostics/` (column `reads_corrected_files` of `06_diagnostics/diagnostics.csv`); they stay tracked.

The values of the manuscript for every analysis are in `05_manuscript/not_adopted/manuscript_values_by_analysis.csv`, written by step 13 in an optional run, with "not run" for the analyses that have no folder. The paper's analysis has no rows there (its values are in `05_manuscript/manuscript_values.csv`). Each statement has two values (decision of the author, 2026-11-20). The main columns (`regenerated`, `differs_from_manuscript`, `status`) are the paper recipe: the paper's own computation (`expression_legacy`) applied to the tables of the analysis, through the template `expression_paper_recipe` of the values map, so that a difference from the paper's value comes only from the option; with `A = "legacy"` every template gives the paper's value (`tests/check_paper_recipe_legacy.R`). The columns ending in `_corrected_recipe` are the corrected recipe: the template `expression_optional`, the former corrected-mode expression with `A` in place of `corrected`, which for some rows is a corrected definition of the quantity or another column (applied to the paper's own outputs it differs from the paper's value in 56 of 263 rows). `05_manuscript/README.md` describes both and the rows that needed more than a change of file name.

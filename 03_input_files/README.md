# 03_input_files

Every input the pipeline reads. Nothing in this folder is produced by the pipeline, and nothing in `04_output/` is needed to rebuild it.

For what the pipeline does with these files, see the [root README](../README.md) and [`01_pipeline/README.md`](../01_pipeline/README.md).

## The rule: GitHub first

A file the pipeline needs lives in this folder if it fits in GitHub. The working ceiling is 50 MB per file, after gzip where the format allows it. R reads the gzipped files directly; none needs to be unpacked by hand.

Only files above that ceiling stay on gannet. They are listed in `gannet_manifest.csv` with URL, size, and md5. `01_pipeline/01_fetch_inputs.R` downloads them into `gannet_cache/` (git-ignored) and verifies the checksums. To put the cache somewhere else, for example on a larger drive, set the environment variable `WGBS_CACHE_DIR` before running.

## Contents

| Path | What it is | Origin | Read by |
|---|---|---|---|
| `metadata/` | Library and animal sample sheets | lab records | all steps |
| `elisa/` | Raw plate absorbance and plate map for the global 5-mC ELISA (N = 32) | plate reader exports, 2019-02-06 | step 03 |
| `qpcr/` | Raw Cq values, nine per-gene replicate tables and one consolidated table; in `legacy/`, the printed output of the legacy analysis that step 12 compares with | qPCR instrument exports; the legacy output is a copy from `08_archive/` | step 12 |
| `manuscript/` | The table that ties each quantitative statement of the draft to an R expression over the pipeline outputs, and the table of files to copy into `05_manuscript/`. The manuscript text itself is not in the repository | written during the audit, 2026-11-19 | step 13 |
| `bismark_reports/` | Bismark alignment, deduplication, and summary reports for the 22 libraries, all from the June 2021 alignment | gannet `/panopea/061021-big/` | steps 02, 03 |
| `genome_features/` | Thirteen feature tracks for the Roslin v1 assembly, gzipped | Venkataraman et al. 2022, OSF doi:10.17605/OSF.IO/YGCTB | steps 04, 07, 08 |
| `annotation/` | blastx hits against Swiss-Prot and the transcript-to-GO map; in `legacy/`, outputs of the legacy analysis that the checks and step 11 compare with (annotation tables, REVIGO exports, three Fisher tables of another analysis) | gannet project directory; the REVIGO exports and the Fisher tables are copies from `08_archive/` | steps 08, 09, 11 |
| `gannet_manifest_alignments.csv` | The 22 sorted, deduplicated alignments and their indexes on gannet (URL and size), read only by the optional variant step `01_pipeline/optional/ct_variant_loci.R`; not downloaded by step 01 | gannet `/panopea/061021-big/` | `01_pipeline/optional/ct_variant_loci.R` |
| `gannet_manifest.csv` | The inputs that stay on gannet, with URL, size and md5: the 22 coverage files (9.8 GB) and the reference genome (0.66 GB), 10.5 GB in total | gannet project directory, `cov_files/`, and `/panopea/Cg-roslin/` | steps 01, 04, 05 |
| `gannet_cache/` | Download cache; created on first run; not tracked | | |

Each subfolder has its own README with a file inventory.

No step reads `08_archive/`: a legacy file that a step needs is copied into the fitting subfolder here, and the README of that subfolder names the file it was copied from.

Step 01 writes the size and md5 of every file in this folder to `04_output/01_fetch_inputs/tables/input_inventory.csv`, README files included. After any edit here, rerun step 01 and commit the table with the edit. `.gitattributes` at the repository root tells git not to convert line endings anywhere in this folder, so the checksums are the same on every operating system.

## The coverage files

The 22 files named `<library>_R1.CpG_report.merged_CpG_evidence.cov` are the starting point of the methylation analysis. Each is the output of Bismark `coverage2cytosine --merge_CpG --zero_based`: one row per CpG dinucleotide with reads from both strands combined, and six tab-separated columns with no header:

| Column | Content |
|---|---|
| 1 | Sequence name (RefSeq accession; NC_047559.1 to NC_047568.1 are the ten linkage groups, NC_001276.1 the mitochondrial genome, NW_ unplaced scaffolds) |
| 2 | Start, zero-based position of the C |
| 3 | End, start + 2 |
| 4 | Percent methylation |
| 5 | Count of methylated reads |
| 6 | Count of unmethylated reads |

They are about 445 MB each (9.8 GB in total), which is why they are not in GitHub. How they were produced is documented in [`07_documentation/upstream_processing/`](../07_documentation/upstream_processing/README.md).

## Inputs and the optional analyses

The paper run and the optional analyses (which the paper did not adopt; off by default) read the same inputs; the switches of `run_pipeline.R` change only which rows and settings the steps use. The one input of an optional analysis alone is `gannet_manifest_alignments.csv`, for the C-to-T variant list. The table of the analyses, their switches and their output folders:

| Analysis | Adopted in the paper? | Switch of `run_pipeline.R` | Setting in `run_config.R` (list `optional`) | Id or token | Outputs | What it changes | Issue |
|---|---|---|---|---|---|---|---|
| The paper | yes | none: `Rscript run_pipeline.R` | every setting at its default | `legacy` | `04_output/<step>/`, `05_manuscript/` | at least 10 reads per CpG in each of the 22 libraries, no upper coverage filter, first line of each coverage file dropped, logistic regression without overdispersion correction, DML at q < 0.01 and a difference of more than 20 points; steps 01 to 09 and 11 to 13 | |
| Coverage filter of the Methods | no | `--coverage-filter` | `coverage_filter_of_methods = TRUE` | `corrected` | `04_output/not_adopted/corrected/` | CpGs at or above the 99.9th percentile of coverage in a library removed; every line of each coverage file read | I-05, I-38 |
| 5x standard | no | `--coverage 5` | `coverage_standard = 5L` | `cov5` | `04_output/not_adopted/cov5/` | 5 reads per CpG in each library (2,297,236 CpGs tested instead of 480,579) | Venkataraman et al. 2022 |
| Overdispersion correction | no | `--overdispersion` | `overdispersion_chisq = TRUE` | `mnchisq` | `04_output/not_adopted/mnchisq/` | methylKit `overdispersion = "MN"`, `test = "Chisq"`; no other test is offered | I-06 |
| C-to-T variant exclusion | no | `--exclude-ct-variants` | `exclude_ct_variants = TRUE` | `ctexcl` | `04_output/not_adopted/ctexcl/` | CpGs listed by `01_pipeline/optional/ct_variant_loci.R` removed from the CpG sets of step 04 and, after the test, from the tested CpGs and the DML (fits, p-values and q-values unchanged); the run refuses to start without the list | I-08 |
| Difference of 50 points | no | `--difference 50` | `dml_difference_standard = 50L` | `d50` | `04_output/not_adopted/d50/` | DML at q < 0.01 and an absolute difference of more than 50 points instead of 20 (methylKit `getMethylDiff`, `difference = 50`); applied after the test, so the fit of the analysis without `d50` is reused (author's decision of 2026-11-21) | Venkataraman et al. 2022 |
| q-value 0.05 | no | `--qvalue 0.05` | `dml_qvalue_standard = 0.05` | `q05` | `04_output/not_adopted/q05/` | DML at q < 0.05 instead of 0.01, with the difference threshold of the analysis (more than 20 points, or more than 50 with `d50`; methylKit `getMethylDiff`, `qvalue = 0.05`); applied after the test, so the fit of the analysis without `q05` is reused (author's request of 2026-11-22) | |
| Any combination of the six | no | the switches together | the settings together | the tokens joined by `-` in the order above | `04_output/not_adopted/<id>/`, for example `corrected-cov5-ctexcl/`, `mnchisq-d50/`, `d50-q05/` | each of its parts | |
| Gene-level GO enrichment | no | `--gene-level-go`, alone (on the paper's DML) or with any of the above | `gene_level_go = TRUE` | | `04_output/not_adopted/<id>/10_go_enrichment_gene_level/` and `.../11_semantic_similarity/` (the paper's DML: `<id>` = `legacy`) | step 10 and the gene-level part of step 11 | I-41 |
| Diagnostics | no | `--diagnostics` (all) or `--diagnostics name1,name2` | `diagnostics = "all"` or names | | `06_diagnostics/` | the checks and controls of the audit, listed in `06_diagnostics/diagnostics.csv` | |
| Side by side | no | made by every optional run; the table of all analyses by `Rscript 01_pipeline/optional/compare_analyses.R` | | | `04_output/not_adopted/comparisons/` (`analysis_comparison.md`), `05_manuscript/not_adopted/` | the tables that list analyses, from every analysis folder on disk, the paper first; the manuscript values of every analysis, computed as the paper computed them and with the corrected recipe of the audit ("not run" where absent) | |

# 02_R_functions

Helper functions shared by the pipeline scripts. Every script in `01_pipeline/` sources the whole folder, so a function defined here is available to all steps and is defined in exactly one place.

## Files

| File | Functions | Purpose |
|---|---|---|
| `packages.R` | `load_packages()` | Loads the packages a step needs; stops, naming every missing package, before any work starts. Never installs. |
| `paths.R` | `input_path()`, `output_dir()`, `output_path()`, `step_cache_dir()`, `gannet_cache_dir()`, `manuscript_path()`, `archive_path()` | Builds every path on top of `here::here()`. The numbered folder names appear only in this file. |
| `logging.R` | `log_msg()`, `start_step()`, `finish_step()`, `write_table()` | Timestamped console messages; `finish_step()` writes `session_info.txt` into the step's output folder. |
| `fetch_gannet.R` | `read_gannet_manifest()`, `fetch_gannet_file()` | Returns the local path of a file listed in `03_input_files/gannet_manifest.csv`, downloading it into the cache when absent and verifying size and md5. |
| `sample_sheet.R` | `read_libraries()`, `sample_sets()` | Reads `wgbs_libraries.csv` in legacy column order and defines the sample sets (all, diploid, triploid, control, heat). |
| `plot_theme.R` | `theme_wgbs()`, `save_figure()` | One ggplot2 theme; saves each figure as PNG and PDF. |
| `coverage_io.R` | `read_cov()`, `first_line_lost()`, `cpg_key()`, `key_to_position()`, `upper_coverage_cutoff()` | Reads Bismark merged CpG coverage files; encodes a CpG position as one number for fast matching; computes the upper coverage cutoff the way `methylKit::filterByCoverage()` does. |
| `methylkit_helpers.R` | `read_filtered_counts()`, `scale_counts_like_methylkit()`, `make_methylraw()`, `dml_table()`, `mean_percent_methylation()`, `dml_eligible()` | Applies methylKit's import, filtering, and normalization rules one library at a time, so that 22 libraries fit in memory; defines the eligible CpGs (mean methylation between half the difference threshold and 100 minus that), a screen that holds all but 2 of the legacy DML and is not a bound (issue I-41). |
| `genomic_overlaps.R` | `feature_track_files()`, `read_feature_track()`, `cpg_ranges()`, `overlaps_track()`, `overlap_pairs()`, `gff_attribute()` | Genome feature tracks as `GRanges` and overlap of CpGs with them; replaces `bedtools intersect`. |
| `annotation.R` | `read_genes()`, `read_transcripts()`, `gene_products()`, `read_go_map()`, `true_transcript_ids()`, `read_legacy_go_terms()`, `read_dml_annotation()` | Gene, transcript, and GO tables; repair of the transcript names of the legacy GO map; the correct and the legacy reader of the DML annotation tables (issue I-40). |
| `legacy_go_db.R` | `legacy_go_db_lib()`, `install_legacy_go_db()`, `use_legacy_go_db()`, `go_source_date()` | Installs GO.db 3.13.0 into a side library in the cache and loads it for step 09, without touching the GO.db of the R installation. |
| `go_tests.R` | `go_gene_level_object()`, `go_gene_level_test()`, `quiet_go()` | The gene-level enrichment test (topGO Fisher and goseq Wallenius with Benjamini-Hochberg q-values), defined once and used by step 10 and by the random-CpG control in `06_diagnostics/`. |

## Conventions

- A function does one thing, takes its inputs as arguments, and returns a value; it does not read global variables other than the configuration list created by `run_config.R`.
- Functions do not change the working directory and do not write files except where writing is their stated purpose.
- Coordinates inside R are 1-based closed intervals (`GenomicRanges`). BED and Bismark coverage inputs are converted on read, and converted back on write.

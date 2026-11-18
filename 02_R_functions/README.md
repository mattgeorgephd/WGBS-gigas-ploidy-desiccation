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

Still to come with the steps that need them: `coverage_io.R` (reading Bismark coverage files), `methylkit_helpers.R` (methylKit objects and coverage filters in legacy and corrected form), `genomic_overlaps.R` (feature tracks and overlap counts), `go_helpers.R` (blastx parsing and GO maps).

## Conventions

- A function does one thing, takes its inputs as arguments, and returns a value; it does not read global variables other than the configuration list created by `run_config.R`.
- Functions do not change the working directory and do not write files except where writing is their stated purpose.
- Coordinates inside R are 1-based closed intervals (`GenomicRanges`). BED and Bismark coverage inputs are converted on read, and converted back on write.

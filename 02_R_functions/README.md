# 02_R_functions

Helper functions shared by the pipeline scripts. Every script in `01_pipeline/` sources the whole folder, so a function defined here is available to all steps and is defined in exactly one place.

## Files

To be added with the pipeline (phase 2 of the plan). Planned contents:

| File | Functions for |
|---|---|
| `packages.R` | Loading the packages a step needs and stopping with a clear message if one is missing |
| `paths.R` | Building input, output, and cache paths with `here::here()`; creating a step's output folder |
| `fetch_gannet.R` | Downloading a file listed in `03_input_files/gannet_manifest.csv` into the cache, with resume and md5 verification |
| `sample_sheet.R` | Reading `wgbs_libraries.csv` and returning sample sets (all, by ploidy, by treatment, by cohort) |
| `coverage_io.R` | Reading a Bismark merged CpG coverage file efficiently |
| `methylkit_helpers.R` | Building methylKit objects one library at a time; coverage filters in legacy and corrected form |
| `genomic_overlaps.R` | Reading GFF and BED tracks into one coordinate convention; overlap counts equivalent to `bedtools intersect -u` and `-wb` |
| `go_helpers.R` | Parsing the blastx table; building transcript and gene to GO maps |
| `plot_theme.R` | The ggplot2 theme and the colour assignments for ploidy, treatment, and DML sets |

## Conventions

- A function does one thing, takes its inputs as arguments, and returns a value; it does not read global variables other than the configuration list created by `run_config.R`.
- Functions do not change the working directory and do not write files except where writing is their stated purpose.
- Coordinates inside R are 1-based closed intervals (`GenomicRanges`). BED and Bismark coverage inputs are converted on read, and converted back on write.

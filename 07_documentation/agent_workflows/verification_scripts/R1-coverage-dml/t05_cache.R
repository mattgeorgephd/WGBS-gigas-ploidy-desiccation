# Does a rerun after changing a setting in run_config.R refit, or silently reuse the cached fit?
suppressPackageStartupMessages({library(here); library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml/t05"
source(file.path(repo, "run_config.R"))
for (f in list.files(file.path(repo, "02_R_functions"), pattern = "\\.R$", full.names = TRUE)) source(f)
heads_dir <- file.path(scratch, "heads"); out_root <- file.path(scratch, "out")     # cache written by the earlier test run (min_coverage 10, 99.9th percentile)
fetch_gannet_file <- function(file, config, ...) file.path(heads_dir, file)
output_dir <- function(step, sub = NULL) { d <- if (is.null(sub)) file.path(out_root, step) else file.path(out_root, step, sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
output_path <- function(step, sub, file) file.path(output_dir(step, sub), file)
step_cache_dir <- function(step) output_dir(step, "cache")
ex <- parse(file.path(repo, "01_pipeline/05_differential_methylation.R"))
for (e in ex) if (is.call(e) && as.character(e[[1]]) == "<-" && as.character(e[[2]]) %in% c("step", "lib", "treat_ploidy", "treat_heat", "fit_mode", "write_mode")) eval(e, globalenv())
config$min_coverage <- 5; config$upper_coverage_percentile <- 99
fit <- fit_mode("corrected")
cat("settings now: min_coverage =", config$min_coverage, ", percentile =", config$upper_coverage_percentile, "\n")
cat("CpGs in the fit returned:", nrow(fit$united), " ; minimum raw-equivalent coverage threshold recorded in fit$filtering: n_cpg_min_coverage of library 1 =", fit$filtering$n_cpg_min_coverage[1],
    "; upper cutoff of library 1 =", fit$filtering$upper_cutoff[1], "\n")

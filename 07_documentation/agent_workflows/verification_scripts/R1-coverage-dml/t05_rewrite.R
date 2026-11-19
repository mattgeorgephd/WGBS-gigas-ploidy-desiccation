# Re-run the committed write stage of step 05 on the cached fits; outputs go to scratch; then compare with the committed tables.
suppressPackageStartupMessages({library(here); library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml/t05"
source(file.path(repo, "run_config.R"))
for (f in list.files(file.path(repo, "02_R_functions"), pattern = "\\.R$", full.names = TRUE)) source(f)
out_root <- file.path(scratch, "rewrite"); unlink(out_root, recursive = TRUE)
fetch_gannet_file <- function(...) stop("must not read coverage files in this test")
output_dir <- function(step, sub = NULL) { d <- if (is.null(sub)) file.path(out_root, step) else file.path(out_root, step, sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
output_path <- function(step, sub, file) file.path(output_dir(step, sub), file)
step_cache_dir <- function(step) file.path(repo, "04_output", step, "cache")      # read only: fit_mode returns the cached fit when the file exists
ex <- parse(file.path(repo, "01_pipeline/05_differential_methylation.R"))
for (e in ex) if (is.call(e) && as.character(e[[1]]) == "<-" && as.character(e[[2]]) %in% c("step", "lib", "treat_ploidy", "treat_heat", "fit_mode", "write_mode")) eval(e, globalenv())
headline <- rbindlist(lapply(config$modes, function(mode) {
  stopifnot(file.exists(file.path(step_cache_dir(step), sprintf("methylkit_%s.rds", mode))))
  fit <- fit_mode(mode)
  cat(sprintf("[%s] cached fit: sample ids %s ; ploidy treatment %s ; heat treatment %s ; methylKit object rows %d\n", mode, paste(fit$united@sample.ids, collapse = ","),
              paste(fit$diffmeth$ploidy@treatment, collapse = ""), paste(fit$diffmeth$heat@treatment, collapse = ""), nrow(fit$united)))
  out <- write_mode(fit); rm(fit); gc(verbose = FALSE); out }))
write_table(headline, output_path(step, "tables", "dml_counts_by_mode.csv"))

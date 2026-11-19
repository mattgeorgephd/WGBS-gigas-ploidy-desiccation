# Equivalence test: step 05's own functions (fit_mode, write_mode, evaluated from the real script)
# against methylKit's native route, on regional extracts of the 22 coverage files.
suppressPackageStartupMessages({library(here); library(data.table); library(methylKit)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml/t05"
stopifnot(normalizePath(here()) == repo)
source(file.path(repo, "run_config.R"))
for (f in list.files(file.path(repo, "02_R_functions"), pattern = "\\.R$", full.names = TRUE)) source(f)
# redirect every path that the step writes or downloads to the scratch directory
heads_dir <- Sys.getenv("T05_HEADS", file.path(scratch, "heads"))
out_root <- Sys.getenv("T05_OUT", file.path(scratch, "out"))
fetch_gannet_file <- function(file, config, ...) file.path(heads_dir, file)
output_dir <- function(step, sub = NULL) { d <- if (is.null(sub)) file.path(out_root, step) else file.path(out_root, step, sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
output_path <- function(step, sub, file) file.path(output_dir(step, sub), file)
step_cache_dir <- function(step) output_dir(step, "cache")
unlink(out_root, recursive = TRUE)
ex <- parse(file.path(repo, "01_pipeline/05_differential_methylation.R"))
keep <- vapply(ex, function(e) is.call(e) && as.character(e[[1]]) == "<-" && as.character(e[[2]]) %in% c("step", "lib", "treat_ploidy", "treat_heat", "fit_mode", "write_mode"), logical(1))
for (e in ex[keep]) eval(e, globalenv())
cat("evaluated from the step script:", paste(vapply(ex[keep], function(e) as.character(e[[2]]), ""), collapse = ", "), "\n")
config$cores <- 2L
files <- as.list(file.path(heads_dir, lib$cov_file))
res <- list()
for (mode in config$modes) {
  cat("\n================ mode", mode, "================\n")
  fit <- fit_mode(mode)
  hl <- write_mode(fit)
  print(hl)
  # methylKit's native route
  if (mode == "legacy") {
    raw <- suppressMessages(methRead(files, sample.id = as.list(as.character(1:22)), assembly = "Roslin", treatment = treat_ploidy,
                                     pipeline = "bismarkCoverage", mincov = 2))
    flt <- filterByCoverage(raw, lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
  } else {
    raw <- suppressMessages(methRead(files, sample.id = as.list(as.character(1:22)), assembly = "Roslin", treatment = treat_ploidy,
                                     pipeline = "bismarkCoverage", mincov = 2, header = FALSE))
    flt <- filterByCoverage(raw, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
  }
  nrm <- normalizeCoverage(flt)
  uni <- suppressMessages(unite(nrm, destrand = FALSE))
  a <- getData(uni); b <- getData(fit$united)
  cat("native united rows:", nrow(a), " pipeline united rows:", nrow(b), "\n")
  cat("identical data frames (all.equal):", isTRUE(all.equal(a, b, check.attributes = FALSE)), "\n")
  cat("identical positions:", identical(a$chr, b$chr) && identical(as.integer(a$start), as.integer(b$start)) && identical(as.integer(a$end), as.integer(b$end)), "\n")
  cat("max abs count difference:", max(abs(as.matrix(a[, -(1:4)]) - as.matrix(b[, -(1:4)]))), "\n")
  cat("sample ids:", identical(uni@sample.ids, fit$united@sample.ids), " treatment:", identical(as.numeric(uni@treatment), as.numeric(fit$united@treatment)), "\n")
  cat("column classes native:", paste(unique(sapply(a, class)), collapse = ","), " pipeline:", paste(unique(sapply(b, class)), collapse = ","), "\n")
  # per-library medians / scaling
  med_native <- sapply(flt, function(x) median(x$coverage))
  cat("scaling factors equal:", isTRUE(all.equal(unname(max(med_native) / med_native), fit$filtering$scaling_factor)), "\n")
  cat("kept per library equal:", identical(as.integer(sapply(flt, nrow)), as.integer(fit$filtering$n_cpg_kept)), "\n")
  # native tests with the legacy covariate data frames
  covariatePLOIDY <- data.frame(PLOIDY = c(rep("D", 5), rep("T", 5), rep("D", 6), rep("T", 6)))
  covariateHEAT <- data.frame(HEAT = c(rep("yes", 10), rep("no", 12)))
  dp <- suppressMessages(calculateDiffMeth(uni, covariates = covariateHEAT, mc.cores = 1))
  uni_h <- suppressMessages(unite(new("methylRawList", nrm@.Data, treatment = treat_heat), destrand = FALSE))
  dh <- suppressMessages(calculateDiffMeth(uni_h, covariates = covariatePLOIDY, mc.cores = 1))
  for (ct in c("ploidy", "heat")) {
    n <- getData(if (ct == "ploidy") dp else dh); p <- getData(fit$diffmeth[[ct]])
    cat(sprintf("[%s] p identical: %s (max abs diff %.3g) ; q identical: %s (max rel diff %.3g) ; meth.diff identical: %s\n", ct,
                identical(n$pvalue, p$pvalue), max(abs(n$pvalue - p$pvalue)), identical(n$qvalue, p$qvalue), max(abs(n$qvalue - p$qvalue) / pmax(n$qvalue, 1e-300)),
                identical(n$meth.diff, p$meth.diff)))
    for (th in c(5, 10, 20)) {
      cat(sprintf("   DML at %d: native %d, pipeline %d\n", th, nrow(getMethylDiff(if (ct == "ploidy") dp else dh, difference = th, qvalue = 0.01)),
                  sum(p$qvalue < 0.01 & abs(p$meth.diff) > th)))
    }
  }
  res[[mode]] <- list(fit = fit, native = a)
}
saveRDS(res, file.path(scratch, "equiv_res.rds"))

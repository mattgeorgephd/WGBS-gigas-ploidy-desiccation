# Equivalence test on one linkage group (NC_047567.1, all 22 libraries):
#   A. the legacy methylKit calls of 08_archive/code/1_coverage_analysis.R, run as written
#   B. the code of 01_pipeline/04_cpg_coverage.R (sections 1, 3 and 5), copied verbatim, on the same files
suppressPackageStartupMessages({library(methylKit); library(data.table)})
args <- commandArgs(trailingOnly = TRUE); sub_dir <- args[1]
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(repo, "02_R_functions", "coverage_io.R"))          # read_cov, cpg_key, first_line_lost, upper_coverage_cutoff
lib <- read.csv(file.path(repo, "03_input_files", "metadata", "wgbs_libraries.csv"), stringsAsFactors = FALSE)
lib <- lib[order(lib$legacy_sample_index), ]; rownames(lib) <- NULL
n_lib <- nrow(lib); files <- file.path(sub_dir, lib$cov_file); stopifnot(all(file.exists(files)))
config <- list(import_min_coverage = 2, upper_coverage_percentile = 99.9, min_coverage = 10, modes = c("legacy", "corrected"))

## ---------------- A. legacy methylKit calls -------------------------------------------------
mk <- function(x, legacy = TRUE, import_mincov = 1) {
  if (legacy) {
    processedFiles <- methylKit::methRead(as.list(files), sample.id = as.list(lib$library_id), assembly = "Roslin",
                                          treatment = as.integer(lib$ploidy == "triploid"), pipeline = "bismarkCoverage", mincov = 1)
    # the misnamed arguments of the legacy script, verbatim
    flt <- methylKit::filterByCoverage(processedFiles, lo.count = x, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
  } else {
    processedFiles <- methylKit::methRead(as.list(files), sample.id = as.list(lib$library_id), assembly = "Roslin", header = FALSE,
                                          treatment = as.integer(lib$ploidy == "triploid"), pipeline = "bismarkCoverage", mincov = import_mincov)
    flt <- methylKit::filterByCoverage(processedFiles, lo.count = x, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
  }
  nrm <- methylKit::normalizeCoverage(flt)
  uni <- methylKit::unite(nrm, destrand = FALSE, mc.cores = 1)
  out <- capture.output(methylKit::getCorrelation(uni, plot = FALSE))          # the call that is commented out in the legacy script
  d <- getData(uni)
  mm <- as.matrix(d[, uni@numCs.index]) / (as.matrix(d[, uni@numCs.index]) + as.matrix(d[, uni@numTs.index]))
  r <- cor(mm, method = "pearson", use = "na.or.complete"); dimnames(r) <- list(lib$library_id, lib$library_id)
  rs <- cor(mm, method = "spearman", use = "na.or.complete")
  pdf(NULL); pc <- methylKit::PCASamples(uni, obj.return = TRUE); dev.off()
  list(n = nrow(uni), r = r, rs = rs, pc = pc, printed = out, rows_per_lib = sapply(processedFiles, nrow),
       key = paste(d$chr, d$start), pm = 100 * mm)
}

## ---------------- B. step 04 code, verbatim -------------------------------------------------
chr_levels <- character(0); per_lib <- vector("list", n_lib); universe <- numeric(0); first_key <- numeric(n_lib); cutoffs <- numeric(n_lib)
for (j in seq_len(n_lib)) {
  d <- read_cov(files[j])
  d[, coverage := nC + nT]
  chr_levels <- union(chr_levels, unique(d$chr))
  key <- cpg_key(d$chr, d$start, chr_levels)
  first_key[j] <- key[1]
  cutoffs[j] <- upper_coverage_cutoff(d$coverage, config$import_min_coverage, config$upper_coverage_percentile)
  per_lib[[j]] <- list(key = key, coverage = d$coverage, nC = d$nC)
  universe <- unique(c(universe, key))
}
universe <- sort(universe); n_cpg <- length(universe)
cov_mat <- matrix(0L, n_cpg, n_lib); nC_mat <- matrix(0L, n_cpg, n_lib)
for (j in seq_len(n_lib)) {
  i <- match(per_lib[[j]]$key, universe)
  cov_mat[i, j] <- per_lib[[j]]$coverage; nC_mat[i, j] <- per_lib[[j]]$nC
}
first_row <- match(first_key, universe)
covered <- function(j, x, mode) {
  cv <- cov_mat[, j]
  ok <- if (mode == "legacy") cv >= x else cv >= x & cv < cutoffs[j]
  if (first_line_lost(mode)) ok[first_row[j]] <- FALSE
  ok
}
percent_methylation <- function(rows, x, mode) {
  pm <- matrix(NA_real_, length(rows), n_lib)
  med <- vapply(seq_len(n_lib), function(j) median(cov_mat[covered(j, x, mode), j]), numeric(1))
  sc <- max(med) / med
  for (j in seq_len(n_lib)) {
    cv <- cov_mat[rows, j]; nc <- nC_mat[rows, j]
    new_cov <- round(sc[j] * cv)
    new_c <- round(new_cov * nc / cv); new_t <- round(new_cov * (cv - nc) / cv)
    pm[, j] <- 100 * new_c / (new_c + new_t)
  }
  colnames(pm) <- lib$library_id
  pm
}
step04 <- function(x, mode) {
  in_all <- rep(TRUE, n_cpg)
  for (j in seq_len(n_lib)) in_all <- in_all & covered(j, x, mode)
  rows <- which(in_all)
  pm <- percent_methylation(rows, x, mode)
  pm <- pm[complete.cases(pm), , drop = FALSE]
  r <- cor(pm)
  sds <- sqrt(rowSums((pm - rowMeans(pm))^2) / (n_lib - 1))
  pc <- prcomp(t(pm[sds > quantile(sds, 0.5), , drop = FALSE]), center = TRUE, scale. = TRUE)
  list(n = nrow(pm), r = r, pc = pc, pm = pm, key = universe[rows])
}

## ---------------- comparison ---------------------------------------------------------------
cmp <- function(a, b, label) {
  off <- upper.tri(a$r)
  ve_a <- 100 * a$pc$sdev^2 / sum(a$pc$sdev^2); ve_b <- 100 * b$pc$sdev^2 / sum(b$pc$sdev^2)
  sdiff <- function(k) min(max(abs(a$pc$x[, k] - b$pc$x[, k])), max(abs(a$pc$x[, k] + b$pc$x[, k])))
  cat(sprintf("%-34s CpGs: methylKit %d, step 04 %d | Pearson range methylKit %.4f..%.4f, step 04 %.4f..%.4f | max |dr| %.2e | PC1 %.4f vs %.4f, PC2 %.4f vs %.4f | max |score diff| PC1 %.3g PC2 %.3g (score SD %.0f) | PCA loci %d vs %d\n",
              label, a$n, b$n, min(a$r[off]), max(a$r[off]), min(b$r[off]), max(b$r[off]), max(abs(a$r - b$r)),
              ve_a[1], ve_b[1], ve_a[2], ve_b[2], sdiff(1), sdiff(2), a$pc$sdev[1], ncol(a$pc$rotation) * 0 + nrow(a$pc$rotation), nrow(b$pc$rotation)))
  # how many CpG x library values differ after normalization
  if (a$n == b$n) {
    d <- abs(a$pm - b$pm) > 1e-9
    cat(sprintf("%-34s percent-methylation values that differ between the two: %d of %d (%.4f%%), largest difference %.3f points; CpGs affected %d\n",
                "", sum(d), length(d), 100 * mean(d), max(abs(a$pm - b$pm)), sum(rowSums(d) > 0)))
  }
}
for (x in c(1, 10)) {
  a <- mk(x, legacy = TRUE); b <- step04(x, "legacy")
  if (x == 1) { cat("rows read by methRead per library minus rows in file:", paste(unique(a$rows_per_lib - sapply(per_lib, function(p) length(p$key))), collapse = ","), "\n")
                cat("getCorrelation printed", length(a$printed), "lines; first:", a$printed[1], "\n") }
  cmp(a, b, sprintf("legacy, %dx", x))
  off <- upper.tri(a$rs); cat(sprintf("%-34s Spearman range on the same methylKit matrix: %.4f..%.4f\n", "", min(a$rs[off]), max(a$rs[off])))
}
# corrected mode: every line read, upper filter applied. Two readings of "as the Methods describe":
a <- mk(10, legacy = FALSE, import_mincov = 2); b <- step04(10, "corrected"); cmp(a, b, "corrected, 10x (import >= 2 reads)")
a <- mk(1, legacy = FALSE, import_mincov = 1);  b <- step04(1, "corrected");  cmp(a, b, "corrected, 1x (import >= 1 read)")
cat("upper cutoffs, percentile over rows >= 2 reads (step 04):", paste(round(cutoffs, 2), collapse = " "), "\n")
cat("upper cutoffs, percentile over rows >= 1 read  (methylKit after mincov = 1):",
    paste(round(sapply(per_lib, function(p) quantile(p$coverage, 0.999)), 2), collapse = " "), "\n")
cat("R", as.character(getRversion()), "methylKit", as.character(packageVersion("methylKit")), "\n")

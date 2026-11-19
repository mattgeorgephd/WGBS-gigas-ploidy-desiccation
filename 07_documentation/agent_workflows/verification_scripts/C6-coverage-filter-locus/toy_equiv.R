# Pipeline functions (sourced read-only from the repository) against methylKit's own route, on toy libraries
# built to hit the edge cases: ties exactly at an integer cutoff, a fractional cutoff, rows with 1 read
# (below the import minimum), and a first line with >= 10 reads.
suppressPackageStartupMessages({library(methylKit); library(data.table)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(repo, "02_R_functions", "coverage_io.R"))
source(file.path(repo, "02_R_functions", "methylkit_helpers.R"))
config <- list(import_min_coverage = 2, min_coverage = 10, upper_coverage_percentile = 99.9)
dir.create("toy", showWarnings = FALSE)
n <- 60000; n_lib <- 6
pos <- data.table(chr = rep(c("chrB", "chrA", "chrC"), each = n / 3), start = rep(seq(50L, by = 9L, length.out = n / 3), 3))
files <- character(n_lib)
set.seed(42)
for (j in seq_len(n_lib)) {
  cov <- pmax(1L, rnbinom(n, size = 3, mu = 9 + j))
  hot <- sample(n, 300); cov[hot] <- cov[hot] + rgeom(300, 0.02)       # heavy upper tail, with ties
  if (j == 2) cov[1] <- 25L                                             # first line with >= 10 reads
  keep <- runif(n) < 0.97; keep[1] <- TRUE                              # some CpGs absent from a library
  nC <- rbinom(n, cov, 0.25)
  d <- data.table(pos, end = pos$start + 2L, pct = sprintf("%.6f", 100 * nC / cov), nC = nC, nT = cov - nC)[keep]
  d <- d[sample(nrow(d))]                                               # unsorted, like the real files
  if (j == 2) d <- rbind(d[chr == "chrB" & start == 50L], d[!(chr == "chrB" & start == 50L)])
  files[j] <- file.path("toy", sprintf("lib%d.cov", j))
  fwrite(d, files[j], sep = "\t", col.names = FALSE, quote = FALSE)
}
treat <- c(0, 0, 0, 1, 1, 1)

run_methylkit <- function(header, upper) {
  raw <- suppressMessages(methRead(as.list(files), sample.id = as.list(paste0("s", 1:n_lib)), assembly = "x",
                                   treatment = treat, pipeline = "bismarkCoverage", mincov = 2, header = header))
  cut <- sapply(raw, function(x) quantile(x$coverage, 99.9 / 100))
  f <- if (upper) filterByCoverage(raw, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
       else filterByCoverage(raw, lo.count = 10, lo.perc = NULL, high.count = NULL, high.perc = 99.9)
  u <- suppressMessages(unite(normalizeCoverage(f), destrand = FALSE))
  list(united = as.data.table(getData(u)), cutoffs = unname(cut), kept = sapply(f, nrow))
}
run_pipeline <- function(mode) {
  kept <- list(); rec <- list(); chr_levels <- character(0); common <- NULL
  for (i in seq_len(n_lib)) {
    r <- read_filtered_counts(files[i], mode, config, chr_levels)
    chr_levels <- r$chr_levels
    common <- if (is.null(common)) r$counts$pos_key else common[common %in% r$counts$pos_key]
    kept[[i]] <- r$counts; rec[[i]] <- r$record
  }
  filtering <- rbindlist(rec); filtering[, scaling_factor := max(median_coverage) / median_coverage]
  p <- key_to_position(common, chr_levels)
  objects <- lapply(seq_len(n_lib), function(i) {
    k <- kept[[i]][match(common, pos_key)]
    make_methylraw(p$chr, p$start, scale_counts_like_methylkit(k$coverage, k$nC, k$nT, filtering$scaling_factor[i]), paste0("s", i))
  })
  u <- suppressMessages(unite(new("methylRawList", objects, treatment = treat), destrand = FALSE))
  list(united = as.data.table(getData(u)), cutoffs = filtering$upper_cutoff, kept = filtering$n_cpg_kept, filtering = filtering)
}
cmp <- function(a, b, label) {
  setkey(a$united, chr, start); setkey(b$united, chr, start)
  cat(sprintf("%-46s methylKit %6d rows, pipeline %6d rows; identical united table: %s; kept per library identical: %s\n",
              label, nrow(a$united), nrow(b$united),
              isTRUE(all.equal(a$united, b$united, check.attributes = FALSE)), identical(as.integer(a$kept), as.integer(b$kept))))
}
mk_leg <- run_methylkit(header = TRUE,  upper = FALSE); pl_leg <- run_pipeline("legacy")
mk_cor <- run_methylkit(header = FALSE, upper = TRUE);  pl_cor <- run_pipeline("corrected")
cmp(mk_leg, pl_leg, "legacy (header = TRUE, high.perc):")
cmp(mk_cor, pl_cor, "corrected (header = FALSE, hi.perc = 99.9):")
cat("methylKit cutoffs (corrected route):", mk_cor$cutoffs, "\n")
cat("pipeline  cutoffs:                  ", pl_cor$cutoffs, "\n")
cat("cutoffs identical:", identical(as.numeric(mk_cor$cutoffs), as.numeric(pl_cor$cutoffs)), "\n")
# the variant a user might run instead: methRead default header (TRUE) with the corrected filter
mk_cor_hdr <- run_methylkit(header = TRUE, upper = TRUE)
cat("methylKit cutoffs if line 1 is also lost (header = TRUE + hi.perc):", mk_cor_hdr$cutoffs, " rows:", nrow(mk_cor_hdr$united), "\n")
# ties: how many rows sit exactly at an integer cutoff, and are they removed by both?
for (i in seq_len(n_lib)) {
  d <- read_cov(files[i]); cv <- d$nC + d$nT; ct <- pl_cor$cutoffs[i]
  cat(sprintf("lib%d cutoff %-8g integer: %-5s rows with coverage == cutoff: %d; rows >= cutoff: %d; removed by pipeline: %d\n",
              i, ct, ct == round(ct), sum(cv == ct), sum(cv >= ct & cv >= 10), pl_cor$filtering$n_removed_by_upper[i]))
}
cat("first line of lib2 (25 reads) present in corrected united:", nrow(pl_cor$united[chr == "chrB" & start == 50L]),
    " in legacy united:", nrow(pl_leg$united[chr == "chrB" & start == 50L]), "\n")

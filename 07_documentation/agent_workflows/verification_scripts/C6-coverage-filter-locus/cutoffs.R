suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
lc <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
first <- fread("first_lines.tsv", header = FALSE, col.names = c("library_id", "chr", "start", "end", "pct", "nC", "nT"))
first[, cov := nC + nT]
# type-7 quantile from a histogram, by hand (no expansion): h = (n - 1) p + 1; value = x[floor(h)] + (h - floor(h)) (x[floor(h) + 1] - x[floor(h)])
q7_hist <- function(value, count, p) {
  o <- order(value); value <- value[o]; count <- count[o]; n <- sum(count); cum <- cumsum(count)
  h <- (n - 1) * p + 1; lo <- floor(h); hi <- ceiling(h)
  xlo <- value[which(cum >= lo)[1]]; xhi <- value[which(cum >= hi)[1]]
  xlo + (h - lo) * (xhi - xlo)
}
rows <- list()
for (lib in lc$library_id) {
  h <- fread(sprintf("pass/hist_%s.tsv", lib), header = FALSE, col.names = c("cov", "n"))
  n_lines <- sum(h$n)
  f <- first[library_id == lib]
  base2 <- h[cov >= 2]
  v <- rep(base2$cov, base2$n)                                   # what methylKit holds as data$coverage after methRead(mincov = 2, header = FALSE)
  q_r  <- unname(quantile(v, 99.9 / 100))                        # methylKit's expression: quantile(data$coverage, hi.perc/100)
  q_hand <- q7_hist(base2$cov, base2$n, 99.9 / 100)
  # legacy import (line 1 lost) + corrected argument name: the cutoff a user would get by fixing only the typo
  hl <- copy(h); hl[cov == f$cov, n := n - 1L]; bl <- hl[cov >= 2]
  q_line1lost <- unname(quantile(rep(bl$cov, bl$n), 99.9 / 100))
  q_base1  <- unname(quantile(rep(h$cov, h$n), 99.9 / 100))                       # methRead(mincov = 1), as in 1_coverage_analysis.R
  b10 <- h[cov >= 10]; q_base10 <- unname(quantile(rep(b10$cov, b10$n), 99.9 / 100))  # methRead default mincov = 10
  types <- sapply(1:9, function(t) unname(quantile(v, 0.999, type = t)))
  rows[[lib]] <- data.table(library_id = lib, n_lines = n_lines, n_ge2 = sum(base2$n), n_ge10 = sum(b10$n),
    median_all = median(rep(h$cov, h$n)), mean_all = sum(as.numeric(h$cov) * h$n) / n_lines, max_cov = max(h$cov),
    cutoff_mine = q_r, cutoff_hand = q_hand, cutoff_line1lost = q_line1lost, cutoff_base1 = q_base1, cutoff_base10 = q_base10,
    n_at_or_above = h[cov >= q_r, sum(n)], n_at_or_above_ge10 = h[cov >= q_r & cov >= 10, sum(n)],
    pct_at_or_above_of_ge2 = 100 * h[cov >= q_r, sum(n)] / sum(base2$n),
    n_equal_cutoff = h[cov == q_r, sum(n)], q_type_range = sprintf("%g-%g", min(types), max(types)),
    q_p0999_literal = unname(quantile(v, 0.999)))
  rm(v); gc(verbose = FALSE)
}
r <- merge(rbindlist(rows), lc[, .(library_id, cohort, ploidy, n_cpg_1x, n_cpg_10x, mean_coverage, median_coverage, upper_cutoff, n_cpg_at_or_above_cutoff)], by = "library_id", sort = FALSE)
r[, `:=`(lines_ok = n_lines == n_cpg_1x, n10_ok = n_ge10 == n_cpg_10x, median_ok = median_all == median_coverage,
         mean_ok = abs(mean_all - mean_coverage) < 1e-9, cutoff_ok = cutoff_mine == upper_cutoff, hand_ok = abs(cutoff_hand - upper_cutoff) < 1e-9,
         n_above_ok = n_at_or_above == n_cpg_at_or_above_cutoff)]
options(width = 250)
print(r[, .(library_id, n_lines, n_ge2, n_ge10, median_all, max_cov, cutoff_mine, cutoff_hand, upper_cutoff, cutoff_line1lost, cutoff_base1, cutoff_base10, n_at_or_above, n_cpg_at_or_above_cutoff, n_equal_cutoff, q_type_range)])
cat("\nall line counts equal n_cpg_1x:", all(r$lines_ok), "\nall n >= 10 equal n_cpg_10x:", all(r$n10_ok), "\nall medians equal:", all(r$median_ok),
    "\nall means equal:", all(r$mean_ok), "\nall cutoffs identical (==) to library_coverage.csv:", all(r$cutoff_ok),
    "\nhand-computed type-7 agrees:", all(r$hand_ok), "\nall counts at or above cutoff equal:", all(r$n_above_ok), "\n")
cat("share of CpGs with >= 2 reads at or above the cutoff, per library: range", sprintf("%.4f%% to %.4f%%", min(r$pct_at_or_above_of_ge2), max(r$pct_at_or_above_of_ge2)), "\n")
cat("cutoff with 0.999 literal vs 99.9/100 identical:", all(r$q_p0999_literal == r$cutoff_mine), "; max abs diff", max(abs(r$q_p0999_literal - r$cutoff_mine)), "\n")
cat("cutoff unchanged if line 1 is lost (ceiling compared):", all(ceiling(r$cutoff_line1lost) == ceiling(r$cutoff_mine)), "\n")
cat("ratio cutoff / median of all covered CpGs: range", range(r$cutoff_mine / r$median_all), "\n")
fwrite(r, "cutoffs_mine.csv")

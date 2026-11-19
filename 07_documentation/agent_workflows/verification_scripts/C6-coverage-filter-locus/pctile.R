suppressPackageStartupMessages(library(data.table))
pl <- fread("locus_per_library.csv")
z <- readRDS("cand_matrix.rds"); cov <- z$cov
out <- rbindlist(lapply(pl$library_id, function(l) {
  h <- fread(sprintf("pass/hist_%s.tsv", l), header = FALSE, col.names = c("cov", "n"))[cov >= 2][order(cov)]
  m <- pl[library_id == l, median_cov]
  data.table(library_id = l, state = ifelse(pl[library_id == l, n_over] > 0, "high", "low"), locus_median_cov = m,
             pct_of_cpgs_ge2_below = round(100 * h[cov < m, sum(n)] / sum(h$n), 2),
             pct_of_tested_cpgs_below = round(100 * mean(cov[, l] < m), 2), median_cov_tested = median(cov[, l]))
}))
print(out[order(state, pct_of_cpgs_ge2_below)])
cat("low libraries: locus median sits at percentile", range(out[state == "low", pct_of_cpgs_ge2_below]), "of CpGs with >= 2 reads, and at percentile",
    range(out[state == "low", pct_of_tested_cpgs_below]), "of the 480,579 tested CpGs\n")
cat("high libraries: percentile", range(out[state == "high", pct_of_cpgs_ge2_below]), "of CpGs with >= 2 reads;", range(out[state == "high", pct_of_tested_cpgs_below]), "of tested CpGs\n")

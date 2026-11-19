# Task 4 in detail: which CpGs make the counts of dml_sensitivity.R (no margin) differ from those of library_methylation_level.R.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
dd <- fread(file.path(S, "out", "relabel_cpgs_where_no_margin_differs_legacy.csv")); a <- fread(file.path(S, "out", "relabel_all_groupings_legacy.csv"))
options(width = 230)
print(dd[, .(cases = .N, distinct_cpgs = uniqueN(paste(chr, start)), groupings = uniqueN(grouping), all_exact_ties = all(exact_tie), in_observed_or_mirror = sum(observed_or_mirror)), by = comparison])
# exactness once more, with integers: 100 * (c_b * n_a - c_a * n_b) == +-20 * n_a * n_b
dd[, tie_check := abs(5 * (c_b * n_a - c_a * n_b)) == n_a * n_b]; cat("integer check of all cases:", all(dd$tie_check), "\n")
cat("\nthe cases in the observed grouping or its mirror image:\n"); print(dd[observed_or_mirror == TRUE, .(comparison = substr(comparison, 1, 32), libraries_in_group_a, chr, start, c_a, n_a, c_b, n_b, d_float)])
for (nm in unique(dd$comparison)) {
  x <- a[comparison == nm][order(n_exact)]
  cat("\n", nm, ": the grouping with the smallest count: exact ", x$n_exact[1], ", without margin ", x$n_S[1], " (", x$libraries_in_group_a[1], ")\n", sep = "")
  print(dd[comparison == nm & grouping == x$grouping[1], .(chr, start, c_a, n_a, c_b, n_b, pct_a = round(100 * c_a / n_a, 4), pct_b = round(100 * c_b / n_b, 4), d_float)])
  xs <- sort(x$n_S); xe <- sort(x$n_exact); n <- length(xs)
  cat("  sorted counts at the two middle positions: exact ", xe[n / 2], " and ", xe[n / 2 + 1], "; without margin ", xs[n / 2], " and ", xs[n / 2 + 1], "\n", sep = "")
  cat("  sorted positions at which the two sets of counts differ: ", sum(xs != xe), " of ", n, "\n", sep = "")
}
# pairs of percentages that make the ties, and what floating point gives
dd[, `:=`(pa = round(100 * c_a / n_a, 6), pb = round(100 * c_b / n_b, 6))]
cat("\nmost frequent pairs of pooled percentages among the cases (group a, group b):\n"); print(dd[, .N, by = .(pa, pb, d_float)][order(-N)][1:12])

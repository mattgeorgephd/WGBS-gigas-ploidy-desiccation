# Exact two-sided p-values of Spearman's test for the 16 rows of level_correlations.csv, from the published table of the
# exact null distribution for n <= 22 (CRAN package pspearman 0.3-1, file R/spearman.data.R, read as data: the file is
# parsed and only accepted if it is one assignment of nested list()/c() calls of numbers; no code of the package is run).
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE); data_file <- args[1]; S <- args[2]
e <- parse(data_file, keep.source = FALSE); stopifnot(length(e) == 1, identical(e[[1]][[1]], as.name("<-")), identical(e[[1]][[2]], as.name("spearman.list")))
only_data <- function(x) if (is.call(x)) (identical(x[[1]], as.name("list")) || identical(x[[1]], as.name("c"))) && all(vapply(as.list(x)[-1], only_data, logical(1))) else is.numeric(x)
stopifnot(only_data(e[[1]][[3]]))
spearman.list <- eval(e[[1]][[3]], baseenv())
exact_p <- function(q, n) {                      # two-sided, as stats::cor.test defines it: twice the smaller tail
  m <- n * (n^2 - 1) / 6 + 1; g <- spearman.list[[n]]; freq <- numeric(m); freq[seq_along(g)] <- g; freq[m + 1 - seq_along(g)] <- g
  stopifnot(abs(sum(freq) / factorial(n) - 1) < 1e-12)
  k <- q / 2 + 1; lower <- sum(freq[1:k]) / sum(freq); upper <- sum(freq[k:m]) / sum(freq)
  min(1, 2 * min(lower, upper))
}
cors <- fread(file.path(S, "out", "level_correlations_recomputed.csv"))
cors[, p_exact_table := mapply(exact_p, sum_d2, n)]
cors[, relative_error_of_committed := p_from_cor.test_exact_TRUE / p_exact_table - 1]
options(width = 220)
print(cors[, .(y, libraries = substr(libraries, 1, 14), n, spearman = signif(spearman, 6), sum_d2, committed = signif(p_from_cor.test_exact_TRUE, 6), exact = signif(p_exact_table, 6),
               exact_own_dp = signif(p_exact_own, 6), committed_over_exact_minus_1 = round(relative_error_of_committed, 4))])
fwrite(cors, file.path(S, "out", "level_correlations_exact_p.csv"))
# my tail computation (python, exact integers) against the table
mine <- c("396" = 4.190884e-05, "418" = 6.455901e-05, "430" = 8.086514e-05, "436" = 9.026788e-05, "486" = 2.124182e-04, "508" = 3.001668e-04)
cat("\nn = 22, my dynamic programme against the table:\n"); for (q in names(mine)) cat(sprintf("  S = %s: %.6e and %.6e\n", q, mine[[q]], exact_p(as.numeric(q), 22)))
cat(sprintf("n = 10, S = 30: %.6e (mine 5.818452e-03); n = 12, S = 50: %.6e (mine 1.591991e-03)\n", exact_p(30, 10), exact_p(50, 12)))

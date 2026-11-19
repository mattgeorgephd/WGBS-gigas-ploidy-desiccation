# cohort_separation.R
#
# Purpose   Do the 22 libraries group by cohort (heat experiment against control
#           experiment) or by ploidy when all CpGs are used, before any locus is
#           tested? Results 3.1 says that correlations show "consistent
#           methylation patterns within treatment groups" and that the PCA "did
#           not reveal discrete separation of samples by ploidy status or heat
#           treatment". This script puts numbers on both statements, from the
#           correlation matrices and PCA scores that step 04 writes.
#             1. Pairwise Pearson correlation: median within a group minus
#                median between groups, with a permutation p-value (group labels
#                shuffled among libraries, group sizes kept).
#             2. Principal components 1 to 4: the largest number of libraries
#                that a single threshold on the component assigns to their
#                group, with a permutation p-value for that number. A threshold
#                chosen to fit the labels always does better than chance, so
#                the count must be compared with its own null distribution.
#                Also the ratio of the spread of the two groups on the
#                component (interquartile ranges), since two groups can differ
#                in spread without differing in location.
#           Both are computed for cohort and for ploidy, at 1x and at 10x, in
#           both modes.
#
# Run       Rscript 06_diagnostics/cohort_separation.R          (about a minute)
#
# Reads     04_output/04_cpg_coverage/tables/library_correlation_{1x,10x}_<mode>.csv
#           04_output/04_cpg_coverage/tables/pca_scores_{1x,10x}_<mode>.csv
#           03_input_files/metadata/wgbs_libraries.csv
# Writes    06_diagnostics/cohort_separation_correlation.csv
#           06_diagnostics/cohort_separation_pca.csv
#
# Notes     - Cohort, heat treatment, sequencing project and read depth coincide
#             in these data (issue I-01). A separation by cohort is a difference
#             between the two experiments and cannot be attributed to heat.
#           - Ploidy is balanced across cohorts (5 + 6 against 5 + 6), so its
#             labels are permuted freely as well; permuting them within cohort
#             gives the same reading and is reported in a second column.
#           - 10,000 permutations; with a fixed seed the output is reproducible.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))

n_perm <- 10000
set.seed(20261119)
lib <- read_libraries()
tab <- function(f) output_path("04_cpg_coverage", "tables", f)
labels <- list(cohort = as.character(lib$treatment), ploidy = as.character(lib$ploidy))
cohort <- as.character(lib$treatment)

# Label permutations: free, and within cohort (for ploidy only; cohort labels
# cannot be permuted within cohort).
permute <- function(x, within = NULL) {
  if (is.null(within)) return(sample(x))
  out <- x
  for (g in unique(within)) out[within == g] <- sample(x[within == g])
  out
}

# ---- 1. Correlations -------------------------------------------------------------------
within_minus_between <- function(r, lab) {
  same <- outer(lab, lab, "==")
  upper <- upper.tri(r)
  median(r[upper & same]) - median(r[upper & !same])
}
cor_rows <- list()
for (mode in config$modes) for (cov in c(1, 10)) {
  m <- fread(tab(sprintf("library_correlation_%dx_%s.csv", cov, mode)))
  stopifnot(identical(m$library_id, lib$library_id))
  r <- as.matrix(m[, -1]); upper <- upper.tri(r)
  for (what in names(labels)) {
    lab <- labels[[what]]
    same <- outer(lab, lab, "==")
    obs <- within_minus_between(r, lab)
    null_free <- replicate(n_perm, within_minus_between(r, permute(lab)))
    null_within <- if (what == "ploidy") replicate(n_perm, within_minus_between(r, permute(lab, cohort))) else NA_real_
    by_pair <- vapply(split(r[upper & same], outer(lab, lab, paste)[upper & same]), median, numeric(1))
    cor_rows[[length(cor_rows) + 1]] <- data.table(
      mode = mode, min_coverage = cov, grouping = what,
      median_r_within = median(r[upper & same]), median_r_between = median(r[upper & !same]),
      within_minus_between = obs,
      p_permutation = (1 + sum(null_free >= obs)) / (1 + n_perm),
      p_permutation_within_cohort = if (what == "ploidy") (1 + sum(null_within >= obs)) / (1 + n_perm) else NA_real_,
      median_r_within_by_group = paste(sprintf("%s: %.4f", sub(" .*", "", names(by_pair)), by_pair), collapse = "; "),
      permutations = n_perm)
  }
}
cor_tab <- rbindlist(cor_rows)
num <- c("median_r_within", "median_r_between", "within_minus_between")
cor_tab[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(cor_tab, here("06_diagnostics", "cohort_separation_correlation.csv"))

# ---- 2. Principal components -----------------------------------------------------------
# Largest number of libraries assigned to their group by one threshold on x.
best_threshold <- function(x, lab) {
  a <- lab == lab[1]
  o <- order(x); a <- a[o]; n <- length(a)
  left_a <- cumsum(a)[-n]; k <- seq_len(n - 1)          # threshold after the k-th smallest value
  correct <- left_a + ((n - sum(a)) - (k - left_a))     # group a below, the other above
  max(correct, n - correct)
}
pca_rows <- list()
for (mode in config$modes) for (cov in c(1, 10)) {
  s <- fread(tab(sprintf("pca_scores_%dx_%s.csv", cov, mode)))
  stopifnot(identical(s$library_id, lib$library_id))
  for (pc in paste0("PC", 1:4)) for (what in names(labels)) {
    x <- s[[pc]]; lab <- labels[[what]]
    obs <- best_threshold(x, lab)
    null_free <- replicate(n_perm, best_threshold(x, permute(lab)))
    null_within <- if (what == "ploidy") replicate(n_perm, best_threshold(x, permute(lab, cohort))) else NA_integer_
    g <- split(x, lab)
    pca_rows[[length(pca_rows) + 1]] <- data.table(
      mode = mode, min_coverage = cov, component = pc, grouping = what,
      libraries_assigned_by_best_threshold = obs, of = length(x),
      null_median = median(null_free), null_95th_percentile = unname(quantile(null_free, 0.95)),
      p_permutation = (1 + sum(null_free >= obs)) / (1 + n_perm),
      p_permutation_within_cohort = if (what == "ploidy") (1 + sum(null_within >= obs)) / (1 + n_perm) else NA_real_,
      group_1 = names(g)[1], range_group_1 = sprintf("%.1f to %.1f", min(g[[1]]), max(g[[1]])),
      group_2 = names(g)[2], range_group_2 = sprintf("%.1f to %.1f", min(g[[2]]), max(g[[2]])),
      group_2_libraries_inside_range_of_group_1 = sum(g[[2]] >= min(g[[1]]) & g[[2]] <= max(g[[1]])),
      group_1_libraries_inside_range_of_group_2 = sum(g[[1]] >= min(g[[2]]) & g[[1]] <= max(g[[2]])),
      iqr_ratio_group_2_over_group_1 = round(IQR(g[[2]]) / IQR(g[[1]]), 2),
      permutations = n_perm)
  }
}
pca_tab <- rbindlist(pca_rows)
fwrite(pca_tab, here("06_diagnostics", "cohort_separation_pca.csv"))

options(width = 220)
cat("Correlations: median within a group minus median between groups\n")
print(cor_tab[, .(mode, min_coverage, grouping, median_r_within, median_r_between, within_minus_between, p_permutation, p_permutation_within_cohort)])
cat("\nPrincipal components, legacy mode: libraries assigned to their group by the best single threshold\n")
print(pca_tab[mode == "legacy", .(min_coverage, component, grouping, assigned = libraries_assigned_by_best_threshold, null_median, null_95th_percentile,
                                  p_permutation, p_permutation_within_cohort, iqr_ratio_group_2_over_group_1)])

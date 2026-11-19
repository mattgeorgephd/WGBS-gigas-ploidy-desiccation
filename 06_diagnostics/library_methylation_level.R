# library_methylation_level.R
#
# Purpose   How much of the difference between groups of libraries is a
#           difference in the overall methylation level of single libraries?
#           The sensitivity analyses (06_diagnostics/dml_sensitivity.R) find a
#           ploidy difference inside the heat cohort, no more CpGs above the
#           threshold inside the control cohort than regroupings give, and a
#           tested ploidy by cohort interaction. This script describes what
#           those contrasts consist of. It tests nothing at single CpGs.
#             1. Level of each library: genome-wide percent CpG methylation
#                (Bismark, step 02), and the mean percent methylation over the
#                tested CpGs and over the eligible ones (mean of the 22
#                libraries between 10% and 90%). The libraries above the
#                largest gap in the sorted genome-wide values are marked
#                "high"; the script stops if the largest gap in the mean over
#                eligible CpGs marks other libraries.
#             2. Group comparisons of the level (exact Wilcoxon rank-sum
#                tests), and its correlation with principal components 1 to 4
#                of step 04, with the means over tested and eligible CpGs,
#                with CHH methylation, and with the ELISA value of the animals
#                that have both measurements. The p-value of Spearman's test
#                is exact only below 10 libraries; above, R returns an
#                Edgeworth series approximation (the column name
#                spearman_p_exact is therefore wrong for 15 of its 16 rows).
#             3. The ploidy contrast inside each cohort across CpGs: how many
#                CpGs differ by more than the DML threshold, in which
#                direction, and whether the two cohorts agree (correlation of
#                the two differences across CpGs).
#             4. Does a CpG's dependence on the level predict the contrast?
#                The slope of each CpG on the genome-wide level is estimated
#                in one cohort and compared with the ploidy difference of the
#                other cohort, so no library is on both sides.
#             5. The same counts with one library left out at a time.
#             6. Relabelling references by complete enumeration: the number of
#                CpGs above the threshold for the observed grouping against
#                the number for every grouping of the same libraries into
#                groups of the same sizes. The reference median, percentiles
#                and minimum are taken over all groupings, the observed one
#                and its mirror image included; largest_other_grouping leaves
#                both out. Six comparisons: ploidy inside each cohort, with
#                all libraries and without the "high" ones; ploidy over both
#                cohorts without the "high" ones, with labels moved within
#                cohort; and the acute heat shock inside the heat cohort
#                (issue I-04), with labels moved within ploidy.
#             7. Similarity of the libraries of one group to each other: the
#                mean pairwise Pearson correlation (step 04, CpGs with 10x
#                in all libraries) among the diploids and among the triploids
#                of each cohort, against every other subset of the same size
#                from that cohort. (The mean, not the median: a subset that
#                holds four of five similar libraries has six similar pairs
#                out of ten, so its median can equal that of the full group.
#                With the median, 6 of the 252 subsets of the heat cohort
#                reach the value of the triploids instead of 1.)
#             8. The CpGs that the sensitivity analyses call with the
#                overdispersion correction and the F test (the lists
#                dml_<analysis>_mn_f_<mode>.csv; the two lists made with the
#                chi-squared test are not read): percent methylation of every
#                library, and whether one of the four groups lies entirely
#                outside the range of the other libraries. A CpG at which one
#                group of animals is uniformly different from all others is
#                what a sequence variant or an allele shared by related
#                animals looks like in bisulfite data (issue I-08).
#
# Run       Rscript 06_diagnostics/library_methylation_level.R   (about 2 minutes, 1.5 GB)
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           04_output/02_alignment_qc/tables/bismark_by_library.csv
#           04_output/03_global_methylation/tables/elisa_percent_5mC.csv
#           04_output/04_cpg_coverage/tables/pca_scores_10x_<mode>.csv
#           04_output/04_cpg_coverage/tables/library_correlation_10x_<mode>.csv
#           04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           04_output/05_differential_methylation/tables/diffmeth_{ploidy,heat}_<mode>.csv.gz
#           06_diagnostics/dml_sensitivity/dml_<analysis>_mn_f_<mode>.csv   (part 8; skipped
#             with a message if 06_diagnostics/dml_sensitivity.R has not been run)
# Writes    06_diagnostics/library_methylation_level/
#             library_level.csv               one row per library
#             level_group_tests.csv           group comparisons of the level
#             level_correlations.csv          level against PCs, CHH, ELISA
#             contrast_direction.csv          CpGs above the threshold, by direction
#             contrast_agreement.csv          correlation of the two within-cohort differences
#             level_slope_prediction.csv      part 4, correlations
#             level_slope_deciles.csv         part 4, by decile of slope
#             leave_one_library_out.csv       part 5
#             relabelling_reference.csv       part 6, one row per comparison
#             relabelling_small_sets.csv      part 6, every grouping of the comparisons
#                                             that have at most 252 groupings
#             pair_similarity.csv             part 7
#             corrected_dml_library_values.csv  part 8, one row per CpG and mode
#             corrected_dml_group_apart.csv     part 8, counts per list
#
# Notes     - "Difference" is the definition methylKit uses for meth.diff: the
#             percent methylation of one group minus that of the other, with
#             the reads of the group's libraries pooled (normalized counts of
#             step 05). The script checks its values against step 05 and stops
#             if they differ.
#           - A CpG counts as above the threshold when the absolute difference
#             exceeds config$dml_difference by more than 1e-9. A difference of
#             exactly 20 points is not above 20; in floating-point arithmetic
#             such a CpG can fall on either side, so ties are counted
#             separately (column n_exactly_at_threshold). Counts here can
#             therefore differ by a few CpGs from counts made without that
#             margin (dml_sensitivity/permutation_reference.csv: 3,350 against
#             3,349 for the control cohort, medians 4,019.5 against 4,018.5
#             for the heat cohort).
#           - Sign: triploid minus diploid for ploidy contrasts, heat cohort
#             minus control cohort for the cohort contrast. The interaction
#             contrast is the ploidy difference of the heat cohort minus that
#             of the control cohort.
#           - Cohort, heat treatment and sequencing project coincide (issue
#             I-01), and each of the four groups is one set of animals from
#             one source. A difference between two groups is a difference
#             between two sets of animals; nothing here attributes it to
#             ploidy or to heat.
#           - "High" libraries are found by a rule that does not use the group
#             labels, but leaving them out in parts 3 and 6 is a choice made
#             after seeing the data. Those rows are exploratory.
#           - Part 7: a group whose libraries resemble each other more than
#             other libraries of the cohort do is not by itself evidence of
#             one cross or one processing batch. The correlation between two
#             libraries falls with sequencing noise and with a difference in
#             level, so the values of the two cohorts cannot be set side by
#             side; at equal read depth the heat-cohort triploids are as
#             similar to each other as ordinary control libraries are
#             (issues.md, I-49, item 6).
#           - The smallest p-value a complete enumeration can give is
#             (groupings as extreme as the observed one) / (all groupings).
#             In part 6, with two groups of equal size every split occurs
#             twice (labels exchanged), so the smallest value is 2/252 for 5
#             against 5; for 2 against 5 it is 1/21. In part 7 the statistic
#             belongs to one subset and the complementary subset has its own
#             value, so the p-value is one-sided (subsets at least as similar
#             as the group) and its smallest value is 1/252 for 5 of 10
#             libraries and 1/924 for 6 of 12.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))
data.table::setDTthreads(1)

out_dir <- here("06_diagnostics", "library_methylation_level")
dir.create(out_dir, showWarnings = FALSE)
cutoff <- config$dml_difference
tie_tolerance <- 1e-9
lib <- read_libraries()
n_lib <- nrow(lib)
grp <- split(seq_len(n_lib), lib$group)
heat_lib <- which(lib$treatment == "heat")
ctrl_lib <- which(lib$treatment == "control")
acute <- which(lib$acute_heat_shock == "yes")
stopifnot(length(acute) == 4, all(acute %in% heat_lib))

above <- function(d) abs(d) > cutoff + tie_tolerance
at_threshold <- function(d) abs(abs(d) - cutoff) <= tie_tolerance

# ---- 1. Level of each library ------------------------------------------------------------
bis <- fread(output_path("02_alignment_qc", "tables", "bismark_by_library.csv"))
stopifnot(identical(bis$library_id, lib$library_id))
elisa <- fread(output_path("03_global_methylation", "tables", "elisa_percent_5mC.csv"))
level <- data.table(
  library_id = lib$library_id, animal_id = lib$animal_id, cohort = lib$cohort,
  ploidy = as.character(lib$ploidy), treatment = as.character(lib$treatment), group = as.character(lib$group),
  acute_heat_shock = lib$acute_heat_shock,
  pct_mCpG_genome = bis$pct_mCpG, pct_mCHG_genome = bis$pct_mCHG, pct_mCHH_genome = bis$pct_mCHH,
  elisa_percent_5mC = elisa$percent_5mC[match(lib$animal_id, elisa$animal_id)])
L <- level$pct_mCpG_genome

# Libraries above the largest gap in the sorted values.
above_largest_gap <- function(x) {
  s <- sort(x); k <- which.max(diff(s))
  list(high = x > s[k], gap = s[k + 1] - s[k], below = s[k], above = s[k + 1],
       second_gap = sort(diff(s), decreasing = TRUE)[2])
}
gap_genome <- above_largest_gap(L)
level[, high_level := gap_genome$high]
high <- which(level$high_level)

# Methylation matrices of one mode: C (methylated reads) and N (all reads) per
# CpG and library, in the library order of the sample sheet.
read_counts <- function(mode) {
  u <- fread(output_path("05_differential_methylation", "tables", sprintf("united_cpg_%s.csv.gz", mode)))
  stopifnot(sum(grepl("^numCs", names(u))) == n_lib)
  C <- as.matrix(u[, paste0("numCs", seq_len(n_lib)), with = FALSE])
  N <- C + as.matrix(u[, paste0("numTs", seq_len(n_lib)), with = FALSE])
  dimnames(C) <- dimnames(N) <- NULL
  list(key = paste(u$chr, u$start), C = C, N = N)
}
pooled <- function(m, idx) 100 * rowSums(m$C[, idx, drop = FALSE]) / rowSums(m$N[, idx, drop = FALSE])
difference <- function(m, a, b) pooled(m, b) - pooled(m, a)           # group b minus group a

# The contrasts of part 3, as functions of the count matrices.
contrasts <- list(
  "ploidy, all 22 libraries (the draft's ploidy contrast)" = function(m) difference(m, which(lib$ploidy == "diploid"), which(lib$ploidy == "triploid")),
  "cohort, all 22 libraries (the draft's heat contrast)"   = function(m) difference(m, ctrl_lib, heat_lib),
  "ploidy inside the heat cohort"                           = function(m) difference(m, grp[["diploid heat"]], grp[["triploid heat"]]),
  "ploidy inside the control cohort"                        = function(m) difference(m, grp[["diploid control"]], grp[["triploid control"]]),
  "interaction (heat-cohort ploidy difference minus control-cohort ploidy difference)" = function(m)
    difference(m, grp[["diploid heat"]], grp[["triploid heat"]]) - difference(m, grp[["diploid control"]], grp[["triploid control"]]),
  "ploidy inside the heat cohort, without the high libraries" = function(m)
    difference(m, setdiff(grp[["diploid heat"]], high), setdiff(grp[["triploid heat"]], high)),
  "ploidy inside the control cohort, without the high libraries" = function(m)
    difference(m, setdiff(grp[["diploid control"]], high), setdiff(grp[["triploid control"]], high)),
  "acute heat shock inside the heat cohort (4 libraries with it minus 6 without)" = function(m)
    difference(m, setdiff(heat_lib, acute), acute))

# Complete enumeration of the regroupings of a set of libraries. `strata` is a
# list with one element per stratum: list(a = libraries of group a, b = of
# group b). Labels move within a stratum only and group sizes are kept. Returns
# the number of CpGs above the threshold for the observed grouping and for
# every grouping, with the libraries of group a in each.
relabel_reference <- function(m, strata, chunk = 100) {
  libs <- unlist(lapply(strata, function(s) c(s$a, s$b)))
  pm <- 100 * m$C[, libs, drop = FALSE] / m$N[, libs, drop = FALSE]
  # A pooled group value lies between the smallest and the largest library
  # value of its group, so at a CpG whose library values span no more than the
  # threshold the difference cannot exceed the threshold by more than rounding
  # error (about 1e-13), far less than the margin of above(). Dropping those
  # CpGs therefore changes no count. This relies on the margin: with a plain
  # '> cutoff', a CpG whose groups sit at exactly 60% and 80% has a span of 20
  # and a computed difference of 20.000000000000007.
  span <- do.call(pmax, as.data.frame(pm)) - do.call(pmin, as.data.frame(pm))
  keep <- which(span > cutoff)
  C <- m$C[keep, , drop = FALSE]; N <- m$N[keep, , drop = FALSE]
  per_stratum <- lapply(strata, function(s) combn(c(s$a, s$b), length(s$a), simplify = FALSE))
  grid <- as.matrix(expand.grid(lapply(per_stratum, seq_along)))
  group_a <- lapply(seq_len(nrow(grid)), function(i)
    sort(unlist(lapply(seq_along(strata), function(k) per_stratum[[k]][[grid[i, k]]]))))
  counts <- integer(nrow(grid))
  for (start in seq(1, nrow(grid), by = chunk)) {
    rows <- start:min(start + chunk - 1, nrow(grid))
    A <- matrix(0, n_lib, length(rows)); B <- A
    for (j in seq_along(rows)) { A[group_a[[rows[j]]], j] <- 1; B[setdiff(libs, group_a[[rows[j]]]), j] <- 1 }
    d <- 100 * ((C %*% B) / (N %*% B) - (C %*% A) / (N %*% A))
    counts[rows] <- colSums(above(d))
  }
  a_obs <- sort(unlist(lapply(strata, `[[`, "a"))); b_obs <- sort(unlist(lapply(strata, `[[`, "b")))
  is_observed <- vapply(group_a, function(a) identical(a, a_obs), logical(1))
  is_mirror <- vapply(group_a, function(a) identical(a, b_obs), logical(1))
  stopifnot(sum(is_observed) == 1)
  list(observed = counts[is_observed], counts = counts, group_a = group_a, is_observed = is_observed, is_mirror = is_mirror,
       n_cpg_screened = length(keep))
}

direction_rows <- list(); agreement_rows <- list(); slope_rows <- list(); decile_rows <- list()
loo_rows <- list(); relabel_rows <- list(); small_rows <- list(); similarity_rows <- list()
value_rows <- list(); apart_rows <- list()
corrected_lists <- c("A1_ploidy_mn_f", "A1_heat_mn_f", "A2_interaction_mn_f", "A4_ploidy_in_control_mn_f", "A4_ploidy_in_heat_mn_f")
sens_dir <- here("06_diagnostics", "dml_sensitivity")

for (mode in config$modes) {
  m <- read_counts(mode)
  pm <- 100 * m$C / m$N
  mean_pm <- rowMeans(pm)
  eligible <- dml_eligible(mean_pm, config)
  level[, paste0("mean_pct_tested_", mode) := colMeans(pm)]
  level[, paste0("mean_pct_eligible_", mode) := colMeans(pm[eligible, , drop = FALSE])]
  gap_mode <- above_largest_gap(level[[paste0("mean_pct_eligible_", mode)]])
  if (!identical(gap_mode$high, level$high_level))
    stop("The largest gap in the mean over eligible CpGs (", mode, ") marks other libraries than the genome-wide level does.")
  pcs <- fread(output_path("04_cpg_coverage", "tables", sprintf("pca_scores_10x_%s.csv", mode)))
  stopifnot(identical(pcs$library_id, lib$library_id))
  for (pc in paste0("PC", 1:4)) level[, paste0(tolower(pc), "_10x_", mode) := pcs[[pc]]]

  # ---- 3. Contrasts across CpGs ----------------------------------------------------------
  d <- lapply(contrasts, function(f) f(m))
  # Check against step 05: the two contrasts of the draft must equal meth.diff.
  for (what in c("ploidy", "heat")) {
    s5 <- fread(output_path("05_differential_methylation", "tables", sprintf("diffmeth_%s_%s.csv.gz", what, mode)))
    mine <- d[[if (what == "ploidy") 1 else 2]]
    stopifnot(identical(paste(s5$chr, s5$start), m$key), max(abs(mine - s5$meth.diff)) < 1e-8)
  }
  for (nm in names(d)) {
    x <- d[[nm]]
    direction_rows[[length(direction_rows) + 1]] <- data.table(
      mode = mode, contrast = nm, n_cpg = length(x), n_eligible = sum(eligible),
      n_above_threshold = sum(above(x)), n_positive = sum(x > cutoff + tie_tolerance), n_negative = sum(x < -cutoff - tie_tolerance),
      n_exactly_at_threshold = sum(at_threshold(x)),
      share_in_majority_direction = max(sum(x > cutoff + tie_tolerance), sum(x < -cutoff - tie_tolerance)) / sum(above(x)),
      mean_difference_eligible = mean(x[eligible]), median_difference_eligible = median(x[eligible]),
      q25_difference_eligible = unname(quantile(x[eligible], 0.25)), q75_difference_eligible = unname(quantile(x[eligible], 0.75)),
      share_eligible_positive = mean(x[eligible] > 0))
  }
  d_heat <- d[["ploidy inside the heat cohort"]]; d_ctrl <- d[["ploidy inside the control cohort"]]
  agreement_sets <- list(
    "all libraries of each cohort" = list(h = d_heat, c = d_ctrl),
    "without the high libraries" = list(h = d[["ploidy inside the heat cohort, without the high libraries"]],
                                        c = d[["ploidy inside the control cohort, without the high libraries"]]))
  for (libs_used in names(agreement_sets)) for (set in c("all tested CpGs", "eligible CpGs")) {
    xh <- agreement_sets[[libs_used]]$h; xc <- agreement_sets[[libs_used]]$c
    k <- if (set == "all tested CpGs") rep(TRUE, length(xh)) else eligible
    both <- above(xh[k]) & above(xc[k])
    agreement_rows[[length(agreement_rows) + 1]] <- data.table(
      mode = mode, libraries = libs_used, cpg_set = set, n_cpg = sum(k),
      pearson = cor(xh[k], xc[k]), spearman = cor(xh[k], xc[k], method = "spearman"),
      n_above_in_heat_cohort = sum(above(xh[k])), n_above_in_control_cohort = sum(above(xc[k])),
      n_above_in_both = sum(both),
      n_above_in_both_same_sign = sum(both & sign(xh[k]) == sign(xc[k])),
      # Expected overlap if the two sets of CpGs were unrelated draws from this CpG set.
      n_above_in_both_expected_if_unrelated = sum(above(xh[k])) * sum(above(xc[k])) / sum(k))
  }

  # ---- 4. Slope on the level, estimated in the other cohort ------------------------------
  slope <- function(idx) {
    l <- L[idx] - mean(L[idx])
    as.vector((pm[, idx] - rowMeans(pm[, idx])) %*% l) / sum(l^2)
  }
  pairs <- list(
    list(from = "control cohort", idx = ctrl_lib, to = "ploidy inside the heat cohort", d = d_heat,
         level_difference = mean(L[grp[["triploid heat"]]]) - mean(L[grp[["diploid heat"]]])),
    list(from = "heat cohort", idx = heat_lib, to = "ploidy inside the control cohort", d = d_ctrl,
         level_difference = mean(L[grp[["triploid control"]]]) - mean(L[grp[["diploid control"]]])))
  for (p in pairs) {
    s <- slope(p$idx)
    for (set in c("all tested CpGs", "eligible CpGs")) {
      k <- if (set == "all tested CpGs") rep(TRUE, length(s)) else eligible
      slope_rows[[length(slope_rows) + 1]] <- data.table(
        mode = mode, slope_estimated_in = p$from, contrast = p$to, cpg_set = set, n_cpg = sum(k),
        level_difference_triploid_minus_diploid = p$level_difference,
        median_slope = median(s[k]), pearson = cor(s[k], p$d[k]), spearman = cor(s[k], p$d[k], method = "spearman"))
    }
    dec <- cut(s[eligible], quantile(s[eligible], 0:10 / 10), include.lowest = TRUE, labels = FALSE)
    decile_rows[[length(decile_rows) + 1]] <- data.table(decile = dec, s = s[eligible], x = p$d[eligible])[
      , .(mode = mode, slope_estimated_in = p$from, contrast = p$to, n_cpg = .N, median_slope = median(s),
          share_above_threshold_triploid_higher = mean(x > cutoff + tie_tolerance),
          share_above_threshold_diploid_higher = mean(x < -cutoff - tie_tolerance),
          median_difference = median(x)), by = decile][order(decile)]
  }

  # ---- 5. One library left out -----------------------------------------------------------
  for (coh in c("heat", "control")) {
    a <- grp[[paste("diploid", coh)]]; b <- grp[[paste("triploid", coh)]]
    full <- sum(above(difference(m, a, b)))
    for (k in c(a, b)) loo_rows[[length(loo_rows) + 1]] <- data.table(
      mode = mode, contrast = sprintf("ploidy inside the %s cohort", coh), library_left_out = lib$library_id[k],
      group = as.character(lib$group[k]), high_level = level$high_level[k],
      n_above_threshold_all_libraries = full,
      n_above_threshold_without_it = sum(above(difference(m, setdiff(a, k), setdiff(b, k)))))
  }

  # ---- 6. Relabelling references ---------------------------------------------------------
  dh <- grp[["diploid heat"]]; th <- grp[["triploid heat"]]; dc <- grp[["diploid control"]]; tc <- grp[["triploid control"]]
  comparisons <- list(
    list(name = "ploidy inside the heat cohort, all 10 libraries", group_a = "diploid", group_b = "triploid",
         labels_moved = "among the 10 libraries of the heat cohort",
         strata = list(list(a = dh, b = th))),
    list(name = "ploidy inside the control cohort, all 12 libraries", group_a = "diploid", group_b = "triploid",
         labels_moved = "among the 12 libraries of the control cohort",
         strata = list(list(a = dc, b = tc))),
    list(name = "ploidy inside the heat cohort, without the high libraries", group_a = "diploid", group_b = "triploid",
         labels_moved = "among the 7 libraries of the heat cohort that are not high",
         strata = list(list(a = setdiff(dh, high), b = setdiff(th, high)))),
    list(name = "ploidy inside the control cohort, without the high libraries", group_a = "diploid", group_b = "triploid",
         labels_moved = "among the 10 libraries of the control cohort that are not high",
         strata = list(list(a = setdiff(dc, high), b = setdiff(tc, high)))),
    list(name = "ploidy, both cohorts, without the high libraries", group_a = "diploid", group_b = "triploid",
         labels_moved = "within cohort, among the 17 libraries that are not high",
         strata = list(list(a = setdiff(dh, high), b = setdiff(th, high)), list(a = setdiff(dc, high), b = setdiff(tc, high)))),
    list(name = "acute heat shock inside the heat cohort", group_a = "with acute heat shock", group_b = "without",
         labels_moved = "within ploidy, among the 10 libraries of the heat cohort",
         strata = list(list(a = intersect(dh, acute), b = setdiff(dh, acute)), list(a = intersect(th, acute), b = setdiff(th, acute)))))
  for (cp in comparisons) {
    r <- relabel_reference(m, cp$strata)
    others <- r$counts[!r$is_observed & !r$is_mirror]
    relabel_rows[[length(relabel_rows) + 1]] <- data.table(
      mode = mode, comparison = cp$name, group_a = cp$group_a,
      n_a = sum(vapply(cp$strata, function(s) length(s$a), integer(1))),
      group_b = cp$group_b, n_b = sum(vapply(cp$strata, function(s) length(s$b), integer(1))),
      labels_moved = cp$labels_moved,
      n_above_threshold_observed = r$observed, groupings = length(r$counts),
      reference_median = median(r$counts), reference_p05 = unname(quantile(r$counts, 0.05)),
      reference_p95 = unname(quantile(r$counts, 0.95)), reference_min = min(r$counts),
      largest_other_grouping = max(others),
      observed_over_reference_median = r$observed / median(r$counts),
      groupings_at_or_above_observed = sum(r$counts >= r$observed),
      p_enumeration = sum(r$counts >= r$observed) / length(r$counts))
    if (length(r$counts) <= 252) small_rows[[length(small_rows) + 1]] <- data.table(
      mode = mode, comparison = cp$name, group_a = cp$group_a,
      libraries_in_group_a = vapply(r$group_a, function(a) paste(lib$library_id[a], collapse = " "), character(1)),
      n_above_threshold = r$counts, observed_grouping = r$is_observed, observed_grouping_with_labels_exchanged = r$is_mirror)[
        order(-n_above_threshold)]
  }

  # ---- 7. Similarity of the libraries of one group ---------------------------------------
  cm <- fread(output_path("04_cpg_coverage", "tables", sprintf("library_correlation_10x_%s.csv", mode)))
  stopifnot(identical(cm$library_id, lib$library_id))
  r <- as.matrix(cm[, -1]); dimnames(r) <- NULL
  mean_within <- function(idx) { s <- r[idx, idx]; mean(s[upper.tri(s)]) }
  for (coh in c("heat", "control")) {
    members <- if (coh == "heat") heat_lib else ctrl_lib
    pair_r <- combn(members, 2, function(p) r[p[1], p[2]])
    pair_ids <- combn(members, 2, simplify = FALSE)
    for (pl in c("diploid", "triploid")) {
      own <- grp[[paste(pl, coh)]]; other <- setdiff(members, own)
      subsets <- combn(members, length(own), simplify = FALSE)
      ref <- vapply(subsets, mean_within, numeric(1))
      obs <- mean_within(own)
      own_pair <- vapply(pair_ids, function(p) all(p %in% own), logical(1))
      top <- order(-pair_r)[seq_len(sum(own_pair))]
      similarity_rows[[length(similarity_rows) + 1]] <- data.table(
        mode = mode, cohort = paste(coh, "cohort"), group = pl, n_libraries = length(own), n_pairs = sum(own_pair),
        mean_r_within_group = obs, min_r_within_group = min(pair_r[own_pair]), max_r_within_group = max(pair_r[own_pair]),
        mean_r_within_other_group = mean_within(other),
        mean_r_between_groups = mean(r[own, other]),
        largest_mean_r_of_another_subset = max(ref[!vapply(subsets, function(s) setequal(s, own), logical(1))]),
        pairs_of_group_among_the_n_pairs_most_similar_in_cohort = sum(own_pair[top]),
        subsets_of_this_size = length(subsets), subsets_at_or_above_observed = sum(ref >= obs - 1e-12),
        p_enumeration = sum(ref >= obs - 1e-12) / length(subsets))
    }
  }

  # ---- 8. Library values at the CpGs called with the overdispersion correction -----------
  files <- file.path(sens_dir, sprintf("dml_%s_%s.csv", corrected_lists, mode))
  present <- file.exists(files)
  if (!any(present)) {
    message("Part 8 skipped for ", mode, " mode: no dml_*_mn_f_", mode, ".csv under 06_diagnostics/dml_sensitivity/")
  } else {
    called <- rbindlist(lapply(which(present), function(i) {
      x <- fread(files[i], select = c("chr", "start")); x[, list_id := corrected_lists[i]]; x
    }))
    called[, key := paste(chr, start)]
    stopifnot(all(called$key %in% m$key))
    per_cpg <- called[, .(lists = paste(list_id, collapse = "; ")), by = .(chr, start, key)]
    rows <- match(per_cpg$key, m$key)
    v <- pm[rows, , drop = FALSE]
    # A group is "apart" when the range of its libraries does not overlap the
    # range of all other libraries.
    apart <- vapply(names(grp), function(g) {
      own <- grp[[g]]; rest <- setdiff(seq_len(n_lib), own)
      lo_own <- do.call(pmin, as.data.frame(v[, own, drop = FALSE])); hi_own <- do.call(pmax, as.data.frame(v[, own, drop = FALSE]))
      lo_rest <- do.call(pmin, as.data.frame(v[, rest, drop = FALSE])); hi_rest <- do.call(pmax, as.data.frame(v[, rest, drop = FALSE]))
      lo_own > hi_rest | hi_own < lo_rest
    }, logical(nrow(v)))
    if (is.null(dim(apart))) apart <- matrix(apart, nrow = 1, dimnames = list(NULL, names(grp)))
    per_cpg[, group_apart_from_all_other_libraries := apply(apart, 1, function(a) if (any(a)) paste(names(grp)[a], collapse = "; ") else "none")]
    for (g in names(grp)) per_cpg[, paste0("range_", gsub(" ", "_", g)) :=
      sprintf("%.0f to %.0f", do.call(pmin, as.data.frame(v[, grp[[g]], drop = FALSE])), do.call(pmax, as.data.frame(v[, grp[[g]], drop = FALSE])))]
    vals <- as.data.table(round(v, 1)); setnames(vals, lib$library_id)
    value_rows[[length(value_rows) + 1]] <- cbind(data.table(mode = mode), per_cpg[, !"key"], vals)
    for (l in corrected_lists[present]) {
      k <- per_cpg$key %in% called[list_id == l, key]
      apart_rows[[length(apart_rows) + 1]] <- data.table(
        mode = mode, list = l, n_cpg = sum(k),
        n_diploid_heat_apart = sum(apart[k, "diploid heat"]), n_triploid_heat_apart = sum(apart[k, "triploid heat"]),
        n_diploid_control_apart = sum(apart[k, "diploid control"]), n_triploid_control_apart = sum(apart[k, "triploid control"]),
        n_no_group_apart = sum(!apply(apart[k, , drop = FALSE], 1, any)))
    }
  }
  rm(m, pm, d); invisible(gc())
}

# ---- 2. Group comparisons and correlations of the level ----------------------------------
measures <- c("pct_mCpG_genome", grep("^mean_pct_", names(level), value = TRUE), "pct_mCHH_genome")
level_comparisons <- list(
  "ploidy inside the heat cohort"    = list(a = grp[["diploid heat"]], b = grp[["triploid heat"]], la = "diploid heat", lb = "triploid heat"),
  "ploidy inside the control cohort" = list(a = grp[["diploid control"]], b = grp[["triploid control"]], la = "diploid control", lb = "triploid control"),
  "cohort"                           = list(a = ctrl_lib, b = heat_lib, la = "control cohort", lb = "heat cohort"),
  "ploidy, all 22 libraries"         = list(a = which(lib$ploidy == "diploid"), b = which(lib$ploidy == "triploid"), la = "diploid", lb = "triploid"),
  "acute heat shock inside the heat cohort" = list(a = setdiff(heat_lib, acute), b = acute, la = "without acute heat shock", lb = "with acute heat shock"))
test_rows <- list()
for (ms in measures) for (nm in names(level_comparisons)) {
  cp <- level_comparisons[[nm]]; x <- level[[ms]][cp$a]; y <- level[[ms]][cp$b]
  ties <- anyDuplicated(c(x, y)) > 0
  test_rows[[length(test_rows) + 1]] <- data.table(
    measure = ms, comparison = nm, group_a = cp$la, n_a = length(x), group_b = cp$lb, n_b = length(y),
    min_a = min(x), median_a = median(x), max_a = max(x), min_b = min(y), median_b = median(y), max_b = max(y),
    difference_of_means_b_minus_a = mean(y) - mean(x),
    a_libraries_above_every_b_library = sum(x > max(y)), b_libraries_above_every_a_library = sum(y > max(x)),
    complete_separation = min(x) > max(y) || min(y) > max(x),
    wilcoxon_p = suppressWarnings(wilcox.test(x, y, exact = !ties)$p.value),
    wilcoxon_p_is_exact = !ties)
}
test_tab <- rbindlist(test_rows)

cor_row <- function(x_name, y_name, x, y, set) {
  ok <- !is.na(x) & !is.na(y)
  data.table(x = x_name, y = y_name, libraries = set, n = sum(ok),
             pearson = cor(x[ok], y[ok]), spearman = cor(x[ok], y[ok], method = "spearman"),
             spearman_p_exact = suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = TRUE)$p.value))
}
cor_rows <- list()
for (mode in config$modes) {
  for (pc in paste0("pc", 1:4)) cor_rows[[length(cor_rows) + 1]] <-
    cor_row("pct_mCpG_genome", sprintf("%s_10x_%s", pc, mode), L, level[[sprintf("%s_10x_%s", pc, mode)]], "all 22")
  for (ms in c("mean_pct_tested_", "mean_pct_eligible_")) cor_rows[[length(cor_rows) + 1]] <-
    cor_row("pct_mCpG_genome", paste0(ms, mode), L, level[[paste0(ms, mode)]], "all 22")
}
cor_rows[[length(cor_rows) + 1]] <- cor_row("pct_mCpG_genome", "pct_mCHH_genome", L, level$pct_mCHH_genome, "all 22")
cor_rows[[length(cor_rows) + 1]] <- cor_row("pct_mCpG_genome", "pct_mCHH_genome", L[heat_lib], level$pct_mCHH_genome[heat_lib], "heat cohort")
cor_rows[[length(cor_rows) + 1]] <- cor_row("pct_mCpG_genome", "pct_mCHH_genome", L[ctrl_lib], level$pct_mCHH_genome[ctrl_lib], "control cohort")
cor_rows[[length(cor_rows) + 1]] <- cor_row("pct_mCpG_genome", "elisa_percent_5mC", L, level$elisa_percent_5mC,
                                            "animals with both a WGBS library and an ELISA value")
cor_tab <- rbindlist(cor_rows)

level[, level_rank := frank(-pct_mCpG_genome)]
level[, gap_rule := sprintf("largest gap in the sorted genome-wide values: %.2f to %.2f (%.2f points; second largest gap %.2f)",
                             gap_genome$below, gap_genome$above, gap_genome$gap, gap_genome$second_gap)]

# ---- Write -------------------------------------------------------------------------------
round_cols <- function(x) {
  num <- names(x)[vapply(x, is.double, logical(1))]
  if (length(num)) x[, (num) := lapply(.SD, signif, 6), .SDcols = num]
  x
}
write_table(round_cols(copy(level)), file.path(out_dir, "library_level.csv"))
write_table(round_cols(test_tab), file.path(out_dir, "level_group_tests.csv"))
write_table(round_cols(cor_tab), file.path(out_dir, "level_correlations.csv"))
write_table(round_cols(rbindlist(direction_rows)), file.path(out_dir, "contrast_direction.csv"))
write_table(round_cols(rbindlist(agreement_rows)), file.path(out_dir, "contrast_agreement.csv"))
write_table(round_cols(rbindlist(slope_rows)), file.path(out_dir, "level_slope_prediction.csv"))
write_table(round_cols(rbindlist(decile_rows)), file.path(out_dir, "level_slope_deciles.csv"))
write_table(rbindlist(loo_rows), file.path(out_dir, "leave_one_library_out.csv"))
relabel_tab <- round_cols(rbindlist(relabel_rows))
write_table(relabel_tab, file.path(out_dir, "relabelling_reference.csv"))
write_table(rbindlist(small_rows), file.path(out_dir, "relabelling_small_sets.csv"))
similarity_tab <- round_cols(rbindlist(similarity_rows))
write_table(similarity_tab, file.path(out_dir, "pair_similarity.csv"))
if (length(value_rows)) {
  write_table(rbindlist(value_rows), file.path(out_dir, "corrected_dml_library_values.csv"))
  write_table(rbindlist(apart_rows), file.path(out_dir, "corrected_dml_group_apart.csv"))
}

options(width = 220)
cat("High libraries (", level$gap_rule[1], "):\n  ", paste(level[high_level == TRUE, sprintf("%s (%s)", library_id, group)], collapse = ", "), "\n\n", sep = "")
print(level[, .(library_id, group, acute_heat_shock, pct_mCpG_genome = round(pct_mCpG_genome, 2),
                mean_pct_eligible_legacy = round(mean_pct_eligible_legacy, 1), high_level, elisa = round(elisa_percent_5mC, 2))])
cat("\nGroup comparisons of the genome-wide level\n")
print(test_tab[measure == "pct_mCpG_genome", .(comparison, median_a = round(median_a, 2), median_b = round(median_b, 2), complete_separation, wilcoxon_p = signif(wilcoxon_p, 3))])
cat("\nCpGs above the threshold, by direction (legacy mode)\n")
print(rbindlist(direction_rows)[mode == "legacy", .(contrast = substr(contrast, 1, 60), n_above_threshold, n_positive, n_negative, n_exactly_at_threshold,
                                                    majority = round(share_in_majority_direction, 3), median_eligible = round(median_difference_eligible, 2))])
cat("\nSimilarity of the libraries of one group (mean pairwise correlation)\n")
print(similarity_tab[, .(mode, cohort, group, within = round(mean_r_within_group, 4), other_group = round(mean_r_within_other_group, 4),
                         between = round(mean_r_between_groups, 4), top_pairs = pairs_of_group_among_the_n_pairs_most_similar_in_cohort, n_pairs,
                         at_or_above = subsets_at_or_above_observed, subsets = subsets_of_this_size, p = signif(p_enumeration, 3))])
if (length(apart_rows)) {
  cat("\nCpGs called with the overdispersion correction: is one group outside the range of all other libraries?\n")
  print(rbindlist(apart_rows))
}
cat("\nRelabelling references\n")
print(relabel_tab[, .(mode, comparison = substr(comparison, 1, 60), observed = n_above_threshold_observed, groupings, reference_median,
                      largest_other_grouping, at_or_above = groupings_at_or_above_observed, p = signif(p_enumeration, 3))])

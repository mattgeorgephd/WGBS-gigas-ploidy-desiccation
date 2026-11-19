# library_methylation_level.R
#
# Purpose   How much of the difference between groups of libraries is a
#           difference in the overall methylation level of single libraries?
#           The sensitivity analyses (06_diagnostics/dml_sensitivity.R) find a
#           ploidy difference inside the heat cohort, none inside the control
#           cohort, and a tested ploidy by cohort interaction. This script
#           describes what those contrasts consist of. It tests nothing at
#           single CpGs.
#             1. Level of each library: genome-wide percent CpG methylation
#                (Bismark, step 02), and the mean percent methylation over the
#                tested CpGs and over the eligible ones (mean of the 22
#                libraries between 10% and 90%). The libraries above the
#                largest gap in the sorted genome-wide values are marked
#                "high"; the script stops if the largest gap in the mean over
#                eligible CpGs marks other libraries.
#             2. Group comparisons of the level (exact Wilcoxon rank-sum
#                tests), and its correlation with principal components 1 to 4
#                of step 04, with non-CpG methylation, and with the ELISA
#                value of the animals that have both measurements.
#             3. The ploidy contrast inside each cohort across CpGs: how many
#                CpGs differ by more than the DML threshold, in which
#                direction, and whether the two cohorts agree (correlation of
#                the two differences across CpGs).
#             4. Does a CpG's dependence on the level predict the contrast?
#                The slope of each CpG on the genome-wide level is estimated
#                in one cohort and compared with the ploidy difference of the
#                other cohort, so no library is on both sides.
#             5. The same counts with one library left out at a time.
#             6. The ploidy contrast among the libraries that are not "high",
#                against every relabelling of those libraries (complete
#                enumeration, labels moved within cohort only).
#
# Run       Rscript 06_diagnostics/library_methylation_level.R   (about 3 minutes, 3 GB)
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           04_output/02_alignment_qc/tables/bismark_by_library.csv
#           04_output/03_global_methylation/tables/elisa_percent_5mC.csv
#           04_output/04_cpg_coverage/tables/pca_scores_10x_<mode>.csv
#           04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           04_output/05_differential_methylation/tables/diffmeth_{ploidy,heat}_<mode>.csv.gz
# Writes    06_diagnostics/library_methylation_level/
#             library_level.csv               one row per library
#             level_group_tests.csv           group comparisons of the level
#             level_correlations.csv          level against PCs, CHH, ELISA
#             contrast_direction.csv          CpGs above the threshold, by direction
#             contrast_agreement.csv          correlation of the two within-cohort differences
#             level_slope_prediction.csv      part 4, correlations
#             level_slope_deciles.csv         part 4, by decile of slope
#             leave_one_library_out.csv       part 5
#             relabelling_reference.csv       part 6
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
#             separately (column n_exactly_at_threshold).
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
#             labels, but leaving them out in part 6 is a choice made after
#             seeing the data. Part 6 is exploratory.
#           - The smallest p-value a complete enumeration can give is
#             (labellings as extreme as the observed one) / (all labellings).
#             With two groups of equal size every split occurs twice (labels
#             exchanged), so the smallest value is 2/252 for 5 against 5; for
#             2 against 5 it is 1/21.

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
    difference(m, grp[["diploid heat"]], grp[["triploid heat"]]) - difference(m, grp[["diploid control"]], grp[["triploid control"]]))

# Complete enumeration of the relabellings of a set of libraries. `sets` is a
# list with one element per cohort: list(a = libraries of group a, b = of
# group b). Labels move within a cohort only and group sizes are kept. Returns
# the number of CpGs above the threshold for the observed labelling and for
# every labelling.
relabel_reference <- function(m, sets, chunk = 100) {
  libs <- unlist(lapply(sets, function(s) c(s$a, s$b)))
  pm <- 100 * m$C[, libs, drop = FALSE] / m$N[, libs, drop = FALSE]
  # A pooled group value lies between the smallest and the largest library
  # value, so a CpG whose library values span no more than the threshold can
  # never be counted. Dropping those CpGs changes no count.
  span <- apply(pm, 1, max) - apply(pm, 1, min)
  keep <- which(span > cutoff)
  C <- m$C[keep, , drop = FALSE]; N <- m$N[keep, , drop = FALSE]
  per_cohort <- lapply(sets, function(s) combn(c(s$a, s$b), length(s$a), simplify = FALSE))
  grid <- as.matrix(expand.grid(lapply(per_cohort, seq_along)))
  counts <- integer(nrow(grid))
  for (start in seq(1, nrow(grid), by = chunk)) {
    rows <- start:min(start + chunk - 1, nrow(grid))
    A <- matrix(0, n_lib, length(rows)); B <- A
    for (j in seq_along(rows)) {
      a <- unlist(lapply(seq_along(sets), function(k) per_cohort[[k]][[grid[rows[j], k]]]))
      A[a, j] <- 1; B[setdiff(libs, a), j] <- 1
    }
    d <- 100 * ((C %*% B) / (N %*% B) - (C %*% A) / (N %*% A))
    counts[rows] <- colSums(above(d))
  }
  a_obs <- unlist(lapply(sets, `[[`, "a")); b_obs <- unlist(lapply(sets, `[[`, "b"))
  observed <- sum(above(100 * (rowSums(C[, b_obs, drop = FALSE]) / rowSums(N[, b_obs, drop = FALSE]) -
                               rowSums(C[, a_obs, drop = FALSE]) / rowSums(N[, a_obs, drop = FALSE]))))
  list(observed = observed, counts = counts, n_cpg_screened = length(keep))
}

direction_rows <- list(); agreement_rows <- list(); slope_rows <- list(); decile_rows <- list()
loo_rows <- list(); relabel_rows <- list()
high <- which(level$high_level)

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
  for (set in c("all tested CpGs", "eligible CpGs")) {
    k <- if (set == "all tested CpGs") rep(TRUE, length(d_heat)) else eligible
    agreement_rows[[length(agreement_rows) + 1]] <- data.table(
      mode = mode, cpg_set = set, n_cpg = sum(k),
      pearson = cor(d_heat[k], d_ctrl[k]), spearman = cor(d_heat[k], d_ctrl[k], method = "spearman"),
      n_above_in_both = sum(above(d_heat[k]) & above(d_ctrl[k])),
      n_above_in_both_same_sign = sum(above(d_heat[k]) & above(d_ctrl[k]) & sign(d_heat[k]) == sign(d_ctrl[k])),
      n_above_in_heat_cohort = sum(above(d_heat[k])), n_above_in_control_cohort = sum(above(d_ctrl[k])))
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

  # ---- 6. Relabelling reference ----------------------------------------------------------
  sets <- list(
    "ploidy inside the heat cohort, all 10 libraries" =
      list(heat = list(a = grp[["diploid heat"]], b = grp[["triploid heat"]])),
    "ploidy inside the control cohort, all 12 libraries" =
      list(control = list(a = grp[["diploid control"]], b = grp[["triploid control"]])),
    "ploidy inside the heat cohort, without the high libraries" =
      list(heat = list(a = setdiff(grp[["diploid heat"]], high), b = setdiff(grp[["triploid heat"]], high))),
    "ploidy inside the control cohort, without the high libraries" =
      list(control = list(a = setdiff(grp[["diploid control"]], high), b = setdiff(grp[["triploid control"]], high))),
    "ploidy, both cohorts, without the high libraries" =
      list(heat = list(a = setdiff(grp[["diploid heat"]], high), b = setdiff(grp[["triploid heat"]], high)),
           control = list(a = setdiff(grp[["diploid control"]], high), b = setdiff(grp[["triploid control"]], high))))
  for (nm in names(sets)) {
    r <- relabel_reference(m, sets[[nm]])
    n_a <- sum(vapply(sets[[nm]], function(s) length(s$a), integer(1)))
    n_b <- sum(vapply(sets[[nm]], function(s) length(s$b), integer(1)))
    relabel_rows[[length(relabel_rows) + 1]] <- data.table(
      mode = mode, comparison = nm, n_diploid = n_a, n_triploid = n_b,
      n_above_threshold_observed = r$observed, labellings = length(r$counts),
      reference_median = median(r$counts), reference_p05 = unname(quantile(r$counts, 0.05)),
      reference_p95 = unname(quantile(r$counts, 0.95)), reference_min = min(r$counts), reference_max = max(r$counts),
      second_largest = sort(r$counts, decreasing = TRUE)[2],
      labellings_at_or_above_observed = sum(r$counts >= r$observed),
      p_enumeration = sum(r$counts >= r$observed) / length(r$counts))
  }
  rm(m, pm, d); invisible(gc())
}

# ---- 2. Group comparisons and correlations of the level ----------------------------------
measures <- c("pct_mCpG_genome", grep("^mean_pct_", names(level), value = TRUE), "pct_mCHH_genome")
comparisons <- list(
  "ploidy inside the heat cohort"    = list(a = grp[["diploid heat"]], b = grp[["triploid heat"]], la = "diploid heat", lb = "triploid heat"),
  "ploidy inside the control cohort" = list(a = grp[["diploid control"]], b = grp[["triploid control"]], la = "diploid control", lb = "triploid control"),
  "cohort"                           = list(a = ctrl_lib, b = heat_lib, la = "control cohort", lb = "heat cohort"),
  "ploidy, all 22 libraries"         = list(a = which(lib$ploidy == "diploid"), b = which(lib$ploidy == "triploid"), la = "diploid", lb = "triploid"))
test_rows <- list()
for (ms in measures) for (nm in names(comparisons)) {
  cp <- comparisons[[nm]]; x <- level[[ms]][cp$a]; y <- level[[ms]][cp$b]
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
  num <- names(x)[vapply(x, function(v) is.double(v), logical(1))]
  x[, (num) := lapply(.SD, signif, 6), .SDcols = num]
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
write_table(round_cols(rbindlist(relabel_rows)), file.path(out_dir, "relabelling_reference.csv"))

options(width = 220)
cat("High libraries (", level$gap_rule[1], "):\n  ", paste(level[high_level == TRUE, sprintf("%s (%s)", library_id, group)], collapse = ", "), "\n\n", sep = "")
print(level[, .(library_id, group, acute_heat_shock, pct_mCpG_genome = round(pct_mCpG_genome, 2),
                mean_pct_eligible_legacy = round(mean_pct_eligible_legacy, 1), high_level, elisa = round(elisa_percent_5mC, 2))])
cat("\nGroup comparisons of the genome-wide level\n")
print(test_tab[measure == "pct_mCpG_genome", .(comparison, median_a = round(median_a, 2), median_b = round(median_b, 2), complete_separation, wilcoxon_p = signif(wilcoxon_p, 3))])
cat("\nCpGs above the threshold, by direction (legacy mode)\n")
print(rbindlist(direction_rows)[mode == "legacy", .(contrast = substr(contrast, 1, 48), n_above_threshold, n_positive, n_negative, n_exactly_at_threshold,
                                                    majority = round(share_in_majority_direction, 3), median_eligible = round(median_difference_eligible, 2))])
cat("\nRelabelling reference\n")
print(rbindlist(relabel_rows)[, .(mode, comparison = substr(comparison, 1, 58), observed = n_above_threshold_observed, labellings, reference_median,
                                  reference_max, at_or_above = labellings_at_or_above_observed, p = signif(p_enumeration, 3))])

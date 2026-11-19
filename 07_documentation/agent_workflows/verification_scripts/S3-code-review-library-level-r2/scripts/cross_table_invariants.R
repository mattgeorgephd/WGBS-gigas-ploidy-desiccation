# Invariants that must hold between the 13 committed tables, and against dml_sensitivity/dml_counts.csv.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; com <- file.path(repo, "06_diagnostics/library_methylation_level")
rd <- function(f) fread(file.path(com, f))
dir <- rd("contrast_direction.csv"); agr <- rd("contrast_agreement.csv"); loo <- rd("leave_one_library_out.csv"); ref <- rd("relabelling_reference.csv")
small <- rd("relabelling_small_sets.csv"); sim <- rd("pair_similarity.csv"); tests <- rd("level_group_tests.csv"); sp <- rd("level_slope_prediction.csv")
dec <- rd("level_slope_deciles.csv"); vals <- rd("corrected_dml_library_values.csv"); ap <- rd("corrected_dml_group_apart.csv"); lev <- rd("library_level.csv"); cors <- rd("level_correlations.csv")
ok <- function(what, cond) cat(sprintf("%-110s %s\n", what, if (isTRUE(cond)) "holds" else "FAILS"))
ok("contrast_direction: n_positive + n_negative = n_above_threshold (16 rows)", all(dir$n_positive + dir$n_negative == dir$n_above_threshold))
ok("contrast_direction: share_in_majority_direction = max(pos, neg) / n_above at 6 digits", all(abs(signif(pmax(dir$n_positive, dir$n_negative) / dir$n_above_threshold, 6) - dir$share_in_majority_direction) < 1e-12))
map <- c("ploidy inside the heat cohort, all 10 libraries" = "ploidy inside the heat cohort", "ploidy inside the control cohort, all 12 libraries" = "ploidy inside the control cohort",
         "ploidy inside the heat cohort, without the high libraries" = "ploidy inside the heat cohort, without the high libraries",
         "ploidy inside the control cohort, without the high libraries" = "ploidy inside the control cohort, without the high libraries",
         "acute heat shock inside the heat cohort" = "acute heat shock inside the heat cohort (4 libraries with it minus 6 without)")
x <- merge(ref[comparison %in% names(map), .(mode, contrast = map[comparison], matrix_path = n_above_threshold_observed)], dir[, .(mode, contrast, rowsums_path = n_above_threshold)], by = c("mode", "contrast"))
ok(sprintf("matrix path (relabelling_reference) = rowSums path (contrast_direction) in %d shared rows", nrow(x)), nrow(x) == 10 && all(x$matrix_path == x$rowsums_path))
x <- merge(unique(loo[, .(mode, contrast, n = n_above_threshold_all_libraries)]), dir[, .(mode, contrast, n_above_threshold)], by = c("mode", "contrast"))
ok("leave_one_library_out: count with all libraries = contrast_direction (4 rows)", nrow(x) == 4 && all(x$n == x$n_above_threshold))
x <- merge(agr[cpg_set == "all tested CpGs" & libraries == "all libraries of each cohort"], dcast(dir[contrast %in% c("ploidy inside the heat cohort", "ploidy inside the control cohort")], mode ~ contrast, value.var = "n_above_threshold"), by = "mode")
ok("contrast_agreement: counts above in each cohort = contrast_direction", all(x$n_above_in_heat_cohort == x[["ploidy inside the heat cohort"]], x$n_above_in_control_cohort == x[["ploidy inside the control cohort"]]))
ok("contrast_agreement: expected overlap = n1 * n2 / n at 6 digits", all(abs(signif(agr$n_above_in_heat_cohort * agr$n_above_in_control_cohort / agr$n_cpg, 6) - agr$n_above_in_both_expected_if_unrelated) < 1e-9))
ok("contrast_agreement: same-sign count <= overlap <= each count", all(agr$n_above_in_both_same_sign <= agr$n_above_in_both, agr$n_above_in_both <= pmin(agr$n_above_in_heat_cohort, agr$n_above_in_control_cohort)))
ok("relabelling_reference: p = groupings at or above / groupings at 6 digits", all(abs(signif(ref$groupings_at_or_above_observed / ref$groupings, 6) - ref$p_enumeration) < 1e-12))
ok("relabelling_reference: groupings are 252, 924, 21, 210, 4410, 100 in both modes", identical(ref$groupings, rep(c(252L, 924L, 21L, 210L, 4410L, 100L), 2)))
ok("relabelling_reference: groupings = product of binomial coefficients of the group sizes", identical(as.numeric(ref$groupings), rep(c(choose(10, 5), choose(12, 6), choose(7, 2), choose(10, 6), choose(7, 2) * choose(10, 6), choose(5, 2)^2), 2)))
s <- small[, .(n = .N, obs = sum(observed_grouping), mir = sum(observed_grouping_with_labels_exchanged), obs_count = n_above_threshold[observed_grouping],
               med = as.numeric(median(n_above_threshold)), mn = min(n_above_threshold), other = max(n_above_threshold[!observed_grouping & !observed_grouping_with_labels_exchanged]),
               ge = sum(n_above_threshold >= n_above_threshold[observed_grouping]), distinct = uniqueN(libraries_in_group_a)), by = .(mode, comparison)]
x <- merge(s, ref, by = c("mode", "comparison"))
ok("relabelling_small_sets: 8 comparisons, rows = groupings, all groupings distinct, one observed row each", nrow(x) == 8 && all(x$n == x$groupings, x$distinct == x$n, x$obs == 1))
ok("relabelling_small_sets: observed count, median, minimum, largest other, at-or-above = relabelling_reference", all(x$obs_count == x$n_above_threshold_observed, x$med == x$reference_median, x$mn == x$reference_min,
                                                                                                                  x$other == x$largest_other_grouping, x$ge == x$groupings_at_or_above_observed))
ok("relabelling_small_sets: a mirror row exists exactly where the two groups have the same size", all((x$mir == 1) == (x$n_a == x$n_b)))
ok("pair_similarity: within-other of one group = within of the other group of the cohort", all(sim[group == "diploid", mean_r_within_other_group] == sim[group == "triploid", mean_r_within_group]))
ok("pair_similarity: p = subsets at or above / subsets at 6 digits", all(abs(signif(sim$subsets_at_or_above_observed / sim$subsets_of_this_size, 6) - sim$p_enumeration) < 1e-12))
ok("level_slope_prediction: level difference = difference of means in level_group_tests", all(sp[contrast == "ploidy inside the heat cohort", level_difference_triploid_minus_diploid] == tests[measure == "pct_mCpG_genome" & comparison == "ploidy inside the heat cohort", difference_of_means_b_minus_a],
    sp[contrast == "ploidy inside the control cohort", level_difference_triploid_minus_diploid] == tests[measure == "pct_mCpG_genome" & comparison == "ploidy inside the control cohort", difference_of_means_b_minus_a]))
x <- merge(dec[, .(n = sum(n_cpg)), by = .(mode, slope_estimated_in)], unique(dir[, .(mode, n_eligible)]), by = "mode")
ok("level_slope_deciles: tenths sum to the number of eligible CpGs", all(x$n == x$n_eligible))
x <- merge(dec[, .(pos = sum(round(share_above_threshold_triploid_higher * n_cpg)), neg = sum(round(share_above_threshold_diploid_higher * n_cpg))), by = .(mode, contrast)],
           agr[cpg_set == "eligible CpGs" & libraries == "all libraries of each cohort", .(mode, h = n_above_in_heat_cohort, c = n_above_in_control_cohort)], by = "mode")
ok("level_slope_deciles: CpGs above the threshold summed over tenths = eligible counts of contrast_agreement", all(ifelse(grepl("heat", x$contrast), x$pos + x$neg == x$h, x$pos + x$neg == x$c)))
x <- vals[, .(n = .N), by = mode]; y <- ap[, .(n = sum(n_cpg)), by = mode]
cat("corrected_dml_library_values rows per mode:", paste(x$mode, x$n, collapse = "; "), "| sum of n_cpg over lists:", paste(y$mode, y$n, collapse = "; "), "(a CpG in two lists is one row)\n")
ok("corrected_dml_group_apart: none-apart count = rows labelled 'none' among the CpGs of each list", {
  all(sapply(seq_len(nrow(ap)), function(i) { v <- vals[mode == ap$mode[i] & grepl(ap$list[i], lists, fixed = TRUE)]; nrow(v) == ap$n_cpg[i] && sum(v$group_apart_from_all_other_libraries == "none") == ap$n_no_group_apart[i] })) })
dc <- fread(file.path(repo, "06_diagnostics/dml_sensitivity/dml_counts.csv"))
x <- merge(unique(dir[, .(mode, n_cpg, n_eligible)]), unique(dc[, .(mode, n_cpg_tested, n_cpg_eligible)]), by = "mode")
ok("n_cpg and n_eligible = dml_sensitivity/dml_counts.csv (480,579 / 65,480 and 454,555 / 61,056)", all(x$n_cpg == x$n_cpg_tested, x$n_eligible == x$n_cpg_eligible))
cmp <- rbind(merge(dir[contrast == "ploidy, all 22 libraries (the draft's ploidy contrast)", .(mode, here = n_above_threshold, pos = n_positive, neg = n_negative)], dc[analysis_id == "A0_ploidy", .(mode, sens = n_effect_above_cutoff, n_dml, n_hyper, n_hypo)], by = "mode"),
             merge(dir[contrast == "cohort, all 22 libraries (the draft's heat contrast)", .(mode, here = n_above_threshold, pos = n_positive, neg = n_negative)], dc[analysis_id == "A0_heat", .(mode, sens = n_effect_above_cutoff, n_dml, n_hyper, n_hypo)], by = "mode"),
             merge(dir[contrast == "ploidy inside the heat cohort", .(mode, here = n_above_threshold, pos = n_positive, neg = n_negative)], dc[analysis_id == "A4_ploidy_in_heat", .(mode, sens = n_effect_above_cutoff, n_dml, n_hyper, n_hypo)], by = "mode"),
             merge(dir[contrast == "ploidy inside the control cohort", .(mode, here = n_above_threshold, pos = n_positive, neg = n_negative)], dc[analysis_id == "A4_ploidy_in_control", .(mode, sens = n_effect_above_cutoff, n_dml, n_hyper, n_hypo)], by = "mode"),
             merge(dir[grepl("^interaction", contrast), .(mode, here = n_above_threshold, pos = n_positive, neg = n_negative)], dc[analysis_id == "A2_interaction", .(mode, sens = n_effect_above_cutoff, n_dml, n_hyper, n_hypo)], by = "mode"))
cat("\ncounts above the threshold here against n_effect_above_cutoff of dml_counts.csv (rows: ploidy 22, heat 22, ploidy in heat, ploidy in control, interaction):\n"); print(cmp)
# is the tie CpG of the control cohort a DML of A4_ploidy_in_control?
a4 <- fread(file.path(repo, "06_diagnostics/dml_sensitivity/dml_A4_ploidy_in_control_legacy.csv"))
cat("\nNC_047563.1:63152425 in dml_A4_ploidy_in_control_legacy.csv:", nrow(a4[chr == "NC_047563.1" & start == 63152425]), "row(s)\n"); print(a4[chr == "NC_047563.1" & start == 63152425])
cat("lists with the correction and the chi-squared test (not in part 8):\n")
for (f in c("dml_A1_ploidy_mn_chisq_legacy.csv", "dml_A1_heat_mn_chisq_legacy.csv")) { z <- fread(file.path(repo, "06_diagnostics/dml_sensitivity", f)); ff <- fread(file.path(repo, "06_diagnostics/dml_sensitivity", sub("chisq", "f", f)))
  cat(" ", f, nrow(z), "CpGs; of the", nrow(ff), "CpGs of the F-test list", nrow(merge(z, ff, by = c("chr", "start"))), "are among them\n") }

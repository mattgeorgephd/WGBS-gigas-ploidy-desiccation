suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t4 <- file.path(repo, "04_output/04_cpg_coverage/tables")
ok <- function(label, cond) cat(sprintf("%-100s %s\n", label, if (isTRUE(cond)) "OK" else "*** FAIL ***"))
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
cc <- fread(file.path(t4, "cpg_counts_by_coverage.csv"))
ok("counts: 50 rows, 2 modes x 5 sets x 5 levels", nrow(cc) == 50 && nrow(unique(cc[, .(mode, sample_set, min_coverage)])) == 50)
ok("counts: all <= any", all(cc$cpg_in_all_libraries <= cc$cpg_in_any_library))
ok("counts: non-increasing with coverage (all and any)", all(cc[order(min_coverage), .(a = all(diff(cpg_in_all_libraries) <= 0), b = all(diff(cpg_in_any_library) <= 0)), by = .(mode, sample_set)][, a & b]))
w <- dcast(cc, mode + min_coverage ~ sample_set, value.var = "cpg_in_all_libraries")
ok("counts: 'all' <= each subset (in all libraries)", all(w$all <= pmin(w$diploid, w$triploid, w$control, w$heat)))
w2 <- dcast(cc, mode + min_coverage ~ sample_set, value.var = "cpg_in_any_library")
ok("counts: 'all' >= each subset (in any library)", all(w2$all >= pmax(w2$diploid, w2$triploid, w2$control, w2$heat)))
m <- merge(cc[mode == "legacy"], cc[mode == "corrected"], by = c("sample_set", "min_coverage"))
ok("counts: corrected <= legacy + 1 (first-line CpG can return) for 'all libraries'", all(m$cpg_in_all_libraries.y <= m$cpg_in_all_libraries.x + 1))
ok("counts: n_libraries", all(cc[, n_libraries == c(all = 22, diploid = 11, triploid = 11, control = 12, heat = 10)[sample_set]]))
ok("counts: pct columns use 13,245,836", isTRUE(all.equal(cc$pct_of_genome_cpg_all, 100 * cc$cpg_in_all_libraries / 13245836)) && isTRUE(all.equal(cc$pct_of_genome_cpg_any, 100 * cc$cpg_in_any_library / 13245836)))
# legacy numbers typed in 08_archive/code/1_coverage_analysis.R comments
leg <- c(all = "5852939,3714898,2297236,480579,29884", diploid = "6583553,4410404,2928292,766054,47820", triploid = "NA,4527090,3023664,861750,54998",
         control = "6337867,3993030,2465533,529068,35027", heat = "7261925,5416388,4194066,1921856,193521")
for (s in names(leg)) { v <- suppressWarnings(as.numeric(strsplit(leg[[s]], ",")[[1]])); g <- cc[mode == "legacy" & sample_set == s][order(min_coverage)]$cpg_in_all_libraries
  ok(sprintf("legacy script comments, set %s (%d values compared)", s, sum(!is.na(v))), all(v == g, na.rm = TRUE)) }
dl <- fread(file.path(t4, "methylation_distribution_by_library.csv")); dg <- fread(file.path(t4, "methylation_distribution_by_group.csv"))
ok("distribution by library: 22 x 2 modes x 2 levels x 10 bins", nrow(dl) == 880)
ok("distribution: percent sums to 100 per library/mode/level", all(abs(dl[, sum(percent_of_cpg), by = .(library_id, mode, min_coverage)]$V1 - 100) < 1e-9))
lcv <- fread(file.path(t4, "library_coverage.csv"))
tot <- dl[, .(n = sum(n_cpg)), by = .(library_id, mode, min_coverage)]; tot <- merge(tot, lcv, by = "library_id")
ok("distribution totals: corrected 1x == n_cpg_1x - n_at_or_above_cutoff", all(tot[mode == "corrected" & min_coverage == 1, n == n_cpg_1x - n_cpg_at_or_above_cutoff]))
ok("distribution totals: legacy 1x == n_cpg_1x - 1", all(tot[mode == "legacy" & min_coverage == 1, n == n_cpg_1x - 1]))
ok("distribution totals: legacy 10x == n_cpg_10x - (first line >= 10)", all(tot[mode == "legacy" & min_coverage == 10, n == n_cpg_10x - (first_line_coverage >= 10)]))
ok("distribution totals: corrected 10x == n_cpg_10x - n_at_or_above_cutoff", all(tot[mode == "corrected" & min_coverage == 10, n == n_cpg_10x - n_cpg_at_or_above_cutoff]))
ok("distribution: ploidy/treatment labels match the sample sheet", all(dl[, ploidy == lib$ploidy[match(library_id, lib$library_id)] & treatment == lib$treatment[match(library_id, lib$library_id)]]))
se <- function(v) sd(v) / sqrt(length(v))
rec <- rbind(dl[, .(grouping = "ploidy", n = .N, mean_percent = mean(percent_of_cpg), se_percent = se(percent_of_cpg)), by = .(mode, min_coverage, group = ploidy, bin)],
             dl[, .(grouping = "treatment", n = .N, mean_percent = mean(percent_of_cpg), se_percent = se(percent_of_cpg)), by = .(mode, min_coverage, group = treatment, bin)],
             dl[, .(grouping = "all", n = .N, mean_percent = mean(percent_of_cpg), se_percent = se(percent_of_cpg)), by = .(mode, min_coverage, group = rep("all libraries", nrow(dl)), bin)])
mm <- merge(dg, rec, by = c("mode", "min_coverage", "grouping", "group", "bin"))
ok("distribution by group: recomputed means, SE and n agree", nrow(mm) == nrow(dg) && all(mm$n.x == mm$n.y) && isTRUE(all.equal(mm$mean_percent.x, mm$mean_percent.y)) && isTRUE(all.equal(mm$se_percent.x, mm$se_percent.y)))
print(dcast(dg[mode == "legacy" & min_coverage == 1 & bin %in% c("0-10", "91-100")], grouping + group + n ~ bin, value.var = "mean_percent"))
# correlation tables symmetric with unit diagonal; summary consistent
cs <- fread(file.path(t4, "library_correlation_summary.csv"))
for (md in c("legacy", "corrected")) for (x in c(1, 10)) {
  r <- fread(file.path(t4, sprintf("library_correlation_%dx_%s.csv", x, md))); M <- as.matrix(r[, -1]); off <- M[upper.tri(M)]
  same <- outer(lib$cohort, lib$cohort, "==")[upper.tri(M)]; s <- cs[mode == md & min_coverage == x]
  ok(sprintf("correlation %dx %s: symmetric, diag 1, ids in sheet order, summary consistent", x, md), isSymmetric(unname(M)) && all(abs(diag(M) - 1) < 1e-12) && identical(r$library_id, lib$library_id) && identical(colnames(M), lib$library_id) &&
       isTRUE(all.equal(c(min(off), max(off), median(off), median(off[same]), median(off[!same])), c(s$r_min, s$r_max, s$r_median, s$r_median_within_cohort, s$r_median_between_cohorts))))
}
pv <- fread(file.path(t4, "pca_variance.csv")); ok("pca variance sums to 100 per mode/level", all(abs(pv[, sum(percent_variance), by = .(mode, min_coverage)]$V1 - 100) < 1e-9))
ok("counts at 10x 'all' equal summary n_cpg and step 05 cpg_analysed", all(cs[min_coverage == 10]$n_cpg == cc[sample_set == "all" & min_coverage == 10][match(cs[min_coverage == 10]$mode, mode)]$cpg_in_all_libraries))

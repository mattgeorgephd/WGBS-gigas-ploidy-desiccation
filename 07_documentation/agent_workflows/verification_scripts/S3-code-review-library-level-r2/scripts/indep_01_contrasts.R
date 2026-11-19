# indep_01_contrasts.R  <mode>
# Independent recomputation of parts 3, 4, 5 and 8 of library_methylation_level.R and of the
# per-library means of part 1. Different code: data.table long format with by-group sums (no
# matrices, no helper of the repository), exact integer arithmetic for every threshold decision.
# Reads only; writes into the scratch folder.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE); mode <- args[1]
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
out <- file.path(S, "out"); committed <- file.path(repo, "06_diagnostics/library_methylation_level")
t0 <- proc.time()

lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE)
stopifnot(identical(lib$legacy_sample_index, 1:22))
bis <- read.csv(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"), stringsAsFactors = FALSE)
stopifnot(identical(bis$library_id, lib$library_id))
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)          # from the read counts, not from the pct column
cat("largest difference between my genome-wide level and the pct_mCpG column:", max(abs(L - bis$pct_mCpG)), "\n")
# high libraries: above the largest gap, found independently
o <- order(L); gaps <- diff(L[o]); kmax <- which(gaps == max(gaps)); stopifnot(length(kmax) == 1)
high <- seq_len(22) %in% o[(kmax + 1):22]
cat("gaps in the sorted genome-wide level, largest three:", paste(round(sort(gaps, decreasing = TRUE)[1:3], 4), collapse = ", "),
    "; position of the largest (number of libraries below it):", kmax, "\n")
cat("high libraries:", paste(lib$library_id[high], collapse = " "), "\n")
lib$high <- high
lib$cell <- paste0(substr(lib$ploidy, 1, 1), substr(lib$treatment, 1, 1), ifelse(lib$high, "H", "L"), ifelse(lib$acute_heat_shock == "yes", "A", "N"))
print(table(lib$cell))

u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
n_cpg <- nrow(u); u[, id := .I]
stopifnot(!anyDuplicated(u[, .(chr, start)]))
# is numCs + numTs the coverage column?
cov_diff <- sapply(1:22, function(j) sum(u[[paste0("coverage", j)]] != u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]]))
cat("CpG x library cells where coverage != numCs + numTs:", sum(cov_diff), "of", 22 * n_cpg, "\n")
cat("smallest numCs + numTs:", min(sapply(1:22, function(j) min(u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]]))),
    "; largest:", max(sapply(1:22, function(j) max(u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]]))), "\n")

long <- rbindlist(lapply(1:22, function(j) data.table(id = u$id, lj = j, C = as.numeric(u[[paste0("numCs", j)]]),
                                                      N = as.numeric(u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]]))))
cell_of_lib <- lib$cell
long[, cell := cell_of_lib[lj]]
long[, pm := 100 * C / N]
setkey(long, id)

# ---- part 1: mean percent methylation per library over tested and eligible CpGs ----------------
mean_pm <- long[, .(m = mean(pm)), keyby = id]$m
eligible <- mean_pm >= 10 & mean_pm <= 90
cat("eligible CpGs:", sum(eligible), "of", n_cpg, "\n")
# how close to the boundaries of the eligible range do means fall? (floating-point sensitivity of the set)
cat("CpGs with mean within 1e-9 of 10 or of 90:", sum(abs(mean_pm - 10) < 1e-9 | abs(mean_pm - 90) < 1e-9),
    "; within 1e-6:", sum(abs(mean_pm - 10) < 1e-6 | abs(mean_pm - 90) < 1e-6), "\n")
long[, elig := eligible[id]]
lib_means <- long[, .(mean_pct_tested = mean(pm), mean_pct_eligible = mean(pm[elig])), keyby = lj]
level_committed <- fread(file.path(committed, "library_level.csv"))
cmp <- function(what, mine, theirs, digits = 6) {
  ok <- isTRUE(all.equal(signif(as.numeric(mine), digits), as.numeric(theirs), tolerance = 1e-12))
  rel <- max(abs(as.numeric(mine) - as.numeric(theirs)) / pmax(abs(as.numeric(theirs)), 1e-300))
  cat(sprintf("  %-62s %s (largest relative difference %.2e)\n", what, if (ok) "identical at 6 significant digits" else "DIFFERENT", rel))
  invisible(ok)
}
cat("library_level.csv:\n")
cmp(paste0("mean_pct_tested_", mode), lib_means$mean_pct_tested, level_committed[[paste0("mean_pct_tested_", mode)]])
cmp(paste0("mean_pct_eligible_", mode), lib_means$mean_pct_eligible, level_committed[[paste0("mean_pct_eligible_", mode)]])
cmp("pct_mCpG_genome", L, level_committed$pct_mCpG_genome)
stopifnot(identical(level_committed$high_level, high))
# the gap rule on the mean over eligible CpGs
oe <- order(lib_means$mean_pct_eligible); ge <- diff(lib_means$mean_pct_eligible[oe])
cat("gaps in the sorted mean over eligible CpGs, largest three:", paste(round(sort(ge, decreasing = TRUE)[1:3], 3), collapse = ", "),
    "; libraries above the largest:", paste(lib$library_id[sort(oe[(which.max(ge) + 1):22])], collapse = " "), "\n")
ot <- order(lib_means$mean_pct_tested); gt <- diff(lib_means$mean_pct_tested[ot])
cat("gaps in the sorted mean over tested CpGs, largest three:", paste(round(sort(gt, decreasing = TRUE)[1:3], 3), collapse = ", "),
    "; libraries above the largest:", paste(lib$library_id[sort(ot[(which.max(gt) + 1):22])], collapse = " "), "\n")

# ---- group sums by cell -------------------------------------------------------------------------
cs <- long[, .(C = sum(C), N = sum(N)), keyby = .(id, cell)]
w <- dcast(cs, id ~ cell, value.var = c("C", "N"))
stopifnot(identical(w$id, seq_len(n_cpg)))
cells <- sort(unique(lib$cell))
gs <- function(sel) {            # pooled counts of the libraries in the cells `sel`
  stopifnot(all(sel %in% cells))
  list(C = Reduce(`+`, lapply(sel, function(cl) w[[paste0("C_", cl)]])), N = Reduce(`+`, lapply(sel, function(cl) w[[paste0("N_", cl)]])))
}
cell_of <- function(keep) sort(unique(lib$cell[keep]))
is <- function(...) Reduce(`&`, list(...))
dh <- lib$ploidy == "diploid" & lib$treatment == "heat"; th <- lib$ploidy == "triploid" & lib$treatment == "heat"
dc <- lib$ploidy == "diploid" & lib$treatment == "control"; tc <- lib$ploidy == "triploid" & lib$treatment == "control"
acute <- lib$acute_heat_shock == "yes"
check_cells <- function(keep) stopifnot(identical(sort(which(keep)), sort(which(lib$cell %in% cell_of(keep)))))   # a set must be a union of cells
two_group <- function(a, b) { check_cells(a); check_cells(b); list(a = gs(cell_of(a)), b = gs(cell_of(b))) }

# exact classification of b minus a against 20 points
exact2 <- function(g) {
  x <- g$b$C * g$a$N - g$a$C * g$b$N
  lhs <- 5 * abs(x); rhs <- g$a$N * g$b$N
  stopifnot(max(lhs) < 2^53, max(rhs) < 2^53, all(rhs > 0))
  list(sign = sign(x), above = lhs > rhs, tie = lhs == rhs,
       d_mk = 100 * (g$b$C / g$b$N - g$a$C / g$a$N),                 # methylKit's arithmetic
       d_script = 100 * g$b$C / g$b$N - 100 * g$a$C / g$a$N,          # the arithmetic of difference() in the script
       na = g$a$N, nb = g$b$N, min_gap_possible = 20 / max(rhs))
}
contrast_sets <- list(
  "ploidy, all 22 libraries (the draft's ploidy contrast)" = two_group(lib$ploidy == "diploid", lib$ploidy == "triploid"),
  "cohort, all 22 libraries (the draft's heat contrast)" = two_group(lib$treatment == "control", lib$treatment == "heat"),
  "ploidy inside the heat cohort" = two_group(dh, th),
  "ploidy inside the control cohort" = two_group(dc, tc),
  "ploidy inside the heat cohort, without the high libraries" = two_group(dh & !high, th & !high),
  "ploidy inside the control cohort, without the high libraries" = two_group(dc & !high, tc & !high),
  "acute heat shock inside the heat cohort (4 libraries with it minus 6 without)" = two_group(lib$treatment == "heat" & !acute, acute))
res <- lapply(contrast_sets, exact2)
# the interaction: four groups; floating point, with the CpGs within 1e-6 of the threshold written out for an exact check
g4 <- list(dh = gs(cell_of(dh)), th = gs(cell_of(th)), dc = gs(cell_of(dc)), tc = gs(cell_of(tc)))
d_int_script <- (100 * g4$th$C / g4$th$N - 100 * g4$dh$C / g4$dh$N) - (100 * g4$tc$C / g4$tc$N - 100 * g4$dc$C / g4$dc$N)
d_int_sens <- (100 * g4$th$C / g4$th$N - 100 * g4$tc$C / g4$tc$N) - (100 * g4$dh$C / g4$dh$N - 100 * g4$dc$C / g4$dc$N)
near <- which(abs(abs(d_int_script) - 20) < 1e-6 | abs(abs(d_int_sens) - 20) < 1e-6)
int_near <- data.table(mode = mode, chr = u$chr[near], start = u$start[near],
                       c_dh = g4$dh$C[near], n_dh = g4$dh$N[near], c_th = g4$th$C[near], n_th = g4$th$N[near],
                       c_dc = g4$dc$C[near], n_dc = g4$dc$N[near], c_tc = g4$tc$C[near], n_tc = g4$tc$N[near],
                       d_script = sprintf("%.17g", d_int_script[near]), d_sens = sprintf("%.17g", d_int_sens[near]))
fwrite(int_near, file.path(out, sprintf("interaction_near_threshold_%s.csv", mode)))
cat("interaction: CpGs within 1e-6 of 20 points:", length(near), "; of these |d| > 20 in the script's arithmetic without margin:",
    sum(abs(d_int_script[near]) > 20), "; in the arithmetic of dml_sensitivity.R:", sum(abs(d_int_sens[near]) > 20), "\n")
res[["interaction (heat-cohort ploidy difference minus control-cohort ploidy difference)"]] <- list(
  sign = sign(d_int_script), above = abs(d_int_script) > 20 + 1e-6, tie = abs(abs(d_int_script) - 20) <= 1e-6,   # ties confirmed exactly in python
  d_mk = d_int_sens, d_script = d_int_script, na = NA, nb = NA, min_gap_possible = NA)

# ---- part 3: contrast_direction.csv ------------------------------------------------------------
dir_committed <- fread(file.path(committed, "contrast_direction.csv"))[mode == get("mode", envir = globalenv())]
mine_dir <- rbindlist(lapply(names(res), function(nm) {
  r <- res[[nm]]; x <- r$d_script
  data.table(contrast = nm, n_cpg = n_cpg, n_eligible = sum(eligible),
             n_above_threshold = sum(r$above), n_positive = sum(r$above & r$sign > 0), n_negative = sum(r$above & r$sign < 0),
             n_exactly_at_threshold = sum(r$tie),
             share_in_majority_direction = max(sum(r$above & r$sign > 0), sum(r$above & r$sign < 0)) / sum(r$above),
             mean_difference_eligible = mean(x[eligible]), median_difference_eligible = median(x[eligible]),
             q25_difference_eligible = unname(quantile(x[eligible], 0.25)), q75_difference_eligible = unname(quantile(x[eligible], 0.75)),
             share_eligible_positive = mean(r$sign[eligible] > 0),
             # what a rule without margin gives in each floating-point form
             n_above_no_margin_methylkit_form = sum(abs(r$d_mk) > 20), n_above_no_margin_script_form = sum(abs(r$d_script) > 20),
             ties_above_20_in_methylkit_form = sum(r$tie & abs(r$d_mk) > 20), ties_above_20_in_script_form = sum(r$tie & abs(r$d_script) > 20),
             ties_at_exactly_20_in_methylkit_form = sum(r$tie & abs(r$d_mk) == 20), ties_below_20_in_methylkit_form = sum(r$tie & abs(r$d_mk) < 20),
             largest_float_error_at_ties = if (any(r$tie)) max(abs(abs(r$d_mk[r$tie]) - 20), abs(abs(r$d_script[r$tie]) - 20)) else NA_real_,
             smallest_distance_of_a_non_tie_from_20 = min(abs(abs(r$d_script[!r$tie]) - 20)),
             largest_pooled_reads_product = if (is.na(r$min_gap_possible)) NA_real_ else 20 / r$min_gap_possible)
}))
fwrite(cbind(mode = mode, mine_dir), file.path(out, sprintf("contrast_direction_recomputed_%s.csv", mode)))
cat("\ncontrast_direction.csv (", mode, "):\n", sep = "")
m <- merge(mine_dir, dir_committed, by = "contrast", suffixes = c("", ".committed"), sort = FALSE)
stopifnot(nrow(m) == 8)
for (col in c("n_cpg", "n_eligible", "n_above_threshold", "n_positive", "n_negative", "n_exactly_at_threshold"))
  cat(sprintf("  %-28s %s\n", col, if (all(m[[col]] == m[[paste0(col, ".committed")]])) "all 8 rows equal" else paste("DIFFERENT in rows", paste(which(m[[col]] != m[[paste0(col, ".committed")]]), collapse = ","))))
for (col in c("share_in_majority_direction", "mean_difference_eligible", "median_difference_eligible", "q25_difference_eligible", "q75_difference_eligible", "share_eligible_positive"))
  cmp(col, m[[col]], m[[paste0(col, ".committed")]])
options(width = 250)
print(mine_dir[, .(contrast = substr(contrast, 1, 45), n_above_threshold, n_positive, n_negative, n_exactly_at_threshold, n_above_no_margin_methylkit_form,
                   n_above_no_margin_script_form, ties_above_20_in_methylkit_form, ties_above_20_in_script_form, ties_at_exactly_20_in_methylkit_form, ties_below_20_in_methylkit_form,
                   largest_float_error_at_ties = signif(largest_float_error_at_ties, 3), smallest_distance_of_a_non_tie_from_20 = signif(smallest_distance_of_a_non_tie_from_20, 3),
                   largest_pooled_reads_product = signif(largest_pooled_reads_product, 3))])

# the exact ties of the two within-cohort contrasts, with methylKit's own value where a cached fit exists
tie_tab <- rbindlist(lapply(c("ploidy inside the heat cohort", "ploidy inside the control cohort"), function(nm) {
  r <- res[[nm]]; g <- contrast_sets[[nm]]; k <- which(r$tie)
  fit_file <- file.path(repo, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_A4_ploidy_in_%s_%s.rds", if (grepl("heat", nm)) "heat" else "control", mode))
  mk <- if (file.exists(fit_file)) readRDS(fit_file)$fit$meth.diff[k] else rep(NA_real_, length(k))
  data.table(mode = mode, contrast = nm, chr = u$chr[k], start = u$start[k], c_diploid = g$a$C[k], n_diploid = g$a$N[k], c_triploid = g$b$C[k], n_triploid = g$b$N[k],
             pct_diploid = 100 * g$a$C[k] / g$a$N[k], pct_triploid = 100 * g$b$C[k] / g$b$N[k],
             d_methylkit_form = sprintf("%.17g", r$d_mk[k]), d_script_form = sprintf("%.17g", r$d_script[k]),
             meth.diff_of_cached_methylkit_fit = sprintf("%.17g", mk),
             counted_by_rule_without_margin = abs(r$d_mk[k]) > 20, counted_by_cached_fit = abs(mk) > 20, counted_by_this_script = FALSE)
}))
fwrite(tie_tab, file.path(out, sprintf("ties_within_cohort_%s.csv", mode)))
cat("\nexact ties of the within-cohort contrasts (", mode, "):\n", sep = ""); print(tie_tab[, .N, by = .(contrast, counted_by_rule_without_margin, counted_by_cached_fit)])

# ---- part 3: contrast_agreement.csv ------------------------------------------------------------
agr_committed <- fread(file.path(committed, "contrast_agreement.csv"))[mode == get("mode", envir = globalenv())]
agr <- rbindlist(lapply(list(c("all libraries of each cohort", "ploidy inside the heat cohort", "ploidy inside the control cohort"),
                             c("without the high libraries", "ploidy inside the heat cohort, without the high libraries", "ploidy inside the control cohort, without the high libraries")),
  function(s) rbindlist(lapply(c("all tested CpGs", "eligible CpGs"), function(set) {
    h <- res[[s[2]]]; cc <- res[[s[3]]]; k <- if (set == "all tested CpGs") rep(TRUE, n_cpg) else eligible
    both <- h$above & cc$above & k
    data.table(libraries = s[1], cpg_set = set, n_cpg = sum(k), pearson = cor(h$d_mk[k], cc$d_mk[k]),
               spearman = cor(frank(h$d_script[k]), frank(cc$d_script[k])),
               spearman_methylkit_form = cor(frank(h$d_mk[k]), frank(cc$d_mk[k])),
               spearman_values_rounded_to_1e_9 = cor(frank(round(h$d_script[k], 9)), frank(round(cc$d_script[k], 9))),
               distinct_values_heat_script_form = uniqueN(h$d_script[k]), distinct_values_heat_rounded = uniqueN(round(h$d_script[k], 9)),
               n_above_in_heat_cohort = sum(h$above & k), n_above_in_control_cohort = sum(cc$above & k), n_above_in_both = sum(both),
               n_above_in_both_same_sign = sum(both & h$sign == cc$sign),
               n_above_in_both_expected_if_unrelated = as.numeric(sum(h$above & k)) * sum(cc$above & k) / sum(k))
  }))))
cat("\ncontrast_agreement.csv (", mode, "):\n", sep = "")
stopifnot(identical(agr$libraries, agr_committed$libraries), identical(agr$cpg_set, agr_committed$cpg_set))
for (col in c("n_cpg", "n_above_in_heat_cohort", "n_above_in_control_cohort", "n_above_in_both", "n_above_in_both_same_sign"))
  cat(sprintf("  %-28s %s\n", col, if (all(agr[[col]] == agr_committed[[col]])) "all 4 rows equal" else "DIFFERENT"))
for (col in c("pearson", "spearman", "n_above_in_both_expected_if_unrelated")) cmp(col, agr[[col]], agr_committed[[col]])
print(agr[, .(libraries, cpg_set, spearman_script_form = sprintf("%.9f", spearman), spearman_methylkit_form = sprintf("%.9f", spearman_methylkit_form),
              spearman_values_rounded_to_1e_9 = sprintf("%.9f", spearman_values_rounded_to_1e_9), distinct_values_heat_script_form, distinct_values_heat_rounded)])
fwrite(cbind(mode = mode, agr), file.path(out, sprintf("contrast_agreement_recomputed_%s.csv", mode)))

# ---- part 5: leave_one_library_out.csv ---------------------------------------------------------
loo_committed <- fread(file.path(committed, "leave_one_library_out.csv"))[mode == get("mode", envir = globalenv())]
loo <- rbindlist(lapply(list(list("heat", dh, th), list("control", dc, tc)), function(s) {
  g <- two_group(s[[2]], s[[3]])
  rbindlist(lapply(which(s[[2]] | s[[3]]), function(k) {
    ck <- as.numeric(u[[paste0("numCs", k)]]); nk <- ck + u[[paste0("numTs", k)]]
    gg <- g
    if (s[[2]][k]) { gg$a$C <- g$a$C - ck; gg$a$N <- g$a$N - nk } else { gg$b$C <- g$b$C - ck; gg$b$N <- g$b$N - nk }
    data.table(contrast = sprintf("ploidy inside the %s cohort", s[[1]]), library_left_out = lib$library_id[k],
               n_above_threshold_all_libraries = sum(exact2(g)$above), n_above_threshold_without_it = sum(exact2(gg)$above), ties_without_it = sum(exact2(gg)$tie))
  }))
}))
cat("\nleave_one_library_out.csv (", mode, "): ", sep = "")
stopifnot(identical(loo$library_left_out, loo_committed$library_left_out))
cat(if (all(loo$n_above_threshold_all_libraries == loo_committed$n_above_threshold_all_libraries) &&
        all(loo$n_above_threshold_without_it == loo_committed$n_above_threshold_without_it)) "all rows equal" else "DIFFERENT", nrow(loo), "rows; exact ties per row:", paste(loo$ties_without_it, collapse = " "), "\n")

# ---- part 4: slopes on the genome-wide level ----------------------------------------------------
slope_committed <- fread(file.path(committed, "level_slope_prediction.csv"))[mode == get("mode", envir = globalenv())]
dec_committed <- fread(file.path(committed, "level_slope_deciles.csv"))[mode == get("mode", envir = globalenv())]
slope_of <- function(keep) {
  lc <- numeric(22); lc[keep] <- L[keep] - mean(L[keep])
  long[lj %in% which(keep), .(s = sum(pm * lc[lj])), keyby = id]$s / sum(lc^2)      # sum(pm * centred level) / sum(centred level^2)
}
type7 <- function(x, p) { x <- sort(x); h <- (length(x) - 1) * p + 1; lo <- floor(h); hi <- ceiling(h); x[lo] + (h - lo) * (x[hi] - x[lo]) }
slope_rows <- list(); dec_rows <- list(); boundary_rows <- list()
for (pr in list(list(from = "control cohort", keep = lib$treatment == "control", to = "ploidy inside the heat cohort", tri = th, di = dh),
                list(from = "heat cohort", keep = lib$treatment == "heat", to = "ploidy inside the control cohort", tri = tc, di = dc))) {
  s <- slope_of(pr$keep); r <- res[[pr$to]]
  # the same slope with the values centred first, as the script does it (exact zero where a CpG is constant in the cohort)
  lc0 <- numeric(22); lc0[pr$keep] <- L[pr$keep] - mean(L[pr$keep])
  s_centred <- long[lj %in% which(pr$keep), .(s = sum((pm - mean(pm)) * lc0[lj])), keyby = id]$s / sum(lc0^2)
  cat(sprintf("slope estimated in the %s: largest difference between the two formulas %.3g; exact zeros %d (uncentred formula) and %d (centred); Spearman with the centred slopes, all tested CpGs %.9f, eligible %.9f\n",
              pr$from, max(abs(s - s_centred)), sum(s == 0), sum(s_centred == 0), cor(frank(s_centred), frank(r$d_script)), cor(frank(s_centred[eligible]), frank(r$d_script[eligible]))))
  for (set in c("all tested CpGs", "eligible CpGs")) {
    k <- if (set == "all tested CpGs") rep(TRUE, n_cpg) else eligible
    slope_rows[[length(slope_rows) + 1]] <- data.table(slope_estimated_in = pr$from, contrast = pr$to, cpg_set = set, n_cpg = sum(k),
      level_difference_triploid_minus_diploid = mean(L[pr$tri]) - mean(L[pr$di]), median_slope = median(s[k]),
      pearson = cor(s[k], r$d_mk[k]), spearman = cor(frank(s[k]), frank(r$d_script[k])),
      spearman_rounded_to_1e_9 = cor(frank(round(s[k], 9)), frank(round(r$d_script[k], 9))), distinct_slopes = uniqueN(s[k]), distinct_slopes_rounded = uniqueN(round(s[k], 9)))
  }
  se <- s[eligible]; br <- sapply(0:10 / 10, function(p) type7(se, p))
  dec <- 1L + rowSums(outer(se, br[2:10], `>`))                      # right-closed intervals: a value equal to a break goes to the lower tenth
  dec_rank <- pmin(10L, 1L + as.integer(floor(10 * (frank(se, ties.method = "first") - 1) / length(se))))   # tenths by rank, for comparison
  dt <- data.table(decile = dec, s = se, above_pos = (r$above & r$sign > 0)[eligible], above_neg = (r$above & r$sign < 0)[eligible], x = r$d_script[eligible])
  dec_rows[[length(dec_rows) + 1]] <- dt[, .(slope_estimated_in = pr$from, contrast = pr$to, n_cpg = .N, median_slope = median(s),
    share_above_threshold_triploid_higher = mean(above_pos), share_above_threshold_diploid_higher = mean(above_neg), median_difference = median(x)), keyby = decile]
  boundary_rows[[length(boundary_rows) + 1]] <- data.table(mode = mode, slope_estimated_in = pr$from, breaks_unique = !anyDuplicated(br),
    slopes_equal_to_an_inner_break = sum(se %in% br[2:10]), slopes_within_1e_9_of_an_inner_break = sum(sapply(br[2:10], function(b) sum(abs(se - b) < 1e-9))),
    cpgs_in_another_tenth_by_rank = sum(dec != dec_rank), eligible_slopes_exactly_zero = sum(se == 0), distinct_slopes = uniqueN(se), n = length(se),
    all_tested_slopes_exactly_zero = sum(s == 0))
}
slope_tab <- rbindlist(slope_rows); dec_tab <- rbindlist(dec_rows)
cat("\nlevel_slope_prediction.csv (", mode, "):\n", sep = "")
stopifnot(identical(slope_tab$slope_estimated_in, slope_committed$slope_estimated_in), identical(slope_tab$cpg_set, slope_committed$cpg_set))
for (col in c("n_cpg", "level_difference_triploid_minus_diploid", "median_slope", "pearson", "spearman")) cmp(col, slope_tab[[col]], slope_committed[[col]])
print(slope_tab[, .(slope_estimated_in, cpg_set, spearman_mine = sprintf("%.9f", spearman), spearman_rounded_to_1e_9 = sprintf("%.9f", spearman_rounded_to_1e_9), distinct_slopes, distinct_slopes_rounded)])
print(cbind(slope_committed[, .(spearman_committed = spearman)]))
fwrite(cbind(mode = mode, slope_tab), file.path(out, sprintf("level_slope_prediction_recomputed_%s.csv", mode)))
cat("level_slope_deciles.csv (", mode, "):\n", sep = "")
stopifnot(all(dec_tab$decile == dec_committed$decile), identical(dec_tab$slope_estimated_in, dec_committed$slope_estimated_in))
cat(sprintf("  %-62s %s\n", "n_cpg", if (all(dec_tab$n_cpg == dec_committed$n_cpg)) "all 20 rows equal" else "DIFFERENT"))
for (col in c("median_slope", "share_above_threshold_triploid_higher", "share_above_threshold_diploid_higher", "median_difference")) cmp(col, dec_tab[[col]], dec_committed[[col]])
print(rbindlist(boundary_rows))
fwrite(rbindlist(boundary_rows), file.path(out, sprintf("decile_boundaries_%s.csv", mode)))

# ---- part 8: library values at the CpGs called with the correction and the F test -----------------
val_committed <- fread(file.path(committed, "corrected_dml_library_values.csv"))[mode == get("mode", envir = globalenv())]
apart_committed <- fread(file.path(committed, "corrected_dml_group_apart.csv"))[mode == get("mode", envir = globalenv())]
sens <- file.path(repo, "06_diagnostics/dml_sensitivity")
lists <- sub(sprintf("^dml_(.*)_%s\\.csv$", mode), "\\1", list.files(sens, pattern = sprintf("^dml_.*_mn_f_%s\\.csv$", mode)))
cat("\npart 8, lists with the correction and the F test found for", mode, ":", paste(lists, collapse = ", "), "\n")
cat("  lists with the correction and the chi-squared test (not read by the script):",
    paste(list.files(sens, pattern = sprintf("^dml_.*_mn_chisq_%s\\.csv$", mode)), collapse = ", "), "\n")
called <- rbindlist(lapply(lists, function(l) { x <- read.csv(file.path(sens, sprintf("dml_%s_%s.csv", l, mode)), stringsAsFactors = FALSE); data.table(chr = x$chr, start = x$start, list = l) }))
called <- merge(called, u[, .(chr, start, id)], by = c("chr", "start"), all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(called$id))
grp4 <- paste(lib$ploidy, lib$treatment)
lv <- long[id %in% called$id]; lv[, g := grp4[lj]]
ap <- rbindlist(lapply(unique(grp4), function(g0) {
  lv[, .(group = g0, apart = min(pm[g == g0]) > max(pm[g != g0]) | max(pm[g == g0]) < min(pm[g != g0]),
         # the same decision by exact cross-multiplication of the counts
         apart_exact = {
           own <- g == g0; lo <- which(own)[which.min(pm[own])]; hi <- which(own)[which.max(pm[own])]
           rlo <- which(!own)[which.min(pm[!own])]; rhi <- which(!own)[which.max(pm[!own])]
           (C[lo] * N[rhi] > C[rhi] * N[lo]) | (C[hi] * N[rlo] < C[rlo] * N[hi]) },
         touching = min(pm[g == g0]) == max(pm[g != g0]) | max(pm[g == g0]) == min(pm[g != g0])), by = id]
}))
stopifnot(identical(ap$apart, ap$apart_exact))
apw <- dcast(ap, id ~ group, value.var = "apart")
per_list <- merge(called, apw, by = "id", sort = FALSE)[, .(n_cpg = .N, n_diploid_heat_apart = sum(`diploid heat`), n_triploid_heat_apart = sum(`triploid heat`),
  n_diploid_control_apart = sum(`diploid control`), n_triploid_control_apart = sum(`triploid control`),
  n_no_group_apart = sum(!(`diploid heat` | `triploid heat` | `diploid control` | `triploid control`)),
  n_two_groups_apart = sum(`diploid heat` + `triploid heat` + `diploid control` + `triploid control` >= 2)), by = list]
cat("corrected_dml_group_apart.csv (", mode, "): ", sep = "")
mm <- merge(per_list, apart_committed, by = "list", suffixes = c("", ".committed"), sort = FALSE)
stopifnot(nrow(mm) == nrow(apart_committed), nrow(mm) == nrow(per_list))
cat(if (all(sapply(c("n_cpg", "n_diploid_heat_apart", "n_triploid_heat_apart", "n_diploid_control_apart", "n_triploid_control_apart", "n_no_group_apart"),
                   function(col) all(mm[[col]] == mm[[paste0(col, ".committed")]])))) "all rows equal" else "DIFFERENT", "\n")
print(per_list)
cat("  groups whose range touches the range of the other libraries without overlapping (strict inequality decides):", sum(ap$touching & !ap$apart), "\n")
# the values table: one row per CpG; values rounded to one decimal; the label column
vw <- dcast(lv, id ~ lj, value.var = "pm"); setnames(vw, as.character(1:22), lib$library_id)
vw <- merge(u[, .(id, chr, start)], vw, by = "id")
vc <- merge(val_committed, vw, by = c("chr", "start"), suffixes = c(".committed", ""), sort = FALSE)
stopifnot(nrow(vc) == nrow(val_committed), nrow(val_committed) == uniqueN(called$id))
worst <- max(sapply(lib$library_id, function(l) max(abs(round(vc[[l]], 1) - vc[[paste0(l, ".committed")]]))))
cat("corrected_dml_library_values.csv (", mode, "): ", nrow(val_committed), " CpGs; largest difference of a library value rounded to one decimal: ", worst, "\n", sep = "")
lab <- merge(val_committed[, .(chr, start, group_apart_from_all_other_libraries, lists)], merge(u[, .(id, chr, start)], apw, by = "id"), by = c("chr", "start"), sort = FALSE)
mylab <- apply(as.matrix(lab[, c("diploid control", "diploid heat", "triploid control", "triploid heat"), with = FALSE]), 1,
               function(a) if (any(a)) paste(c("diploid control", "diploid heat", "triploid control", "triploid heat")[a], collapse = "; ") else "none")
cat("  label column group_apart_from_all_other_libraries:", if (identical(unname(mylab), lab$group_apart_from_all_other_libraries)) "all rows equal" else "DIFFERENT", "\n")
mylists <- called[, .(lists = paste(list, collapse = "; ")), by = .(chr, start)]
ll <- merge(val_committed[, .(chr, start, lists)], mylists, by = c("chr", "start"), suffixes = c(".committed", ""), sort = FALSE)
same_sets <- mapply(function(a, b) setequal(strsplit(a, "; ")[[1]], strsplit(b, "; ")[[1]]), ll$lists, ll$lists.committed)
cat("  column lists: the same set of lists in", sum(same_sets), "of", length(same_sets), "rows; CpGs in more than one list:", sum(grepl(";", ll$lists.committed)), "\n")
# range strings
rng <- lv[, .(lo = min(pm), hi = max(pm)), by = .(id, g)]
rng[, txt := sprintf("%.0f to %.0f", lo, hi)]
rw <- merge(u[, .(id, chr, start)], dcast(rng, id ~ g, value.var = "txt"), by = "id")
rc <- merge(val_committed, rw, by = c("chr", "start"), sort = FALSE)
cat("  range columns:", if (all(rc$range_diploid_control == rc$`diploid control`, rc$range_diploid_heat == rc$`diploid heat`,
                                rc$range_triploid_control == rc$`triploid control`, rc$range_triploid_heat == rc$`triploid heat`)) "all rows equal" else "DIFFERENT", "\n")
cat("\nelapsed seconds:", round((proc.time() - t0)[3]), "\n"); print(gc())

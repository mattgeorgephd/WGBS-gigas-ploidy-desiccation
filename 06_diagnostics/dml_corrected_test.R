# dml_corrected_test.R
#
# Purpose   What the lists of methylKit's DML model with the overdispersion
#           correction can be read as. The correction is fitted by
#           06_diagnostics/dml_sensitivity.R (analyses *_mn_chisq: overdispersion
#           "MN", chi-squared test, the setting of Venkataraman et al. 2022). This
#           script reads those fits and asks of them what issues.md, I-06 (items
#           2, 4, 6 and 7), I-01, I-02 and I-07 need:
#             1. Level under relabelling, all 22 libraries. The share of p-values
#                below 0.05, 0.01 and 0.001 with the ploidy labels reassigned
#                within cohort, for the ploidy test (heat as covariate) and the
#                tested interaction, without and with the correction.
#             2. Level under regrouping within one cohort: the ploidy test among
#                the ten heat-treated libraries and among the twelve controls.
#             3. Length of the corrected lists under relabelling: how many CpGs
#                above the difference cutoff a relabelled grouping calls at the
#                p-value threshold of the observed list, and how many it calls
#                whatever the other CpGs do. With the same labellings, the share
#                of the CpGs above the cutoff with a corrected p-value below
#                0.05, 0.01, 0.001 and 0.0001.
#             4. Where the correction is not applied. methylKit sets a dispersion
#                estimate that is not above 1 to 1, so that the corrected p-value
#                is the uncorrected one there. Share of the eligible CpGs and of
#                each corrected list at which that happens.
#             5. Stability of the corrected lists of the draft's two contrasts:
#                corrected mode (coverage filter), the Benjamini-Hochberg
#                adjustment in place of SLIM, and looser q-value thresholds.
#             6. The corrected ploidy CpGs inside each cohort.
#             7. The p-values at the loci that the manuscript values map names
#                (03_input_files/manuscript/manuscript_values_map.csv).
#           A diagnostics script: nothing in 01_pipeline/ reads its outputs.
#
# Run       Rscript 06_diagnostics/dml_corrected_test.R
#             after 06_diagnostics/dml_sensitivity.R (full run; this script reads
#             its cached fits and its tables) and
#             06_diagnostics/library_methylation_level.R (library_level.csv names
#             the libraries with a raised overall methylation level, issue I-49).
#             About 25 minutes on one core, under 2 GB of memory.
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           06_diagnostics/dml_sensitivity/cache/fit_<analysis>_<mode>.rds   (git-ignored; the script
#             checks that each fit was made from the united table it reads, with the stated settings)
#           06_diagnostics/dml_sensitivity/dml_<analysis>_<mode>.csv   the lists
#           06_diagnostics/library_methylation_level/library_level.csv   column high_level
#
# Writes    06_diagnostics/dml_corrected_test/
#             relabel_level_pooled.csv          part 1, one row per sample and labelling
#             relabel_level_within_cohort.csv   part 2, one row per cohort and grouping
#             relabel_level_summary.csv         parts 1 and 2 summarised: per design, test,
#                                               kind of relabelling and p level
#             relabel_corrected_lists.csv       part 3, one row per design and labelling
#             relabel_corrected_lists_summary.csv   observed_n_above_cutoff: CpGs above the
#                                               difference cutoff with the observed labels;
#                                               observed_n_in_shares: the CpGs over which the
#                                               observed_share_* columns are computed (those
#                                               above the cutoff, except in the row of the heat
#                                               contrast: the draft's heat DML, q < 0.01 and
#                                               above the cutoff, one CpG fewer)
#             correction_not_applied.csv        part 4
#             corrected_lists_stability.csv     part 5, one row per analysis, mode and measure
#             corrected_ploidy_dml_by_cohort.csv  part 6, one row per corrected ploidy CpG
#             named_loci.csv                    part 7
#             session_info.txt
#
# Notes     - The test is methylKit's logReg (methylKit 1.28.0), written out in
#             base R so that it can be repeated for many labellings: binomial
#             glm of the full model and of the model without the tested term;
#             deviance difference referred to chi-squared with 1 degree of
#             freedom (uncorrected), or divided first by the dispersion
#             estimate (corrected). The dispersion estimate is the Pearson
#             chi-squared statistic of the full model over its residual degrees
#             of freedom, set to 1 where it is not above 1; a p-value that
#             cannot be computed is set to 1, as methylKit does. Where the
#             estimate itself cannot be computed (NaN), methylKit does not set
#             it to 1: ifelse(phi > 1, phi, 1) returns NA, the p-value is NA and
#             is replaced by 1. The script does the same, and such a CpG is not
#             counted as one at which the correction is not applied (it has
#             p = 1 and is never called). For the
#             observed labels the script checks its p-values against the cached
#             fits of dml_sensitivity.R and stops if they differ by more than
#             `equal_tolerance` (relative).
#           - Relabelling keeps the group sizes within each cohort, as the
#             permutation reference of dml_sensitivity.R does: the heat label
#             cannot be relabelled, because it does not vary within a cohort
#             (issue I-01). Within one cohort a grouping and its mirror image
#             give the same p-values, so one of each pair is used (the first
#             library of the cohort is never in the group labelled 1).
#           - Part 1 uses two independent seeded samples of 3,000 eligible CpGs
#             (eligible: mean methylation between 10% and 90%, as in
#             dml_sensitivity.R) and 150 relabellings for each. The seeds are
#             those of the check of 2026-11-19 (issues.md, I-06), so that the
#             uncorrected and corrected columns repeat that check. A share in a
#             sample is compared with the share that the observed labels give
#             at the same CpGs ("at or above observed"), and also with the
#             share over all eligible CpGs of the full data (the value of
#             dml_counts.csv; "at or above reference").
#           - Part 3 counts, among the CpGs above the difference cutoff, those
#             with a corrected p-value at or below the largest p-value that has
#             q < config$dml_qvalue in the observed fit (a fixed threshold, so
#             that no q-values are computed for the relabellings), and those
#             with a p-value below config$dml_qvalue / (number of CpGs tested),
#             which have q below the cutoff whatever the other CpGs do. Where
#             the groups do not differ and the test holds its level, the second
#             count is about 0.01 per labelling.
#           - Kinds of relabelling (part 1): the number of exchanges of two
#             libraries that turn the real grouping of the heat cohort into the
#             relabelled one (0, 1 or 2; a grouping and its mirror image count
#             as one), and how the five libraries with a raised overall
#             methylation level (issue I-49; library_level.csv, column
#             high_level) are split between the two labels (3:2, or 4:1 and
#             5:0). Part 2 splits the regroupings of the controls by whether the
#             two such control libraries are in one group, and those of the heat
#             cohort by whether its three such libraries are.
#             The mirror rule of column exchanges_in_heat_cohort is that of a
#             test inside one cohort, where a grouping and its mirror image give
#             the same p-values. In the design of all 22 libraries (parts 1 and
#             3) the mirror of the heat cohort's grouping, with the controls
#             unchanged, is a different labelling and a different test. There
#             the column is a measure of how close a grouping is to the real
#             heat-cohort grouping, counting a grouping and its mirror as the
#             same: "kept or changed by one exchange" (0 or 1) includes
#             labellings that are 4 or 5 direct exchanges away, and "two
#             exchanges" those 2 or 3 direct exchanges away. Column
#             direct_exchanges_in_heat_cohort gives the number of direct
#             exchanges (0 to 5) without the mirror rule.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))
if (!requireNamespace("R.utils", quietly = TRUE)) stop("Missing R package: R.utils (needed to read .csv.gz files).", call. = FALSE)

t0 <- start_step("dml_corrected_test", "the DML model with the overdispersion correction: level, lists under relabelling, stability")

# ---- Settings -------------------------------------------------------------------------
n_sample          <- 3000            # eligible CpGs per sample (part 1, and part 2 on the first sample)
sample_seeds      <- c(606, 607)     # the two samples of part 1; the second is drawn from the CpGs not in the first
labelling_seeds   <- c(707, 708)     # relabellings of part 1, one seed per sample
n_relabel         <- 150             # relabellings per sample (part 1)
control_seed      <- 802             # regroupings of the controls (parts 2 and 3)
n_control         <- 150             # regroupings of the controls besides the observed one
list_seed         <- 20261119        # labellings of the pooled design in part 3: those of the permutation reference of dml_sensitivity.R
n_list_pooled     <- 200             # of which the first n are used besides the observed one
p_levels          <- c(0.05, 0.01, 0.001)
p_levels_lists    <- c(0.05, 0.01, 0.001, 1e-4)
equal_tolerance   <- 1e-6            # largest relative difference accepted against the cached fits

out_dir  <- here("06_diagnostics", "dml_corrected_test")
sens_dir <- here("06_diagnostics", "dml_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, name)

# ---- Inputs ---------------------------------------------------------------------------------
lib <- read_libraries()
trip <- as.integer(lib$ploidy == "triploid")
heat <- as.integer(lib$treatment == "heat")
cohort_rows <- split(seq_len(nrow(lib)), lib$cohort)
heat_rows <- which(lib$treatment == "heat"); control_rows <- which(lib$treatment == "control")
if (!setequal(heat_rows, cohort_rows[[unique(lib$cohort[heat_rows])]]) || length(unique(lib$cohort)) != 2) {
  stop("Expected two cohorts, one of them holding every heat-treated library.", call. = FALSE)
}
level_file <- here("06_diagnostics", "library_methylation_level", "library_level.csv")
if (!file.exists(level_file)) stop("Missing ", level_file, ". Run 06_diagnostics/library_methylation_level.R first.", call. = FALSE)
level <- fread(level_file)
if (!identical(level$library_id, lib$library_id)) stop("library_level.csv does not list the libraries of the library sheet in order.", call. = FALSE)
high <- level$high_level
log_msg("libraries with a raised overall level (library_level.csv): %s", paste(lib$library_id[high], collapse = ", "))

united_path <- function(mode) {
  f <- output_path("05_differential_methylation", "tables", sprintf("united_cpg_%s.csv.gz", mode))
  if (!file.exists(f)) stop("Missing ", f, ". Run 01_pipeline/05_differential_methylation.R first.", call. = FALSE)
  f
}
read_counts <- function(mode) {
  u <- fread(united_path(mode))
  cs <- as.matrix(u[, paste0("numCs", seq_len(nrow(lib))), with = FALSE]) + 0
  ns <- cs + as.matrix(u[, paste0("numTs", seq_len(nrow(lib))), with = FALSE])
  list(pos = u[, .(chr, start)], cs = cs, ns = ns, eligible = which(dml_eligible(mean_percent_methylation(u), config)),
       md5 = unname(tools::md5sum(united_path(mode))))
}

# A cached fit of dml_sensitivity.R, checked against the united table and the settings.
read_fit <- function(analysis_id, mode, counts, overdispersion) {
  f <- file.path(sens_dir, "cache", sprintf("fit_%s_%s.rds", analysis_id, mode))
  if (!file.exists(f)) stop("Missing ", f, ". Run 06_diagnostics/dml_sensitivity.R (full run) first.", call. = FALSE)
  x <- readRDS(f)
  s <- x$signature
  if (!identical(s[["united_md5"]], counts$md5) || !identical(s[["cpg_set"]], "all") || !identical(s[["overdispersion"]], overdispersion) ||
      !identical(s[["test"]], "Chisq") || nrow(x$fit) != nrow(counts$pos)) {
    stop(basename(f), " was not made from the current united table with overdispersion '", overdispersion,
         "' and the chi-squared test. Run 06_diagnostics/dml_sensitivity.R again.", call. = FALSE)
  }
  x$fit
}
# Share of the eligible CpGs of the full data with a p-value below `level` in a
# legacy fit (the value of dml_counts.csv, columns share_eligible_p_below_*).
full_share <- function(id, level) mean(fits[[id]]$pvalue[legacy$eligible] < level)

# ---- The test --------------------------------------------------------------------------------
# glm.fit warns at some CpGs (a fitted probability numerically 0 or 1, no
# convergence); those warnings are expected, any other stops the run.
quiet_glm <- function(X, y, w) {
  withCallingHandlers(stats::glm.fit(X, y / w, weights = w, family = stats::binomial()),
                      warning = function(cond) {
                        msg <- conditionMessage(cond)
                        if (!grepl("fitted probabilities numerically 0 or 1|algorithm did not converge", msg)) {
                          stop("Unexpected warning: ", msg, call. = FALSE)
                        }
                        invokeRestart("muffleWarning")
                      })
}
# One CpG: p-values without and with the correction, and the dispersion
# estimate before it is set to 1. `X_red` NULL means the intercept-only model.
one_test <- function(y, w, X_full, X_red, dev_red = NULL) {
  full <- quiet_glm(X_full, y, w)
  if (is.null(dev_red)) dev_red <- if (is.null(X_red)) full$null.deviance else quiet_glm(X_red, y, w)$deviance
  mu <- full$fitted.values
  phi <- sum((y - w * mu)^2 / (mu * (w - w * mu))) / (length(w) - ncol(X_full))
  dev <- dev_red - full$deviance
  p_none <- stats::pchisq(dev, 1, lower.tail = FALSE)
  p_mn <- stats::pchisq(dev / ifelse(phi > 1, phi, 1), 1, lower.tail = FALSE)
  c(p_uncorrected = if (is.na(p_none)) 1 else p_none, p_corrected = if (is.na(p_mn)) 1 else p_mn, dispersion = phi)
}
# Rows `rows` of count matrices `cs`, `ns` (columns: the libraries of the design).
run_tests <- function(cs, ns, rows, X_full, X_red) {
  if (length(rows) == 0) return(matrix(numeric(0), 0, 3, dimnames = list(NULL, c("p_uncorrected", "p_corrected", "dispersion"))))
  t(vapply(rows, function(i) one_test(cs[i, ], ns[i, ], X_full, X_red), numeric(3)))
}
# The correction is not applied where methylKit sets the dispersion estimate
# to 1: where it is not above 1. A NaN estimate is not set to 1 by methylKit
# (ifelse(NaN > 1, ...) is NA, so the p-value is NA and replaced by 1), and
# one_test() does the same; such a CpG has p = 1, is never a call, and is not
# counted here.
not_corrected <- function(dispersion) !is.na(dispersion) & dispersion <= 1
largest_relative_difference <- function(x, y) max(abs(x - y) / pmax(abs(y), .Machine$double.xmin))
check_against <- function(mine, cached, what) {
  worst <- largest_relative_difference(mine, cached)
  if (!is.finite(worst) || worst > equal_tolerance) {
    stop(what, ": the p-values computed here differ from the cached fit (largest relative difference ", signif(worst, 3), ").", call. = FALSE)
  }
  log_msg("  %s: p-values agree with the cached fit at %d CpGs (largest relative difference %.1e)", what, length(mine), worst)
}
shares <- function(p, levels) vapply(levels, function(l) mean(p < l), numeric(1))
level_names <- function(levels) gsub(".", "_", format(levels, scientific = FALSE, drop0trailing = TRUE), fixed = TRUE)

# Groupings of one cohort: one of each mirror pair (the first library of the
# cohort is never labelled 1), observed first; all of them, or the observed one
# and `n_random` others drawn with `seed`.
cohort_groupings <- function(rows, n_random = Inf, seed = NULL) {
  obs <- trip[rows]
  if (obs[1] != 0) stop("The first library of a cohort must be diploid for the mirror rule.", call. = FALSE)
  sets <- utils::combn(length(rows), sum(obs))
  sets <- sets[, sets[1, ] != 1, drop = FALSE]
  labs <- apply(sets, 2, function(s) { z <- integer(length(rows)); z[s] <- 1L; z })
  is_obs <- colSums(labs != obs) == 0
  if (sum(is_obs) != 1) stop("The observed grouping was not found among the groupings of a cohort.", call. = FALSE)
  others <- which(!is_obs)
  if (length(others) > n_random) { set.seed(seed); others <- sort(sample(others, n_random)) }
  list(labels = cbind(labs[, is_obs], labs[, others, drop = FALSE]), n_distinct = ncol(labs))
}
# Labellings of all 22 libraries, ploidy permuted within each cohort (observed first).
pooled_labellings <- function(n, seed) {
  set.seed(seed)
  cbind(trip, vapply(seq_len(n), function(draw) {
    z <- trip; for (b in cohort_rows) z[b] <- trip[b][sample.int(length(b))]; z
  }, integer(length(trip))))
}
# Exchanges of two libraries between the groups of the heat cohort, counting a
# grouping and its mirror image as the same (the rule of a test inside the
# cohort; in the pooled design the mirror is a different labelling, see the
# header), and without that rule.
exchanges_in_heat_cohort <- function(labels) {
  moved <- colSums(labels[heat_rows, , drop = FALSE] != trip[heat_rows])
  pmin(moved, length(heat_rows) - moved) / 2
}
direct_exchanges_in_heat_cohort <- function(labels) colSums(labels[heat_rows, , drop = FALSE] != trip[heat_rows]) / 2
high_split <- function(labels) {
  ones <- colSums(labels[high, , drop = FALSE])
  ifelse(abs(2 * ones - sum(high)) == 1, "3:2", "4:1 or 5:0")
}

legacy <- read_counts("legacy")
n_tested <- nrow(legacy$pos)
p_single <- config$dml_qvalue / n_tested
fits <- list(
  A0_ploidy = read_fit("A0_ploidy", "legacy", legacy, "none"),
  A0_heat = read_fit("A0_heat", "legacy", legacy, "none"),
  A1_ploidy_mn_chisq = read_fit("A1_ploidy_mn_chisq", "legacy", legacy, "MN"),
  A1_heat_mn_chisq = read_fit("A1_heat_mn_chisq", "legacy", legacy, "MN"),
  A2_interaction = read_fit("A2_interaction", "legacy", legacy, "none"),
  A2_interaction_mn_chisq = read_fit("A2_interaction_mn_chisq", "legacy", legacy, "MN"),
  A4_ploidy_in_control = read_fit("A4_ploidy_in_control", "legacy", legacy, "none"),
  A4_ploidy_in_control_mn_chisq = read_fit("A4_ploidy_in_control_mn_chisq", "legacy", legacy, "MN"),
  A4_ploidy_in_heat = read_fit("A4_ploidy_in_heat", "legacy", legacy, "none"),
  A4_ploidy_in_heat_mn_chisq = read_fit("A4_ploidy_in_heat_mn_chisq", "legacy", legacy, "MN"))

# ---- 1. Level under relabelling, all 22 libraries ------------------------------------------------
log_msg("part 1: level under relabelling, all 22 libraries (%d samples of %d eligible CpGs, %d relabellings each)",
        length(sample_seeds), n_sample, n_relabel)
samples <- list()
for (k in seq_along(sample_seeds)) {
  set.seed(sample_seeds[k])
  pool <- setdiff(legacy$eligible, unlist(samples))
  samples[[k]] <- sort(sample(pool, n_sample))
}
X_red_ploidy <- cbind(1, heat)
pooled_level <- rbindlist(lapply(seq_along(samples), function(k) {
  rows <- samples[[k]]
  cs <- legacy$cs[rows, , drop = FALSE]; ns <- legacy$ns[rows, , drop = FALSE]
  dev_red <- vapply(seq_along(rows), function(i) quiet_glm(X_red_ploidy, cs[i, ], ns[i, ])$deviance, numeric(1))
  labels <- pooled_labellings(n_relabel, labelling_seeds[k])
  started <- Sys.time()
  res <- rbindlist(lapply(seq_len(ncol(labels)), function(j) {
    z <- labels[, j]
    X_add <- cbind(1, z, heat); X_fac <- cbind(1, z, heat, z * heat)
    tests <- t(vapply(seq_along(rows), function(i) {
      y <- cs[i, ]; w <- ns[i, ]
      c(one_test(y, w, X_add, NULL, dev_red[i])[1:2], one_test(y, w, X_fac, X_add)[1:2])
    }, numeric(4)))
    if (j == 1) {
      check_against(tests[, 2], fits$A1_ploidy_mn_chisq$pvalue[rows], sprintf("part 1, sample %d, ploidy corrected", k))
      check_against(tests[, 1], fits$A0_ploidy$pvalue[rows], sprintf("part 1, sample %d, ploidy uncorrected", k))
      check_against(tests[, 4], fits$A2_interaction_mn_chisq$pvalue[rows], sprintf("part 1, sample %d, interaction corrected", k))
      check_against(tests[, 3], fits$A2_interaction$pvalue[rows], sprintf("part 1, sample %d, interaction uncorrected", k))
    }
    s <- c(shares(tests[, 1], p_levels), shares(tests[, 2], p_levels), shares(tests[, 3], p_levels), shares(tests[, 4], p_levels))
    names(s) <- paste0(rep(c("ploidy_uncorrected", "ploidy_corrected", "interaction_uncorrected", "interaction_corrected"), each = 3),
                       "_p_below_", level_names(p_levels))
    as.data.table(as.list(s))
  }))
  log_msg("  sample %d: %d labellings in %.1f min", k, ncol(labels), as.numeric(difftime(Sys.time(), started, units = "mins")))
  cbind(data.table(sample = k, sample_seed = sample_seeds[k], labelling_seed = labelling_seeds[k], labelling = seq_len(ncol(labels)) - 1L,
                   observed = seq_len(ncol(labels)) == 1L, exchanges_in_heat_cohort = exchanges_in_heat_cohort(labels),
                   direct_exchanges_in_heat_cohort = direct_exchanges_in_heat_cohort(labels),
                   high_level_libraries_split = high_split(labels)), res)
}))
write_table(pooled_level, out_file("relabel_level_pooled.csv"))

# ---- 2. Level under regrouping within one cohort -------------------------------------------------
log_msg("part 2: level under regrouping within one cohort (first sample of part 1)")
within_level <- rbindlist(lapply(c("heat", "control"), function(cohort_i) {
  rows_lib <- if (cohort_i == "heat") heat_rows else control_rows
  g <- if (cohort_i == "heat") cohort_groupings(rows_lib) else cohort_groupings(rows_lib, n_control, control_seed)
  rows <- samples[[1]]
  cs <- legacy$cs[rows, rows_lib, drop = FALSE]; ns <- legacy$ns[rows, rows_lib, drop = FALSE]
  analysis <- if (cohort_i == "heat") "A4_ploidy_in_heat" else "A4_ploidy_in_control"
  high_here <- high[rows_lib]
  res <- rbindlist(lapply(seq_len(ncol(g$labels)), function(j) {
    z <- g$labels[, j]
    tests <- run_tests(cs, ns, seq_along(rows), cbind(1, z), NULL)
    if (j == 1) {
      check_against(tests[, "p_corrected"], fits[[paste0(analysis, "_mn_chisq")]]$pvalue[rows], paste("part 2,", analysis, "corrected"))
      check_against(tests[, "p_uncorrected"], fits[[analysis]]$pvalue[rows], paste("part 2,", analysis, "uncorrected"))
    }
    s <- c(shares(tests[, "p_uncorrected"], p_levels), shares(tests[, "p_corrected"], p_levels))
    names(s) <- paste0(rep(c("uncorrected", "corrected"), each = 3), "_p_below_", level_names(p_levels))
    moved <- sum(z != trip[rows_lib])
    cbind(data.table(cohort = cohort_i, analysis = analysis, grouping = j - 1L, observed = j == 1L,
                     n_groupings_distinct = g$n_distinct, exchanges = min(moved, length(z) - moved) / 2,
                     high_level_libraries_in_one_group = length(unique(z[high_here])) == 1), as.data.table(as.list(s)))
  }))
  log_msg("  %s cohort: %d groupings of %d distinct", cohort_i, nrow(res), g$n_distinct)
  res
}))
write_table(within_level, out_file("relabel_level_within_cohort.csv"))

# Summary of parts 1 and 2: one row per design, test, kind of relabelling and p level.
summarise_level <- function(d, design, test_column_stem, reference_analysis, subset_name, keep, by_sample) {
  rbindlist(lapply(seq_along(p_levels), function(l) {
    col <- paste0(test_column_stem, "_p_below_", level_names(p_levels)[l])
    obs <- d[observed == TRUE]
    rel <- d[observed == FALSE][keep[d$observed == FALSE]]
    observed_of_row <- if (by_sample) obs[[col]][match(rel$sample, obs$sample)] else rep(obs[[col]], nrow(rel))
    ref <- if (is.na(reference_analysis)) NA_real_ else full_share(reference_analysis, p_levels[l])
    data.table(design = design, test = test_column_stem, relabellings = subset_name, p_level = p_levels[l], n_relabellings = nrow(rel),
               observed_share_in_sample = paste(signif(obs[[col]], 4), collapse = "; "),
               mean_share = mean(rel[[col]]), median_share = stats::median(rel[[col]]),
               q05_share = unname(stats::quantile(rel[[col]], 0.05)), q95_share = unname(stats::quantile(rel[[col]], 0.95)),
               max_share = max(rel[[col]]), n_at_or_above_observed = sum(rel[[col]] >= observed_of_row - 1e-12),
               reference_analysis = reference_analysis, reference_share_full_data = ref,
               n_at_or_above_reference = if (is.na(ref)) NA_integer_ else sum(rel[[col]] >= ref - 1e-12))
  }))
}
p_rel <- pooled_level[observed == FALSE]
pooled_subsets <- list("all" = rep(TRUE, nrow(p_rel)),
                       "heat-cohort grouping kept or changed by one exchange" = p_rel$exchanges_in_heat_cohort <= 1,
                       "two exchanges in the heat cohort" = p_rel$exchanges_in_heat_cohort == 2,
                       "high-level libraries split 3:2" = p_rel$high_level_libraries_split == "3:2",
                       "high-level libraries split 4:1 or 5:0" = p_rel$high_level_libraries_split == "4:1 or 5:0")
expand <- function(keep_rel, d) { k <- rep(FALSE, nrow(d)); k[d$observed == FALSE] <- keep_rel; k }
level_summary <- rbindlist(c(
  unlist(lapply(names(pooled_subsets), function(nm) {
    keep <- expand(pooled_subsets[[nm]], pooled_level)
    list(summarise_level(pooled_level, "all 22 libraries, ploidy (heat as covariate)", "ploidy_uncorrected", "A0_ploidy", nm, keep, TRUE),
         summarise_level(pooled_level, "all 22 libraries, ploidy (heat as covariate)", "ploidy_corrected", "A1_ploidy_mn_chisq", nm, keep, TRUE),
         summarise_level(pooled_level, "all 22 libraries, ploidy by heat interaction", "interaction_uncorrected", "A2_interaction", nm, keep, TRUE),
         summarise_level(pooled_level, "all 22 libraries, ploidy by heat interaction", "interaction_corrected", "A2_interaction_mn_chisq", nm, keep, TRUE))
  }), recursive = FALSE),
  # The cohort contrast cannot be relabelled; its share is set against the relabelled ploidy contrast.
  list(summarise_level(pooled_level, "all 22 libraries, relabelled ploidy against the heat contrast", "ploidy_corrected", "A1_heat_mn_chisq", "all",
                       expand(rep(TRUE, nrow(p_rel)), pooled_level), TRUE)),
  unlist(lapply(c("heat", "control"), function(cohort_i) {
    d <- within_level[cohort == cohort_i]
    rel <- d[observed == FALSE]
    subsets <- list("all" = rep(TRUE, nrow(rel)),
                    "high-level libraries of the cohort in one group" = rel$high_level_libraries_in_one_group,
                    "high-level libraries of the cohort split" = !rel$high_level_libraries_in_one_group)
    an <- d$analysis[1]
    unlist(lapply(names(subsets), function(nm) {
      keep <- expand(subsets[[nm]], d)
      list(summarise_level(d, paste("ploidy among the", if (cohort_i == "heat") "heat-treated libraries" else "controls"), "uncorrected", an, nm, keep, FALSE),
           summarise_level(d, paste("ploidy among the", if (cohort_i == "heat") "heat-treated libraries" else "controls"), "corrected", paste0(an, "_mn_chisq"), nm, keep, FALSE))
    }), recursive = FALSE)
  }), recursive = FALSE)))
write_table(level_summary, out_file("relabel_level_summary.csv"))

# ---- 3. Length of the corrected lists under relabelling -----------------------------------------
log_msg("part 3: corrected lists under relabelling (all %d tested CpGs)", n_tested)
list_designs <- list(
  pooled = list(analysis = "A1_ploidy_mn_chisq", rows_lib = seq_len(nrow(lib)), covariate = heat,
                labels = pooled_labellings(999L, list_seed)[, seq_len(n_list_pooled + 1L)]),
  heat_cohort = list(analysis = "A4_ploidy_in_heat_mn_chisq", rows_lib = heat_rows, covariate = NULL,
                     labels = cohort_groupings(heat_rows)$labels),
  control_cohort = list(analysis = "A4_ploidy_in_control_mn_chisq", rows_lib = control_rows, covariate = NULL,
                        labels = cohort_groupings(control_rows, n_control, control_seed)$labels))
list_rows <- rbindlist(lapply(names(list_designs), function(nm) {
  d <- list_designs[[nm]]
  fit <- fits[[d$analysis]]
  p_star <- max(fit$pvalue[fit$qvalue < config$dml_qvalue])
  cs <- legacy$cs[, d$rows_lib, drop = FALSE]; ns <- legacy$ns[, d$rows_lib, drop = FALSE]
  rs_c <- rowSums(cs); rs_n <- rowSums(ns)
  X_red <- if (is.null(d$covariate)) NULL else cbind(1, d$covariate)
  started <- Sys.time()
  res <- rbindlist(lapply(seq_len(ncol(d$labels)), function(j) {
    z <- d$labels[, j]
    c1 <- drop(cs %*% z); n1 <- drop(ns %*% z)
    effect <- 100 * (c1 / n1 - (rs_c - c1) / (rs_n - n1))
    above <- which(abs(effect) > config$dml_difference)
    X_full <- if (is.null(d$covariate)) cbind(1, z) else cbind(1, z, d$covariate)
    tests <- run_tests(cs, ns, above, X_full, X_red)
    called <- tests[, "p_corrected"] <= p_star
    if (j == 1) {
      check_against(tests[, "p_corrected"], fit$pvalue[above], paste("part 3,", d$analysis))
      n_dml <- sum(fit$qvalue < config$dml_qvalue & abs(fit$meth.diff) > config$dml_difference)
      if (sum(called) != n_dml) stop("Part 3, ", d$analysis, ": ", sum(called), " CpGs called at the p threshold, ", n_dml, " DML in the fit.", call. = FALSE)
    }
    ps <- vapply(p_levels_lists, function(l) sum(tests[, "p_corrected"] < l), numeric(1))
    names(ps) <- paste0("n_above_cutoff_corrected_p_below_", level_names(p_levels_lists))
    moved <- sum(z != trip[d$rows_lib])
    cbind(data.table(design = nm, analysis = d$analysis, labelling = j - 1L, observed = j == 1L,
                     exchanges_in_heat_cohort = if (nm == "pooled") exchanges_in_heat_cohort(matrix(z))
                                                else if (nm == "heat_cohort") min(moved, length(z) - moved) / 2 else NA_real_,
                     direct_exchanges_in_heat_cohort = if (nm == "pooled") direct_exchanges_in_heat_cohort(matrix(z)) else NA_real_,
                     p_threshold_of_observed_list = p_star, n_above_cutoff = length(above), n_called = sum(called),
                     n_called_correction_not_applied = sum(called & not_corrected(tests[, "dispersion"])),
                     n_called_on_their_own = sum(tests[, "p_corrected"] < p_single),
                     n_called_on_their_own_correction_not_applied = sum(tests[, "p_corrected"] < p_single & not_corrected(tests[, "dispersion"])),
                     median_dispersion_above_cutoff = stats::median(tests[, "dispersion"], na.rm = TRUE)),
          as.data.table(as.list(ps)))
  }))
  log_msg("  %s (%s): %d labellings in %.1f min; observed list %d, relabelled median %g", nm, d$analysis, nrow(res),
          as.numeric(difftime(Sys.time(), started, units = "mins")), res$n_called[1], stats::median(res$n_called[-1]))
  res
}))
write_table(list_rows, out_file("relabel_corrected_lists.csv"))

list_summary <- rbindlist(lapply(split(list_rows, list_rows$design), function(d) {
  o <- d[observed == TRUE]; r <- d[observed == FALSE]
  q <- function(x, p) unname(stats::quantile(x, p))
  share_cols <- grep("^n_above_cutoff_corrected_p_below_", names(d), value = TRUE)
  base <- data.table(design = o$design, analysis = o$analysis, n_relabellings = nrow(r), p_threshold_of_observed_list = o$p_threshold_of_observed_list,
                     observed_n_above_cutoff = o$n_above_cutoff, observed_n_in_shares = o$n_above_cutoff, observed_n_called = o$n_called,
                     observed_n_called_correction_not_applied = o$n_called_correction_not_applied,
                     observed_n_called_on_their_own = o$n_called_on_their_own,
                     relabelled_n_called_median = stats::median(r$n_called), relabelled_n_called_q05 = q(r$n_called, 0.05),
                     relabelled_n_called_q95 = q(r$n_called, 0.95), relabelled_n_called_max = max(r$n_called),
                     relabellings_n_called_at_or_above_observed = sum(r$n_called >= o$n_called),
                     relabelled_share_of_calls_correction_not_applied = if (sum(r$n_called) > 0) sum(r$n_called_correction_not_applied) / sum(r$n_called) else NA_real_,
                     relabelled_n_called_on_their_own_mean = mean(r$n_called_on_their_own),
                     relabelled_n_called_on_their_own_median = stats::median(r$n_called_on_their_own),
                     relabelled_n_called_on_their_own_max = max(r$n_called_on_their_own),
                     relabellings_on_their_own_at_or_above_observed = sum(r$n_called_on_their_own >= o$n_called_on_their_own),
                     expected_on_their_own_per_labelling_if_level_held = config$dml_qvalue,
                     observed_median_dispersion_above_cutoff = o$median_dispersion_above_cutoff,
                     relabelled_median_dispersion_above_cutoff = stats::median(r$median_dispersion_above_cutoff))
  for (cc in share_cols) {
    obs_share <- o[[cc]] / o$n_above_cutoff; rel_share <- r[[cc]] / r$n_above_cutoff
    base[, (sub("^n_above_cutoff_", "observed_share_", cc)) := obs_share]
    base[, (sub("^n_above_cutoff_", "relabelled_mean_share_", cc)) := mean(rel_share)]
    base[, (sub("^n_above_cutoff_", "relabellings_share_at_or_above_observed_", cc)) := sum(rel_share >= obs_share - 1e-12)]
  }
  base
}), fill = TRUE)
# The heat contrast cannot be relabelled. Its DML (step 05, the draft's list) with
# a corrected p-value below each level, set against the relabelled ploidy lists.
# The shares of this row are over the DML (observed_n_in_shares), which are all
# but one of the CpGs above the cutoff (observed_n_above_cutoff; the other one
# has q >= 0.01, issues.md, I-05).
heat_dml <- which(fits$A0_heat$qvalue < config$dml_qvalue & abs(fits$A0_heat$meth.diff) > config$dml_difference)
pooled_rel <- list_rows[design == "pooled" & observed == FALSE]
heat_row <- data.table(design = "heat contrast, the draft's heat DML, against the relabelled ploidy lists of the pooled design",
                       analysis = "A1_heat_mn_chisq", n_relabellings = nrow(pooled_rel),
                       observed_n_above_cutoff = sum(abs(fits$A0_heat$meth.diff) > config$dml_difference),
                       observed_n_in_shares = length(heat_dml))
for (l in seq_along(p_levels_lists)) {
  nm <- level_names(p_levels_lists)[l]
  obs_share <- mean(fits$A1_heat_mn_chisq$pvalue[heat_dml] < p_levels_lists[l])
  rel_share <- pooled_rel[[paste0("n_above_cutoff_corrected_p_below_", nm)]] / pooled_rel$n_above_cutoff
  heat_row[, (paste0("observed_share_corrected_p_below_", nm)) := obs_share]
  heat_row[, (paste0("relabelled_mean_share_corrected_p_below_", nm)) := mean(rel_share)]
  heat_row[, (paste0("relabellings_share_at_or_above_observed_corrected_p_below_", nm)) := sum(rel_share >= obs_share - 1e-12)]
}
list_summary <- rbind(list_summary, heat_row, fill = TRUE)
write_table(list_summary, out_file("relabel_corrected_lists_summary.csv"))

# ---- 4. Where the correction is not applied -------------------------------------------------------
log_msg("part 4: CpGs at which the dispersion estimate is not above 1")
corrected <- read_counts("corrected")
fits_corrected <- list(
  A0_ploidy = read_fit("A0_ploidy", "corrected", corrected, "none"),
  A0_heat = read_fit("A0_heat", "corrected", corrected, "none"),
  A1_ploidy_mn_chisq = read_fit("A1_ploidy_mn_chisq", "corrected", corrected, "MN"),
  A1_heat_mn_chisq = read_fit("A1_heat_mn_chisq", "corrected", corrected, "MN"),
  A2_interaction = read_fit("A2_interaction", "corrected", corrected, "none"),
  A2_interaction_mn_chisq = read_fit("A2_interaction_mn_chisq", "corrected", corrected, "MN"))
correction_specs <- list(
  A1_ploidy_mn_chisq = list(rows_lib = seq_len(nrow(lib)), X_full = cbind(1, trip, heat), X_red = cbind(1, heat), modes = c("legacy", "corrected")),
  A1_heat_mn_chisq = list(rows_lib = seq_len(nrow(lib)), X_full = cbind(1, heat, trip), X_red = cbind(1, trip), modes = c("legacy", "corrected")),
  A2_interaction_mn_chisq = list(rows_lib = seq_len(nrow(lib)), X_full = cbind(1, trip * heat, trip, heat), X_red = cbind(1, trip, heat),
                                 modes = c("legacy", "corrected")),
  A4_ploidy_in_control_mn_chisq = list(rows_lib = control_rows, X_full = cbind(1, trip[control_rows]), X_red = NULL, modes = "legacy"),
  A4_ploidy_in_heat_mn_chisq = list(rows_lib = heat_rows, X_full = cbind(1, trip[heat_rows]), X_red = NULL, modes = "legacy"))
dispersion_cache <- list()
not_applied <- rbindlist(unlist(lapply(names(correction_specs), function(id) {
  sp <- correction_specs[[id]]
  lapply(sp$modes, function(mode_i) {
    cnt <- if (mode_i == "legacy") legacy else corrected
    fit <- if (mode_i == "legacy") fits[[id]] else fits_corrected[[id]]
    cs <- cnt$cs[, sp$rows_lib, drop = FALSE]; ns <- cnt$ns[, sp$rows_lib, drop = FALSE]
    # The dispersion depends on the full model only: A1_ploidy and A1_heat share it.
    design_key <- paste(mode_i, paste(sort(apply(sp$X_full, 2, paste, collapse = "")), collapse = "|"))
    if (is.null(dispersion_cache[[design_key]])) {
      dispersion_cache[[design_key]] <<- run_tests(cs, ns, cnt$eligible, sp$X_full, sp$X_red)[, "dispersion"]
    }
    phi_eligible <- dispersion_cache[[design_key]]
    listed <- fread(file.path(sens_dir, sprintf("dml_%s_%s.csv", id, mode_i)), select = c("chr", "start"))
    rows <- match(paste(listed$chr, listed$start), paste(cnt$pos$chr, cnt$pos$start))
    if (anyNA(rows)) stop("dml_", id, "_", mode_i, ".csv holds CpGs that are not in the united table.", call. = FALSE)
    tests <- run_tests(cs, ns, rows, sp$X_full, sp$X_red)
    if (length(rows) > 0) check_against(tests[, "p_corrected"], fit$pvalue[rows], paste("part 4,", id, mode_i))
    data.table(analysis_id = id, mode = mode_i, n_libraries = length(sp$rows_lib), residual_df = length(sp$rows_lib) - ncol(sp$X_full),
               n_cpg_eligible = length(cnt$eligible), n_eligible_dispersion_not_above_1 = sum(not_corrected(phi_eligible)),
               share_eligible_dispersion_not_above_1 = mean(not_corrected(phi_eligible)),
               n_dml = length(rows), n_dml_dispersion_not_above_1 = sum(not_corrected(tests[, "dispersion"])),
               median_dispersion_dml = if (length(rows) > 0) stats::median(tests[, "dispersion"]) else NA_real_)
  })
}), recursive = FALSE))
write_table(not_applied, out_file("correction_not_applied.csv"))

# ---- 5. Stability of the corrected lists of the draft's two contrasts ----------------------------------
log_msg("part 5: stability of the corrected lists of A1")
key <- function(cnt) paste(cnt$pos$chr, cnt$pos$start)
legacy_key <- key(legacy); corrected_key <- key(corrected)
stability <- rbindlist(lapply(c("A1_ploidy_mn_chisq", "A1_heat_mn_chisq"), function(id) {
  m <- function(mode_i, measure, value) data.table(analysis_id = id, mode = mode_i, measure = measure, value = as.character(value))
  per_mode <- rbindlist(lapply(c("legacy", "corrected"), function(mode_i) {
    fit <- if (mode_i == "legacy") fits[[id]] else fits_corrected[[id]]
    above <- abs(fit$meth.diff) > config$dml_difference
    bh <- stats::p.adjust(fit$pvalue, method = "BH")
    rbind(m(mode_i, "CpGs tested", nrow(fit)),
          m(mode_i, "CpGs with q < 0.01", sum(fit$qvalue < config$dml_qvalue)),
          m(mode_i, "DML (q < 0.01, above the difference cutoff)", sum(fit$qvalue < config$dml_qvalue & above)),
          m(mode_i, "DML with q < 0.05 in place of 0.01", sum(fit$qvalue < 0.05 & above)),
          m(mode_i, "DML with q < 0.1 in place of 0.01", sum(fit$qvalue < 0.1 & above)),
          m(mode_i, "CpGs with Benjamini-Hochberg adjusted p < 0.01", sum(bh < config$dml_qvalue)),
          m(mode_i, "DML with the Benjamini-Hochberg adjustment in place of SLIM (adjusted p < 0.01)", sum(bh < config$dml_qvalue & above)))
  }))
  fl <- fits[[id]]; fc <- fits_corrected[[id]]
  dml_l <- which(fl$qvalue < config$dml_qvalue & abs(fl$meth.diff) > config$dml_difference)
  dml_c_keys <- corrected_key[fc$qvalue < config$dml_qvalue & abs(fc$meth.diff) > config$dml_difference]
  dropped <- dml_l[!legacy_key[dml_l] %in% dml_c_keys]
  in_c <- match(legacy_key[dropped], corrected_key)
  removed_by_filter <- dropped[is.na(in_c)]
  kept_tested <- dropped[!is.na(in_c)]
  q_c <- fc$qvalue[in_c[!is.na(in_c)]]
  sig_l <- which(fl$qvalue < config$dml_qvalue)
  sig_removed <- sig_l[!legacy_key[sig_l] %in% corrected_key]
  rng <- function(x) if (length(x) > 0) sprintf("%.2g to %.2g", min(x), max(x)) else "none"
  rbind(per_mode,
        m("both", "legacy DML not in the corrected-mode list", length(dropped)),
        m("both", "of these, not tested in corrected mode (removed by the coverage filter)", length(removed_by_filter)),
        m("both", "of these, tested in corrected mode with q >= 0.01", sum(q_c >= config$dml_qvalue)),
        m("both", "q-values of the latter in legacy mode", rng(fl$qvalue[kept_tested[q_c >= config$dml_qvalue]])),
        m("both", "q-values of the latter in corrected mode", rng(q_c[q_c >= config$dml_qvalue])),
        m("both", "corrected-mode DML not in the legacy list", sum(!dml_c_keys %in% legacy_key[dml_l])),
        m("legacy", paste0("CpGs with q < 0.01 on the mitochondrial genome (", config$mito_sequence, ")"), sum(legacy$pos$chr[sig_l] == config$mito_sequence)),
        m("both", "legacy CpGs with q < 0.01 that the coverage filter removes", length(sig_removed)),
        m("both", paste0("of these, on the mitochondrial genome (", config$mito_sequence, ")"), sum(legacy$pos$chr[sig_removed] == config$mito_sequence)),
        m("both", "of these, above the difference cutoff", sum(abs(fl$meth.diff[sig_removed]) > config$dml_difference)))
}))
write_table(stability, out_file("corrected_lists_stability.csv"))

# ---- 6. The corrected ploidy CpGs inside each cohort ----------------------------------------------
log_msg("part 6: the corrected ploidy CpGs inside each cohort")
is_dml <- function(fit) fit$qvalue < config$dml_qvalue & abs(fit$meth.diff) > config$dml_difference
p_rows <- which(is_dml(fits$A1_ploidy_mn_chisq))
shared_rows <- which(is_dml(fits$A0_ploidy) & is_dml(fits$A0_heat))
by_cohort <- data.table(chr = legacy$pos$chr[p_rows], start = legacy$pos$start[p_rows],
                        effect_size_all = fits$A1_ploidy_mn_chisq$meth.diff[p_rows],
                        pvalue_legacy = fits$A1_ploidy_mn_chisq$pvalue[p_rows],
                        qvalue_legacy = fits$A1_ploidy_mn_chisq$qvalue[p_rows],
                        qvalue_corrected_mode = fits_corrected$A1_ploidy_mn_chisq$qvalue[match(legacy_key[p_rows], corrected_key)],
                        effect_size_heat_cohort = fits$A4_ploidy_in_heat_mn_chisq$meth.diff[p_rows],
                        qvalue_heat_cohort = fits$A4_ploidy_in_heat_mn_chisq$qvalue[p_rows],
                        called_heat_cohort = is_dml(fits$A4_ploidy_in_heat_mn_chisq)[p_rows],
                        effect_size_control_cohort = fits$A4_ploidy_in_control_mn_chisq$meth.diff[p_rows],
                        qvalue_control_cohort = fits$A4_ploidy_in_control_mn_chisq$qvalue[p_rows],
                        called_control_cohort = is_dml(fits$A4_ploidy_in_control_mn_chisq)[p_rows],
                        in_corrected_heat_list = is_dml(fits$A1_heat_mn_chisq)[p_rows],
                        pvalue_corrected_heat = fits$A1_heat_mn_chisq$pvalue[p_rows],
                        in_draft_ploidy_list = is_dml(fits$A0_ploidy)[p_rows],
                        in_draft_shared_loci = p_rows %in% shared_rows)
write_table(by_cohort, out_file("corrected_ploidy_dml_by_cohort.csv"))

# ---- 7. Named loci of the manuscript values map ----------------------------------------------------
named <- data.table(
  label = c("CAD protein", "CAD protein", "nesprin-1", "adenylate kinase 9", "adenylate kinase 9"),
  chr = c("NC_047564.1", "NC_047564.1", "NC_047561.1", "NC_047563.1", "NC_047563.1"),
  start = c(8425110L, 8425110L, 28607016L, 5569229L, 5569232L),
  analysis_uncorrected = c("A4_ploidy_in_heat", "A2_interaction", "A2_interaction", "A2_interaction", "A2_interaction"),
  values_map_rows = c("D4.3-016", "D4.9-008", "D4.9-010", "D4.9-014", "D4.9-014"))
named[, analysis_corrected := paste0(analysis_uncorrected, "_mn_chisq")]
named[, row := match(paste(chr, start), legacy_key)]
if (anyNA(named$row)) stop("A named locus is not among the tested CpGs.", call. = FALSE)
named[, p_uncorrected := vapply(seq_len(.N), function(k) fits[[analysis_uncorrected[k]]]$pvalue[row[k]], numeric(1))]
named[, p_corrected := vapply(seq_len(.N), function(k) fits[[analysis_corrected[k]]]$pvalue[row[k]], numeric(1))]
write_table(named[, !"row"], out_file("named_loci.csv"))

writeLines(capture.output(print(utils::sessionInfo())), out_file("session_info.txt"))
log_msg("dml_corrected_test finished in %.1f min; outputs in %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir)

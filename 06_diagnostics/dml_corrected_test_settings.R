# dml_corrected_test_settings.R
#
# Purpose   The relabelling check of 06_diagnostics/dml_corrected_test.R and
#           dml_corrected_test_5x.R for every corrected control-only and
#           interaction list of the optional analyses (author's request of
#           2026-11-24, "the other settings": upper filter, C-to-T exclusion, 50
#           points, q < 0.05). One run checks one fit (argument or DML_FIT), one
#           of the eight corrected fits with option ctrl or inter:
#             mnchisq-ctrl, corrected-mnchisq-ctrl, cov5-mnchisq-ctrl, corrected-cov5-mnchisq-ctrl,
#             mnchisq-inter, corrected-mnchisq-inter, cov5-mnchisq-inter, corrected-cov5-mnchisq-inter
#           and with it the eight analyses that share that fit (ctexcl, d50 and q05,
#           each on or off), because those settings act after the test: the
#           exclusion removes CpGs (p- and q-values unchanged), d50 raises the
#           difference threshold to 50 points, q05 the q-value cutoff to 0.05.
#           For each relabelling the test is made once, at the CpGs whose
#           difference exceeds 20 points under that labelling, and each variant is
#           counted from it:
#             called        corrected p at or below the threshold of the real list
#                           (with a relative tolerance of 1e-9, because the CpG that sets
#                           the threshold is itself recomputed to about 1e-14)
#                           of that variant (the largest p-value with q below its
#                           cutoff in the real fit), difference above its threshold,
#                           not excluded;
#             on_their_own  corrected p below cutoff / (CpGs tested), which passes
#                           the cutoff whatever the other CpGs do.
#           Labellings: those of the 10x and 5x checks (control-only: 150
#           regroupings of the controls, seed 802; interaction: the first 200
#           labellings of the permutation reference of dml_sensitivity.R, seed
#           20261119, ploidy reassigned within cohort).
#           Checks: with the real labels, the p-values and differences equal the
#           pipeline's (relative difference at most 1e-6), and each variant calls
#           exactly the DML count of analysis_comparison.csv. The run stops before
#           the relabellings if either fails.
#           A diagnostics script: nothing in 01_pipeline/ reads its outputs; not in
#           diagnostics.csv (it reads git-ignored tables of optional analyses).
#
# Run       WGBS_CORES=4 Rscript 06_diagnostics/dml_corrected_test_settings.R <fit id>
#           (DML5X_TRIAL=<n>: n relabellings, nothing written)
#
# Reads     04_output/not_adopted/<fit>/05_differential_methylation/tables/
#             united_cpg_<fit>.csv.gz, diffmeth_<ploidy|interaction>_<fit>.csv.gz   (git-ignored)
#           04_output/not_adopted/<fit>-ctexcl/05_differential_methylation/tables/ct_excluded_cpgs_<fit>-ctexcl.csv.gz
#           04_output/not_adopted/comparisons/analysis_comparison.csv
#
# Writes    06_diagnostics/dml_corrected_test_settings/
#             relabel_<fit>.csv     one row per labelling and variant
#             summary_<fit>.csv     one row per variant
#             summary.csv           the rows of every fit run so far (rebuilt from the summary_<fit>.csv files)

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))
if (!requireNamespace("R.utils", quietly = TRUE)) stop("Missing R package: R.utils (needed to read .csv.gz files).", call. = FALSE)

fits_allowed <- c("mnchisq-ctrl", "corrected-mnchisq-ctrl", "cov5-mnchisq-ctrl", "corrected-cov5-mnchisq-ctrl",
                  "mnchisq-inter", "corrected-mnchisq-inter", "cov5-mnchisq-inter", "corrected-cov5-mnchisq-inter")
args <- commandArgs(trailingOnly = TRUE)
fit_id <- if (length(args)) args[1] else Sys.getenv("DML_FIT", "")
if (!fit_id %in% fits_allowed) stop("Give one of: ", paste(fits_allowed, collapse = ", "), call. = FALSE)
design <- if (grepl("-ctrl$", fit_id)) "ctrl" else "inter"
t0 <- start_step("dml_corrected_test_settings", paste("relabelling check of the corrected lists of fit", fit_id))

control_seed <- 802; n_control <- 150; list_seed <- 20261119; n_list_pooled <- 200
equal_tolerance <- 1e-6
cores <- max(1L, as.integer(Sys.getenv("WGBS_CORES", "1")))
trial <- suppressWarnings(as.integer(Sys.getenv("DML5X_TRIAL", "")))
writing <- is.na(trial)
if (!identical(config$dml_qvalue, 0.01) || !identical(as.numeric(config$dml_difference), 20)) {
  stop("Run with the paper's thresholds: unset WGBS_ANALYSIS and the optional settings.", call. = FALSE)
}
out_dir <- here("06_diagnostics", "dml_corrected_test_settings")
if (writing) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tab <- function(id, file) here("04_output", "not_adopted", id, "05_differential_methylation", "tables", file)
need <- function(f) { if (!file.exists(f)) stop("Missing ", f, ". Run that optional analysis (or its step 05) first.", call. = FALSE); f }

# ---- Inputs ---------------------------------------------------------------------------------
lib <- read_libraries()
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
cohort_rows <- split(seq_len(nrow(lib)), lib$cohort)
control_rows <- which(lib$treatment == "control")

u <- fread(need(tab(fit_id, sprintf("united_cpg_%s.csv.gz", fit_id))))
cs <- as.matrix(u[, paste0("numCs", seq_len(nrow(lib))), with = FALSE]) + 0
ns <- cs + as.matrix(u[, paste0("numTs", seq_len(nrow(lib))), with = FALSE])
pos <- u[, .(chr, start)]; n_tested <- nrow(u); rm(u); invisible(gc())
contrast <- if (design == "ctrl") "ploidy" else "interaction"
fit <- fread(need(tab(fit_id, sprintf("diffmeth_%s_%s.csv.gz", contrast, fit_id))), select = c("chr", "start", "pvalue", "qvalue", "meth.diff", "united_row"))
if (nrow(fit) != n_tested || !identical(fit$united_row, seq_len(n_tested)) || !identical(fit$chr, pos$chr) || !identical(fit$start, pos$start)) {
  stop("The ", contrast, " table does not follow the united table.", call. = FALSE)
}
excl_tab <- fread(need(tab(paste0(fit_id, "-ctexcl"), sprintf("ct_excluded_cpgs_%s-ctexcl.csv.gz", fit_id))), select = c("chr", "start"))
excluded <- paste(pos$chr, pos$start) %in% paste(excl_tab$chr, excl_tab$start)
log_msg("fit %s (%s): %d CpGs tested, %d of them on the C-to-T list", fit_id, design, n_tested, sum(excluded))

# The eight variants and the pipeline's DML counts.
cmp <- fread(need(here("04_output", "not_adopted", "comparisons", "analysis_comparison.csv")))
variants <- CJ(ctexcl = c(FALSE, TRUE), d50 = c(FALSE, TRUE), q05 = c(FALSE, TRUE))
variants[, id := paste0(fit_id, ifelse(ctexcl, "-ctexcl", ""), ifelse(d50, "-d50", ""), ifelse(q05, "-q05", ""))]
variants[, `:=`(difference = ifelse(d50, 50, 20), cutoff = ifelse(q05, 0.05, 0.01))]
count_col <- if (design == "ctrl") "ploidy_dml" else "interaction_dml"
variants[, pipeline_dml := cmp[[count_col]][match(id, cmp$analysis)]]
if (anyNA(variants$pipeline_dml)) stop("analysis_comparison.csv lacks: ", paste(variants[is.na(pipeline_dml), id], collapse = ", "), call. = FALSE)
variants[, p_star := vapply(cutoff, function(q) { p <- fit$pvalue[fit$qvalue < q]; if (length(p)) max(p) else 0 }, numeric(1))]
variants[, p_single := cutoff / n_tested]

# ---- The test (copied from dml_corrected_test_5x.R) -------------------------------------------
quiet_glm <- function(X, y, w) {
  withCallingHandlers(stats::glm.fit(X, y / w, weights = w, family = stats::binomial()),
                      warning = function(cond) {
                        msg <- conditionMessage(cond)
                        if (!grepl("fitted probabilities numerically 0 or 1|algorithm did not converge", msg)) stop("Unexpected warning: ", msg, call. = FALSE)
                        invokeRestart("muffleWarning")
                      })
}
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
run_tests <- function(cs, ns, rows, X_full, X_red) {
  if (length(rows) == 0) return(matrix(numeric(0), 0, 3, dimnames = list(NULL, c("p_uncorrected", "p_corrected", "dispersion"))))
  t(vapply(rows, function(i) one_test(cs[i, ], ns[i, ], X_full, X_red), numeric(3)))
}
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
pooled_labellings <- function(n, seed) {
  set.seed(seed)
  cbind(trip, vapply(seq_len(n), function(draw) {
    z <- trip; for (b in cohort_rows) z[b] <- trip[b][sample.int(length(b))]; z
  }, integer(length(trip))))
}

# ---- The design ----------------------------------------------------------------------------
if (design == "ctrl") {
  rows_lib <- control_rows
  labels <- cohort_groupings(control_rows, n_control, control_seed)$labels
  X_full <- function(z) cbind(1, z); X_red <- function(z) NULL
  effect <- function(z) {
    c_ <- cs[, control_rows, drop = FALSE]; n_ <- ns[, control_rows, drop = FALSE]
    c1 <- drop(c_ %*% z); n1 <- drop(n_ %*% z)
    100 * (c1 / n1 - (rowSums(c_) - c1) / (rowSums(n_) - n1))
  }
} else {
  rows_lib <- seq_len(nrow(lib))
  labels <- pooled_labellings(999L, list_seed)[, seq_len(n_list_pooled + 1L)]
  X_full <- function(z) cbind(1, z * heat, z, heat); X_red <- function(z) cbind(1, z, heat)
  effect <- function(z) {
    pct <- function(g) { j <- which(g); 100 * rowSums(cs[, j, drop = FALSE]) / rowSums(ns[, j, drop = FALSE]) }
    (pct(z == 1 & heat == 1) - pct(z == 1 & heat == 0)) - (pct(z == 0 & heat == 1) - pct(z == 0 & heat == 0))
  }
}
if (!is.na(trial)) labels <- labels[, seq_len(min(ncol(labels), trial + 1L)), drop = FALSE]
c_lib <- cs[, rows_lib, drop = FALSE]; n_lib <- ns[, rows_lib, drop = FALSE]

evaluate <- function(k) {
  z <- labels[, k]
  eff <- effect(z)
  above <- which(abs(eff) > 20)
  tests <- run_tests(c_lib, n_lib, above, X_full(z), X_red(z))
  p <- tests[, "p_corrected"]; a_eff <- abs(eff[above]); ex <- excluded[above]
  lab_no <- k - 1L; is_obs <- k == 1L
  out <- variants[, .(id, labelling = lab_no, observed = is_obs,
    n_called = vapply(seq_len(.N), function(v) sum(p <= p_star[v] * (1 + 1e-9) & a_eff > difference[v] & !(ctexcl[v] & ex)), numeric(1)),
    n_called_on_their_own = vapply(seq_len(.N), function(v) sum(p < p_single[v] & a_eff > difference[v] & !(ctexcl[v] & ex)), numeric(1)),
    n_above_cutoff = vapply(seq_len(.N), function(v) sum(a_eff > difference[v] & !(ctexcl[v] & ex)), numeric(1)))]
  if (k == 1) attr(out, "check") <- list(p = max(abs(p - fit$pvalue[above]) / pmax(abs(fit$pvalue[above]), .Machine$double.xmin)),
                                         d = max(abs(eff - fit$meth.diff)))
  out
}

# Real labels first: stop before the relabellings if anything disagrees with the pipeline.
started <- Sys.time()
obs <- evaluate(1L); chk <- attr(obs, "check")
if (!is.finite(chk$p) || chk$p > equal_tolerance) stop("Real-label p-values differ from the pipeline's (", signif(chk$p, 3), ").", call. = FALSE)
if (chk$d > 1e-8) stop("Real-label differences differ from the pipeline's (", signif(chk$d, 3), ").", call. = FALSE)
mism <- merge(obs, variants[, .(id, pipeline_dml)], by = "id")[n_called != pipeline_dml]
if (nrow(mism)) stop("Real labels do not reproduce the pipeline's DML counts: ", paste(sprintf("%s %d vs %d", mism$id, mism$n_called, mism$pipeline_dml), collapse = "; "), call. = FALSE)
log_msg("real labels: p-values (%.1e) and differences agree with the pipeline; DML counts reproduced for all 8 variants: %s",
        chk$p, paste(sprintf("%s=%d", sub(paste0("^", fit_id), "", obs$id), obs$n_called), collapse = " "))
res <- parallel::mclapply(2:ncol(labels), evaluate, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
if (any(bad)) stop("Labelling ", which(bad)[1] + 1L, " failed: ", as.character(res[[which(bad)[1]]]), call. = FALSE)
rel <- rbindlist(c(list(obs), res))
log_msg("%d labellings in %.1f min", ncol(labels), as.numeric(difftime(Sys.time(), started, units = "mins")))

q_ <- function(x, p) unname(stats::quantile(x, p))
summ <- rel[, {
  o <- .SD[observed == TRUE]; r <- .SD[observed == FALSE]
  .(fit = fit_id, design = design, n_relabellings = nrow(r), observed_n_called = o$n_called,
    relabelled_n_called_median = stats::median(r$n_called), relabelled_n_called_q05 = q_(r$n_called, 0.05),
    relabelled_n_called_q95 = q_(r$n_called, 0.95), relabelled_n_called_max = max(r$n_called),
    relabellings_at_or_above_observed = sum(r$n_called >= o$n_called),
    observed_n_called_on_their_own = o$n_called_on_their_own,
    relabelled_n_called_on_their_own_median = stats::median(r$n_called_on_their_own),
    relabellings_on_their_own_at_or_above_observed = sum(r$n_called_on_their_own >= o$n_called_on_their_own))
}, by = id]
summ <- merge(variants[, .(id, ctexcl, d50, q05, difference, cutoff, p_star, pipeline_dml)], summ, by = "id")
print(summ[, .(id, observed_n_called, relabelled_n_called_median, relabelled_n_called_q95, relabelled_n_called_max, relabellings_at_or_above_observed, n_relabellings)])

if (writing) {
  write_table(rel, file.path(out_dir, sprintf("relabel_%s.csv", fit_id)))
  write_table(summ, file.path(out_dir, sprintf("summary_%s.csv", fit_id)))
  all <- rbindlist(lapply(list.files(out_dir, pattern = "^summary_.+\\.csv$", full.names = TRUE), fread), fill = TRUE)
  setorder(all, design, fit, id)
  write_table(all, file.path(out_dir, "summary.csv"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
}
log_msg("dml_corrected_test_settings (%s) finished in %.1f min", fit_id, as.numeric(difftime(Sys.time(), t0, units = "mins")))

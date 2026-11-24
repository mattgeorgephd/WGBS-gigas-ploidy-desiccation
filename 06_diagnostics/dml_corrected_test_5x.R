# dml_corrected_test_5x.R
#
# Purpose   The relabelling check of 06_diagnostics/dml_corrected_test.R, made for
#           the two corrected 5x analyses with the new contrasts (author's request
#           of 2026-11-24): cov5-mnchisq-ctrl (the ploidy contrast among the 12
#           control libraries; 120 DML) and cov5-mnchisq-inter (the tested ploidy
#           by heat interaction, all 22 libraries; 1,718 DML). These are the
#           paper's settings with the coverage standard changed to 5x, as the 10x
#           check was made at the paper's settings. Both use overdispersion "MN"
#           with the chi-squared test. Issues: I-06 (item 2), I-07, I-02.
#             A. Level under relabelling. The share of p-values below 0.05, 0.01
#                and 0.001 at two seeded samples of 3,000 eligible CpGs, with the
#                real labels and with relabelled libraries: the controls
#                regrouped (control-only contrast) and the ploidy labels
#                reassigned within each cohort (interaction).
#             B. Length of the lists under relabelling. Among the CpGs whose
#                difference passes the 20-point threshold under a labelling,
#                those with a corrected p-value at or below the largest p-value
#                that has q < 0.01 in the real fit (the list's own threshold),
#                and those with p below 0.01 / (CpGs tested), which pass q < 0.01
#                whatever the other CpGs do.
#           A diagnostics script: nothing in 01_pipeline/ reads its outputs. Not in
#           diagnostics.csv, because it reads tables of optional analyses that are
#           git-ignored (rerun those analyses first; see Reads).
#
# Run       WGBS_CORES=2 Rscript 06_diagnostics/dml_corrected_test_5x.R
#           Optional: DML5X_TRIAL=<n> runs n relabellings of each design and writes
#           nothing (for timing).
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           04_output/not_adopted/cov5-mnchisq-ctrl/05_differential_methylation/tables/
#             united_cpg_cov5-mnchisq-ctrl.csv.gz, diffmeth_ploidy_cov5-mnchisq-ctrl.csv.gz
#           04_output/not_adopted/cov5-mnchisq-inter/05_differential_methylation/tables/
#             united_cpg_cov5-mnchisq-inter.csv.gz, diffmeth_interaction_cov5-mnchisq-inter.csv.gz
#             (all git-ignored: run the two analyses, or step 05 for them, first)
#           04_output/not_adopted/comparisons/05_differential_methylation/tables/dml_counts_by_mode.csv
#           06_diagnostics/library_methylation_level/library_level.csv   column high_level
#
# Writes    06_diagnostics/dml_corrected_test_5x/
#             relabel_level.csv           part A, one row per design, sample and labelling
#             relabel_level_summary.csv   part A summarised per design and p level
#             relabel_lists.csv           part B, one row per design and labelling
#             relabel_lists_summary.csv   part B summarised per design
#             session_info.txt
#
# Notes     - The test is that of dml_corrected_test.R (methylKit 1.28.0 logReg
#             written out in base R; the functions quiet_glm(), one_test() and
#             cohort_groupings(), pooled_labellings() are copied from it
#             unchanged). For the real labels the script checks its p-values
#             against the pipeline's tables (relative difference at most
#             `equal_tolerance`) and that the list's own threshold calls exactly
#             the pipeline's DML; it stops otherwise.
#           - The labellings are those of the 10x check: the same 150 regroupings
#             of the controls (seed 802) and, for the interaction, the first 200
#             labellings of the permutation reference of dml_sensitivity.R (seed
#             20261119) in part B and the 150 of each sample (seeds 707, 708) in
#             part A. Ploidy is reassigned within each cohort; heat cannot be
#             relabelled (I-01).
#           - The interaction difference of a labelling is computed as the
#             pipeline computes it (interaction_groups() of
#             01_pipeline/05_differential_methylation.R): (triploid heat minus
#             triploid control) minus (diploid heat minus diploid control), each
#             group's percentage pooled over its libraries' reads.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))
if (!requireNamespace("R.utils", quietly = TRUE)) stop("Missing R package: R.utils (needed to read .csv.gz files).", call. = FALSE)

t0 <- start_step("dml_corrected_test_5x", "relabelling check of the corrected 5x control-only and interaction lists")

# ---- Settings -------------------------------------------------------------------------
n_sample        <- 3000
sample_seeds    <- c(606, 607)
labelling_seeds <- c(707, 708)
n_relabel       <- 150
control_seed    <- 802
n_control       <- 150
list_seed       <- 20261119
n_list_pooled   <- 200
p_levels        <- c(0.05, 0.01, 0.001)
p_levels_lists  <- c(0.05, 0.01, 0.001, 1e-4)
equal_tolerance <- 1e-6
cores <- max(1L, as.integer(Sys.getenv("WGBS_CORES", "1")))
trial <- suppressWarnings(as.integer(Sys.getenv("DML5X_TRIAL", "")))
trial <- if (is.na(trial)) NA_integer_ else trial

if (!identical(config$dml_qvalue, 0.01) || !identical(as.numeric(config$dml_difference), 20)) {
  stop("Run with the paper's thresholds (q < 0.01, more than 20 points): unset WGBS_ANALYSIS and the optional settings.", call. = FALSE)
}

out_dir <- here("06_diagnostics", "dml_corrected_test_5x")
if (is.na(trial)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, name)
tab <- function(id, file) here("04_output", "not_adopted", id, "05_differential_methylation", "tables", file)
need <- function(f) { if (!file.exists(f)) stop("Missing ", f, ". Run the optional analysis (or its step 05) first.", call. = FALSE); f }

# ---- Inputs ---------------------------------------------------------------------------------
lib <- read_libraries()
trip <- as.integer(lib$ploidy == "triploid")
heat <- as.integer(lib$treatment == "heat")
cohort_rows <- split(seq_len(nrow(lib)), lib$cohort)
heat_rows <- which(lib$treatment == "heat"); control_rows <- which(lib$treatment == "control")
if (!setequal(heat_rows, cohort_rows[[unique(lib$cohort[heat_rows])]]) || length(unique(lib$cohort)) != 2) {
  stop("Expected two cohorts, one of them holding every heat-treated library.", call. = FALSE)
}
level <- fread(need(here("06_diagnostics", "library_methylation_level", "library_level.csv")))
if (!identical(level$library_id, lib$library_id)) stop("library_level.csv does not list the libraries of the library sheet in order.", call. = FALSE)
high <- level$high_level

f_united_ctrl  <- need(tab("cov5-mnchisq-ctrl", "united_cpg_cov5-mnchisq-ctrl.csv.gz"))
f_united_inter <- need(tab("cov5-mnchisq-inter", "united_cpg_cov5-mnchisq-inter.csv.gz"))
md5 <- unname(tools::md5sum(c(f_united_ctrl, f_united_inter)))
if (md5[1] != md5[2]) stop("The united tables of the two analyses differ; both should be the 5x table of their base fit cov5-mnchisq.", call. = FALSE)
log_msg("united table (md5 %s): reading", md5[1])
u <- fread(f_united_ctrl)
cs <- as.matrix(u[, paste0("numCs", seq_len(nrow(lib))), with = FALSE]) + 0
ns <- cs + as.matrix(u[, paste0("numTs", seq_len(nrow(lib))), with = FALSE])
pos <- u[, .(chr, start)]
eligible <- which(dml_eligible(mean_percent_methylation(u), config))
rm(u); invisible(gc())
n_tested <- nrow(pos)
p_single <- config$dml_qvalue / n_tested
log_msg("%d CpGs tested, %d eligible", n_tested, length(eligible))

read_fit <- function(f) {
  x <- fread(need(f), select = c("chr", "start", "pvalue", "qvalue", "meth.diff", "united_row"))
  if (nrow(x) != n_tested || !identical(x$united_row, seq_len(n_tested)) || !identical(x$chr, pos$chr) || !identical(x$start, pos$start)) {
    stop(basename(f), " does not follow the rows of the united table.", call. = FALSE)
  }
  x
}
fits <- list(ctrl = read_fit(tab("cov5-mnchisq-ctrl", "diffmeth_ploidy_cov5-mnchisq-ctrl.csv.gz")),
             inter = read_fit(tab("cov5-mnchisq-inter", "diffmeth_interaction_cov5-mnchisq-inter.csv.gz")))
counts <- fread(need(here("04_output", "not_adopted", "comparisons", "05_differential_methylation", "tables", "dml_counts_by_mode.csv")))
n_dml <- c(ctrl = counts[mode == "cov5-mnchisq-ctrl", ploidy_dml], inter = counts[mode == "cov5-mnchisq-inter", interaction_dml])
for (k in names(fits)) {
  own <- with(fits[[k]], sum(qvalue < config$dml_qvalue & abs(meth.diff) > config$dml_difference))
  if (length(n_dml[[k]]) != 1 || own != n_dml[[k]]) stop("DML count of ", k, " differs from dml_counts_by_mode.csv.", call. = FALSE)
}
log_msg("DML of the real labels: control-only ploidy %d, interaction %d", n_dml[["ctrl"]], n_dml[["inter"]])

# ---- The test (copied from dml_corrected_test.R) ------------------------------------------------
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
not_corrected <- function(dispersion) !is.na(dispersion) & dispersion <= 1
check_against <- function(mine, cached, what) {
  worst <- max(abs(mine - cached) / pmax(abs(cached), .Machine$double.xmin))
  if (!is.finite(worst) || worst > equal_tolerance) {
    stop(what, ": the p-values computed here differ from the pipeline's (largest relative difference ", signif(worst, 3), ").", call. = FALSE)
  }
  log_msg("  %s: p-values agree with the pipeline's at %d CpGs (largest relative difference %.1e)", what, length(mine), worst)
}
shares <- function(p, levels) vapply(levels, function(l) mean(p < l), numeric(1))
level_names <- function(levels) gsub(".", "_", format(levels, scientific = FALSE, drop0trailing = TRUE), fixed = TRUE)
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
high_split_controls <- function(z) {   # the two raised control libraries (I-49): together or apart
  h <- which(high[control_rows]); if (length(h) != 2) return(NA_character_)
  if (z[h[1]] == z[h[2]]) "together" else "apart"
}
high_split_pooled <- function(z) { ones <- sum(z[high]); if (abs(2 * ones - sum(high)) == 1) "3:2" else "4:1 or 5:0" }

# Run fun(j) for every labelling j, in parallel; any failure stops the run.
over_labellings <- function(n, fun) {
  res <- parallel::mclapply(seq_len(n), fun, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("Labelling ", which(bad)[1], " failed: ", as.character(res[[which(bad)[1]]]), call. = FALSE)
  rbindlist(res, fill = TRUE)
}
limit <- function(labels) if (is.na(trial)) labels else labels[, seq_len(min(ncol(labels), trial + 1L)), drop = FALSE]

# The two designs. Effect of a labelling at every CpG, in percentage points.
ctrl_groupings <- cohort_groupings(control_rows, n_control, control_seed)
designs <- list(
  ctrl = list(name = "control-only ploidy (cov5-mnchisq-ctrl)", rows_lib = control_rows,
              X_full = function(z) cbind(1, z), X_red = function(z) NULL,
              effect = function(z) {
                c_ <- cs[, control_rows, drop = FALSE]; n_ <- ns[, control_rows, drop = FALSE]
                c1 <- drop(c_ %*% z); n1 <- drop(n_ %*% z)
                100 * (c1 / n1 - (rowSums(c_) - c1) / (rowSums(n_) - n1))
              },
              split = high_split_controls),
  inter = list(name = "ploidy by heat interaction (cov5-mnchisq-inter)", rows_lib = seq_len(nrow(lib)),
               X_full = function(z) cbind(1, z * heat, z, heat), X_red = function(z) cbind(1, z, heat),
               effect = function(z) {
                 pct <- function(g) { j <- which(g); 100 * rowSums(cs[, j, drop = FALSE]) / rowSums(ns[, j, drop = FALSE]) }
                 (pct(z == 1 & heat == 1) - pct(z == 1 & heat == 0)) - (pct(z == 0 & heat == 1) - pct(z == 0 & heat == 0))
               },
               split = high_split_pooled))

# ---- A. Level under relabelling ---------------------------------------------------------------
samples <- list()
for (k in seq_along(sample_seeds)) {
  set.seed(sample_seeds[k])
  pool <- setdiff(eligible, unlist(samples))
  samples[[k]] <- sort(sample(pool, n_sample))
}
log_msg("part A: level under relabelling (%d samples of %d eligible CpGs)", length(samples), n_sample)
level_rows <- rbindlist(lapply(names(designs), function(dn) {
  d <- designs[[dn]]
  rbindlist(lapply(seq_along(samples), function(k) {
    rows <- samples[[k]]
    labels <- limit(if (dn == "ctrl") ctrl_groupings$labels else pooled_labellings(n_relabel, labelling_seeds[k]))
    c_ <- cs[rows, d$rows_lib, drop = FALSE]; n_ <- ns[rows, d$rows_lib, drop = FALSE]
    started <- Sys.time()
    res <- over_labellings(ncol(labels), function(j) {
      z <- labels[, j]
      tests <- run_tests(c_, n_, seq_along(rows), d$X_full(z), d$X_red(z))
      s <- c(shares(tests[, "p_uncorrected"], p_levels), shares(tests[, "p_corrected"], p_levels))
      names(s) <- paste0(rep(c("uncorrected", "corrected"), each = length(p_levels)), "_p_below_", level_names(p_levels))
      out <- cbind(data.table(design = dn, sample = k, labelling = j - 1L, observed = j == 1L, high_level_split = d$split(z)), as.data.table(as.list(s)))
      if (j == 1) out[, max_rel_diff := max(abs(tests[, "p_corrected"] - fits[[dn]]$pvalue[rows]) / pmax(abs(fits[[dn]]$pvalue[rows]), .Machine$double.xmin))]
      out
    })
    worst <- res[observed == TRUE, max_rel_diff]
    if (!is.finite(worst) || worst > equal_tolerance) stop("part A, ", dn, ", sample ", k, ": real-label p-values differ from the pipeline's (", signif(worst, 3), ").", call. = FALSE)
    log_msg("  %s, sample %d: %d labellings in %.1f min; real labels agree with the pipeline (largest relative difference %.1e)",
            dn, k, nrow(res), as.numeric(difftime(Sys.time(), started, units = "mins")), worst)
    res[, max_rel_diff := NULL]
  }), fill = TRUE)
}), fill = TRUE)

level_summary <- rbindlist(lapply(split(level_rows, by = c("design", "sample")), function(d) {
  o <- d[observed == TRUE]; r <- d[observed == FALSE]
  rbindlist(lapply(level_names(p_levels), function(l) {
    col <- paste0("corrected_p_below_", l); colu <- paste0("uncorrected_p_below_", l)
    data.table(design = o$design, sample = o$sample, p_level = as.numeric(sub("_", ".", sub("^0_", "0.", l))), n_relabellings = nrow(r),
               observed_share_corrected = o[[col]], relabelled_median_share_corrected = stats::median(r[[col]]),
               relabelled_mean_share_corrected = mean(r[[col]]),
               relabellings_at_or_above_observed_corrected = sum(r[[col]] >= o[[col]] - 1e-12),
               observed_share_uncorrected = o[[colu]], relabelled_median_share_uncorrected = stats::median(r[[colu]]))
  }))
}))

# ---- B. Length of the lists under relabelling ----------------------------------------------------
log_msg("part B: lists under relabelling (all %d tested CpGs)", n_tested)
list_rows <- rbindlist(lapply(names(designs), function(dn) {
  d <- designs[[dn]]
  fit <- fits[[dn]]
  p_star <- max(fit$pvalue[fit$qvalue < config$dml_qvalue])
  labels <- limit(if (dn == "ctrl") ctrl_groupings$labels else pooled_labellings(999L, list_seed)[, seq_len(n_list_pooled + 1L)])
  c_ <- cs[, d$rows_lib, drop = FALSE]; n_ <- ns[, d$rows_lib, drop = FALSE]
  started <- Sys.time()
  res <- over_labellings(ncol(labels), function(j) {
    z <- labels[, j]
    effect <- d$effect(z)
    above <- which(abs(effect) > config$dml_difference)
    tests <- run_tests(c_, n_, above, d$X_full(z), d$X_red(z))
    called <- tests[, "p_corrected"] <= p_star
    ps <- vapply(p_levels_lists, function(l) sum(tests[, "p_corrected"] < l), numeric(1))
    names(ps) <- paste0("n_above_cutoff_corrected_p_below_", level_names(p_levels_lists))
    out <- cbind(data.table(design = dn, labelling = j - 1L, observed = j == 1L, high_level_split = d$split(z),
                            p_threshold_of_observed_list = p_star, n_above_cutoff = length(above), n_called = sum(called),
                            n_called_hyper = sum(called & effect[above] > 0),
                            n_called_correction_not_applied = sum(called & not_corrected(tests[, "dispersion"])),
                            n_called_on_their_own = sum(tests[, "p_corrected"] < p_single),
                            median_dispersion_above_cutoff = stats::median(tests[, "dispersion"], na.rm = TRUE)),
                 as.data.table(as.list(ps)))
    if (j == 1) {
      out[, max_rel_diff := max(abs(tests[, "p_corrected"] - fit$pvalue[above]) / pmax(abs(fit$pvalue[above]), .Machine$double.xmin))]
      out[, max_effect_diff := max(abs(effect - fit$meth.diff))]
    }
    out
  })
  o <- res[observed == TRUE]
  if (!is.finite(o$max_rel_diff) || o$max_rel_diff > equal_tolerance) stop("part B, ", dn, ": real-label p-values differ from the pipeline's (", signif(o$max_rel_diff, 3), ").", call. = FALSE)
  if (o$max_effect_diff > 1e-8) stop("part B, ", dn, ": the real-label differences differ from the pipeline's meth.diff (", signif(o$max_effect_diff, 3), ").", call. = FALSE)
  if (o$n_called != n_dml[[dn]]) stop("part B, ", dn, ": ", o$n_called, " CpGs called at the list's threshold, ", n_dml[[dn]], " DML in the pipeline.", call. = FALSE)
  log_msg("  %s: %d labellings in %.1f min; real labels: %d called (= the pipeline's DML), p and differences agree (%.1e); relabelled median %g",
          dn, nrow(res), as.numeric(difftime(Sys.time(), started, units = "mins")), o$n_called, o$max_rel_diff, stats::median(res$n_called[-1]))
  res[, c("max_rel_diff", "max_effect_diff") := NULL]
}), fill = TRUE)

q_ <- function(x, p) unname(stats::quantile(x, p))
list_summary <- rbindlist(lapply(split(list_rows, by = "design"), function(d) {
  o <- d[observed == TRUE]; r <- d[observed == FALSE]
  base <- data.table(design = o$design, analysis = designs[[o$design]]$name, n_relabellings = nrow(r),
                     p_threshold_of_observed_list = o$p_threshold_of_observed_list,
                     observed_n_above_cutoff = o$n_above_cutoff, relabelled_n_above_cutoff_median = stats::median(r$n_above_cutoff),
                     observed_n_called = o$n_called, observed_n_called_hyper = o$n_called_hyper,
                     relabelled_n_called_median = stats::median(r$n_called), relabelled_n_called_q05 = q_(r$n_called, 0.05),
                     relabelled_n_called_q95 = q_(r$n_called, 0.95), relabelled_n_called_max = max(r$n_called),
                     relabellings_n_called_at_or_above_observed = sum(r$n_called >= o$n_called),
                     observed_n_called_on_their_own = o$n_called_on_their_own,
                     relabelled_n_called_on_their_own_median = stats::median(r$n_called_on_their_own),
                     relabelled_n_called_on_their_own_max = max(r$n_called_on_their_own),
                     relabellings_on_their_own_at_or_above_observed = sum(r$n_called_on_their_own >= o$n_called_on_their_own),
                     observed_median_dispersion_above_cutoff = o$median_dispersion_above_cutoff,
                     relabelled_median_dispersion_above_cutoff = stats::median(r$median_dispersion_above_cutoff))
  for (cc in grep("^n_above_cutoff_corrected_p_below_", names(d), value = TRUE)) {
    obs_share <- o[[cc]] / o$n_above_cutoff; rel_share <- r[[cc]] / r$n_above_cutoff
    base[, (sub("^n_above_cutoff_", "observed_share_", cc)) := obs_share]
    base[, (sub("^n_above_cutoff_", "relabelled_mean_share_", cc)) := mean(rel_share)]
    base[, (sub("^n_above_cutoff_", "relabellings_share_at_or_above_observed_", cc)) := sum(rel_share >= obs_share - 1e-12)]
  }
  # by how the libraries with a raised overall level are split (I-49)
  for (s in sort(unique(r$high_level_split))) {
    rs <- r[high_level_split == s]
    base[, (paste0("relabelled_n_called_median_high_", gsub("[^A-Za-z0-9]+", "_", s))) := stats::median(rs$n_called)]
    base[, (paste0("relabellings_high_", gsub("[^A-Za-z0-9]+", "_", s))) := nrow(rs)]
    base[, (paste0("relabellings_at_or_above_observed_high_", gsub("[^A-Za-z0-9]+", "_", s))) := sum(rs$n_called >= o$n_called)]
  }
  base[, observed_high_level_split := o$high_level_split]
  base
}), fill = TRUE)

print(list_summary[, .(design, n_relabellings, observed_n_called, relabelled_n_called_median, relabelled_n_called_q95,
                       relabelled_n_called_max, relabellings_n_called_at_or_above_observed)])
print(level_summary[, .(design, sample, p_level, observed_share_corrected, relabelled_median_share_corrected, relabellings_at_or_above_observed_corrected)])

if (is.na(trial)) {
  write_table(level_rows, out_file("relabel_level.csv"))
  write_table(level_summary, out_file("relabel_level_summary.csv"))
  write_table(list_rows, out_file("relabel_lists.csv"))
  write_table(list_summary, out_file("relabel_lists_summary.csv"))
  writeLines(capture.output(sessionInfo()), out_file("session_info.txt"))
} else {
  log_msg("trial run (%d relabellings per design): nothing written", trial)
}
log_msg("dml_corrected_test_5x finished in %.1f min; outputs in %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir)

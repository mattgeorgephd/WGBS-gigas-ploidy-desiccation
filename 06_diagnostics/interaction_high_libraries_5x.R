# interaction_high_libraries_5x.R
#
# Purpose   Do the five libraries with a raised overall methylation level (issue
#           I-49: three heat-cohort diploids, zr3534_3 to 5, and two control
#           triploids, zr3644_14 and 15) carry the corrected 5x interaction list
#           of the optional analysis cov5-mnchisq-inter (1,718 DML; issues I-02
#           and I-49)? Author's request of 2026-11-24. Three parts:
#             1. At the 1,718 DML. The interaction test is repeated without the
#                five, without each single library, and without each of the 150
#                sets of five libraries with the same make-up (three of the five
#                heat-cohort diploids, two of the six control triploids; the five
#                are one of them). A DML is kept when its corrected p-value is at
#                or below the threshold of the 22-library list and its interaction
#                difference exceeds 20 points with the same sign.
#             2. Without the five, genome-wide. The 17 other libraries, real
#                labels, against regroupings of the same libraries within cohort
#                (2 diploids and 5 triploids in the heat cohort, 6 and 4 among the
#                controls: 21 x 210 = 4,410 groupings, the real one included): the
#                number of CpGs whose interaction difference exceeds 20 points
#                for all 4,410 (the measure of I-02, update of 2026-11-19, at
#                10x), and the corrected calls for the real grouping and 200 of
#                the others drawn with a seed.
#             3. Without each of the 150 sets of part 1, genome-wide: the corrected
#                calls of the real labels, so that the loss of the five can be set
#                against the loss of any five libraries of the same make-up.
#           Calls are counted at the p-value threshold of the 22-library list
#           (8.52e-6, the largest p-value with q < 0.01) and at p below 0.01 /
#           (CpGs tested); no q-values are computed for the reduced designs. The
#           CpGs are those of the 22-library analysis (5x in every library).
#           Exploratory: the rule that finds the five ignores the labels, but
#           setting them aside was decided after seeing the data (I-49, item 5).
#           Not in diagnostics.csv (reads git-ignored tables of an optional
#           analysis). Nothing in 01_pipeline/ reads its outputs.
#
# Run       DML5X_PARTS=1,2 WGBS_CORES=4 Rscript 06_diagnostics/interaction_high_libraries_5x.R
#           DML5X_PARTS=3   WGBS_CORES=4 Rscript 06_diagnostics/interaction_high_libraries_5x.R
#           (default: all three parts; DML5X_TRIAL=<n> limits each part to n
#           evaluations and writes nothing)
#
# Reads     04_output/not_adopted/cov5-mnchisq-inter/05_differential_methylation/tables/
#             united_cpg_cov5-mnchisq-inter.csv.gz, diffmeth_interaction_cov5-mnchisq-inter.csv.gz
#           06_diagnostics/library_methylation_level/library_level.csv   column high_level
#
# Writes    06_diagnostics/interaction_high_libraries_5x/
#             part1_at_dml.csv          one row per set of libraries left out
#             part2_groupings_above20.csv   one row per grouping of the 17 (4,410)
#             part2_corrected_calls.csv one row per grouping tested (201)
#             part3_drop_sets.csv       one row per set of five left out (150)
#             summary.csv               the headline numbers of each part
#             session_info.txt
#
# Notes     The test is that of 06_diagnostics/dml_corrected_test_5x.R (methylKit
#           1.28.0 logReg in base R, overdispersion MN, chi-squared); its
#           functions are copied unchanged. With all 22 libraries the script
#           checks that it calls exactly the pipeline's 1,718 DML.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))
if (!requireNamespace("R.utils", quietly = TRUE)) stop("Missing R package: R.utils (needed to read .csv.gz files).", call. = FALSE)
t0 <- start_step("interaction_high_libraries_5x", "do the five high-level libraries carry the corrected 5x interaction list?")

group_seed <- 20261124
n_groupings_tested <- 200
cores <- max(1L, as.integer(Sys.getenv("WGBS_CORES", "1")))
parts <- trimws(strsplit(Sys.getenv("DML5X_PARTS", "1,2,3"), ",")[[1]])
trial <- suppressWarnings(as.integer(Sys.getenv("DML5X_TRIAL", "")))
if (!identical(config$dml_qvalue, 0.01) || !identical(as.numeric(config$dml_difference), 20)) {
  stop("Run with the paper's thresholds: unset WGBS_ANALYSIS and the optional settings.", call. = FALSE)
}
out_dir <- here("06_diagnostics", "interaction_high_libraries_5x")
writing <- is.na(trial)
if (writing) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, name)
tab <- function(file) here("04_output", "not_adopted", "cov5-mnchisq-inter", "05_differential_methylation", "tables", file)
need <- function(f) { if (!file.exists(f)) stop("Missing ", f, ". Run the optional analysis cov5-mnchisq-inter (or its step 05) first.", call. = FALSE); f }
limit_n <- function(n) if (is.na(trial)) n else min(n, trial)

# ---- Inputs ---------------------------------------------------------------------------------
lib <- read_libraries()
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
level <- fread(need(here("06_diagnostics", "library_methylation_level", "library_level.csv")))
if (!identical(level$library_id, lib$library_id)) stop("library_level.csv does not list the libraries of the library sheet in order.", call. = FALSE)
high <- which(level$high_level)
if (length(high) != 5 || sum(heat[high] == 1 & trip[high] == 0) != 3 || sum(heat[high] == 0 & trip[high] == 1) != 2) {
  stop("Expected the five high libraries of I-49: three heat-cohort diploids and two control triploids.", call. = FALSE)
}
log_msg("the five: %s", paste(lib$library_id[high], collapse = ", "))

u <- fread(need(tab("united_cpg_cov5-mnchisq-inter.csv.gz")))
cs <- as.matrix(u[, paste0("numCs", seq_len(nrow(lib))), with = FALSE]) + 0
ns <- cs + as.matrix(u[, paste0("numTs", seq_len(nrow(lib))), with = FALSE])
n_tested <- nrow(u); rm(u); invisible(gc())
fit <- fread(need(tab("diffmeth_interaction_cov5-mnchisq-inter.csv.gz")), select = c("pvalue", "qvalue", "meth.diff", "united_row"))
if (nrow(fit) != n_tested || !identical(fit$united_row, seq_len(n_tested))) stop("The interaction table does not follow the united table.", call. = FALSE)
p_star <- max(fit$pvalue[fit$qvalue < config$dml_qvalue])
p_single <- config$dml_qvalue / n_tested
dml <- which(fit$qvalue < config$dml_qvalue & abs(fit$meth.diff) > config$dml_difference)
dml_sign <- sign(fit$meth.diff[dml])
log_msg("%d CpGs tested; %d interaction DML; list threshold p <= %.3g", n_tested, length(dml), p_star)

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
over <- function(n, fun) {
  res <- parallel::mclapply(seq_len(n), fun, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("Evaluation ", which(bad)[1], " failed: ", as.character(res[[which(bad)[1]]]), call. = FALSE)
  rbindlist(res, fill = TRUE)
}

# Interaction difference and design for the libraries `keep` with ploidy labels `z` (one per kept library).
interaction_effect <- function(keep, z, rows = NULL) {
  h <- heat[keep]
  c_ <- if (is.null(rows)) cs[, keep, drop = FALSE] else cs[rows, keep, drop = FALSE]
  n_ <- if (is.null(rows)) ns[, keep, drop = FALSE] else ns[rows, keep, drop = FALSE]
  pct <- function(g) { j <- which(g); 100 * rowSums(c_[, j, drop = FALSE]) / rowSums(n_[, j, drop = FALSE]) }
  (pct(z == 1 & h == 1) - pct(z == 1 & h == 0)) - (pct(z == 0 & h == 1) - pct(z == 0 & h == 0))
}
test_at <- function(keep, z, rows) {
  h <- heat[keep]
  run_tests(cs[, keep, drop = FALSE], ns[, keep, drop = FALSE], rows, cbind(1, z * h, z, h), cbind(1, z, h))
}
# Genome-wide corrected calls for libraries `keep`, labels `z`.
genome_calls <- function(keep, z) {
  eff <- interaction_effect(keep, z)
  above <- which(abs(eff) > config$dml_difference)
  tests <- test_at(keep, z, above)
  called <- tests[, "p_corrected"] <= p_star
  list(n_above20 = length(above), n_called = sum(called), n_called_negative = sum(called & eff[above] < 0),
       n_called_on_their_own = sum(tests[, "p_corrected"] < p_single),
       n_called_in_dml22 = sum(called & above %in% dml), eff = eff, above = above, tests = tests)
}

# The 150 sets of five with the make-up of the high five (the high five among them).
heat_dip <- which(heat == 1 & trip == 0); ctrl_trip <- which(heat == 0 & trip == 1)
drop_sets <- unlist(lapply(combn(heat_dip, 3, simplify = FALSE), function(a)
  lapply(combn(ctrl_trip, 2, simplify = FALSE), function(b) sort(c(a, b)))), recursive = FALSE)
is_high_set <- vapply(drop_sets, function(s) identical(s, sort(high)), logical(1))
if (length(drop_sets) != 150 || sum(is_high_set) != 1) stop("Expected 150 sets of five, one of them the high five.", call. = FALSE)
set_label <- function(s) paste(lib$library_id[s], collapse = ";")
high_in <- function(s) length(intersect(s, high))

summary_rows <- list()

# ---- 1. At the 1,718 DML ------------------------------------------------------------------------
if ("1" %in% parts) {
  log_msg("part 1: the %d DML without the five, without each library, and without each of the 150 sets", length(dml))
  sets <- c(list(integer(0), sort(high)), as.list(seq_len(nrow(lib))), drop_sets[!is_high_set])
  kind <- c("none", "the five high libraries", rep("one library", nrow(lib)), rep("five of the same make-up", sum(!is_high_set)))
  n_eval <- limit_n(length(sets))
  started <- Sys.time()
  p1 <- over(n_eval, function(k) {
    drop <- sets[[k]]; keep <- setdiff(seq_len(nrow(lib)), drop)
    eff <- interaction_effect(keep, trip[keep], rows = dml)
    tests <- test_at(keep, trip[keep], dml)
    kept <- tests[, "p_corrected"] <= p_star & abs(eff) > config$dml_difference & sign(eff) == dml_sign
    data.table(kind = kind[k], left_out = if (length(drop)) set_label(drop) else "", n_left_out = length(drop),
               high_left_out = high_in(drop), n_dml = length(dml), n_kept = sum(kept),
               n_same_sign = sum(sign(eff) == dml_sign), n_p_at_threshold = sum(tests[, "p_corrected"] <= p_star),
               median_abs_difference = stats::median(abs(eff)), median_p_corrected = stats::median(tests[, "p_corrected"]))
  })
  log_msg("  %d evaluations in %.1f min", nrow(p1), as.numeric(difftime(Sys.time(), started, units = "mins")))
  if (p1[kind == "none", n_kept] != length(dml)) stop("Part 1: with all 22 libraries the test keeps ", p1[kind == "none", n_kept], " of ", length(dml), " DML.", call. = FALSE)
  if (writing) write_table(p1, out_file("part1_at_dml.csv"))
  others <- p1[kind == "five of the same make-up"]
  if (nrow(others)) summary_rows$p1 <- data.table(part = 1, measure = "DML kept (of 1,718) without the set left out",
    observed = p1[kind == "the five high libraries", n_kept], reference = "the other 149 sets of five of the same make-up",
    reference_median = stats::median(others$n_kept), reference_min = min(others$n_kept), reference_max = max(others$n_kept),
    reference_at_or_below_observed = sum(others$n_kept <= p1[kind == "the five high libraries", n_kept]), n_reference = nrow(others))
  print(p1[kind %in% c("none", "the five high libraries", "one library"), .(kind, left_out, n_kept, n_same_sign, median_abs_difference)])
}

# ---- 2. Without the five, genome-wide --------------------------------------------------------------
if ("2" %in% parts) {
  keep <- setdiff(seq_len(nrow(lib)), high)
  hk <- which(heat[keep] == 1); ck <- which(heat[keep] == 0)
  g_heat <- combn(hk, sum(trip[keep][hk] == 0), simplify = FALSE)     # which heat libraries are labelled diploid
  g_ctrl <- combn(ck, sum(trip[keep][ck] == 1), simplify = FALSE)     # which control libraries are labelled triploid
  grid <- CJ(i = seq_along(g_heat), j = seq_along(g_ctrl))
  lab <- function(i, j) { z <- integer(length(keep)); z[hk] <- 1L; z[g_heat[[i]]] <- 0L; z[g_ctrl[[j]]] <- 1L; z }
  real <- trip[keep]
  grid[, is_real := vapply(seq_len(.N), function(r) identical(lab(i[r], j[r]), real), logical(1))]
  if (nrow(grid) != 4410 || sum(grid$is_real) != 1) stop("Expected 4,410 groupings of the 17 libraries, one of them the real one.", call. = FALSE)
  n_eval <- limit_n(nrow(grid))
  log_msg("part 2: without the five; CpGs above 20 points for %d groupings", n_eval)
  started <- Sys.time()
  ord <- c(which(grid$is_real), setdiff(seq_len(nrow(grid)), which(grid$is_real)))
  p2a <- over(n_eval, function(k) {
    r <- ord[k]; z <- lab(grid$i[r], grid$j[r]); eff <- interaction_effect(keep, z)
    data.table(grouping = r, real = grid$is_real[r], heat_diploids = paste(lib$library_id[keep][g_heat[[grid$i[r]]]], collapse = ";"),
               n_above20 = sum(abs(eff) > config$dml_difference), n_above20_negative = sum(eff < -config$dml_difference))
  })
  log_msg("  %d groupings in %.1f min", nrow(p2a), as.numeric(difftime(Sys.time(), started, units = "mins")))
  if (writing) write_table(p2a, out_file("part2_groupings_above20.csv"))
  set.seed(group_seed)
  tested <- c(which(grid$is_real), sort(sample(which(!grid$is_real), n_groupings_tested)))
  n_eval <- limit_n(length(tested))
  log_msg("part 2: corrected calls for the real grouping and %d others", n_eval - 1L)
  started <- Sys.time()
  p2b <- over(n_eval, function(k) {
    r <- tested[k]; z <- lab(grid$i[r], grid$j[r]); g <- genome_calls(keep, z)
    data.table(grouping = r, real = grid$is_real[r], heat_diploids = paste(lib$library_id[keep][g_heat[[grid$i[r]]]], collapse = ";"),
               n_above20 = g$n_above20, n_called = g$n_called, n_called_negative = g$n_called_negative,
               n_called_on_their_own = g$n_called_on_their_own, n_called_in_dml22 = g$n_called_in_dml22)
  })
  log_msg("  %d groupings in %.1f min", nrow(p2b), as.numeric(difftime(Sys.time(), started, units = "mins")))
  if (writing) write_table(p2b, out_file("part2_corrected_calls.csv"))
  o <- p2a[real == TRUE]; r_ <- p2a[real == FALSE]
  summary_rows$p2a <- data.table(part = 2, measure = "CpGs with an interaction difference above 20 points, 17 libraries",
    observed = o$n_above20, reference = "all other groupings of the 17 within cohort", reference_median = stats::median(r_$n_above20),
    reference_min = min(r_$n_above20), reference_max = max(r_$n_above20), reference_at_or_above_observed = sum(r_$n_above20 >= o$n_above20), n_reference = nrow(r_))
  o <- p2b[real == TRUE]; r_ <- p2b[real == FALSE]
  if (nrow(r_)) summary_rows$p2b <- data.table(part = 2, measure = "corrected calls at the 22-library threshold, 17 libraries",
    observed = o$n_called, reference = "groupings of the 17 within cohort drawn with the seed", reference_median = stats::median(r_$n_called),
    reference_min = min(r_$n_called), reference_max = max(r_$n_called), reference_at_or_above_observed = sum(r_$n_called >= o$n_called), n_reference = nrow(r_))
  print(p2b[real == TRUE])
}

# ---- 3. Without each set of five, genome-wide ------------------------------------------------------
if ("3" %in% parts) {
  sets <- c(list(sort(high)), drop_sets[!is_high_set])
  n_eval <- limit_n(length(sets))
  log_msg("part 3: corrected calls without each of %d sets of five", n_eval)
  started <- Sys.time()
  p3 <- over(n_eval, function(k) {
    drop <- sets[[k]]; keep <- setdiff(seq_len(nrow(lib)), drop); g <- genome_calls(keep, trip[keep])
    data.table(left_out = set_label(drop), the_five = k == 1L, high_left_out = high_in(drop), n_above20 = g$n_above20,
               n_called = g$n_called, n_called_negative = g$n_called_negative, n_called_on_their_own = g$n_called_on_their_own,
               n_called_in_dml22 = g$n_called_in_dml22)
  })
  log_msg("  %d sets in %.1f min", nrow(p3), as.numeric(difftime(Sys.time(), started, units = "mins")))
  if (writing) write_table(p3, out_file("part3_drop_sets.csv"))
  o <- p3[the_five == TRUE]; r_ <- p3[the_five == FALSE]
  if (nrow(r_)) summary_rows$p3 <- data.table(part = 3, measure = "corrected calls of the real labels without the set left out",
    observed = o$n_called, reference = "the other 149 sets of five of the same make-up", reference_median = stats::median(r_$n_called),
    reference_min = min(r_$n_called), reference_max = max(r_$n_called), reference_at_or_below_observed = sum(r_$n_called <= o$n_called), n_reference = nrow(r_))
  print(p3[order(n_called)][1:min(.N, 10)])
}

if (length(summary_rows)) {
  s <- rbindlist(summary_rows, fill = TRUE); print(s)
  if (writing) {
    f <- out_file("summary.csv")
    if (file.exists(f)) { old <- fread(f); s <- rbind(old[!part %in% s$part], s, fill = TRUE); setorder(s, part) }
    write_table(s, f)
  }
}
if (writing) writeLines(capture.output(sessionInfo()), out_file("session_info.txt"))
log_msg("interaction_high_libraries_5x finished in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

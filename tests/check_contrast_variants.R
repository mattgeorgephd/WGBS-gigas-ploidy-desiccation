# tests/check_contrast_variants.R
#
# Unit tests of the options ctrl and inter of step 05 (cached_variant_fit(),
# variant_design(), interaction_groups(), write_mode()) and of the padding of
# rebuild_comparison(), on a small simulated methylKit object: 400 CpGs on two
# linkage groups and two scaffolds, the 22 libraries of the library sheet. Fits
# nothing of the real data, reads no coverage file and writes only into a
# temporary folder (the cache folder of step 05 and not_adopted_path() are
# pointed there). Seconds, well under 1 GB. From the repository root:
#   Rscript tests/check_contrast_variants.R
# Checks:
#  1. ctrl, by sequence, equals the call of analysis A4_ploidy_in_control of
#     06_diagnostics/dml_sensitivity.R on the whole table (reorganize() to the 12
#     controls, calculateDiffMeth(overdispersion "none", test "Chisq")): p-values,
#     SLIM q-values over all CpGs and differences identical; and with mnchisq,
#     the call of A4_ploidy_in_control_mn_chisq;
#  2. ctrl: meth.diff is triploid minus diploid among the controls, from pooled
#     counts; the heat contrast is the base fit's, unchanged;
#  3. inter, by sequence, equals the call of A2_interaction on the whole table
#     (p, q identical); its p-values equal the likelihood ratio test of ploidy:heat
#     of a binomial glm written as ploidy * heat; its difference equals
#     (triploid heat - triploid control) - (diploid heat - diploid control) from
#     pooled counts; the ploidy and heat contrasts are the base fit's; the same
#     with mnchisq against A2_interaction_mn_chisq;
#  4. the cache: a second call reuses methylkit_<fit id>.rds and fits nothing; a
#     base fit with another fingerprint refits the new contrast; parts saved by a
#     stopped run are reused;
#  5. write_mode(): inter writes dml_interaction_<id>.csv and .bed by the DML rule
#     of the analysis, and the counts interaction_dml, _hyper, _hypo and
#     interaction_in_shared after the 12 columns of the other analyses; an
#     analysis without inter keeps exactly those 12 columns; ctexcl removes the
#     listed CpGs from the interaction list too;
#  6. rebuild_comparison() pads the rows of tables without the interaction
#     columns, keeps the paper's row bytes, and writes what it wrote before when
#     no table has extra columns.
# Exit status 1 when a check fails.

suppressPackageStartupMessages(library(here))
Sys.unsetenv(c("WGBS_ANALYSIS", "WGBS_GENE_LEVEL_GO", "WGBS_DIAGNOSTICS"))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))
for (e in parse(here("01_pipeline", "05_differential_methylation.R"), keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      ((is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) || as.character(e[[2]]) %in% c("step", "lib"))) {
    eval(e, globalenv())
  }
}
config$cores <- 1L
log_msg <- function(...) invisible(NULL)
fails <- character(0); n_checks <- 0L
check <- function(ok, what) {
  n_checks <<- n_checks + 1L
  cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n")
  if (!isTRUE(ok)) fails <<- c(fails, what)
}
tmp <- tempfile("contrast_variants"); dir.create(tmp)
step_cache_dir <- function(step) { d <- file.path(tmp, "cache"); dir.create(d, showWarnings = FALSE); d }
not_adopted_path <- function(...) file.path(tmp, "not_adopted", ...)
paper_path <- function(step, sub, file) file.path(tmp, "paper", step, sub, file)

# ---- Simulated united table -----------------------------------------------------------
set.seed(20261122)
n <- 400
chr <- rep(c(config$linkage_groups[1:2], "NW_000001.1", "NW_000002.1"), times = c(150, 150, 60, 40))
start <- unlist(lapply(table(factor(chr, levels = unique(chr))), function(k) sort(sample.int(1e6, k)) * 2L))
base_p <- stats::plogis(rnorm(n, 0.5, 1.5))
eff <- cbind(ploidy = rnorm(n, 0, 0.6) * (runif(n) < 0.3), heat = rnorm(n, 0, 0.6) * (runif(n) < 0.3),
             inter = rnorm(n, 0, 1) * (runif(n) < 0.3))
objects <- lapply(seq_len(nrow(lib)), function(i) {
  trip <- lib$ploidy[i] == "triploid"; heat <- lib$treatment[i] == "heat"
  p <- stats::plogis(stats::qlogis(base_p) + eff[, "ploidy"] * trip + eff[, "heat"] * heat + eff[, "inter"] * trip * heat)
  coverage <- sample(10:60, n, replace = TRUE)
  numCs <- stats::rbinom(n, coverage, p)
  make_methylraw(chr, start, list(coverage = coverage, numCs = numCs, numTs = coverage - numCs), as.character(lib$legacy_sample_index[i]))
})
treat_ploidy <- as.integer(lib$ploidy == "triploid"); treat_heat <- as.integer(lib$treatment == "heat")
united <- suppressMessages(unite(new("methylRawList", objects, treatment = treat_ploidy), destrand = FALSE))
stopifnot(nrow(united) == n)
base_fit <- function(overdispersion) {
  od <- if (overdispersion == "MN") list(overdispersion = "MN", test = "Chisq") else list()
  dm <- function(obj, cov) suppressMessages(suppressWarnings(do.call(calculateDiffMeth, c(list(obj, covariates = cov, mc.cores = 1), od))))
  list(mode = "base", filtering = data.frame(library_id = lib$library_id), united = united,
       diffmeth = list(ploidy = dm(united, data.frame(heat = lib$treatment)),
                       heat = dm(reorganize(united, sample.ids = united@sample.ids, treatment = treat_heat), data.frame(ploidy = lib$ploidy))),
       fingerprint = list(stand_in = paste("base fit", overdispersion)))
}
bases <- list(none = base_fit("none"), MN = base_fit("MN"))
calls <- character(0)
fit_base_of <- function(a) { calls <<- c(calls, paste("base", a$fit_id)); bases[[a$overdispersion]] }
counted_parts <- function(a, base, fingerprint) { calls <<- c(calls, "parts"); suppressWarnings(fit_variant_by_sequence(a, base, fingerprint)) }
variant <- function(id, ...) { calls <<- character(0); cached_variant_fit(analysis_spec(id), fit_base = fit_base_of, fit_parts = counted_parts, ...) }
dd <- function(x) as.data.table(getData(x))
# The calls of analyses A4 and A2 of 06_diagnostics/dml_sensitivity.R (fit_analysis()) on the whole table.
diagnostic_call <- function(keep, tested, covariates, overdispersion) {
  object <- reorganize(united, sample.ids = united@sample.ids[keep], treatment = tested)
  dd(suppressMessages(suppressWarnings(calculateDiffMeth(object, covariates = covariates, overdispersion = overdispersion,
                                                         test = "Chisq", mc.cores = 1))))
}
pooled <- function(j) { x <- dd(united); m <- rowSums(as.matrix(x[, paste0("numCs", j), with = FALSE]))
  100 * m / (m + rowSums(as.matrix(x[, paste0("numTs", j), with = FALSE]))) }
group <- function(p, t) which(lib$ploidy == p & lib$treatment == t)

# ---- 1, 2: ctrl -----------------------------------------------------------------------------
controls <- lib$treatment == "control"
for (od in c("none", "MN")) {
  id <- if (od == "MN") "mnchisq-ctrl" else "ctrl"
  f <- variant(id)
  a4 <- diagnostic_call(controls, treat_ploidy[controls], NULL, od)
  got <- dd(f$diffmeth$ploidy)
  check(identical(got$pvalue, a4$pvalue) && identical(got$qvalue, a4$qvalue) && identical(got$meth.diff, a4$meth.diff) &&
          identical(got$chr, a4$chr) && identical(got$start, a4$start),
        sprintf("%s by sequence: p, q and meth.diff identical to the call of A4_ploidy_in_control%s on the whole table", id,
                if (od == "MN") "_mn_chisq" else ""))
  check(identical(f$diffmeth$heat, bases[[od]]$diffmeth$heat) && identical(f$united, united) && is.null(f$diffmeth$interaction) &&
          identical(f$diffmeth$ploidy@sample.ids, united@sample.ids[controls]) && identical(f$diffmeth$ploidy@treatment, treat_ploidy[controls]),
        sprintf("%s: heat contrast and united table are the base fit's; the ploidy contrast is on the 12 controls", id))
  check(identical(calls, c(paste("base", analysis_spec(id)$base_fit_id), "parts")), sprintf("%s: the base fit is %s; the new contrast is fitted", id,
                                                                                              analysis_spec(id)$base_fit_id))
  if (od == "none") {
    trip_minus_dip <- pooled(c(group("triploid", "control"))) - pooled(group("diploid", "control"))
    check(max(abs(got$meth.diff - trip_minus_dip)) < 1e-9 && sum(controls) == 12,
          "ctrl: meth.diff is triploid minus diploid among the 12 controls (pooled counts)")
  }
}

# ---- 3: inter ---------------------------------------------------------------------------------
tested_inter <- as.integer(lib$ploidy == "triploid" & lib$treatment == "heat")
for (od in c("none", "MN")) {
  id <- if (od == "MN") "mnchisq-inter" else "inter"
  f <- variant(id)
  a2 <- diagnostic_call(rep(TRUE, 22), tested_inter, data.frame(ploidy = lib$ploidy, heat = lib$treatment), od)
  got <- dd(f$diffmeth$interaction)
  check(identical(got$pvalue, a2$pvalue) && identical(got$qvalue, a2$qvalue),
        sprintf("%s by sequence: p and q identical to the call of A2_interaction%s on the whole table", id, if (od == "MN") "_mn_chisq" else ""))
  check(identical(f$diffmeth$ploidy, bases[[od]]$diffmeth$ploidy) && identical(f$diffmeth$heat, bases[[od]]$diffmeth$heat),
        sprintf("%s: the ploidy and heat contrasts are the base fit's", id))
  inter_diff <- (pooled(group("triploid", "heat")) - pooled(group("triploid", "control"))) -
    (pooled(group("diploid", "heat")) - pooled(group("diploid", "control")))
  check(max(abs(got$meth.diff - inter_diff)) < 1e-9 && identical(f$interaction_groups$interaction_diff, got$meth.diff) &&
          max(abs(a2$meth.diff - (pooled(group("triploid", "heat")) - pooled(which(!(lib$ploidy == "triploid" & lib$treatment == "heat")))))) < 1e-9,
        sprintf("%s: the difference is (triploid heat - triploid control) - (diploid heat - diploid control), not methylKit's meth.diff (triploid heat minus the rest)", id))
  if (od == "none") {
    x <- dd(united)
    frame <- data.frame(ploidy = lib$ploidy, heat = lib$treatment)
    additive <- stats::model.matrix(~ ploidy + heat, frame); factorial <- stats::model.matrix(~ ploidy * heat, frame)
    p_glm <- vapply(seq_len(n), function(i) {
      y <- unlist(x[i, paste0("numCs", 1:22), with = FALSE]); m <- y + unlist(x[i, paste0("numTs", 1:22), with = FALSE])
      d0 <- suppressWarnings(stats::glm.fit(additive, y / m, weights = m, family = stats::binomial()))$deviance
      d1 <- suppressWarnings(stats::glm.fit(factorial, y / m, weights = m, family = stats::binomial()))$deviance
      stats::pchisq(d0 - d1, 1, lower.tail = FALSE)
    }, numeric(1))
    check(max(abs(got$pvalue - p_glm) / pmax(p_glm, .Machine$double.xmin)) < 1e-6,
          "inter: p-values equal the likelihood ratio test of ploidy:heat of the glm written as ploidy * heat (relative difference below 1e-6)")
  }
}
check(inherits(tryCatch(variant_design(list(id = "x", contrast_variant = "none")), error = function(e) e), "error"),
      "variant_design() stops for an analysis without ctrl or inter")

# ---- 4: cache ---------------------------------------------------------------------------------
f1 <- variant("ctrl")
check(identical(calls, "base legacy") && file.exists(file.path(tmp, "cache", "methylkit_ctrl.rds")) &&
        !dir.exists(file.path(tmp, "cache", "parts_ctrl")), "ctrl again: methylkit_ctrl.rds reused, nothing fitted, no parts left")
cached_file <- readRDS(file.path(tmp, "cache", "methylkit_ctrl.rds"))
check(is.null(cached_file$united) && length(cached_file$pvalue) == n && identical(cached_file$fingerprint$base_fit_id, "legacy"),
      "methylkit_ctrl.rds holds the new contrast only (no copy of the united table) and names its base fit")
saved_base <- bases$none
bases$none$fingerprint <- list(stand_in = "another base fit")
f2 <- variant("ctrl")
check(identical(calls, c("base legacy", "parts")), "a base fit with another fingerprint refits the new contrast")
bases$none <- saved_base
invisible(variant("ctrl"))                                         # back to the first base fit (refit)
# a run stopped after two parts: the saved parts are reused, the rest is fitted
a_ctrl <- analysis_spec("ctrl")
fp <- variant_fingerprint(a_ctrl, bases$none$fingerprint)
unlink(file.path(tmp, "cache", "methylkit_ctrl.rds"))
# (calculateDiffMeth() is masked in this environment to count the calls; it is not part of the fingerprint)
n_fitted <- 0L; stop_after <- 2L
calculateDiffMeth <- function(...) {
  n_fitted <<- n_fitted + 1L
  if (n_fitted > stop_after) stop("stopped by the test")
  methylKit::calculateDiffMeth(...)
}
invisible(tryCatch(fit_variant_by_sequence(a_ctrl, bases$none, fp), error = function(e) NULL))
saved_parts <- list.files(file.path(tmp, "cache", "parts_ctrl"), pattern = "_ploidy_.*\\.rds$")
n_fitted <- 0L; stop_after <- Inf
f3 <- variant("ctrl")
rm(calculateDiffMeth)
check(length(saved_parts) == 2 && n_fitted == length(sequence_parts(united$chr)) - 2L &&
        identical(dd(f3$diffmeth$ploidy), dd(f1$diffmeth$ploidy)),
      sprintf("a stopped fit: %d saved parts reused, %d fitted, result identical", length(saved_parts), n_fitted))

# ---- 5: write_mode ------------------------------------------------------------------------------
config$run_kind <- "optional"
out_cols <- c("mode", "cpg_analysed", "ploidy_dml", "ploidy_hyper", "ploidy_hypo", "heat_dml", "heat_hyper", "heat_hypo",
              "shared", "shared_both_hyper", "shared_both_hypo", "shared_opposite")
w <- function(id, fit) suppressWarnings(write_mode(fit, analysis_spec(id)))
h_inter <- w("inter-q05", variant("inter"))
t05 <- function(id, file) not_adopted_path(id, step, "tables", file)
di <- fread(t05("inter-q05", "dml_interaction_inter-q05.csv"))
all_i <- fread(t05("inter-q05", "diffmeth_interaction_inter-q05.csv.gz"))
sh <- fread(t05("inter-q05", "dml_shared_inter-q05.csv"))
bed <- fread(not_adopted_path("inter-q05", step, "bed", "dml_interaction_inter-q05.bed"), header = FALSE)
check(identical(names(h_inter), c(out_cols, "interaction_dml", "interaction_hyper", "interaction_hypo", "interaction_in_shared")) &&
        h_inter$interaction_dml == nrow(di) && nrow(di) == all_i[qvalue < 0.05 & abs(meth.diff) > 20, .N] && nrow(di) > 0 &&
        h_inter$interaction_hyper == sum(di$meth.diff > 0) && h_inter$interaction_hypo == sum(di$meth.diff < 0) &&
        h_inter$interaction_in_shared == sum(di$united_row %in% sh$united_row) && nrow(bed) == nrow(di) &&
        all(c("pct_diploid_control", "pct_triploid_heat") %in% names(di)),
      sprintf("inter-q05: dml_interaction (%d at q < 0.05 and more than 20 points; %d in shared), its .bed, and the 4 interaction columns after the 12",
              nrow(di), h_inter$interaction_in_shared))
h_ctrl <- w("ctrl", variant("ctrl"))
check(identical(names(h_ctrl), out_cols) && !file.exists(t05("ctrl", "dml_interaction_ctrl.csv")),
      "ctrl: dml_counts_by_mode has the 12 columns of the paper's table, no interaction list")
# ctexcl with inter: list a few interaction DML as variants
ct_listed <- di[seq_len(min(3, nrow(di))), .(chr, start)]
ct_variant_table <- function() { v <- copy(ct_listed); setkeyv(v, c("chr", "start")); v }
ct_variant_assessability <- function() { v <- dd(united)[, .(chr, start)]; v[, assessable := TRUE]; setkeyv(v, c("chr", "start")); v }
h_ct <- w("inter-ctexcl-q05", variant("inter"))
di_ct <- fread(t05("inter-ctexcl-q05", "dml_interaction_inter-ctexcl-q05.csv"))
summ <- fread(t05("inter-ctexcl-q05", "ct_exclusion_summary_inter-ctexcl-q05.csv"))
check(nrow(di_ct) == nrow(di) - nrow(ct_listed) && h_ct$interaction_dml == nrow(di_ct) &&
        nrow(merge(di_ct, ct_listed, by = c("chr", "start"))) == 0 && identical(summ$contrast, c("ploidy", "heat", "interaction")) &&
        summ[contrast == "interaction", dml_removed] == nrow(ct_listed),
      "inter-ctexcl-q05: the listed CpGs leave the interaction list too (summary row 'interaction')")

# ---- 6: rebuild_comparison() padding -----------------------------------------------------------
dir.create(dirname(paper_path(step, "tables", "x")), recursive = TRUE, showWarnings = FALSE)
write_table(h_ctrl[, mode := "legacy"], paper_path(step, "tables", "dml_counts_by_mode.csv"))
for (id in c("ctrl", "inter-q05")) {
  fwrite(data.table(id = id, filter = "legacy", min_coverage = "10", overdispersion = "none", test = "not passed",
                    contrast_variant = analysis_spec(id)$contrast_variant, ct_excluded = "FALSE", dml_difference = "20",
                    dml_qvalue = format(analysis_spec(id)$dml_qvalue), fit_id = analysis_spec(id)$fit_id, fit_fingerprint_md5 = "",
                    ct_variant_list_md5 = "", git_commit = "", written_utc = "", steps_done = "04;05"),
         not_adopted_path(id, "analysis.csv"))
}
write_table(h_ctrl[, mode := "ctrl"], t05("ctrl", "dml_counts_by_mode.csv"))
write_table(h_inter, t05("inter-q05", "dml_counts_by_mode.csv"))
unlink(not_adopted_path("inter-ctexcl-q05"), recursive = TRUE)
cmp_file <- not_adopted_path("comparisons", step, "tables", "dml_counts_by_mode.csv")
rebuild_comparison(step, "dml_counts_by_mode.csv")
cmp <- readLines(cmp_file)
paper_line <- readLines(paper_path(step, "tables", "dml_counts_by_mode.csv"))[2]
inter_line <- readLines(t05("inter-q05", "dml_counts_by_mode.csv"))[2]
ctrl_line <- readLines(t05("ctrl", "dml_counts_by_mode.csv"))[2]
check(length(cmp) == 4 && cmp[1] == readLines(t05("inter-q05", "dml_counts_by_mode.csv"))[1] &&
        cmp[2] == paste0(paper_line, ",,,,") && cmp[3] == paste0(ctrl_line, ",,,,") && cmp[4] == inter_line &&
        nrow(fread(cmp_file)) == 3 && is.na(fread(cmp_file)[mode == "ctrl", interaction_dml]),
      "rebuild_comparison(): the longest header; rows without the interaction columns padded with 4 empty fields; the inter row as written")
unlink(not_adopted_path("inter-q05"), recursive = TRUE)
rebuild_comparison(step, "dml_counts_by_mode.csv")
check(identical(readLines(cmp_file), c(readLines(paper_path(step, "tables", "dml_counts_by_mode.csv")), ctrl_line)),
      "rebuild_comparison() without an inter table: header and rows exactly as in the tables (as before the padding rule)")
writeLines(c("mode,other", "x,1"), t05("ctrl", "dml_counts_by_mode.csv"))
check(inherits(tryCatch(rebuild_comparison(step, "dml_counts_by_mode.csv"), error = function(e) e), "error"),
      "rebuild_comparison() still stops on columns that differ other than by columns added at the end")

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d of %d checks pass\n", n_checks - length(fails), n_checks))
quit(status = if (length(fails)) 1 else 0)

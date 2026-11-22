# tests/check_contrast_variants_real.R
#
# Cross-check of an analysis with option ctrl or inter, after its step 05, against
# 06_diagnostics/dml_sensitivity.R, which fits the same contrasts (A4_ploidy_in_control
# and A2_interaction, with _mn_chisq for mnchisq) in the modes legacy and corrected:
#   - the new contrast: p-values, q-values and differences of
#     diffmeth_<contrast>_<id>.csv.gz against the diagnostic's cached fit
#     (06_diagnostics/dml_sensitivity/cache/fit_<A>_<mode>.rds, git-ignored) and its
#     DML list against dml_<A>_<mode>.csv (tracked);
#   - the contrasts taken from the base fit: the DML lists of the analysis (heat for
#     ctrl; ploidy, heat and shared for inter) against those of the analysis without
#     the token (the paper's folder for legacy).
# Reads only tables; seconds, under 1 GB. From the repository root, for an analysis
# without ctexcl, d50 and q05 (the diagnostic applies q < 0.01 and 20 points):
#   Rscript tests/check_contrast_variants_real.R ctrl        (also mnchisq-ctrl, inter, mnchisq-inter,
#                                                             corrected-inter, corrected-mnchisq-inter)
# Exit status 1 when a check fails.

suppressPackageStartupMessages(library(here))
Sys.unsetenv(c("WGBS_ANALYSIS", "WGBS_GENE_LEVEL_GO", "WGBS_DIAGNOSTICS"))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
suppressPackageStartupMessages(library(data.table))
id <- commandArgs(trailingOnly = TRUE)[1]
a <- analysis_spec(id)
if (a$contrast_variant == "none" || a$ct_excluded || a$dml_difference != 20 || a$dml_qvalue != 0.01 || a$min_coverage != 10) {
  stop("Give an analysis with ctrl or inter at 10x, without ctexcl, d50 and q05 (the diagnostic's settings).", call. = FALSE)
}
mode <- a$filter                                     # the diagnostic's mode: legacy or corrected
diag_id <- paste0(if (a$contrast_variant == "ctrl") "A4_ploidy_in_control" else "A2_interaction", if (a$overdispersion == "MN") "_mn_chisq" else "")
contrast <- if (a$contrast_variant == "ctrl") "ploidy" else "interaction"
fails <- 0L
check <- function(ok, what) { cat(if (isTRUE(ok)) "ok  " else "FAIL", what, "\n"); if (!isTRUE(ok)) fails <<- fails + 1L }
t05 <- function(A, file) analysis_input(A, "05_differential_methylation", "tables", file)
rel <- function(x, y) max(abs(x - y) / pmax(abs(y), .Machine$double.xmin))

# The new contrast against the diagnostic's fit and list.
mine <- fread(t05(id, sprintf("diffmeth_%s_%s.csv.gz", contrast, id)))
fit <- readRDS(here("06_diagnostics", "dml_sensitivity", "cache", sprintf("fit_%s_%s.rds", diag_id, mode)))$fit
check(nrow(mine) == nrow(fit), sprintf("%s: %d CpGs tested, as in %s (%s)", id, nrow(mine), diag_id, mode))
check(rel(mine$pvalue, fit$pvalue) < 1e-12 && rel(mine$qvalue, fit$qvalue) < 1e-12,
      sprintf("p and q equal those of %s (largest relative differences %.1e, %.1e; the tables carry 15 digits)", diag_id,
              rel(mine$pvalue, fit$pvalue), rel(mine$qvalue, fit$qvalue)))
if (contrast == "ploidy") check(max(abs(mine$meth.diff - fit$meth.diff)) < 1e-10, "meth.diff equals that of the diagnostic (triploid minus diploid among the controls)")
d_mine <- fread(t05(id, sprintf("dml_%s_%s.csv", contrast, id)))
d_diag <- fread(here("06_diagnostics", "dml_sensitivity", sprintf("dml_%s_%s.csv", diag_id, mode)))
check(nrow(d_mine) == nrow(d_diag) && identical(paste(d_mine$chr, d_mine$start), paste(d_diag$chr, d_diag$start)) &&
        identical(d_mine$direction, d_diag$direction) && max(abs(d_mine$meth.diff - d_diag$effect_size), 0) < 1e-9,
      sprintf("the %s DML list equals dml_%s_%s.csv: %d DML (%d hyper, %d hypo), same loci, directions and effect sizes",
              contrast, diag_id, mode, nrow(d_mine), sum(d_mine$direction == "hyper"), sum(d_mine$direction == "hypo")))
if (contrast == "interaction") {
  g <- grep("^pct_", names(d_diag), value = TRUE)
  check(length(g) == 4 && max(abs(as.matrix(d_mine[, g, with = FALSE]) - as.matrix(d_diag[, g, with = FALSE])), 0) < 1e-9,
        "the four group percentages of the interaction DML equal those of the diagnostic")
}

# The contrasts taken from the base fit.
base <- a$base_fit_id
reused <- if (a$contrast_variant == "ctrl") "heat" else c("ploidy", "heat", "shared")
for (set in reused) {
  x <- fread(t05(id, sprintf("dml_%s_%s.csv", set, id))); y <- fread(t05(base, sprintf("dml_%s_%s.csv", set, base)))
  check(isTRUE(all.equal(x, y, check.attributes = FALSE)), sprintf("dml_%s of %s equals that of %s (%d loci)", set, id, base, nrow(x)))
}
k <- fread(t05(id, "dml_counts_by_mode.csv"))
if (contrast == "interaction") {
  sh <- fread(t05(id, sprintf("dml_shared_%s.csv", id)))
  check(k$interaction_dml == nrow(d_mine) && k$interaction_in_shared == sum(d_mine$united_row %in% sh$united_row),
        sprintf("dml_counts_by_mode: interaction_dml %d, interaction_in_shared %d", k$interaction_dml, k$interaction_in_shared))
} else {
  check(identical(names(k), names(fread(t05(base, "dml_counts_by_mode.csv")))) && k$ploidy_dml == nrow(d_mine),
        sprintf("dml_counts_by_mode: the columns of %s, ploidy_dml %d, shared %d", base, k$ploidy_dml, k$shared))
}
cat(if (fails) sprintf("\n%d check(s) fail\n", fails) else "\nall checks pass\n")
quit(status = if (fails) 1 else 0)

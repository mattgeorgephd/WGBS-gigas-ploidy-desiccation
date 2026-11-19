# The within-cohort ploidy differences (computed as the script's difference() does) against methylKit's meth.diff of the
# cached fits A4 (legacy mode), the interaction against the effect size of the A2 list, and the largest difference to step 05.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); N <- C + as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
pooled <- function(idx) 100 * rowSums(C[, idx, drop = FALSE]) / rowSums(N[, idx, drop = FALSE])
d_heat <- pooled(6:10) - pooled(1:5); d_ctrl <- pooled(17:22) - pooled(11:16)
for (a in c("heat", "control")) {
  f <- readRDS(file.path(repo, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_A4_ploidy_in_%s_legacy.rds", a)))$fit
  d <- if (a == "heat") d_heat else d_ctrl
  cat(sprintf("ploidy inside the %s cohort: %d CpGs; largest |difference() - meth.diff of methylKit| = %.3g; |meth.diff| > 20 in %d CpGs; above() counts %d; CpGs where the two decisions differ: %d\n",
              a, length(d), max(abs(d - f$meth.diff)), sum(abs(f$meth.diff) > 20), sum(abs(d) > 20 + 1e-9), sum((abs(f$meth.diff) > 20) != (abs(d) > 20 + 1e-9))))
}
for (w in c("ploidy", "heat")) {
  s5 <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/diffmeth_%s_legacy.csv.gz", w)))
  d <- if (w == "ploidy") pooled(c(6:10, 17:22)) - pooled(c(1:5, 11:16)) else pooled(1:10) - pooled(11:22)
  cat(sprintf("step 05, %s: largest |difference() - meth.diff| = %.3g (the script accepts up to 1e-8); decisions differ at %d CpGs\n", w, max(abs(d - s5$meth.diff)), sum((abs(s5$meth.diff) > 20) != (abs(d) > 20 + 1e-9))))
}

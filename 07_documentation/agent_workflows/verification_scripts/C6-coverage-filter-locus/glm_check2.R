suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
u <- fread(t5("united_cpg_corrected.csv.gz"))
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); T <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); W <- C + T
cat("entries where coverage != numCs + numTs after methylKit rounding:", sum(W != as.matrix(u[, paste0("coverage", 1:22), with = FALSE])), "of", length(W), "\n")
wdiff <- function(g) 100 * (rowSums(C[, g]) / rowSums(W[, g]) - rowSums(C[, !g]) / rowSums(W[, !g]))
for (contrast in c("ploidy", "heat")) {
  d <- fread(t5(sprintf("diffmeth_%s_corrected.csv.gz", contrast)))
  g <- if (contrast == "ploidy") trip else heat; cv <- if (contrast == "ploidy") heat else trip
  md <- wdiff(g)
  cat(sprintf("[%s] meth.diff recomputed for all %d CpGs: max abs diff %.3g; |diff| > 20: mine %d, table %d\n", contrast, nrow(d), max(abs(md - d$meth.diff)), sum(abs(md) > 20), sum(abs(d$meth.diff) > 20)))
  nd <- which(abs(md) > 20 & !(d$qvalue < 0.01))
  cat(sprintf("   CpGs with |diff| > 20 that are not DML (q >= 0.01): %d; their q range %.3g to %.3g; p range %.3g to %.3g\n", length(nd), min(d$qvalue[nd]), max(d$qvalue[nd]), min(d$pvalue[nd]), max(d$pvalue[nd])))
  p <- vapply(nd, function(i) { f1 <- suppressWarnings(glm(cbind(C[i, ], T[i, ]) ~ g + cv, family = binomial)); f0 <- suppressWarnings(glm(cbind(C[i, ], T[i, ]) ~ cv, family = binomial))
    pchisq(f0$deviance - f1$deviance, 1, lower.tail = FALSE) }, numeric(1))
  cat(sprintf("   refit of those: max relative diff of p %.3g; smallest refit p %.3g (largest p among DML in the table: %.3g)\n",
              max(abs(p - d$pvalue[nd]) / pmax(d$pvalue[nd], 1e-300)), min(p), max(d$pvalue[d$qvalue < 0.01 & abs(d$meth.diff) > 20])))
  cat(sprintf("   q is monotone in p across all CpGs: %s; largest p with q < 0.01: %.4g; smallest p with q >= 0.01: %.4g\n",
              !is.unsorted(d$qvalue[order(d$pvalue)]), max(d$pvalue[d$qvalue < 0.01]), min(d$pvalue[d$qvalue >= 0.01])))
}

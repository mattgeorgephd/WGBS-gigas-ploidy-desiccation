suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- file.path(repo, "04_output/05_differential_methylation/tables")
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
tp <- as.integer(lib$ploidy == "triploid"); th <- as.integer(lib$treatment == "heat")
lrt <- function(y, w, tr, cv) {
  full <- suppressWarnings(glm.fit(cbind(1, tr, cv), y / w, weights = w, family = binomial())); red <- suppressWarnings(glm.fit(cbind(1, cv), y / w, weights = w, family = binomial()))
  pchisq(red$deviance - full$deviance, 1, lower.tail = FALSE)
}
set.seed(11)
for (md in c("legacy", "corrected")) {
  u <- fread(file.path(t5, sprintf("united_cpg_%s.csv.gz", md)))
  nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
  for (ct in c("ploidy", "heat")) {
    d <- fread(file.path(t5, sprintf("diffmeth_%s_%s.csv.gz", ct, md)))
    dml <- which(d$qvalue < 0.01 & abs(d$meth.diff) > 20)
    idx <- c(sample(nrow(d), 1500), sample(dml, 500))
    tr <- if (ct == "ploidy") tp else th; cv <- if (ct == "ploidy") th else tp
    p <- vapply(idx, function(i) lrt(nc[i, ], nc[i, ] + nt[i, ], tr, cv), numeric(1))
    p[is.na(p)] <- 1
    rel <- abs(p - d$pvalue[idx]) / pmax(d$pvalue[idx], 1e-300)
    # wrong-model controls: no covariate, and swapped roles
    p_nocov <- vapply(idx[1:300], function(i) { y <- nc[i, ]; w <- nc[i, ] + nt[i, ]; f <- suppressWarnings(glm.fit(cbind(1, tr), y / w, weights = w, family = binomial())); pchisq(f$null.deviance - f$deviance, 1, lower.tail = FALSE) }, numeric(1))
    cat(sprintf("[%s %s] %d CpGs: max abs p diff %.3g ; rows with rel diff > 1e-6: %d ; median rel diff %.2g | control without covariate: rows with rel diff > 1e-6: %d of 300\n", md, ct, length(idx),
                max(abs(p - d$pvalue[idx])), sum(rel > 1e-6 & abs(p - d$pvalue[idx]) > 1e-9), median(rel), sum(abs(p_nocov - d$pvalue[idx[1:300]]) / pmax(d$pvalue[idx[1:300]], 1e-300) > 1e-6)))
  }
  rm(u, nc, nt); gc(verbose = FALSE)
}

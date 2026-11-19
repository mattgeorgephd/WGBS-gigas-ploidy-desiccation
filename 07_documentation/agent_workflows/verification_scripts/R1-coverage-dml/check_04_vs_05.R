suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- file.path(repo, "04_output/05_differential_methylation/tables"); t4 <- file.path(repo, "04_output/04_cpg_coverage/tables")
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
pv <- fread(file.path(t4, "pca_variance.csv")); cs <- fread(file.path(t4, "library_correlation_summary.csv"))
for (md in c("legacy", "corrected")) {
  u <- fread(file.path(t5, sprintf("united_cpg_%s.csv.gz", md)))
  nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
  mm <- nc / (nc + nt); colnames(mm) <- lib$library_id
  cat(sprintf("[%s] united rows %d ; step 04 n_cpg at 10x %d ; NaN rows %d\n", md, nrow(u), cs[mode == md & min_coverage == 10]$n_cpg, sum(!complete.cases(mm))))
  r5 <- cor(mm)
  r4f <- fread(file.path(t4, sprintf("library_correlation_10x_%s.csv", md))); r4 <- as.matrix(r4f[, -1]); rownames(r4) <- r4f$library_id
  cat(sprintf("   correlation: identical dimnames %s ; max |r(step 04) - r(from step 05 united)| = %.3g\n", identical(colnames(r4), colnames(r5)) && identical(rownames(r4), rownames(r5)), max(abs(r4 - r5))))
  off <- r5[upper.tri(r5)]; cat(sprintf("   r range from united: %.6f to %.6f ; step 04 summary: %.6f to %.6f\n", min(off), max(off), cs[mode == md & min_coverage == 10]$r_min, cs[mode == md & min_coverage == 10]$r_max))
  sds <- sqrt(rowSums((mm - rowMeans(mm))^2) / 21)
  pc <- prcomp(t(mm[sds > quantile(sds, 0.5), ]), center = TRUE, scale. = TRUE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  got <- pv[mode == md & min_coverage == 10][order(pc)]$percent_variance
  cat(sprintf("   PCA: PC1 %.4f%% PC2 %.4f%% from united ; step 04: %.4f%% %.4f%% ; max abs difference over 22 PCs %.3g\n", ve[1], ve[2], got[1], got[2], max(abs(ve - got))))
  sc <- fread(file.path(t4, sprintf("pca_scores_10x_%s.csv", md)))
  cat(sprintf("   PC1 score: max | |step 04| - |united| | = %.3g (score sd %.1f) ; PC2: %.3g\n", max(abs(abs(sc$PC1) - abs(pc$x[, 1]))), sd(pc$x[, 1]), max(abs(abs(sc$PC2) - abs(pc$x[, 2])))))
  rm(u, nc, nt, mm); gc(verbose = FALSE)
}

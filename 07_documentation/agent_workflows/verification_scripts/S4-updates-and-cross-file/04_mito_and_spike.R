# S4 check, part 4.
# (i) Mitochondrial CpGs (first bullet of the I-01 update): ratio between cohorts, coverage, methylated reads per CpG.
# (ii) "It does not replicate. That rules out a ploidy difference common to both experiments" (fourth bullet):
#      how large a common ploidy difference would the sign count of contrast_agreement.csv show?
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
out <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S4-updates-and-cross-file"
load(file.path(out, "united_legacy_matrices.RData"))
DH <- 1:5; TH <- 6:10; DC <- 11:16; TC <- 17:22; high <- c(3, 4, 5, 18, 19)
tol <- 1e-9
mt <- chr == "NC_001276.1"
mp <- 100 * colSums(C[mt, ]) / colSums(N[mt, ])
cat(sprintf("mitochondrial CpGs: mean of the ten heat library values %.3f%% ; of the twelve control values %.3f%% ; ratio %.2f ; medians %.3f and %.3f, ratio %.2f\n",
            mean(mp[1:10]), mean(mp[11:22]), mean(mp[11:22]) / mean(mp[1:10]), median(mp[1:10]), median(mp[11:22]), median(mp[11:22]) / median(mp[1:10])))
cat(sprintf("reads pooled over the libraries of a cohort: heat %.3f%% ; control %.3f%% ; ratio %.2f\n",
            100 * sum(C[mt, 1:10]) / sum(N[mt, 1:10]), 100 * sum(C[mt, 11:22]) / sum(N[mt, 11:22]), (sum(C[mt, 11:22]) / sum(N[mt, 11:22])) / (sum(C[mt, 1:10]) / sum(N[mt, 1:10]))))
cat(sprintf("ratio of the range ends: %.2f (lower ends) and %.2f (upper ends)\n", min(mp[11:22]) / min(mp[1:10]), max(mp[11:22]) / max(mp[1:10])))
medcov <- apply(N, 2, median)
mtcov <- colMeans(N[mt, ]); mtC <- colMeans(C[mt, ])
tab <- data.table(library = lib$library_id, pct = round(mp, 3), mean_cov_mito = round(mtcov, 1), median_cov_all = medcov, rel_cov = round(mtcov / medcov, 2), methylated_reads_per_cpg = round(mtC, 2))
print(tab)
cat(sprintf("relative coverage of the mitochondrial CpGs (mean coverage over library median): heat %.1f to %.1f ; control %.1f to %.1f\n", min(tab$rel_cov[1:10]), max(tab$rel_cov[1:10]), min(tab$rel_cov[11:22]), max(tab$rel_cov[11:22])))
cat(sprintf("methylated reads per mitochondrial CpG (normalized counts): heat %.2f to %.2f ; control %.2f to %.2f\n", min(mtC[1:10]), max(mtC[1:10]), min(mtC[11:22]), max(mtC[11:22])))
cat("Spearman, across the 22 libraries, of mitochondrial percent methylation with relative mitochondrial coverage:", round(cor(mp, mtcov / medcov, method = "spearman"), 3), "; within heat", round(cor(mp[1:10], (mtcov / medcov)[1:10], method = "spearman"), 3), "; within control", round(cor(mp[11:22], (mtcov / medcov)[11:22], method = "spearman"), 3), "\n")
# is the excess spread over the mitochondrial CpGs or concentrated at a few?
pm_h <- 100 * rowSums(C[mt, 1:10]) / rowSums(N[mt, 1:10]); pm_c <- 100 * rowSums(C[mt, 11:22]) / rowSums(N[mt, 11:22])
cat(sprintf("per CpG, control minus heat (pooled): quartiles %s ; CpGs with control > heat: %d of %d ; largest %.1f\n", paste(round(quantile(pm_c - pm_h, c(.25, .5, .75)), 2), collapse = " "), sum(pm_c > pm_h), sum(mt), max(pm_c - pm_h)))

cat("\n== spike-in: a ploidy difference common to both cohorts, added to K eligible CpGs, 17 libraries (five high ones set aside)\n")
pct <- function(Cm, S) 100 * rowSums(Cm[, S, drop = FALSE]) / rowSums(N[, S, drop = FALSE])
DHn <- setdiff(DH, high); TCn <- setdiff(TC, high)
stat <- function(Cm) { dh <- pct(Cm, TH) - pct(Cm, DHn); dc <- pct(Cm, TCn) - pct(Cm, DC); b <- abs(dh) > 20 + tol & abs(dc) > 20 + tol
  c(both = sum(b), same = sum(b & sign(dh) == sign(dc)), r = cor(dh, dc)) }
base <- stat(C); print(base)
m22 <- rowMeans(P); pool <- which(elig & m22 > 30 & m22 < 60)
set.seed(4402)
for (shift in c(25, 40)) for (K in c(50, 100, 200, 500)) {
  rs <- t(vapply(1:20, function(i) {
    idx <- sample(pool, K); Cm <- C
    add <- round(N[idx, c(TH, TC)] * shift / 100)
    Cm[idx, c(TH, TC)] <- pmin(N[idx, c(TH, TC)], C[idx, c(TH, TC)] + add)
    s <- stat(Cm)
    dh <- pct(Cm, TH) - pct(Cm, DHn); dc <- pct(Cm, TCn) - pct(Cm, DC)
    c(s, spiked_both_same = sum(abs(dh[idx]) > 20 + tol & abs(dc[idx]) > 20 + tol & sign(dh[idx]) == sign(dc[idx])))
  }, numeric(4)))
  cat(sprintf("shift +%d points in all triploid libraries at K = %3d CpGs: CpGs above 20 in both %.0f, same sign %.0f (%.1f%%; without spike %.1f%%) ; spiked CpGs above 20 in both with the same sign: %.0f of %d ; Pearson %.4f\n",
              shift, K, mean(rs[, "both"]), mean(rs[, "same"]), 100 * mean(rs[, "same"] / rs[, "both"]), 100 * base["same"] / base["both"], mean(rs[, "spiked_both_same"]), K, mean(rs[, "r"])))
}
cat("done\n")

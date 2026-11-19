# 09_summaries.R
# Summary numbers quoted in the report: similarity structure at equal coverage,
# level against PC1 and the eligible mean among the 17 libraries that are not
# high, the enumeration over 4,410 regroupings, and the reverse slope direction.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); id <- lib$library_id
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
dip <- which(lib$ploidy == "diploid"); tri <- which(lib$ploidy == "triploid")
dh <- intersect(dip, heat); th <- intersect(tri, heat); dc <- intersect(dip, ctrl); tc <- intersect(tri, ctrl)
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"))
nh <- setdiff(1:22, high)

cat("==== similarity at equal coverage (mean of 5 thinning seeds) and noise-corrected ====\n")
Mt <- as.matrix(fread(file.path(scratch, "04_thinned_correlation_22.csv"))[, -1]); Rd <- readRDS(file.path(scratch, "03_correlations.rds"))$Rd
st <- as.matrix(fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_correlation_10x_legacy.csv"))[, -1])
within <- function(M, a) { s <- M[a, a]; s[upper.tri(s)] }
core7 <- which(id %in% paste0("zr3644_", c(3, 4, 5, 6, 13, 16, 18)))
c10 <- setdiff(ctrl, high)
rep_line <- function(lab, v) cat(sprintf("  %-62s mean %.4f  range %.4f to %.4f  (n = %d)\n", lab, mean(v), min(v), max(v), length(v)))
for (nm in c("step 04 (observed)", "thinned to equal coverage", "noise-corrected")) {
  M <- switch(nm, "step 04 (observed)" = st, "thinned to equal coverage" = Mt, "noise-corrected" = Rd)
  cat(nm, "\n")
  rep_line("heat triploids among themselves", within(M, th))
  rep_line("ten control libraries that are not high, among themselves", within(M, c10))
  rep_line("seven control libraries 3, 4, 5, 6, 13, 16, 18 among themselves", within(M, core7))
  rep_line("heat triploids with the ten control libraries not high", as.vector(M[th, c10]))
  rep_line("heat triploids with the seven control libraries", as.vector(M[th, core7]))
  rep_line("D11 and D12 with the five heat triploids", as.vector(M[1:2, th]))
  rep_line("D11 and D12 with the ten control libraries not high", as.vector(M[1:2, c10]))
  rep_line("D11 with D12", M[1, 2])
  rep_line("control diploids among themselves", within(M, dc)); rep_line("control triploids among themselves", within(M, tc))
}
cat("mean correlation of each library that is not high with the other 16 that are not high (thinned to equal coverage), sorted:\n")
typ <- sapply(nh, function(j) mean(Mt[j, setdiff(nh, j)])); names(typ) <- id[nh]; print(round(sort(typ, decreasing = TRUE), 4))
typ_d <- sapply(nh, function(j) mean(Rd[j, setdiff(nh, j)])); names(typ_d) <- id[nh]; cat("the same, noise-corrected:\n"); print(round(sort(typ_d, decreasing = TRUE), 4))

cat("\n==== level, PC1 and eligible mean among the 17 libraries that are not high ====\n")
pcs <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/pca_scores_10x_legacy.csv"))
lv <- fread(file.path(repo, "06_diagnostics/library_methylation_level/library_level.csv"))
me <- lv$mean_pct_eligible_legacy; mt <- lv$mean_pct_tested_legacy
f <- function(x, y, idx, lab) cat(sprintf("  %-55s Pearson %.3f  Spearman %.3f (p %.3f)\n", lab, cor(x[idx], y[idx]), cor(x[idx], y[idx], method = "spearman"),
                                         suppressWarnings(cor.test(x[idx], y[idx], method = "spearman", exact = TRUE)$p.value)))
f(L, pcs$PC1, 1:22, "level vs PC1, all 22"); f(L, pcs$PC1, nh, "level vs PC1, the 17 not high"); f(L, pcs$PC1, setdiff(nh, which(id == "zr3644_1")), "level vs PC1, the 17 without zr3644_1")
f(L, me, 1:22, "level vs eligible mean, all 22"); f(L, me, nh, "level vs eligible mean, the 17 not high"); f(L, me, setdiff(nh, which(id == "zr3644_1")), "level vs eligible mean, the 17 without zr3644_1")
f(L, mt, nh, "level vs mean over all tested CpGs, the 17 not high")
f(L, me, heat, "level vs eligible mean, heat cohort"); f(L, me, ctrl, "level vs eligible mean, control cohort")
cat("eligible mean, 17 libraries: exact Wilcoxon heat 2 D vs 5 T p =", round(wilcox.test(me[setdiff(dh, high)], me[th], exact = TRUE)$p.value, 4),
    "; control 6 D vs 4 T p =", round(wilcox.test(me[dc], me[setdiff(tc, high)], exact = TRUE)$p.value, 4), "\n")

cat("\n==== enumeration: both cohorts without the five, labels moved within cohort (4,410 regroupings) ====\n")
m <- readRDS(file.path(scratch, "counts_legacy.rds")); C <- m$C; N <- m$N
h7 <- setdiff(heat, high); c10 <- setdiff(ctrl, high)
pm <- 100 * C[, c(h7, c10)] / N[, c(h7, c10)]
span <- apply(pm, 1, max) - apply(pm, 1, min); keep <- which(span > 20)
Ck <- C[keep, ]; Nk <- N[keep, ]; storage.mode(Ck) <- "double"; storage.mode(Nk) <- "double"
ch <- combn(h7, 2, simplify = FALSE); cc <- combn(c10, 6, simplify = FALSE)
grid <- expand.grid(i = seq_along(ch), j = seq_along(cc))
cnt <- integer(nrow(grid))
for (s0 in seq(1, nrow(grid), by = 150)) {
  rows <- s0:min(s0 + 149, nrow(grid))
  A <- matrix(0, 22, length(rows)); B <- A
  for (k in seq_along(rows)) { a <- c(ch[[grid$i[rows[k]]]], cc[[grid$j[rows[k]]]]); A[a, k] <- 1; B[setdiff(c(h7, c10), a), k] <- 1 }
  d <- 100 * ((Ck %*% B) / (Nk %*% B) - (Ck %*% A) / (Nk %*% A))
  cnt[rows] <- colSums(abs(d) > 20)
}
obs <- cnt[sapply(seq_len(nrow(grid)), function(r) setequal(c(ch[[grid$i[r]]], cc[[grid$j[r]]]), c(setdiff(dh, high), dc)))]
cat("observed", obs, "; median", median(cnt), "; largest other", max(cnt[cnt != obs]), "; at or above observed", sum(cnt >= obs), "of", length(cnt), "\n")

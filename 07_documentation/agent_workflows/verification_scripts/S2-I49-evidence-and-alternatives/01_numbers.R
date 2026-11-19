# 01_numbers.R
# Task (a): recompute the key numbers of I-49 with independent code.
# Nothing is sourced from the repository; only the input tables are read.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 200, digits = 7)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"

lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
stopifnot(identical(lib$library_id, bis$library_id))
id <- lib$library_id
dh <- which(lib$ploidy == "diploid" & lib$treatment == "heat")
th <- which(lib$ploidy == "triploid" & lib$treatment == "heat")
dc <- which(lib$ploidy == "diploid" & lib$treatment == "control")
tc <- which(lib$ploidy == "triploid" & lib$treatment == "control")
stopifnot(length(dh) == 5, length(th) == 5, length(dc) == 6, length(tc) == 6)

cat("==== 1. Genome-wide level, from the Bismark counts ====\n")
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)
cat("max |my level - pct_mCpG column|:", max(abs(L - bis$pct_mCpG)), "\n")
names(L) <- id
s <- sort(L)
print(round(s, 3))
g <- diff(s)
o <- order(-g)
cat("gaps, largest first:\n")
print(data.table(below = names(s)[o], below_value = round(s[o], 3), above = names(s)[o + 1], above_value = round(s[o + 1], 3), gap = round(g[o], 3))[1:6])
high <- which(L > s[which.max(g)])
cat("libraries above the largest gap:", id[high], "\n")
cat("range of the five:", round(range(L[high]), 3), " range of the other 17:", round(range(L[-high]), 3), "\n")

m <- readRDS(file.path(scratch, "counts_legacy.rds"))
C <- m$C; N <- m$N
pm <- 100 * C / N
mean_pm <- rowMeans(pm)
elig <- mean_pm > 10 & mean_pm < 90
elig_incl <- mean_pm >= 10 & mean_pm <= 90
cat("\ntested CpGs:", nrow(pm), " eligible (strict):", sum(elig), " eligible (inclusive bounds):", sum(elig_incl), "\n")
me <- colMeans(pm[elig, ]); names(me) <- id
mt <- colMeans(pm); names(mt) <- id
cat("mean over eligible CpGs, sorted:\n"); print(round(sort(me), 2))
ge <- diff(sort(me)); cat("largest gaps in the mean over eligible:", round(sort(ge, decreasing = TRUE)[1:3], 2), "\n")
cat("libraries above the largest gap of the eligible mean:", names(sort(me))[(which.max(ge) + 1):22], "\n")
cat("five:", round(range(me[high]), 1), " other 17:", round(range(me[-high]), 1), "\n")
cat("cor(level, mean eligible) Pearson:", round(cor(L, me), 4), " Spearman:", round(cor(L, me, method = "spearman"), 4), "\n")
cat("cor(level, mean tested) Pearson:", round(cor(L, mt), 4), "\n")

# PCA of the tested CpGs: own computation (centred, not scaled, and scaled), to compare with 0.93 and 9.5%
pcs <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/pca_scores_10x_legacy.csv"))
cat("cor(level, PC1 of step 04):", round(cor(L, pcs$PC1), 4), "\n")
sdv <- apply(pm, 1, sd)
keep <- sdv > 0
x <- t(pm[keep, ])
p1 <- prcomp(x, center = TRUE, scale. = FALSE)
p2 <- prcomp(x, center = TRUE, scale. = TRUE)
cat("own PCA (centred, unscaled): PC1 variance share", round(100 * p1$sdev[1]^2 / sum(p1$sdev^2), 2), "; |cor| with level", round(abs(cor(L, p1$x[, 1])), 4), "\n")
cat("own PCA (centred, scaled):   PC1 variance share", round(100 * p2$sdev[1]^2 / sum(p2$sdev^2), 2), "; |cor| with level", round(abs(cor(L, p2$x[, 1])), 4), "\n")
rm(x, p1, p2); invisible(gc())

cat("\n==== 2. Non-CpG methylation ====\n")
chh <- 100 * bis$meth_CHH / (bis$meth_CHH + bis$unmeth_CHH)
chg <- 100 * bis$meth_CHG / (bis$meth_CHG + bis$unmeth_CHG)
nonCpG <- 100 * (bis$meth_CHH + bis$meth_CHG) / (bis$meth_CHH + bis$unmeth_CHH + bis$meth_CHG + bis$unmeth_CHG)
print(data.table(id, group = paste(lib$ploidy, lib$treatment), high = seq_along(id) %in% high, CpG = round(L, 3), CHG = round(chg, 4), CHH = round(chh, 4), nonCpG = round(nonCpG, 4)))
cat("CHH range:", round(range(chh), 4), " CHG range:", round(range(chg), 4), " CHG+CHH range:", round(range(nonCpG), 4), "\n")
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
for (nm in c("chh", "chg", "nonCpG")) {
  v <- get(nm)
  cat(sprintf("%s vs CpG level: Spearman all %.4f, heat %.4f (p %.4f), control %.4f (p %.4f); Pearson heat %.3f, control %.3f\n", nm,
              cor(L, v, method = "spearman"), cor(L[heat], v[heat], method = "spearman"),
              cor.test(L[heat], v[heat], method = "spearman")$p.value,
              cor(L[ctrl], v[ctrl], method = "spearman"), cor.test(L[ctrl], v[ctrl], method = "spearman")$p.value,
              cor(L[heat], v[heat]), cor(L[ctrl], v[ctrl])))
}
# heat cohort without the two acute-heat-shock triploids, whose non-CpG values are the two lowest
h8 <- setdiff(heat, which(id %in% c("zr3534_9", "zr3534_10")))
cat("heat cohort without zr3534_9 and zr3534_10: Spearman CHH vs level", round(cor(L[h8], chh[h8], method = "spearman"), 4),
    " p", round(cor.test(L[h8], chh[h8], method = "spearman")$p.value, 4), "\n")
cat("difference in mean CHH, five high minus 17 others:", round(mean(chh[high]) - mean(chh[-high]), 4), "\n")
cat("difference in mean CpG level, five high minus 17 others:", round(mean(L[high]) - mean(L[-high]), 3), "\n")

cat("\n==== 3. Wilcoxon tests of the level ====\n")
w <- function(a, b) wilcox.test(L[a], L[b], exact = TRUE)$p.value
cat("heat cohort D vs T:", round(w(dh, th), 5), " min D - max T:", round(min(L[dh]) - max(L[th]), 4), "\n")
cat("control cohort D vs T:", round(w(dc, tc), 5), " triploids above every diploid:", sum(L[tc] > max(L[dc])), "\n")
cat("heat vs control cohort:", round(w(heat, ctrl), 4), "\n")
cat("diploid vs triploid, all 22:", round(w(c(dh, dc), c(th, tc)), 4), "\n")
cat("group means:", round(c(DH = mean(L[dh]), TH = mean(L[th]), DC = mean(L[dc]), TC = mean(L[tc])), 3), "\n")
# Without the five high libraries
cat("-- without the five --\n")
dh2 <- setdiff(dh, high); tc2 <- setdiff(tc, high)
cat("heat cohort, 2 D vs 5 T: exact Wilcoxon p", round(w(dh2, th), 4), " D values", round(L[dh2], 3), " T values", round(L[th], 3), "\n")
cat("control cohort, 6 D vs 4 T: exact Wilcoxon p", round(w(dc, tc2), 4), " D values", round(L[dc], 3), " T values", round(L[tc2], 3), "\n")
d22 <- data.table(L = L, ploidy = factor(lib$ploidy), treatment = factor(lib$treatment))
a22 <- anova(lm(L ~ ploidy * treatment, d22)); cat("2x2 ANOVA, 22 libraries, interaction F and p:", round(a22$`F value`[3], 2), signif(a22$`Pr(>F)`[3], 3), "\n")
d17 <- d22[-high]
f17 <- lm(L ~ ploidy * treatment, d17)
a17 <- anova(f17); cat("2x2 ANOVA, 17 libraries (five set aside), interaction F(1,", a17$Df[4], ") and p:", round(a17$`F value`[3], 2), signif(a17$`Pr(>F)`[3], 3), "\n")
print(round(coef(summary(f17)), 4))
gm17 <- d17[, .(n = .N, mean = round(mean(L), 3)), by = .(ploidy, treatment)]; print(gm17)
# the same with the mean over eligible CpGs
d22e <- data.table(L = me, ploidy = factor(lib$ploidy), treatment = factor(lib$treatment))
a <- anova(lm(L ~ ploidy * treatment, d22e)); cat("eligible mean: 22 libraries interaction F and p:", round(a$`F value`[3], 2), signif(a$`Pr(>F)`[3], 3), "\n")
a <- anova(lm(L ~ ploidy * treatment, d22e[-high])); cat("eligible mean: 17 libraries interaction F and p:", round(a$`F value`[3], 2), signif(a$`Pr(>F)`[3], 3), "\n")
print(d22e[-high, .(n = .N, mean = round(mean(L), 2)), by = .(ploidy, treatment)])
# permutation check for the 17-library interaction: labels moved within cohort (ploidy), all regroupings
int_stat <- function(v, a1, b1, a2, b2) (mean(v[b1]) - mean(v[a1])) - (mean(v[b2]) - mean(v[a2]))
obs <- int_stat(L, dh2, th, dc, tc2)
h7 <- c(dh2, th); c10 <- c(dc, tc2)
ch <- combn(h7, 2, simplify = FALSE); cc <- combn(c10, 6, simplify = FALSE)
ref <- unlist(lapply(ch, function(a1) vapply(cc, function(a2) int_stat(L, a1, setdiff(h7, a1), a2, setdiff(c10, a2)), numeric(1))))
cat("interaction contrast in the level, 17 libraries:", round(obs, 3), "; regroupings (within cohort) with |contrast| >= observed:",
    sum(abs(ref) >= abs(obs) - 1e-12), "of", length(ref), "\n")
obs22 <- int_stat(L, dh, th, dc, tc)
cat("interaction contrast in the level, 22 libraries:", round(obs22, 3), "\n")

cat("\n==== 4. Contrasts across CpGs (pooled reads, > 20 points) ====\n")
pool <- function(idx) 100 * rowSums(C[, idx, drop = FALSE]) / rowSums(N[, idx, drop = FALSE])
dHeat <- pool(th) - pool(dh)
dCtrl <- pool(tc) - pool(dc)
dAll <- pool(c(th, tc)) - pool(c(dh, dc))
dCoh <- pool(heat) - pool(ctrl)
cnt <- function(d, lab) cat(sprintf("%-45s above 20: %6d  positive %6d  negative %6d  majority %.3f  median eligible %.2f  (|d| >= 20: %d)\n",
                                    lab, sum(abs(d) > 20), sum(d > 20), sum(d < -20), max(sum(d > 20), sum(d < -20)) / sum(abs(d) > 20), median(d[elig]), sum(abs(d) >= 20)))
cnt(dAll, "ploidy, all 22"); cnt(dCoh, "cohort, all 22"); cnt(dHeat, "ploidy inside heat cohort"); cnt(dCtrl, "ploidy inside control cohort")
cnt(dHeat - dCtrl, "interaction")
# check against step 05
s5 <- fread(file.path(repo, "04_output/05_differential_methylation/tables/diffmeth_ploidy_legacy.csv.gz"))
cat("max |my all-22 ploidy difference - step 05 meth.diff|:", max(abs(dAll - s5$meth.diff)), "\n")
cat("correlation of the two within-cohort differences: Pearson all", round(cor(dHeat, dCtrl), 4), " Spearman all", round(cor(dHeat, dCtrl, method = "spearman"), 4),
    " Pearson eligible", round(cor(dHeat[elig], dCtrl[elig]), 4), " Spearman eligible", round(cor(dHeat[elig], dCtrl[elig], method = "spearman"), 4), "\n")
both <- abs(dHeat) > 20 & abs(dCtrl) > 20
cat("above 20 in both:", sum(both), " same sign:", sum(both & sign(dHeat) == sign(dCtrl)), "\n")
dHeat2 <- pool(th) - pool(dh2); dCtrl2 <- pool(tc2) - pool(dc)
cnt(dHeat2, "heat cohort without its 3 high"); cnt(dCtrl2, "control cohort without its 2 high")
cat("without the five: Pearson all", round(cor(dHeat2, dCtrl2), 4), " Spearman all", round(cor(dHeat2, dCtrl2, method = "spearman"), 4),
    " Pearson eligible", round(cor(dHeat2[elig], dCtrl2[elig]), 4), "\n")
both2 <- abs(dHeat2) > 20 & abs(dCtrl2) > 20
cat("above 20 in both (without the five):", sum(both2), " same sign:", sum(both2 & sign(dHeat2) == sign(dCtrl2)), "\n")
saveRDS(list(L = L, me = me, elig = elig, high = high, dHeat = dHeat, dCtrl = dCtrl), file.path(scratch, "derived_legacy.rds"))

cat("\n==== 5. Enumerations ====\n")
enum <- function(members, k) {
  cmb <- combn(members, k, simplify = FALSE)
  n <- vapply(cmb, function(a) { b <- setdiff(members, a); sum(abs(pool(b) - pool(a)) > 20) }, numeric(1))
  list(cmb = cmb, n = n)
}
e7 <- enum(c(dh2, th), 2)
tab7 <- data.table(group_of_two = vapply(e7$cmb, function(a) paste(id[a], collapse = " "), ""), n = e7$n)[order(-n)]
print(tab7)
cat("observed (zr3534_1 zr3534_2):", tab7[group_of_two == "zr3534_1 zr3534_2", n], " median of 21:", median(tab7$n), " largest other:", max(tab7[group_of_two != "zr3534_1 zr3534_2", n]),
    " at or above observed:", sum(tab7$n >= tab7[group_of_two == "zr3534_1 zr3534_2", n]), "\n")
e10 <- enum(heat, 5)
obs10 <- e10$n[vapply(e10$cmb, function(a) setequal(a, dh), TRUE)]
cat("heat cohort 5 vs 5: observed", obs10, " median", median(e10$n), " at or above:", sum(e10$n >= obs10), "of", length(e10$n),
    " largest other (excluding the mirror):", max(e10$n[!vapply(e10$cmb, function(a) setequal(a, dh) || setequal(a, th), TRUE)]), "\n")
e12 <- enum(ctrl, 6)
obs12 <- e12$n[vapply(e12$cmb, function(a) setequal(a, dc), TRUE)]
cat("control cohort 6 vs 6: observed", obs12, " median", median(e12$n), " at or above:", sum(e12$n >= obs12), "of", length(e12$n),
    " largest other:", max(e12$n[!vapply(e12$cmb, function(a) setequal(a, dc) || setequal(a, tc), TRUE)]), "\n")
e10c <- enum(c(dc, tc2), 6)
obs10c <- e10c$n[vapply(e10c$cmb, function(a) setequal(a, dc), TRUE)]
cat("control cohort without its 2 high, 6 vs 4: observed", obs10c, " median", median(e10c$n), " at or above:", sum(e10c$n >= obs10c), "of", length(e10c$n), " largest other:", max(e10c$n[e10c$n != obs10c | !vapply(e10c$cmb, function(a) setequal(a, dc), TRUE)]), "\n")
# direction of the 2-vs-5 difference
cat("2 vs 5, direction: diploid higher", sum(dHeat2 < -20), "of", sum(abs(dHeat2) > 20), "=", round(sum(dHeat2 < -20) / sum(abs(dHeat2) > 20), 3), "\n")

cat("\n==== 6. Pairwise correlations, heat cohort ====\n")
cm <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_correlation_10x_legacy.csv"))
R <- as.matrix(cm[, -1]); rownames(R) <- colnames(R)
stopifnot(identical(colnames(R), id))
# own correlation from the united table
Rown <- cor(pm)
cat("max |own correlation of the united table - step 04 matrix|:", signif(max(abs(Rown - R)), 3), "\n")
mw <- function(M, idx) { s <- M[idx, idx]; mean(s[upper.tri(s)]) }
sub5 <- combn(heat, 5, simplify = FALSE)
v <- vapply(sub5, function(a) mw(R, a), numeric(1))
cat("mean pairwise r, heat triploids:", round(mw(R, th), 4), " rank among 252:", rank(-v)[vapply(sub5, function(a) setequal(a, th), TRUE)],
    " next best:", round(sort(v, decreasing = TRUE)[2], 4), "(", id[sub5[[order(-v)[2]]]], ")\n")
cat("mean pairwise r, heat diploids:", round(mw(R, dh), 4), " between:", round(mean(R[dh, th]), 4), "\n")
pr <- t(combn(heat, 2)); prd <- data.table(a = id[pr[, 1]], b = id[pr[, 2]], r = R[pr], both_T = pr[, 1] %in% th & pr[, 2] %in% th)[order(-r)]
print(prd)
cat("range of the T-T pairs:", round(range(prd[both_T == TRUE, r]), 4), " range of the other 35:", round(range(prd[both_T == FALSE, r]), 4), "\n")
cat("control: among diploids", round(mw(R, dc), 4), " among triploids", round(mw(R, tc), 4), " between", round(mean(R[dc, tc]), 4), "\n")

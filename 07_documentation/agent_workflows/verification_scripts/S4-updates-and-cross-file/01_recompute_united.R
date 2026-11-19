# S4 check, part 1: numbers of the updates of I-01, I-02, I-04, I-07, I-08 recomputed
# from the legacy united table and the cached fits. Read-only on the repository.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
root <- "/home/user/WGBS-gigas-ploidy-desiccation"
out  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S4-updates-and-cross-file"
u <- fread(file.path(root, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
lib <- fread(file.path(root, "03_input_files/metadata/wgbs_libraries.csv"))
stopifnot(nrow(u) == 480579, nrow(lib) == 22)
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); storage.mode(C) <- "double"
Tm <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); storage.mode(Tm) <- "double"
COV <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE]); storage.mode(COV) <- "double"
N <- C + Tm
cat("rows where coverage != numCs + numTs in at least one library:", sum(rowSums(COV != N) > 0), "; cells:", sum(COV != N), "; largest absolute difference:", max(abs(COV - N)), "\n")
rm(Tm, COV)
chr <- u$chr; pos <- u$start
rm(u); invisible(gc())
P <- 100 * C / N
DH <- 1:5; TH <- 6:10; DC <- 11:16; TC <- 17:22
high <- c(3, 4, 5, 18, 19); acute <- c(4, 5, 9, 10)
stopifnot(all(lib$ploidy[DH] == "diploid"), all(lib$ploidy[TH] == "triploid"), all(lib$treatment[1:10] == "heat"),
          all(lib$ploidy[DC] == "diploid"), all(lib$ploidy[TC] == "triploid"), all(lib$acute_heat_shock[acute] == "yes"),
          identical(lib$library_id[high], c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15")))
pct <- function(S) 100 * rowSums(C[, S, drop = FALSE]) / rowSums(N[, S, drop = FALSE])
tol <- 1e-9
ab  <- function(d) abs(d) > 20 + tol     # rule of library_methylation_level.R
ab0 <- function(d) abs(d) > 20           # rule of dml_sensitivity.R
cat("== (a) mitochondrial CpGs\n")
mt <- chr == "NC_001276.1"
cat("tested mitochondrial CpGs:", sum(mt), "\n")
mp <- 100 * colSums(C[mt, ]) / colSums(N[mt, ])
print(round(setNames(mp, lib$library_id), 3))
cat(sprintf("heat libraries: %.2f to %.2f ; control libraries: %.2f to %.2f ; ratio of means %.2f ; ratio of medians %.2f\n",
            min(mp[1:10]), max(mp[1:10]), min(mp[11:22]), max(mp[11:22]), mean(mp[11:22]) / mean(mp[1:10]), median(mp[11:22]) / median(mp[1:10])))
cat("mean of per-CpG percentages (unweighted), by library:\n"); print(round(colMeans(P[mt, ]), 3))

cat("\n== (b) eligible CpGs and contrasts\n")
m22 <- rowMeans(P)
elig <- m22 >= 10 & m22 <= 90
cat("eligible (10 <= mean <= 90):", sum(elig), " strict:", sum(m22 > 10 & m22 < 90), "\n")
d_all <- pct(c(TH, TC)) - pct(c(DH, DC))
d_coh <- pct(1:10) - pct(11:22)
d_h <- pct(TH) - pct(DH); d_c <- pct(TC) - pct(DC); d_int <- d_h - d_c
cnt <- function(d) c(above_tol = sum(ab(d)), above_strict = sum(ab0(d)), pos = sum(d > 20 + tol), neg = sum(d < -20 - tol), at = sum(abs(abs(d) - 20) <= tol))
print(rbind(ploidy_all = cnt(d_all), cohort = cnt(d_coh), ploidy_heat = cnt(d_h), ploidy_control = cnt(d_c), interaction = cnt(d_int)))
cat(sprintf("interaction: share negative %.4f ; median over eligible %.2f\n", sum(d_int < -20 - tol) / sum(ab(d_int)), median(d_int[elig])))
cat(sprintf("pearson d_h vs d_c, all tested %.4f ; eligible %.4f ; spearman all %.4f\n", cor(d_h, d_c), cor(d_h[elig], d_c[elig]), cor(d_h, d_c, method = "spearman")))
both <- ab(d_h) & ab(d_c)
cat("above 20 in both cohorts:", sum(both), " same sign:", sum(both & sign(d_h) == sign(d_c)), "\n")
print(table(heat = sign(d_h[both]), control = sign(d_c[both])))
ph <- mean(d_h[both] > 0); pc <- mean(d_c[both] > 0)
cat(sprintf("  expected same-sign share from the margins of these CpGs: %.3f ; observed %.3f\n", ph * pc + (1 - ph) * (1 - pc), mean(sign(d_h[both]) == sign(d_c[both]))))

cat("\n== (c) the draft's 1,038 against the within-cohort differences (20-point rule only)\n")
L <- ab0(d_all)
cat("pooled above 20:", sum(L), "; of these above 20 in heat cohort:", sum(L & ab0(d_h)), "; in control cohort:", sum(L & ab0(d_c)),
    "; in both:", sum(L & ab0(d_h) & ab0(d_c)), "; in neither:", sum(L & !ab0(d_h) & !ab0(d_c)), "\n")
cat("  same sign as pooled and above 20: heat", sum(L & ab0(d_h) & sign(d_h) == sign(d_all)), "; control", sum(L & ab0(d_c) & sign(d_c) == sign(d_all)), "\n")

cat("\n== (d) without the five high libraries\n")
DHn <- setdiff(DH, high); TCn <- setdiff(TC, high)
d_h2 <- pct(TH) - pct(DHn); d_c2 <- pct(TCn) - pct(DC); d_int2 <- d_h2 - d_c2
d_all2 <- pct(c(TH, TCn)) - pct(c(DHn, DC))
print(rbind(ploidy_17 = cnt(d_all2), ploidy_heat_7 = cnt(d_h2), ploidy_control_10 = cnt(d_c2), interaction_17 = cnt(d_int2)))
cat(sprintf("interaction without the five: share negative %.4f ; median over eligible %.2f\n", sum(d_int2 < -20 - tol) / sum(ab(d_int2)), median(d_int2[elig])))
cat(sprintf("pearson d_h2 vs d_c2 all %.4f ; eligible %.4f\n", cor(d_h2, d_c2), cor(d_h2[elig], d_c2[elig])))
both2 <- ab(d_h2) & ab(d_c2)
cat("above 20 in both cohorts (17 libraries):", sum(both2), " same sign:", sum(both2 & sign(d_h2) == sign(d_c2)), "\n")
print(table(heat = sign(d_h2[both2]), control = sign(d_c2[both2])))
ph <- mean(d_h2[both2] > 0); pc <- mean(d_c2[both2] > 0)
cat(sprintf("  margins: heat positive %.3f, control positive %.3f ; expected same-sign share if signs unrelated %.3f ; observed %.3f\n",
            ph, pc, ph * pc + (1 - ph) * (1 - pc), mean(sign(d_h2[both2]) == sign(d_c2[both2]))))
# clustering: how many 1 kb clusters do the 923 CpGs form?
o <- order(chr[both2], pos[both2]); cc <- chr[both2][o]; pp <- pos[both2][o]
newc <- c(TRUE, cc[-1] != cc[-length(cc)] | diff(pp) > 1000)
cat("  the CpGs above 20 in both (17 libraries) form", sum(newc), "clusters at 1 kb\n")
cl <- cumsum(newc); same <- (sign(d_h2[both2]) == sign(d_c2[both2]))[o]
cls <- tapply(same, cl, mean)
cat(sprintf("  per cluster: share of clusters with majority same sign %.3f (n clusters %d)\n", mean(cls > 0.5), length(cls)))

save(C, N, P, chr, pos, elig, lib, file = file.path(out, "united_legacy_matrices.RData"), compress = FALSE)

cat("\n== (e) cached fits\n")
rd <- function(a, m = "legacy") readRDS(file.path(root, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_%s.rds", a, m)))$fit
h0 <- rd("A0_heat"); hf <- rd("A1_heat_mn_f"); p0 <- rd("A0_ploidy"); pf <- rd("A1_ploidy_mn_f")
dmlh <- h0$qvalue < 0.01 & abs(h0$meth.diff) > 20
cat("heat DML (A0):", sum(dmlh), "; with corrected p < 0.05:", sum(hf$pvalue[dmlh] < 0.05), sprintf("(%.1f%%)", 100 * mean(hf$pvalue[dmlh] < 0.05)), "\n")
dmlp <- p0$qvalue < 0.01 & abs(p0$meth.diff) > 20
cat("ploidy DML (A0):", sum(dmlp), "; with corrected p < 0.05:", sum(pf$pvalue[dmlp] < 0.05), sprintf("(%.1f%%)", 100 * mean(pf$pvalue[dmlp] < 0.05)), "\n")
sigh <- hf$qvalue < 0.01
cat("corrected heat: q < 0.01 at", sum(sigh), "; mitochondrial among them:", sum(sigh & mt), "; DML:", sum(sigh & abs(hf$meth.diff) > 20), "\n")
cat("  eligible CpGs: corrected heat p < 0.05 share", round(mean(hf$pvalue[elig] < 0.05), 4), "; ploidy", round(mean(pf$pvalue[elig] < 0.05), 4), "\n")
cat("  mitochondrial CpGs with corrected heat q < 0.01:", sum(sigh & mt), "of", sum(mt), "; max |diff| among them", round(max(abs(hf$meth.diff[sigh & mt])), 2), "\n")
cat("  share of eligible CpGs that are mitochondrial:", sum(elig & mt), "\n")
cp <- which(pf$qvalue < 0.01 & abs(pf$meth.diff) > 20); ch <- which(hf$qvalue < 0.01 & abs(hf$meth.diff) > 20)
cat("corrected ploidy DML rows:", paste(chr[cp], pos[cp]), "\n"); cat("corrected heat DML:", length(ch), "; in both corrected lists:", paste(chr[intersect(cp, ch)], pos[intersect(cp, ch)]), "\n")
cat("  smallest corrected p, ploidy list:", paste(chr[cp][which.min(pf$pvalue[cp])], pos[cp][which.min(pf$pvalue[cp])]), "; heat list:", paste(chr[ch][which.min(hf$pvalue[ch])], pos[ch][which.min(hf$pvalue[ch])]), "\n")
i <- which(chr == "NC_047565.1" & pos == 35030701)
cat("NC_047565.1:35030701 library percentages:\n"); print(round(setNames(P[i, ], lib$library_id), 2))
nb <- which(chr == "NC_047565.1" & abs(pos - 35030701) <= 50)
for (j in nb) { cat(pos[j], ":", paste(round(P[j, ], 1), collapse = " "), "\n") }
# interaction, corrected test: is the corrected p the uncorrected one at the two calls?
for (m in c("legacy", "corrected")) {
  a <- rd("A2_interaction", m); b <- rd("A2_interaction_mn_f", m)
  dl <- fread(file.path(root, "06_diagnostics/dml_sensitivity", sprintf("dml_A2_interaction_mn_f_%s.csv", m)))
  cat("mode", m, ": corrected interaction list\n"); print(dl)
  uu <- if (m == "legacy") data.table(chr = chr, start = pos) else fread(file.path(root, "04_output/05_differential_methylation/tables/united_cpg_corrected.csv.gz"), select = c("chr", "start"))
  k <- match(paste(dl$chr, dl$start), paste(uu$chr, uu$start))
  print(data.table(chr = dl$chr, start = dl$start, p_uncorrected = a$pvalue[k], p_corrected = b$pvalue[k], equal = abs(a$pvalue[k] / b$pvalue[k] - 1) < 1e-8))
  cat("  q < 0.01 without the difference rule (uncorrected):", sum(a$qvalue < 0.01), "; corrected:", sum(b$qvalue < 0.01), "\n")
}
# A4 lists
a4h <- rd("A4_ploidy_in_heat"); a4c <- rd("A4_ploidy_in_control"); a4hf <- rd("A4_ploidy_in_heat_mn_f"); a4cf <- rd("A4_ploidy_in_control_mn_f")
lh <- a4h$qvalue < 0.01 & abs(a4h$meth.diff) > 20; lc <- a4c$qvalue < 0.01 & abs(a4c$meth.diff) > 20
cat("A4 heat-only DML", sum(lh), "; control-only", sum(lc), "; of the 1,038 in heat-only", sum(dmlp & lh), "; in control-only", sum(dmlp & lc), "; in both", sum(dmlp & lh & lc), "; neither", sum(dmlp & !lh & !lc), "\n")
cat("  control-only DML that are draft DML:", sum(lc & dmlp), "of", sum(lc), "; heat-only DML that are draft DML:", sum(lh & dmlp), "of", sum(lh), "\n")
cat("  A4 corrected: heat", sum(a4hf$qvalue < 0.01 & abs(a4hf$meth.diff) > 20), "; control", sum(a4cf$qvalue < 0.01 & abs(a4cf$meth.diff) > 20), "\n")
a3p <- rd("A3_ploidy_no_acute"); a3h <- rd("A3_heat_no_acute")
l3p <- a3p$qvalue < 0.01 & abs(a3p$meth.diff) > 20; l3h <- a3h$qvalue < 0.01 & abs(a3h$meth.diff) > 20
cat("A3: ploidy", sum(l3p), "of which draft DML", sum(l3p & dmlp), sprintf("(%.1f%% of 1,038)", 100 * sum(l3p & dmlp) / sum(dmlp)), "; heat", sum(l3h), "of which draft DML", sum(l3h & dmlh), sprintf("(%.1f%% of 2,351)", 100 * sum(l3h & dmlh) / sum(dmlh)), "\n")
cat("  A3 ploidy: same direction among the 591:", sum(l3p & dmlp & sign(a3p$meth.diff) == sign(p0$meth.diff)), "; hyper", sum(l3p & a3p$meth.diff > 0), "hypo", sum(l3p & a3p$meth.diff < 0), "\n")
cat("done\n")

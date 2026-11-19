# 07_alt678_elisa_hypothesis.R
# Tasks (b)6, (b)7, (b)8 and remaining number checks:
#  - ELISA group means, rank tests, values capped at the top standard, the six
#    animals with both measurements, size of the difference in the two assays
#  - frequency of "high" libraries in the two cohorts and in project zr3644
#  - acute heat shock by ploidy
#  - the draft's DML against the level slope (bullet 2 of "Why it matters")
#  - direction shares of the regroupings (is 59% "no preferred direction"?)
#  - neighbours of NC_047565.1:35,030,701
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
id <- lib$library_id
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG); names(L) <- id
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
dh <- which(lib$ploidy == "diploid" & lib$treatment == "heat"); th <- which(lib$ploidy == "triploid" & lib$treatment == "heat")
dc <- which(lib$ploidy == "diploid" & lib$treatment == "control"); tc <- which(lib$ploidy == "triploid" & lib$treatment == "control")
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"))

cat("==== (b)8 ELISA ====\n")
el <- fread(file.path(repo, "04_output/03_global_methylation/tables/elisa_percent_5mC.csv"))
el[, group := paste(ploidy, treatment)]
print(el[, .(n = .N, mean = round(mean(percent_5mC), 3), median = round(median(percent_5mC), 3), min = round(min(percent_5mC), 2), max = round(max(percent_5mC), 2),
             above_top_standard = sum(above_top_standard)), by = group])
gs <- fread(file.path(repo, "04_output/03_global_methylation/tables/group_summary.csv")); print(gs[, .(assay, ploidy, treatment, n, mean = round(mean, 3))])
sc <- fread(file.path(repo, "04_output/03_global_methylation/tables/elisa_standard_curve.csv"))
top_abs <- sc[in_fit == TRUE, max(mean_abs)]
fit <- lm(percent_5mC ~ mean_abs, el)       # percent is linear in absorbance for a fixed DNA input
cat("absorbance of the highest standard used in the fit:", round(top_abs, 4), "; percent 5-mC is linear in absorbance (R2", round(summary(fit)$r.squared, 6), "); value at that absorbance:", round(predict(fit, data.table(mean_abs = top_abs)), 3), "\n")
cat("samples above it:", el[mean_abs > top_abs, paste(animal_id, round(percent_5mC, 2), collapse = "; ")], "\n")
cap <- predict(fit, data.table(mean_abs = top_abs))
el[, capped := pmin(percent_5mC, cap)]
print(el[, .(mean_as_reported = round(mean(percent_5mC), 2), mean_capped_at_top_standard = round(mean(capped), 2)), by = group])
g <- function(p, t) el[ploidy == p & treatment == t, percent_5mC]
wt <- function(a, b) wilcox.test(a, b, exact = TRUE)$p.value
cat("exact Wilcoxon (rank tests do not depend on extrapolation): diploid heat vs triploid heat p =", round(wt(g("diploid", "heat"), g("triploid", "heat")), 4),
    "; diploid heat vs diploid control p =", round(wt(g("diploid", "heat"), g("diploid", "control")), 4),
    "; triploid heat vs triploid control p =", round(wt(g("triploid", "heat"), g("triploid", "control")), 4), "\n")
cat("diploid heat animals above every triploid heat animal:", sum(g("diploid", "heat") > max(g("triploid", "heat"))), "of 8\n")
cat("ratio of group means, heat diploid over heat triploid: ELISA", round(mean(g("diploid", "heat")) / mean(g("triploid", "heat")), 2),
    "(capped", round(el[group == "diploid heat", mean(capped)] / el[group == "triploid heat", mean(capped)], 2), "); WGBS genome-wide", round(mean(L[dh]) / mean(L[th]), 3), "\n")
cat("ratio heat diploid over control diploid: ELISA", round(mean(g("diploid", "heat")) / mean(g("diploid", "control")), 2), "; WGBS (other experiment's diploids)", round(mean(L[dh]) / mean(L[dc]), 3), "\n")
both <- data.table(animal = lib$animal_id, library = id, wgbs = L, elisa = el$percent_5mC[match(lib$animal_id, el$animal_id)], flagged = el$above_top_standard[match(lib$animal_id, el$animal_id)])[!is.na(elisa)]
print(both)
ct <- function(x, y, m) { r <- suppressWarnings(cor.test(x, y, method = m, exact = TRUE)); c(estimate = unname(r$estimate), p = r$p.value) }
cat("six animals: Spearman", round(ct(both$wgbs, both$elisa, "spearman"), 4), "; Pearson", round(ct(both$wgbs, both$elisa, "pearson"), 4), "\n")
b5 <- both[animal != "D13"]
cat("without D13 (five animals): Spearman", round(ct(b5$wgbs, b5$elisa, "spearman"), 4), "; Pearson", round(ct(b5$wgbs, b5$elisa, "pearson"), 4), "\n")
cat("six animals, means: ELISA D", round(both[grepl("^D", animal), mean(elisa)], 2), "T", round(both[grepl("^T", animal), mean(elisa)], 2),
    "; WGBS D", round(both[grepl("^D", animal), mean(wgbs)], 2), "T", round(both[grepl("^T", animal), mean(wgbs)], 2), "\n")
cat("animals with a WGBS library but no ELISA value:", paste(setdiff(lib$animal_id[heat], el$animal_id), collapse = " "), "\n")

cat("\n==== frequency of high libraries ====\n")
ft <- function(a, b, c, d) fisher.test(matrix(c(a, b, c, d), 2))$p.value
cat("heat cohort 3 of 10 against control cohort 2 of 12: Fisher exact p =", round(ft(3, 7, 2, 10), 3), "\n")
cat("heat cohort 3 of 10 against project zr3644 2 of 24 (levels above 11%): Fisher exact p =", round(ft(3, 7, 2, 22), 3), "\n")
cat("heat diploids 3 of 5 against all untreated libraries of project zr3644 2 of 24: p =", round(ft(3, 2, 2, 22), 3), "; against the 12 control libraries of the study 2 of 12: p =", round(ft(3, 2, 2, 10), 3), "\n")
cat("within the heat cohort, high by ploidy 3 of 5 against 0 of 5: p =", round(ft(3, 2, 0, 5), 3), "; within the control cohort 0 of 6 against 2 of 6: p =", round(ft(0, 6, 2, 4), 3), "\n")

cat("\n==== acute heat shock by ploidy (heat cohort) ====\n")
d <- data.table(L = L[heat], ploidy = lib$ploidy[heat], acute = lib$acute_heat_shock[heat], chh = (100 * bis$meth_CHH / (bis$meth_CHH + bis$unmeth_CHH))[heat])
print(d[, .(n = .N, mean_level = round(mean(L), 3), values = paste(round(L, 2), collapse = " "), chh = paste(round(chh, 3), collapse = " ")), by = .(ploidy, acute)])
print(round(coef(summary(lm(L ~ ploidy * acute, d))), 4))

m <- readRDS(file.path(scratch, "counts_legacy.rds")); C <- m$C; N <- m$N; chr <- m$chr; start <- m$start
pm <- 100 * C / N
mp <- rowMeans(pm); elig <- mp > 10 & mp < 90
pool <- function(idx) 100 * rowSums(C[, idx, drop = FALSE]) / rowSums(N[, idx, drop = FALSE])

cat("\n==== the draft's DML against the slope on the level ====\n")
# slope on the level using all 22 libraries, with values and level centred within cohort
lc <- L; lc[heat] <- L[heat] - mean(L[heat]); lc[ctrl] <- L[ctrl] - mean(L[ctrl])
pc <- pm; pc[, heat] <- pm[, heat] - rowMeans(pm[, heat]); pc[, ctrl] <- pm[, ctrl] - rowMeans(pm[, ctrl])
s_all <- as.vector(pc %*% lc) / sum(lc^2); rm(pc)
dec <- rep(NA_integer_, length(s_all)); dec[elig] <- cut(s_all[elig], quantile(s_all[elig], 0:10 / 10), include.lowest = TRUE, labels = FALSE)
dP <- fread(file.path(repo, "04_output/05_differential_methylation/tables/diffmeth_ploidy_legacy.csv.gz"))
dH <- fread(file.path(repo, "04_output/05_differential_methylation/tables/diffmeth_heat_legacy.csv.gz"))
stopifnot(identical(dP$start, start), identical(dH$start, start))
dmlP <- dP$qvalue < 0.01 & abs(dP$meth.diff) > 20; dmlH <- dH$qvalue < 0.01 & abs(dH$meth.diff) > 20
cat("ploidy DML:", sum(dmlP), "(eligible:", sum(dmlP & elig), "); heat DML:", sum(dmlH), "(eligible:", sum(dmlH & elig), ")\n")
dHeat <- pool(th) - pool(dh); dCtrl <- pool(tc) - pool(dc)
sets <- list("draft ploidy DML" = dmlP & elig, "draft heat DML" = dmlH & elig, "shared loci (in both lists)" = dmlP & dmlH & elig,
             "heat cohort ploidy, above 20 points" = abs(dHeat) > 20 & elig, "control cohort ploidy, above 20 points" = abs(dCtrl) > 20 & elig)
print(rbindlist(lapply(names(sets), function(nm) { k <- sets[[nm]]; data.table(set = nm, n = sum(k), share_in_top_tenth_of_slope = round(mean(dec[k] == 10), 3),
  share_in_top_three_tenths = round(mean(dec[k] >= 8), 3), share_in_bottom_tenth = round(mean(dec[k] == 1), 3)) })))
# reference: lists made by the 20-point rule under random within-cohort relabelling of ploidy
set.seed(5); ref <- t(sapply(1:200, function(b) {
  a <- c(sample(heat, 5), sample(ctrl, 6)); bb <- setdiff(1:22, a); d <- pool(bb) - pool(a); k <- abs(d) > 20 & elig
  c(n = sum(k), top = mean(dec[k] == 10), top3 = mean(dec[k] >= 8), n_high_a = sum(a %in% high)) }))
cat("200 random relabellings of ploidy within cohort: share of the list in the top tenth of slope: median", round(median(ref[, "top"]), 3), " 5th to 95th percentile", round(quantile(ref[, "top"], c(.05, .95)), 3), "\n")
cat("  by imbalance of the five high libraries between the two groups (|n in group a - 2.5|):\n")
print(as.data.table(ref)[, .(relabellings = .N, median_n = median(n), median_share_top_tenth = round(median(top), 3)), by = .(high_in_group_a = n_high_a)][order(high_in_group_a)])
cat("level difference between the groups of the draft's contrasts: ploidy (T - D)", round(mean(L[c(th, tc)]) - mean(L[c(dh, dc)]), 3), "; cohort (heat - control)", round(mean(L[heat]) - mean(L[ctrl]), 3), "\n")
cat("share of the eligible CpGs above 20 points with diploids higher in the heat cohort that lie in the top tenth / top two tenths of the control-cohort slope:\n")
lcc <- L[ctrl] - mean(L[ctrl]); sC <- as.vector((pm[, ctrl] - rowMeans(pm[, ctrl])) %*% lcc) / sum(lcc^2)
decC <- rep(NA_integer_, length(sC)); decC[elig] <- cut(sC[elig], quantile(sC[elig], 0:10 / 10), include.lowest = TRUE, labels = FALSE)
k <- elig & dHeat < -20; cat("  n =", sum(k), "; top tenth", round(mean(decC[k] == 10), 3), "; top two tenths", round(mean(decC[k] >= 9), 3), "; top three", round(mean(decC[k] >= 8), 3), "\n")

cat("\n==== direction shares of the regroupings ====\n")
maj <- function(a, b) { d <- pool(b) - pool(a); c(n = sum(abs(d) > 20), share_a_higher = sum(d < -20) / sum(abs(d) > 20)) }
seven <- setdiff(heat, high)
r7 <- t(sapply(combn(seven, 2, simplify = FALSE), function(a) maj(a, setdiff(seven, a))))
r7 <- data.table(pair = sapply(combn(seven, 2, simplify = FALSE), function(a) paste(id[a], collapse = " ")), r7)
r7[, majority := pmax(share_a_higher, 1 - share_a_higher)]
print(r7[order(-n)][, .(pair, n, share_pair_higher = round(share_a_higher, 3))])
cat("observed pair (D11, D12): share with the pair higher", round(r7[pair == "zr3534_1 zr3534_2", share_a_higher], 3), "; majority share of the other 20 regroupings: median", round(median(r7[pair != "zr3534_1 zr3534_2", majority]), 3),
    " range", round(range(r7[pair != "zr3534_1 zr3534_2", majority]), 3), "\n")
r10 <- t(sapply(combn(heat, 5, simplify = FALSE), function(a) maj(a, setdiff(heat, a))))
mj10 <- pmax(r10[, 2], 1 - r10[, 2]); obs <- which(sapply(combn(heat, 5, simplify = FALSE), function(a) setequal(a, dh)))
cat("heat cohort 5 vs 5: observed majority share", round(mj10[obs], 3), "; other regroupings: median", round(median(mj10[-obs]), 3), " 95th percentile", round(quantile(mj10[-obs], 0.95), 3), " max", round(max(mj10[-obs]), 3),
    "; regroupings with a majority share at least the observed:", sum(mj10 >= mj10[obs] - 1e-12), "of 252\n")
nh <- sapply(combn(heat, 5, simplify = FALSE), function(a) sum(a %in% high))
cat("  majority share by the number of high libraries in the first group:\n"); print(data.table(n_high_in_group = nh, majority = mj10)[, .(regroupings = .N, median_majority = round(median(majority), 3), max = round(max(majority), 3)), by = n_high_in_group][order(n_high_in_group)])

cat("\n==== NC_047565.1:35,030,701 and its neighbours among the tested CpGs ====\n")
i0 <- which(chr == "NC_047565.1" & start == 35030701)
cat("row", i0, "; percent methylation:\n"); print(round(setNames(pm[i0, ], id), 2)); cat("reads (methylated/all):", paste0(C[i0, ], "/", N[i0, ]), "\n")
nb <- which(chr == "NC_047565.1" & abs(start - 35030701) <= 2000)
x <- data.table(start = start[nb], distance = start[nb] - 35030701, DH = round(pool(dh)[nb], 1), TH = round(pool(th)[nb], 1), DC = round(pool(dc)[nb], 1), TC = round(pool(tc)[nb], 1),
                TH_range = apply(pm[nb, th, drop = FALSE], 1, function(v) paste(round(range(v)), collapse = "-")), others_range = apply(pm[nb, -th, drop = FALSE], 1, function(v) paste(round(range(v)), collapse = "-")))
print(x)

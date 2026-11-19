# 03_alt23_similarity.R
# Task (b)2 and (b)3: is the mutual similarity of the five heat-cohort triploids
# explained by sampling noise (coverage) or by the level?
#  - thinning: every library of the heat cohort is subsampled (without replacement,
#    hypergeometric) to the minimum raw coverage of the ten libraries at each CpG;
#    pairwise Pearson correlations and the rank of the triploid subset are recomputed.
#  - noise-corrected (disattenuated) correlations, no random draw.
#  - pairwise r against the coverage of the two libraries.
#  - (b)3: the seven libraries that are not high; correlations of residuals after
#    regressing every CpG on the level.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
id <- lib$library_id
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
dh <- which(lib$ploidy == "diploid" & lib$treatment == "heat"); th <- which(lib$ploidy == "triploid" & lib$treatment == "heat")
dc <- which(lib$ploidy == "diploid" & lib$treatment == "control"); tc <- which(lib$ploidy == "triploid" & lib$treatment == "control")
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"))
raw <- readRDS(file.path(scratch, "raw_counts_legacy.rds"))
CR <- raw$CR; NR <- raw$NR; chr <- raw$chr
pm <- 100 * CR / NR
mean_pm <- rowMeans(pm)
elig <- mean_pm > 10 & mean_pm < 90
R0 <- cor(pm); dimnames(R0) <- list(id, id)
step04 <- as.matrix(fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_correlation_10x_legacy.csv"))[, -1])
cat("max |correlation from reconstructed raw counts - step 04|:", signif(max(abs(R0 - step04)), 3), "\n")

mw <- function(M, idx) { s <- M[idx, idx]; mean(s[upper.tri(s)]) }
summarise <- function(M, members, own, label) {
  subs <- combn(members, length(own), simplify = FALSE)
  v <- vapply(subs, function(a) mw(M, a), numeric(1))
  is_own <- vapply(subs, function(a) setequal(a, own), TRUE)
  pr <- t(combn(members, 2)); r <- M[pr]
  own_pair <- pr[, 1] %in% own & pr[, 2] %in% own
  top <- order(-r)[seq_len(sum(own_pair))]
  data.table(what = label, mean_own = round(mw(M, own), 4), mean_other_group = round(mw(M, setdiff(members, own)), 4),
             mean_between = round(mean(M[own, setdiff(members, own)]), 4),
             rank_of_own_subset = rank(-v)[is_own], n_subsets = length(v), best_other_subset = round(max(v[!is_own]), 4),
             own_pairs_among_top = sum(own_pair[top]), n_own_pairs = sum(own_pair),
             min_own_pair = round(min(r[own_pair]), 4), max_other_pair = round(max(r[!own_pair]), 4))
}

cat("\n==== observed ====\n")
print(summarise(R0, heat, th, "heat cohort, triploids, observed"))

cat("\n==== coverage of the heat libraries at the tested CpGs (reconstructed raw counts) ====\n")
covt <- colMeans(NR); cove <- colMeans(NR[elig, ]); medt <- apply(NR, 2, median)
print(data.table(id = id[heat], ploidy = lib$ploidy[heat], mean_cov_tested = round(covt[heat], 2), median_cov_tested = medt[heat], mean_cov_eligible = round(cove[heat], 2),
                 mean_inverse_cov = round(colMeans(1 / NR[, heat]), 5)))
cat("mean raw coverage at tested CpGs: diploids", round(mean(covt[dh]), 2), " triploids", round(mean(covt[th]), 2), "\n")
cat("mean raw coverage at eligible CpGs: diploids", round(mean(cove[dh]), 2), " triploids", round(mean(cove[th]), 2), "\n")
pr <- t(combn(heat, 2))
pair <- data.table(a = id[pr[, 1]], b = id[pr[, 2]], r = R0[pr], TT = pr[, 1] %in% th & pr[, 2] %in% th, DD = pr[, 1] %in% dh & pr[, 2] %in% dh,
                   mean_cov = (covt[pr[, 1]] + covt[pr[, 2]]) / 2, harm_cov = 2 / (1 / covt[pr[, 1]] + 1 / covt[pr[, 2]]),
                   mean_inv = (colMeans(1 / NR)[pr[, 1]] + colMeans(1 / NR)[pr[, 2]]) / 2,
                   level_gap = abs(L[pr[, 1]] - L[pr[, 2]]))
cat("45 pairs of the heat cohort: Spearman of r with mean coverage", round(cor(pair$r, pair$mean_cov, method = "spearman"), 3),
    "(p", signif(cor.test(pair$r, pair$mean_cov, method = "spearman", exact = FALSE)$p.value, 2), "); with harmonic-mean coverage", round(cor(pair$r, pair$harm_cov, method = "spearman"), 3),
    "; with mean inverse coverage", round(cor(pair$r, pair$mean_inv, method = "spearman"), 3),
    "; with the difference in level", round(cor(pair$r, pair$level_gap, method = "spearman"), 3), "\n")
cat("among the 35 pairs that are not triploid-triploid: Spearman r vs mean coverage", round(pair[TT == FALSE, cor(r, mean_cov, method = "spearman")], 3),
    "; r vs level difference", round(pair[TT == FALSE, cor(r, level_gap, method = "spearman")], 3), "\n")
cat("mean coverage of pairs: T-T", round(pair[TT == TRUE, mean(mean_cov)], 2), " D-D", round(pair[DD == TRUE, mean(mean_cov)], 2), " D-T", round(pair[!TT & !DD, mean(mean_cov)], 2), "\n")
print(summary(lm(r ~ mean_inv + level_gap + TT, pair))$coefficients)

cat("\n==== noise-corrected correlations (no random draw) ====\n")
# var(observed percent) = var(true) + mean binomial sampling variance; cov is unaffected by independent sampling.
# Sampling variance of 100*c/n estimated without bias by 100^2 * p(1-p)/(n-1).
disatt <- function(cols, rows = rep(TRUE, nrow(pm))) {
  x <- pm[rows, cols]; n <- NR[rows, cols]
  p <- x / 100
  noise <- colMeans(1e4 * p * (1 - p) / (n - 1))
  S <- cov(x)
  vt <- diag(S) - noise
  M <- S / sqrt(outer(vt, vt)); diag(M) <- 1
  list(M = M, share_noise = noise / diag(S))
}
da <- disatt(1:22)
Rd <- da$M; dimnames(Rd) <- list(id, id)
cat("share of each library's variance across CpGs that is sampling noise (heat cohort):\n"); print(round(setNames(da$share_noise[heat], id[heat]), 4))
print(rbind(summarise(Rd, heat, th, "heat cohort, triploids, noise-corrected, all tested CpGs"),
            summarise(Rd, heat, dh, "heat cohort, diploids, noise-corrected, all tested CpGs")))
cat("noise-corrected correlations, heat cohort:\n"); print(round(Rd[heat, heat], 4))
dae <- disatt(1:22, elig); Rde <- dae$M; dimnames(Rde) <- list(id, id)
R0e <- cor(pm[elig, ]); dimnames(R0e) <- list(id, id)
cat("share of variance that is sampling noise, eligible CpGs only (heat cohort):\n"); print(round(setNames(dae$share_noise[heat], id[heat]), 4))
print(rbind(summarise(R0e, heat, th, "heat cohort, triploids, eligible CpGs, observed"),
            summarise(Rde, heat, th, "heat cohort, triploids, eligible CpGs, noise-corrected")))

cat("\n==== thinning to the minimum coverage of the ten heat libraries at each CpG ====\n")
nmin <- do.call(pmin, as.data.frame(NR[, heat]))
cat("common coverage: median", median(nmin), " mean", round(mean(nmin), 2), " range", range(nmin), "\n")
thin_rows <- list(); thin_mats <- list()
for (seed in 1:8) {
  set.seed(1000 + seed)
  X <- sapply(heat, function(j) 100 * rhyper(nrow(NR), CR[, j], NR[, j] - CR[, j], nmin) / nmin)
  M <- matrix(NA_real_, 22, 22); M[heat, heat] <- cor(X)
  Me <- matrix(NA_real_, 22, 22); Me[heat, heat] <- cor(X[elig, ])
  thin_rows[[length(thin_rows) + 1]] <- cbind(seed = seed, rbind(summarise(M, heat, th, "thinned, all tested CpGs"), summarise(Me, heat, th, "thinned, eligible CpGs")))
  thin_mats[[seed]] <- M[heat, heat]
}
thin_tab <- rbindlist(thin_rows); print(thin_tab)
fwrite(thin_tab, file.path(scratch, "03_thinning_heat.csv"))
Mavg <- Reduce(`+`, thin_mats) / length(thin_mats); dimnames(Mavg) <- list(id[heat], id[heat])
cat("mean thinned correlation matrix over the 8 seeds:\n"); print(round(Mavg, 4))
# thinning to a fixed coverage of 10 reads everywhere (strongest equalization)
cat("\n-- thinning every heat library to exactly 10 reads at every CpG --\n")
rows10 <- list()
for (seed in 1:4) {
  set.seed(2000 + seed)
  X <- sapply(heat, function(j) 100 * rhyper(nrow(NR), CR[, j], NR[, j] - CR[, j], 10L) / 10)
  M <- matrix(NA_real_, 22, 22); M[heat, heat] <- cor(X)
  rows10[[seed]] <- cbind(seed = seed, summarise(M, heat, th, "thinned to 10 reads, all tested CpGs"))
}
print(rbindlist(rows10))

cat("\n==== (b)3 the seven heat libraries that are not high ====\n")
seven <- setdiff(heat, high)
pr7 <- t(combn(seven, 2))
p7 <- data.table(a = id[pr7[, 1]], b = id[pr7[, 2]], r_observed = round(R0[pr7], 4), r_noise_corrected = round(Rd[pr7], 4),
                 type = ifelse(pr7[, 1] %in% th & pr7[, 2] %in% th, "T-T", ifelse(pr7[, 1] %in% dh & pr7[, 2] %in% dh, "D-D", "D-T")))[order(-r_observed)]
p7[, rank := seq_len(.N)]; print(p7)
print(rbind(summarise(R0, seven, th, "seven not high, triploids, observed"), summarise(Rd, seven, th, "seven not high, triploids, noise-corrected")))
cat("D11-D12 correlation:", round(R0[1, 2], 4), "; D11 with the five triploids:", round(R0[1, th], 4), "; D12 with the five triploids:", round(R0[2, th], 4), "\n")
cat("rank of the D11-D12 pair among the 21 pairs of the seven:", p7[a == "zr3534_1" & b == "zr3534_2", rank], "; among the 45 pairs of the cohort:", which(pair[order(-r), paste(a, b)] == "zr3534_1 zr3534_2"), "\n")
cat("T11 (zr3534_6) with the other four triploids:", round(R0[6, setdiff(th, 6)], 4), "\n")

cat("\n==== (b)3 correlations after removing the level from every CpG (heat cohort) ====\n")
resid_cor <- function(cols, lev, rows = rep(TRUE, nrow(pm))) {
  x <- pm[rows, cols]
  l <- lev - mean(lev)
  xc <- x - rowMeans(x)
  b <- as.vector(xc %*% l) / sum(l^2)
  res <- xc - outer(b, l)
  cor(res)
}
Rr <- matrix(NA_real_, 22, 22); Rr[heat, heat] <- resid_cor(heat, L[heat])
cat("(residual correlations are centred near -1/(n-2); only their ordering is read)\n")
print(rbind(summarise(Rr, heat, th, "heat cohort, residuals on the genome-wide level, all tested CpGs"),
            summarise(Rr, heat, dh, "heat cohort diploids, residuals on the genome-wide level")))
me <- colMeans(pm[elig, ])
Rr2 <- matrix(NA_real_, 22, 22); Rr2[heat, heat] <- resid_cor(heat, me[heat])
print(summarise(Rr2, heat, th, "heat cohort, residuals on the eligible-CpG mean, all tested CpGs"))
x <- Rr[heat, heat]; dimnames(x) <- list(id[heat], id[heat]); print(round(x, 3))
# centring only (no level): correlations of deviations from the cohort mean
Rc <- matrix(NA_real_, 22, 22); Rc[heat, heat] <- cor(pm[, heat] - rowMeans(pm[, heat]))
print(summarise(Rc, heat, th, "heat cohort, deviations from the CpG mean of the cohort (no level removed)"))

cat("\n==== by chromosome: is the triploid similarity everywhere? ====\n")
big <- names(which(table(chr) > 5000))
bychr <- rbindlist(lapply(big, function(cc) { M <- matrix(NA_real_, 22, 22); M[heat, heat] <- cor(pm[chr == cc, heat]); cbind(chr = cc, n_cpg = sum(chr == cc), summarise(M, heat, th, "")) }))
print(bychr[, .(chr, n_cpg, mean_own, mean_other_group, mean_between, rank_of_own_subset, own_pairs_among_top)])

cat("\n==== control cohort for comparison ====\n")
print(rbind(summarise(R0, ctrl, tc, "control cohort, triploids, observed"), summarise(Rd, ctrl, tc, "control cohort, triploids, noise-corrected"),
            summarise(R0, ctrl, dc, "control cohort, diploids, observed"), summarise(Rd, ctrl, dc, "control cohort, diploids, noise-corrected")))
ten <- setdiff(ctrl, high)
print(rbind(summarise(R0, ten, setdiff(tc, high), "control cohort without its two high, triploids (4), observed"),
            summarise(Rd, ten, setdiff(tc, high), "control cohort without its two high, triploids (4), noise-corrected")))
cat("noise-corrected correlations, control cohort:\n"); print(round(Rd[ctrl, ctrl], 3))
cat("share of variance that is sampling noise (control cohort):\n"); print(round(setNames(da$share_noise[ctrl], id[ctrl]), 4))
saveRDS(list(R0 = R0, Rd = Rd), file.path(scratch, "03_correlations.rds"))

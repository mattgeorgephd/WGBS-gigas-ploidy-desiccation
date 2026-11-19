# 04_cross_cohort_similarity.R
# Follow-up to (b)2/(b)3: are the five heat-cohort triploids similar only to each
# other, or also to the libraries of the control cohort?  If a cross or a
# processing batch made them similar, they should not be closer to the control
# libraries than the heat-cohort diploids are.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 230)
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
raw <- readRDS(file.path(scratch, "raw_counts_legacy.rds")); CR <- raw$CR; NR <- raw$NR
pm <- 100 * CR / NR
step04 <- as.matrix(fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_correlation_10x_legacy.csv"))[, -1]); dimnames(step04) <- list(id, id)
cs <- readRDS(file.path(scratch, "03_correlations.rds")); Rd <- cs$Rd

short <- sub("zr3534_", "H", sub("zr3644_", "C", id))
lab <- paste0(short, ifelse(lib$ploidy == "diploid", "d", "t"), ifelse(seq_along(id) %in% high, "*", ""))

cat("==== step 04 correlations (observed): each heat library with the control libraries ====\n")
x <- step04[heat, ctrl]; dimnames(x) <- list(lab[heat], lab[ctrl]); print(round(x, 3))
ctrl_nh <- setdiff(ctrl, high)
tab <- data.table(library = id[heat], ploidy = lib$ploidy[heat], high = heat %in% high, level = round(L[heat], 2),
                  mean_r_with_12_control = round(rowMeans(step04[heat, ctrl]), 4),
                  mean_r_with_10_control_not_high = round(rowMeans(step04[heat, ctrl_nh]), 4),
                  noise_corrected_with_10_control_not_high = round(rowMeans(Rd[heat, ctrl_nh]), 4),
                  mean_r_with_other_heat = round(sapply(heat, function(j) mean(step04[j, setdiff(heat, j)])), 4))
print(tab)
cat("for how many of the 10 control libraries that are not high is each heat triploid more similar than D11 and than D12 (step 04 values)?\n")
for (t in th) cat(" ", id[t], ": above D11 at", sum(step04[t, ctrl_nh] > step04[1, ctrl_nh]), "of 10; above D12 at", sum(step04[t, ctrl_nh] > step04[2, ctrl_nh]), "of 10\n")
cat("and for the two high control libraries (zr3644_14, zr3644_15)?\n")
for (t in th) cat(" ", id[t], ": above D11 at", sum(step04[t, intersect(ctrl, high)] > step04[1, intersect(ctrl, high)]), "of 2; above D12 at", sum(step04[t, intersect(ctrl, high)] > step04[2, intersect(ctrl, high)]), "of 2\n")

cat("\n==== thinning all 22 libraries to the minimum coverage of the 22 at each CpG ====\n")
nmin <- do.call(pmin, as.data.frame(NR))
cat("common coverage: median", median(nmin), " mean", round(mean(nmin), 2), "\n")
acc <- matrix(0, 22, 22); nseed <- 5
for (seed in 1:nseed) {
  set.seed(3000 + seed)
  X <- sapply(1:22, function(j) 100 * rhyper(nrow(NR), CR[, j], NR[, j] - CR[, j], nmin) / nmin)
  M <- cor(X); acc <- acc + M
  cat(sprintf("seed %d: mean r heat T-T %.4f; heat T with 10 control not high %.4f; D11,D12 with 10 control not high %.4f; D11-D12 %.4f; D11/D12 with heat T %.4f; among 10 control not high %.4f\n", seed,
              mean(M[th, th][upper.tri(M[th, th])]), mean(M[th, ctrl_nh]), mean(M[1:2, ctrl_nh]), M[1, 2], mean(M[1:2, th]), mean(M[ctrl_nh, ctrl_nh][upper.tri(M[ctrl_nh, ctrl_nh])])))
}
Mt <- acc / nseed; dimnames(Mt) <- list(lab, lab)
cat("mean thinned correlation matrix (5 seeds), all 22 libraries (H = heat cohort, C = control cohort, d/t = ploidy, * = high):\n")
print(round(Mt, 3))
fwrite(data.table(library = lab, round(Mt, 4)), file.path(scratch, "04_thinned_correlation_22.csv"))
tab2 <- data.table(library = lab, cohort = lib$cohort,
                   mean_r_with_heat_triploids = round(sapply(1:22, function(j) mean(Mt[j, setdiff(th, j)])), 4),
                   mean_r_with_control_not_high = round(sapply(1:22, function(j) mean(Mt[j, setdiff(ctrl_nh, j)])), 4),
                   mean_r_with_all_other_21 = round(sapply(1:22, function(j) mean(Mt[j, -j])), 4))
print(tab2[order(-mean_r_with_all_other_21)])
cat("heat triploids more similar (thinned) to each control library that is not high than D11 / D12 are? counts of 10:\n")
for (t in th) cat(" ", id[t], ":", sum(Mt[t, ctrl_nh] > Mt[1, ctrl_nh]), "/", sum(Mt[t, ctrl_nh] > Mt[2, ctrl_nh]), "\n")

cat("\n==== average-linkage clustering of the 22 libraries on 1 - r (thinned to equal coverage) ====\n")
hc <- hclust(as.dist(1 - Mt), method = "average")
for (k in 2:7) { cl <- cutree(hc, k); cat("k =", k, ":", paste(tapply(lab, cl, paste, collapse = " "), collapse = "  |  "), "\n") }
cat("merge order (heights):\n"); print(data.table(step = seq_len(nrow(hc$merge)), height = round(hc$height, 4),
  merged = apply(hc$merge, 1, function(m) paste(ifelse(m < 0, lab[abs(m)], paste0("[", m, "]")), collapse = " + "))))

cat("\n==== the same on the noise-corrected correlations (no thinning) ====\n")
dimnames(Rd) <- list(lab, lab)
hc2 <- hclust(as.dist(1 - Rd), method = "average")
for (k in 2:7) { cl <- cutree(hc2, k); cat("k =", k, ":", paste(tapply(lab, cl, paste, collapse = " "), collapse = "  |  "), "\n") }
cat("noise-corrected correlations, heat libraries (rows) with control libraries (columns):\n"); print(round(Rd[heat, ctrl], 3))
cat("mean noise-corrected r: heat T-T", round(mean(Rd[th, th][upper.tri(Rd[th, th])]), 4), "; heat T with the 10 control not high", round(mean(Rd[th, ctrl_nh]), 4),
    "; D11 and D12 with the 10 control not high", round(mean(Rd[1:2, ctrl_nh]), 4), "; D11 and D12 with heat T", round(mean(Rd[1:2, th]), 4), "; D11-D12", round(Rd[1, 2], 4),
    "; among the 10 control not high", round(mean(Rd[ctrl_nh, ctrl_nh][upper.tri(Rd[ctrl_nh, ctrl_nh])]), 4), "\n")

cat("\n==== (b)3 level removed with slopes estimated in the control cohort ====\n")
# slope of each CpG on the genome-wide level among the 12 control libraries; subtract slope * (level - mean level of
# the heat cohort) from the heat libraries; correlations of the adjusted values.
lc <- L[ctrl] - mean(L[ctrl])
s_ctrl <- as.vector((pm[, ctrl] - rowMeans(pm[, ctrl])) %*% lc) / sum(lc^2)
adj <- pm[, heat] - outer(s_ctrl, L[heat] - mean(L[heat]))
Ma <- cor(adj); dimnames(Ma) <- list(lab[heat], lab[heat])
mw <- function(M, idx) { s <- M[idx, idx]; mean(s[upper.tri(s)]) }
subs <- combn(10, 5, simplify = FALSE); v <- vapply(subs, function(a) mw(Ma, a), numeric(1))
cat("adjusted with control-cohort slopes: mean r T-T", round(mw(Ma, 6:10), 4), "; D-D", round(mw(Ma, 1:5), 4), "; between", round(mean(Ma[1:5, 6:10]), 4),
    "; rank of the triploid subset among 252:", rank(-v)[vapply(subs, function(a) setequal(a, 6:10), TRUE)], "\n")
print(round(Ma, 3))
pr <- t(combn(10, 2)); r <- Ma[pr]; tt <- pr[, 1] > 5 & pr[, 2] > 5
cat("T-T pairs among the 10 most similar after adjustment:", sum(tt[order(-r)[1:10]]), "\n")

# 06_alt5_slope_permutation.R
# Task (b)5: the slope analysis of item 4.  Can the association between a CpG's
# slope on the level (estimated in one cohort) and the ploidy difference of the
# other cohort arise without a level effect, for example because CpGs with a
# high variance among animals have both large slopes and large differences?
# Test: permute the level values among the libraries of the cohort in which the
# slope is estimated and recompute the decile table.
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
m <- readRDS(file.path(scratch, "counts_legacy.rds")); C <- m$C; N <- m$N
pm <- 100 * C / N
elig <- { mp <- rowMeans(pm); mp > 10 & mp < 90 }
pool <- function(idx) 100 * rowSums(C[, idx, drop = FALSE]) / rowSums(N[, idx, drop = FALSE])
dHeat <- (pool(th) - pool(dh))[elig]; dCtrl <- (pool(tc) - pool(dc))[elig]
dHeat2 <- (pool(th) - pool(setdiff(dh, high)))[elig]        # D11, D12 against the five triploids
dCtrl2 <- (pool(setdiff(tc, high)) - pool(dc))[elig]
pe <- pm[elig, ]
rm(C, N, pm, m); invisible(gc())

slope <- function(idx, lev) { l <- lev - mean(lev); as.vector((pe[, idx] - rowMeans(pe[, idx])) %*% l) / sum(l^2) }
dec_tab <- function(s, d) {
  dec <- cut(s, quantile(s, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  data.table(decile = dec, d = d)[, .(n = .N, triploid_higher = mean(d > 20), diploid_higher = mean(d < -20), median_d = median(d)), by = decile][order(decile)]
}
top_bottom <- function(s, d) {
  q <- quantile(s, c(0.1, 0.9)); top <- s >= q[2]; bot <- s <= q[1]
  c(top_dip = mean(d[top] < -20), top_trip = mean(d[top] > 20), bot_dip = mean(d[bot] < -20), bot_trip = mean(d[bot] > 20), spearman = cor(s, d, method = "spearman"))
}

cat("==== observed: slope in the control cohort, ploidy difference of the heat cohort ====\n")
sC <- slope(ctrl, L[ctrl])
print(dec_tab(sC, dHeat)[, .(decile, n, triploid_higher = round(triploid_higher, 4), diploid_higher = round(diploid_higher, 4), median_d = round(median_d, 2))])
obsC <- top_bottom(sC, dHeat); print(round(obsC, 4))
cat("==== observed: slope in the heat cohort, ploidy difference of the control cohort ====\n")
sH <- slope(heat, L[heat])
print(dec_tab(sH, dCtrl)[, .(decile, n, triploid_higher = round(triploid_higher, 4), diploid_higher = round(diploid_higher, 4), median_d = round(median_d, 2))])
obsH <- top_bottom(sH, dCtrl); print(round(obsH, 4))

cat("\n==== permutation of the level among the 12 control libraries (500 random permutations) ====\n")
set.seed(20261119)
B <- 500
permC <- t(sapply(1:B, function(b) top_bottom(slope(ctrl, sample(L[ctrl])), dHeat)))
sm <- function(x) c(min = min(x), q05 = unname(quantile(x, 0.05)), median = median(x), q95 = unname(quantile(x, 0.95)), max = max(x))
cat("top tenth, share with diploids higher by more than 20 points (observed", round(obsC["top_dip"], 4), "):\n"); print(round(sm(permC[, "top_dip"]), 4))
cat("top tenth, share with triploids higher (observed", round(obsC["top_trip"], 4), "):\n"); print(round(sm(permC[, "top_trip"]), 4))
cat("bottom tenth, share with diploids higher (observed", round(obsC["bot_dip"], 4), "):\n"); print(round(sm(permC[, "bot_dip"]), 4))
cat("bottom tenth, share with triploids higher (observed", round(obsC["bot_trip"], 4), "):\n"); print(round(sm(permC[, "bot_trip"]), 4))
cat("Spearman slope vs difference (observed", round(obsC["spearman"], 4), "):\n"); print(round(sm(permC[, "spearman"]), 4))
cat("permutations with a top-tenth diploid-higher share at least the observed:", sum(permC[, "top_dip"] >= obsC["top_dip"]), "of", B, "\n")
cat("permutations with |Spearman| at least the observed:", sum(abs(permC[, "spearman"]) >= abs(obsC["spearman"])), "of", B, "\n")
cat("asymmetry top minus bottom tenth in the diploid-higher share: observed", round(obsC["top_dip"] - obsC["bot_dip"], 4), "; permutations:\n")
print(round(sm(permC[, "top_dip"] - permC[, "bot_dip"]), 4))
cat("the larger of the two extreme tenths (variance effect alone): diploid-higher share, permutations:\n"); print(round(sm(pmax(permC[, "top_dip"], permC[, "bot_dip"])), 4))
cat("diploid-higher share among all eligible CpGs (no selection):", round(mean(dHeat < -20), 4), "; triploid-higher:", round(mean(dHeat > 20), 4), "\n")
fwrite(as.data.table(permC), file.path(scratch, "06_perm_control_slope.csv"))

cat("\n==== exhaustive: every pair of control libraries declared 'high' (slope = pair mean minus mean of the other ten) ====\n")
prs <- combn(ctrl, 2, simplify = FALSE)
ex <- rbindlist(lapply(prs, function(p) { s <- rowMeans(pe[, p]) - rowMeans(pe[, setdiff(ctrl, p)]); as.data.table(c(list(pair = paste(id[p], collapse = " ")), as.list(top_bottom(s, dHeat)))) }))
setorder(ex, -top_dip); print(ex[, .(pair, top_dip = round(top_dip, 4), top_trip = round(top_trip, 4), bot_dip = round(bot_dip, 4), bot_trip = round(bot_trip, 4), spearman = round(spearman, 3))][c(1:8, 60:66)])
cat("rank of the pair zr3644_14 zr3644_15 by top-tenth diploid-higher share:", which(ex$pair == "zr3644_14 zr3644_15"), "of 66\n")
cat("pairs without zr3644_14 and zr3644_15 (45 pairs): top-tenth diploid-higher share range", round(range(ex[!grepl("zr3644_14|zr3644_15", pair), top_dip]), 4),
    "; median", round(median(ex[!grepl("zr3644_14|zr3644_15", pair), top_dip]), 4), "\n")

cat("\n==== permutation of the level among the 10 heat libraries (500), predicting the control-cohort difference ====\n")
permH <- t(sapply(1:B, function(b) top_bottom(slope(heat, sample(L[heat])), dCtrl)))
cat("top tenth, share with triploids higher (observed", round(obsH["top_trip"], 4), "):\n"); print(round(sm(permH[, "top_trip"]), 4))
cat("top tenth, share with diploids higher (observed", round(obsH["top_dip"], 4), "):\n"); print(round(sm(permH[, "top_dip"]), 4))
cat("bottom tenth, share with triploids higher (observed", round(obsH["bot_trip"], 4), "):\n"); print(round(sm(permH[, "bot_trip"]), 4))
cat("Spearman (observed", round(obsH["spearman"], 4), "):\n"); print(round(sm(permH[, "spearman"]), 4))
cat("permutations with a top-tenth triploid-higher share at least the observed:", sum(permH[, "top_trip"] >= obsH["top_trip"]), "of", B, "\n")
cat("permutations with |Spearman| at least the observed:", sum(abs(permH[, "spearman"]) >= abs(obsH["spearman"])), "of", B, "\n")

cat("\n==== is the level a continuous axis?  slope among the libraries that are NOT high ====\n")
c10 <- setdiff(ctrl, high); h7 <- setdiff(heat, high)
sC10 <- slope(c10, L[c10])
cat("slope among the 10 control libraries that are not high (level 9.46 to 10.41) against the heat-cohort difference (all ten heat libraries):\n")
print(round(top_bottom(sC10, dHeat), 4))
cat("  correlation of this slope with the slope from all 12 control libraries:", round(cor(sC10, sC), 3), "\n")
set.seed(7); p10 <- t(sapply(1:300, function(b) top_bottom(slope(c10, sample(L[c10])), dHeat)))
cat("  permutation reference (300) for the top-tenth diploid-higher share:", round(sm(p10[, "top_dip"]), 4), "\n")
cat("  permutation reference for Spearman:", round(sm(p10[, "spearman"]), 4), "\n")
sH7 <- slope(h7, L[h7])
cat("slope among the 7 heat libraries that are not high against the control-cohort difference:\n"); print(round(top_bottom(sH7, dCtrl), 4))

cat("\n==== does the control-cohort slope also predict the difference between D11, D12 and the five triploids? ====\n")
print(dec_tab(sC, dHeat2)[, .(decile, n, triploid_higher = round(triploid_higher, 4), diploid_higher = round(diploid_higher, 4), median_d = round(median_d, 2))])
obs2 <- top_bottom(sC, dHeat2); print(round(obs2, 4))
set.seed(11); p2 <- t(sapply(1:300, function(b) top_bottom(slope(ctrl, sample(L[ctrl])), dHeat2)))
cat("permutation reference (300): top-tenth diploid-higher share", round(sm(p2[, "top_dip"]), 4), "\n   Spearman", round(sm(p2[, "spearman"]), 4), "\n")
cat("level difference D11, D12 minus five triploids:", round(mean(L[setdiff(dh, high)]) - mean(L[th]), 3), "points\n")
cat("\n==== and the heat-cohort slope against the control difference without the two high control libraries ====\n")
print(round(top_bottom(sH, dCtrl2), 4))
cat("level difference four control triploids minus six control diploids:", round(mean(L[setdiff(tc, high)]) - mean(L[dc]), 3), "points\n")

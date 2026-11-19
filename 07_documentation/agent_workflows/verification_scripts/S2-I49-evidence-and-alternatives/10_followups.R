# 10_followups.R
# (1) Permutation of the level among the control libraries, classified by where the largest level values land.
# (2) Figure 1B interaction without the five, also with the values before deduplication.
# (3) Corrected mode: the two new observations (cross-cohort similarity; draft lists without the five) are not recomputed here.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); id <- lib$library_id
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)
Lpre <- 100 * bis$meth_CpG_predup / (bis$meth_CpG_predup + bis$unmeth_CpG_predup)
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
dh <- which(lib$ploidy == "diploid" & lib$treatment == "heat"); th <- which(lib$ploidy == "triploid" & lib$treatment == "heat")
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"))

cat("==== Figure 1B interaction with and without the five ====\n")
for (nm in c("after deduplication", "before deduplication")) {
  v <- if (nm == "after deduplication") L else Lpre
  d <- data.table(v = v, ploidy = factor(lib$ploidy), treatment = factor(lib$treatment))
  a <- anova(lm(v ~ ploidy * treatment, d)); b <- anova(lm(v ~ ploidy * treatment, d[-high]))
  est <- function(x) coef(lm(v ~ ploidy * treatment, x))[4]
  cat(sprintf("%s: 22 libraries F(1,%d) = %.2f, p = %.4f, contrast %.3f; 17 libraries F(1,%d) = %.2f, p = %.4f, contrast %.3f\n", nm, a$Df[4], a$`F value`[3], a$`Pr(>F)`[3], est(d),
              b$Df[4], b$`F value`[3], b$`Pr(>F)`[3], est(d[-high])))
}
# robustness: leave out, instead of the five, the single highest library of each group, or the highest and lowest of each group
d <- data.table(v = L, g = paste(lib$ploidy, lib$treatment), ploidy = factor(lib$ploidy), treatment = factor(lib$treatment), id = id)
trim <- d[, .SD[order(v)][-c(1, .N)], by = g]
a <- anova(lm(v ~ ploidy * treatment, trim)); cat("each group without its highest and its lowest library (14 libraries): F(1,", a$Df[4], ") =", round(a$`F value`[3], 2), ", p =", signif(a$`Pr(>F)`[3], 3), "\n")
med <- d[, .(median = round(median(v), 3)), by = g]; print(med)
cat("interaction contrast in group medians (22 libraries):", round((med[g == "triploid heat", median] - med[g == "diploid heat", median]) - (med[g == "triploid control", median] - med[g == "diploid control", median]), 3), "\n")

cat("\n==== permutation of the level among the 12 control libraries, by where the two largest values land ====\n")
m <- readRDS(file.path(scratch, "counts_legacy.rds")); C <- m$C; N <- m$N
pm <- 100 * C / N; elig <- { mp <- rowMeans(pm); mp > 10 & mp < 90 }
pool <- function(idx) 100 * rowSums(C[, idx, drop = FALSE]) / rowSums(N[, idx, drop = FALSE])
dHeat <- (pool(th) - pool(dh))[elig]; pe <- pm[elig, ctrl]; rm(C, N, pm, m); invisible(gc())
Lc <- L[ctrl]; i14 <- which(id[ctrl] == "zr3644_14"); i15 <- which(id[ctrl] == "zr3644_15")
set.seed(20261119); B <- 500
res <- rbindlist(lapply(1:B, function(b) {
  lp <- sample(Lc); l <- lp - mean(lp); s <- as.vector((pe - rowMeans(pe)) %*% l) / sum(l^2)
  q <- quantile(s, c(0.1, 0.9))
  data.table(top_dip = mean(dHeat[s >= q[2]] < -20), top_trip = mean(dHeat[s >= q[2]] > 20), bot_dip = mean(dHeat[s <= q[1]] < -20),
             rank14 = rank(-lp)[i14], rank15 = rank(-lp)[i15]) }))
res[, class := fifelse(rank14 <= 2, "zr3644_14 has one of the two largest values", fifelse(rank14 >= 11, "zr3644_14 has one of the two smallest values", "zr3644_14 in between"))]
print(res[, .(permutations = .N, top_tenth_diploid_higher_median = round(median(top_dip), 3), min = round(min(top_dip), 3), max = round(max(top_dip), 3),
              bottom_tenth_diploid_higher_median = round(median(bot_dip), 3), top_tenth_triploid_higher_median = round(median(top_trip), 3)), by = class])
cat("observed: top tenth 0.465 diploids higher, 0.011 triploids higher; bottom tenth 0.102 diploids higher\n")
cat("permutations with top-tenth diploid-higher share above 0.30:", res[top_dip > 0.30, .N], "; of these with zr3644_14 among the two largest values:", res[top_dip > 0.30 & rank14 <= 2, .N],
    "; with zr3644_14 or zr3644_15 the largest:", res[top_dip > 0.30 & (rank14 == 1 | rank15 == 1), .N], "\n")
cat("Spearman between the level rank given to zr3644_14 and the top-tenth share:", round(cor(res$rank14, res$top_dip, method = "spearman"), 3), "\n")

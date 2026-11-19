# 02_alt1_artifacts.R
# Task (b)1: is the "level" an artifact of depth, duplication, mapping
# efficiency or coverage at the tested CpGs?  Also: mitochondrial CpGs as an
# internal control of non-conversion in CpG context, and the profile of the
# level across classes of CpGs (proportional or confined to intermediate CpGs).
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
covtab <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
filt <- fread(file.path(repo, "04_output/05_differential_methylation/tables/library_filtering_legacy.csv"))
stopifnot(identical(lib$library_id, bis$library_id), identical(lib$library_id, covtab$library_id), identical(lib$library_id, filt$library_id))
id <- lib$library_id
heat <- which(lib$treatment == "heat"); ctrl <- which(lib$treatment == "control")
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)
Lpre <- 100 * bis$meth_CpG_predup / (bis$meth_CpG_predup + bis$unmeth_CpG_predup)
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"))

# ---- raw counts from the united (normalized) table --------------------------------------
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cv <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE]); storage.mode(cv) <- "integer"
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); storage.mode(C) <- "integer"
Tm <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); storage.mode(Tm) <- "integer"
chr <- u$chr; start <- u$start
rm(u); invisible(gc())
sf <- filt$scaling_factor
covR <- matrix(0L, nrow(cv), 22); CR <- covR
bad <- integer(22)
for (j in 1:22) {
  cr <- as.integer(round(cv[, j] / sf[j]))
  # numCs = round(new_cov * nC / cov)  =>  nC is the integer for which this holds; start from the proportional guess
  c0 <- as.integer(round(C[, j] * cr / cv[, j]))
  chk_cov <- round(sf[j] * cr)
  chk_c <- round(chk_cov * (c0 / cr)); chk_t <- round(chk_cov * ((cr - c0) / cr))
  bad[j] <- sum(chk_cov != cv[, j] | chk_c != C[, j] | chk_t != Tm[, j])
  covR[, j] <- cr; CR[, j] <- c0
}
cat("cells where the reconstructed raw counts do not reproduce the normalized counts, per library:\n"); print(setNames(bad, id))
cat("minimum reconstructed raw coverage:", min(covR), "\n")
saveRDS(list(chr = chr, start = start, CR = CR, NR = covR), file.path(scratch, "raw_counts_legacy.rds"), compress = FALSE)

N <- C + Tm
pm <- 100 * C / N
mean_pm <- rowMeans(pm)
elig <- mean_pm > 10 & mean_pm < 90
me <- colMeans(pm[elig, ])

cat("\n==== (b)1 library-level quantities ====\n")
cov_tested <- colMeans(covR); cov_elig <- colMeans(covR[elig, ])
med_tested <- apply(covR, 2, median)
methylated <- mean_pm >= 90; unmeth <- mean_pm <= 10
cov_ratio <- colMeans(covR[methylated, ]) / colMeans(covR[unmeth, ])
cov_ratio_elig <- colMeans(covR[elig, ]) / colMeans(covR[unmeth, ])
q <- data.table(id, cohort = lib$cohort, high = seq_along(id) %in% high, level = round(L, 3), level_predup = round(Lpre, 3), mean_elig = round(me, 2),
                read_pairs_M = round(bis$read_pairs / 1e6, 1), unique_pairs_M = round(bis$unique_pairs / 1e6, 1), map_eff = round(bis$mapping_efficiency_pct, 2),
                dup_pct = round(bis$duplication_pct, 2), total_C_M = round(bis$total_C / 1e6), n10x_M = round(covtab$n_cpg_10x / 1e6, 2),
                mean_cov_genome = round(covtab$mean_coverage, 2), cov_tested = round(cov_tested, 2), cov_elig = round(cov_elig, 2),
                cov_meth_over_unmeth = round(cov_ratio, 3), cov_elig_over_unmeth = round(cov_ratio_elig, 3),
                ambig_pct = round(100 * bis$ambiguous_pairs / bis$read_pairs, 2), unaligned_pct = round(100 * bis$unaligned_pairs / bis$read_pairs, 2))
print(q)
vars <- c("level_predup", "read_pairs_M", "unique_pairs_M", "map_eff", "dup_pct", "total_C_M", "n10x_M", "mean_cov_genome", "cov_tested", "cov_elig",
          "cov_meth_over_unmeth", "cov_elig_over_unmeth", "ambig_pct", "unaligned_pct")
full <- data.table(level = L, level_predup = Lpre, mean_elig = me, read_pairs_M = bis$read_pairs, unique_pairs_M = bis$unique_pairs, map_eff = bis$mapping_efficiency_pct,
                   dup_pct = bis$duplication_pct, total_C_M = bis$total_C, n10x_M = covtab$n_cpg_10x, mean_cov_genome = covtab$mean_coverage,
                   cov_tested = cov_tested, cov_elig = cov_elig, cov_meth_over_unmeth = cov_ratio, cov_elig_over_unmeth = cov_ratio_elig,
                   ambig_pct = bis$ambiguous_pairs / bis$read_pairs, unaligned_pct = bis$unaligned_pairs / bis$read_pairs)
res <- rbindlist(lapply(c("level", "mean_elig"), function(y) rbindlist(lapply(vars, function(v) {
  f <- function(idx, m) cor(full[[y]][idx], full[[v]][idx], method = m)
  p <- function(idx) suppressWarnings(cor.test(full[[y]][idx], full[[v]][idx], method = "spearman", exact = TRUE)$p.value)
  data.table(y = y, x = v, spearman_heat = round(f(heat, "spearman"), 2), p_heat = signif(p(heat), 2), spearman_control = round(f(ctrl, "spearman"), 2), p_control = signif(p(ctrl), 2),
             pearson_heat = round(f(heat, "pearson"), 2), pearson_control = round(f(ctrl, "pearson"), 2),
             pearson_17_not_high_heat = round(f(setdiff(heat, high), "pearson"), 2), pearson_not_high_control = round(f(setdiff(ctrl, high), "pearson"), 2))
}))))
print(res)
fwrite(res, file.path(scratch, "02_level_vs_qc.csv")); fwrite(q, file.path(scratch, "02_library_qc.csv"))
cat("level before minus after deduplication: range", round(range(Lpre - L), 3), "; cor(level, level_predup) =", round(cor(L, Lpre), 4),
    "; the five highest before deduplication:", id[order(-Lpre)][1:5], "\n")
cat("high vs not high, Wilcoxon exact p within cohort, for each QC quantity:\n")
for (v in vars[-1]) cat(sprintf("  %-22s heat (3 vs 7) p = %.3f   control (2 vs 10) p = %.3f\n", v,
                                wilcox.test(full[[v]][intersect(heat, high)], full[[v]][setdiff(heat, high)], exact = TRUE)$p.value,
                                wilcox.test(full[[v]][intersect(ctrl, high)], full[[v]][setdiff(ctrl, high)], exact = TRUE)$p.value))

cat("\n==== mitochondrial CpGs (NC_001276.1) among the tested CpGs ====\n")
mt <- chr == "NC_001276.1"
cat("tested CpGs on the mitochondrial genome:", sum(mt), "\n")
mito <- 100 * colSums(CR[mt, ]) / colSums(covR[mt, ])
print(data.table(id, high = seq_along(id) %in% high, mito_pct_mCpG = round(mito, 4), mito_reads_M = round(colSums(covR[mt, ]) / 1e6, 3), level = round(L, 2)))
cat("mito CpG methylation: range", round(range(mito), 3), "; five high:", round(mito[high], 3), "\n")
cat("Spearman mito vs level: heat", round(cor(mito[heat], L[heat], method = "spearman"), 2), " control", round(cor(mito[ctrl], L[ctrl], method = "spearman"), 2), "\n")
cat("mean mito, high minus not high: heat", round(mean(mito[intersect(heat, high)]) - mean(mito[setdiff(heat, high)]), 4),
    " control", round(mean(mito[intersect(ctrl, high)]) - mean(mito[setdiff(ctrl, high)]), 4), "\n")

cat("\n==== profile of the level across classes of CpGs ====\n")
# Classes are defined in the OTHER cohort (mean percent methylation of its libraries that are not high), so that
# the libraries compared do not define the classes.
brk <- c(-0.001, 1, 10, 30, 50, 70, 90, 99, 100)
prof <- function(target, definer, label) {
  cls <- cut(rowMeans(pm[, definer]), brk)
  tab <- data.table(class = levels(cls), n_cpg = as.vector(table(cls)))
  hi <- intersect(target, high); lo <- setdiff(target, high)
  tab[, mean_high := sapply(levels(cls), function(l) mean(pm[cls == l, hi]))]
  tab[, mean_not_high := sapply(levels(cls), function(l) mean(pm[cls == l, lo]))]
  tab[, difference := mean_high - mean_not_high]
  tab[, ratio := mean_high / mean_not_high]
  tab[, range_high := sapply(levels(cls), function(l) paste(round(range(colMeans(pm[cls == l, hi, drop = FALSE])), 2), collapse = " to "))]
  tab[, range_not_high := sapply(levels(cls), function(l) paste(round(range(colMeans(pm[cls == l, lo, drop = FALSE])), 2), collapse = " to "))]
  cat(label, "\n"); print(tab[, .(class, n_cpg, mean_high = round(mean_high, 2), mean_not_high = round(mean_not_high, 2), difference = round(difference, 2), ratio = round(ratio, 3), range_high, range_not_high)])
}
prof(heat, setdiff(ctrl, high), "heat cohort (3 high vs 7 others); classes by the mean of the 10 control libraries that are not high")
prof(ctrl, setdiff(heat, high), "control cohort (2 high vs 10 others); classes by the mean of the 7 heat libraries that are not high")

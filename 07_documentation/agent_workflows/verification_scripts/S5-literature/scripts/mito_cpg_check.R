# Mitochondrial CpGs of the united table (legacy mode): what the apparent methylation looks like.
# Purpose: the literature attributes apparent mtDNA methylation in bisulfite data to incomplete
# conversion of intact (circular, supercoiled) mtDNA. That predicts an excess spread over all
# mitochondrial CpGs; reads of another origin, a variant or site-specific methylation predict an
# excess at particular positions. This script describes which of the two the 363 tested CpGs show.
# Input : mito/united_legacy_mito.csv (header + the NC_001276.1 rows of united_cpg_legacy.csv.gz,
#         extracted with zcat | grep; normalized counts as methylKit tests them)
#         03_input_files/metadata/wgbs_libraries.csv (run from the repository root)
# Output: printed tables; mito/mito_per_library.csv, mito/mito_per_cpg.csv
# Read-only on the repository. One thread.
suppressMessages(library(data.table)); setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
scratch <- args[1]
u <- fread(file.path(scratch, "mito", "united_legacy_mito.csv"))
meta <- fread("03_input_files/metadata/wgbs_libraries.csv")
stopifnot(nrow(u) == 363, nrow(meta) == 22)
cov <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE])
nC  <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE])
heat <- meta$cohort == "heat_experiment"
lib <- data.table(library = meta$library_id, cohort = meta$cohort, ploidy = meta$ploidy,
                  reads_at_mito_cpgs = colSums(cov), methylated_reads = colSums(nC),
                  pooled_pct = 100 * colSums(nC) / colSums(cov),
                  mean_of_cpg_pct = colMeans(100 * nC / cov),
                  median_of_cpg_pct = apply(100 * nC / cov, 2, median),
                  cpgs_with_zero_methylated_reads = colSums(nC == 0),
                  max_cpg_pct = apply(100 * nC / cov, 2, max),
                  mean_cov = colMeans(cov))
print(lib[, .(library, cohort, ploidy, reads_at_mito_cpgs, pooled_pct = round(pooled_pct, 3),
              mean_of_cpg_pct = round(mean_of_cpg_pct, 3), median_of_cpg_pct = round(median_of_cpg_pct, 3),
              zero = cpgs_with_zero_methylated_reads, max_cpg_pct = round(max_cpg_pct, 1), mean_cov = round(mean_cov))])
cat("\npooled percent, range by cohort:\n")
print(lib[, .(min = round(min(pooled_pct), 2), max = round(max(pooled_pct), 2),
              mean_cov_min = round(min(mean_cov)), mean_cov_max = round(max(mean_cov))), by = cohort])
cat("ratio of cohort means of the pooled percent (control / heat):",
    round(mean(lib[cohort != "heat_experiment", pooled_pct]) / mean(lib[cohort == "heat_experiment", pooled_pct]), 2), "\n")
cat("ratio of cohort medians:",
    round(median(lib[cohort != "heat_experiment", pooled_pct]) / median(lib[cohort == "heat_experiment", pooled_pct]), 2), "\n")

# per CpG, pooled over the libraries of each cohort
cp <- data.table(pos = u$start,
                 heat_cov = rowSums(cov[, heat]), heat_pct = 100 * rowSums(nC[, heat]) / rowSums(cov[, heat]),
                 ctrl_cov = rowSums(cov[, !heat]), ctrl_pct = 100 * rowSums(nC[, !heat]) / rowSums(cov[, !heat]))
cp[, diff := ctrl_pct - heat_pct]
cat("\nper-CpG pooled percent by cohort (363 CpGs):\n")
print(rbind(heat = quantile(cp$heat_pct, c(0, .05, .25, .5, .75, .95, 1)),
            control = quantile(cp$ctrl_pct, c(0, .05, .25, .5, .75, .95, 1))), digits = 3)
cat("CpGs with a higher pooled percent in the control cohort:", sum(cp$diff > 0), "of", nrow(cp), "\n")
cat("CpGs with pooled percent above 10 in the control cohort:", sum(cp$ctrl_pct > 10), "; in the heat cohort:", sum(cp$heat_pct > 10), "\n")
cat("CpGs with pooled percent above 5 in the control cohort:", sum(cp$ctrl_pct > 5), "; in the heat cohort:", sum(cp$heat_pct > 5), "\n")
# how much of the cohort difference in methylated reads comes from the top CpGs
exc <- cp$ctrl_pct - cp$heat_pct
o <- order(-exc)
cat("share of the summed per-CpG excess (control minus heat) carried by the 10 CpGs with the largest excess:",
    round(100 * sum(exc[o][1:10]) / sum(exc), 1), "% ; by the 36 largest (10% of CpGs):",
    round(100 * sum(exc[o][1:36]) / sum(exc), 1), "%\n")
cat("pooled percent without the 36 CpGs of largest excess: heat",
    round(100 * sum(nC[-o[1:36], heat]) / sum(cov[-o[1:36], heat]), 3), " control",
    round(100 * sum(nC[-o[1:36], !heat]) / sum(cov[-o[1:36], !heat]), 3), "\n")
# dispersion of per-CpG values against binomial sampling at a common rate, per cohort
for (g in c("heat", "control")) {
  idx <- if (g == "heat") heat else !heat
  k <- rowSums(nC[, idx]); n <- rowSums(cov[, idx]); p <- sum(k) / sum(n)
  x2 <- sum((k - n * p)^2 / (n * p * (1 - p)))
  cat(g, ": Pearson chi-squared of per-CpG counts against one common rate =", round(x2), "on", length(k) - 1,
      "df ; ratio", round(x2 / (length(k) - 1), 2), "\n")
}
cat("correlation across CpGs of the per-CpG pooled percent, heat against control: Pearson",
    round(cor(cp$heat_pct, cp$ctrl_pct), 3), " Spearman", round(cor(cp$heat_pct, cp$ctrl_pct, method = "spearman"), 3), "\n")
cat("correlation across CpGs of pooled percent with coverage (Mechta et al. report an inverse relation): heat",
    round(cor(cp$heat_pct, cp$heat_cov, method = "spearman"), 3), " control",
    round(cor(cp$ctrl_pct, cp$ctrl_cov, method = "spearman"), 3), "\n")
cat("\nthe 12 CpGs with the largest control-cohort percent:\n")
print(cp[order(-ctrl_pct)][1:12, .(pos, heat_cov = round(heat_cov), heat_pct = round(heat_pct, 2),
                                    ctrl_cov = round(ctrl_cov), ctrl_pct = round(ctrl_pct, 2))])
# along the genome: thirds of the 363 CpGs by position
cp[, third := cut(seq_len(.N), 3, labels = c("first third", "middle third", "last third"))]
print(cp[, .(n = .N, from = min(pos), to = max(pos),
             heat_pct = round(100 * sum(heat_pct * heat_cov / 100) / sum(heat_cov), 3),
             ctrl_pct = round(100 * sum(ctrl_pct * ctrl_cov / 100) / sum(ctrl_cov), 3)), by = third])
# coverage of the mitochondrial CpGs relative to the scaled median library coverage is not available here;
# report the mean normalized coverage at mitochondrial CpGs per cohort instead
cat("\nmean normalized coverage per mitochondrial CpG and library: heat",
    round(mean(cov[, heat])), " control", round(mean(cov[, !heat])), "\n")
fwrite(lib, file.path(scratch, "mito", "mito_per_library.csv"))
fwrite(cp, file.path(scratch, "mito", "mito_per_cpg.csv"))

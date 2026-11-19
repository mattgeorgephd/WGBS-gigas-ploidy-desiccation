suppressPackageStartupMessages(library(data.table))
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
e <- cpg[elig_p == TRUE]
e[, covq := cut(mean_cov, quantile(mean_cov, 0:4 / 4), include.lowest = TRUE)]
cat("eligible CpGs:", nrow(e), "\n")
show <- function(d, lab) {
  cat("\n", lab, " n =", nrow(d), "\n")
  print(d[, .(n = .N, median_meth = round(median(mean_ct), 1), p_ploidy_pct = round(100 * mean(dml_ploidy), 2), p_heat_pct = round(100 * mean(dml_heat), 2)), keyby = covq])
}
show(e[in_CDS == TRUE & in_TE == FALSE], "CDS, not TE")
show(e[in_CDS == TRUE & in_TE == FALSE & mean_ct >= 40 & mean_ct <= 70], "CDS, not TE, methylation 40-70%")
show(e[in_intron == TRUE & in_CDS == FALSE & in_TE == FALSE & mean_ct >= 40 & mean_ct <= 70], "intron, not TE, methylation 40-70%")
show(e[in_gene == FALSE & mean_ct >= 40 & mean_ct <= 70], "outside genes, methylation 40-70%")
show(e[in_TE == TRUE & mean_ct >= 40 & mean_ct <= 70], "TE, methylation 40-70%")
cat("\nlocus NC_047563.1:39,633,323-39,635,520: mean coverage of its tested CpGs (median, range):\n")
l <- cpg[chr == "NC_047563.1" & start >= 39633322 & start <= 39635520]
print(l[, .(n = .N, median_mean_cov = median(mean_cov), min = min(mean_cov), max = max(mean_cov), median_meth = median(mean_ct))])
cat("genome-wide tested CpGs: quantiles of mean coverage\n"); print(quantile(cpg$mean_cov, c(0, .25, .5, .75, .9, .99, .999, 1)))
cat("\ncomposition of the eligible CpGs outside introns (the comparison group of the intron test):\n")
x <- e[in_intron == FALSE]
cat("n", nrow(x), " CDS", round(100 * mean(x$in_CDS), 1), "%  exon_UTR", round(100 * mean(x$in_exon_UTR), 1), "%  intergenic", round(100 * mean(x$in_intergenic), 1),
    "%  flanks", round(100 * mean(x$in_flanks), 1), "%\n")
cat("composition of the eligible CpGs outside CDS: n", nrow(e[in_CDS == FALSE]), " intron", round(100 * mean(e[in_CDS == FALSE]$in_intron), 1), "%\n")
cat("\nDML neighbours: share of DML with another DML of the same list within 1 kb\n")
for (s in c("dml_ploidy", "dml_heat")) {
  d <- cpg[get(s) == TRUE][order(chr, start)]
  gap_prev <- c(Inf, diff(d$start)); gap_prev[c(TRUE, d$chr[-1] != d$chr[-nrow(d)])] <- Inf
  gap_next <- c(gap_prev[-1], Inf)
  reg <- cumsum(gap_prev > 1000)
  cat(s, ": n", nrow(d), " with neighbour within 1 kb:", sum(pmin(gap_prev, gap_next) <= 1000), sprintf("(%.1f%%)", 100 * mean(pmin(gap_prev, gap_next) <= 1000)),
      " regions:", max(reg), " largest region:", max(table(reg)), "\n")
}

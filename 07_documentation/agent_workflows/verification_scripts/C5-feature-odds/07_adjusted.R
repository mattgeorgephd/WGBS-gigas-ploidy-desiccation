# Is the binary "eligible" background enough? P(DML) by methylation level and coverage,
# adjusted odds ratios, and cluster-aware intervals for the stratified estimates.
suppressPackageStartupMessages({library(data.table); library(splines); library(sandwich)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
set.seed(4107)
cpg[, elig := elig_p]                       # the pipeline's eligible set (65,480)
loc <- cpg$chr == "NC_047563.1" & cpg$start >= 39633322 & cpg$start <= 39635520
dsets <- c(ploidy = "dml_ploidy", heat = "dml_heat")
feats <- c("gene", "intron", "CDS", "TE", "intergenic")

cat("== P(DML) among eligible CpGs by quartile of mean coverage ==\n")
e <- cpg[elig | dml_ploidy | dml_heat]
e[, covq := cut(mean_cov, quantile(mean_cov, 0:4 / 4), include.lowest = TRUE)]
print(e[, .(n = .N, median_cov = median(mean_cov), p_ploidy_pct = round(100 * mean(dml_ploidy), 2), p_heat_pct = round(100 * mean(dml_heat), 2),
            median_meth = round(median(mean_ct), 1), pct_gene = round(100 * mean(in_gene), 1), pct_CDS = round(100 * mean(in_CDS), 1),
            pct_TE = round(100 * mean(in_TE), 1), pct_intergenic = round(100 * mean(in_intergenic), 1)), keyby = covq])
cat("\nsame, within methylation 40-60% only (holds the methylation level roughly fixed):\n")
print(e[mean_ct >= 40 & mean_ct <= 60, .(n = .N, p_ploidy_pct = round(100 * mean(dml_ploidy), 2), p_heat_pct = round(100 * mean(dml_heat), 2),
            pct_gene = round(100 * mean(in_gene), 1), pct_TE = round(100 * mean(in_TE), 1)), keyby = covq])

cat("\n== spread of P(DML) inside the eligible band (10-point bins of mean methylation) ==\n")
e[, m10 := cut(mean_ct, seq(0, 100, 10), include.lowest = TRUE)]
print(e[elig == TRUE, .(n = .N, p_ploidy_pct = round(100 * mean(dml_ploidy), 2), p_heat_pct = round(100 * mean(dml_heat), 2),
            pct_gene = round(100 * mean(in_gene), 1), pct_intron = round(100 * mean(in_intron), 1), pct_CDS = round(100 * mean(in_CDS), 1),
            pct_TE = round(100 * mean(in_TE), 1), pct_intergenic = round(100 * mean(in_intergenic), 1)), keyby = m10])

cat("\n== logistic models among eligible CpGs; 95% intervals: naive / clustered on 100 kb windows ==\n")
e[, win := paste(chr, start %/% 100000L)]
e[, inloc := chr == "NC_047563.1" & start >= 39633322 & start <= 39635520]
fitrow <- function(d, form, lab) {
  m <- glm(form, binomial, d)
  b <- coef(m)["ff"]; sn <- sqrt(diag(vcov(m))["ff"]); sc <- sqrt(diag(vcovCL(m, cluster = d$win, type = "HC0"))["ff"])
  sprintf("%s %.2f (%.2f-%.2f / %.2f-%.2f)", lab, exp(b), exp(b - 1.96 * sn), exp(b + 1.96 * sn), exp(b - 1.96 * sc), exp(b + 1.96 * sc))
}
for (s in names(dsets)) for (f in feats) {
  d <- e[elig | get(dsets[[s]])]
  d[, yy := as.integer(get(dsets[[s]]))][, ff := as.integer(get(paste0("in_", f)))]
  cat(sprintf("%-6s %-10s | ", s, f),
      fitrow(d, yy ~ ff, "crude"), "|",
      fitrow(d, yy ~ ff + ns(mean_ct, 6), "+meth"), "|",
      fitrow(d, yy ~ ff + ns(mean_ct, 6) + log(mean_cov), "+meth+cov"), "|",
      fitrow(d[inloc == FALSE], yy ~ ff + ns(mean_ct, 6), "+meth, no locus"), "|",
      fitrow(d[inloc == FALSE], yy ~ ff + ns(mean_ct, 6) + log(mean_cov), "+meth+cov, no locus"), "\n")
}

cat("\n== mutually adjusted model: DML ~ CDS + exon_UTR + intron + TE + intergenic-or-flank + ns(meth, 6); clustered 95% intervals ==\n")
for (s in names(dsets)) {
  d <- e[elig | get(dsets[[s]])]
  d[, yy := as.integer(get(dsets[[s]]))]
  d[, nongenic := as.integer(!in_gene)]
  for (drop_locus in c(FALSE, TRUE)) {
    dd <- if (drop_locus) d[inloc == FALSE] else d
    m <- glm(yy ~ in_CDS + in_exon_UTR + nongenic + in_TE + ns(mean_ct, 6), binomial, dd)
    V <- vcovCL(m, cluster = dd$win, type = "HC0")
    k <- c("in_CDSTRUE", "in_exon_UTRTRUE", "nongenic", "in_TETRUE")
    cat(s, if (drop_locus) "(without the NC_047563.1 locus)" else "(all)", "; reference = genic, not CDS, not UTR (mostly intron)\n")
    for (j in k) cat(sprintf("   %-16s OR %.2f (%.2f-%.2f)\n", j, exp(coef(m)[j]), exp(coef(m)[j] - 1.96 * sqrt(V[j, j])), exp(coef(m)[j] + 1.96 * sqrt(V[j, j]))))
  }
}

cat("\n== Mantel-Haenszel OR stratified by 5-point methylation bins, 100 kb block bootstrap (1000 resamples) ==\n")
mh_boot <- function(d, y, fc, B = 1000) {
  d <- d[, .(y = get(y), f = get(fc), k = as.integer(cut(mean_ct, seq(0, 100, 5), include.lowest = TRUE)), win)]
  blocks <- unique(d$win); nb <- length(blocks); d[, bi := match(win, blocks)]
  W <- matrix(0L, B, nb)
  for (i in seq_len(B)) W[i, ] <- tabulate(sample.int(nb, nb, replace = TRUE), nb)
  num <- den <- numeric(B); num0 <- den0 <- 0
  for (kk in sort(unique(d$k))) {
    x <- d[k == kk, .(a = sum(y & f), b = sum(y & !f), c = sum(!y & f), d = sum(!y & !f)), by = bi]
    M <- matrix(0, nb, 4); M[x$bi, ] <- as.matrix(x[, .(a, b, c, d)])
    t0 <- colSums(M); n0 <- sum(t0)
    if (n0 > 0) { num0 <- num0 + t0[1] * t0[4] / n0; den0 <- den0 + t0[2] * t0[3] / n0 }
    Tm <- W %*% M; n <- rowSums(Tm); ok <- n > 0
    num[ok] <- num[ok] + Tm[ok, 1] * Tm[ok, 4] / n[ok]; den[ok] <- den[ok] + Tm[ok, 2] * Tm[ok, 3] / n[ok]
  }
  c(or = num0 / den0, quantile(num / den, c(0.025, 0.975), na.rm = TRUE))
}
for (s in names(dsets)) for (f in feats) {
  d <- e[elig | get(dsets[[s]])]
  r <- mh_boot(d, dsets[[s]], paste0("in_", f))
  r2 <- mh_boot(d[inloc == FALSE], dsets[[s]], paste0("in_", f))
  cat(sprintf("%-6s %-10s MH OR %.2f (bootstrap %.2f-%.2f) ; without locus %.2f (%.2f-%.2f)\n", s, f, r[1], r[2], r[3], r2[1], r2[2], r2[3]))
}

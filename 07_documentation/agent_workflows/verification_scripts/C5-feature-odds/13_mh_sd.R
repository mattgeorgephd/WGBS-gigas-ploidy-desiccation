suppressPackageStartupMessages(library(data.table))
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
pm <- 100 * nc / (nc + nt); rm(u, nc, nt); invisible(gc())
grp <- list(1:5, 6:10, 11:16, 17:22)
rv <- function(m) { k <- ncol(m); (rowSums(m^2) - k * rowMeans(m)^2) / (k - 1) }
v <- sapply(grp, function(i) rv(pm[, i, drop = FALSE])); df <- sapply(grp, length) - 1
cpg[, sd_within := sqrt(as.vector(v %*% df) / sum(df))]
set.seed(77)
e <- cpg[elig_p == TRUE | dml_ploidy | dml_heat]
e[, m5 := cut(mean_ct, seq(0, 100, 5), include.lowest = TRUE)]
e[, sd10 := cut(sd_within, quantile(sd_within, 0:10 / 10), include.lowest = TRUE)]
e[, win := paste(chr, start %/% 100000L)]
mh_boot <- function(d, y, fc, B = 1000) {
  d <- d[, .(y = get(y), f = get(fc), k = as.integer(interaction(m5, sd10, drop = TRUE)), win)]
  blocks <- unique(d$win); nb <- length(blocks); d[, bi := match(win, blocks)]
  W <- matrix(0L, B, nb); for (i in seq_len(B)) W[i, ] <- tabulate(sample.int(nb, nb, replace = TRUE), nb)
  num <- den <- numeric(B); num0 <- den0 <- 0
  for (kk in sort(unique(d$k))) {
    x <- d[k == kk, .(a = sum(y & f), b = sum(y & !f), c = sum(!y & f), d = sum(!y & !f)), by = bi]
    M <- matrix(0, nb, 4); M[x$bi, ] <- as.matrix(x[, .(a, b, c, d)])
    t0 <- colSums(M); n0 <- sum(t0); num0 <- num0 + t0[1] * t0[4] / n0; den0 <- den0 + t0[2] * t0[3] / n0
    Tm <- W %*% M; n <- rowSums(Tm); ok <- n > 0
    num[ok] <- num[ok] + Tm[ok, 1] * Tm[ok, 4] / n[ok]; den[ok] <- den[ok] + Tm[ok, 2] * Tm[ok, 3] / n[ok]
  }
  c(num0 / den0, quantile(num / den, c(0.025, 0.975), na.rm = TRUE))
}
cat("Mantel-Haenszel OR among eligible CpGs, strata = 5-point methylation bin x decile of within-group SD; 100 kb block bootstrap\n")
for (s in c("dml_ploidy", "dml_heat")) for (f in c("gene", "intron", "CDS", "TE", "intergenic")) {
  d <- e[elig_p == TRUE | get(s)]
  r <- mh_boot(d, s, paste0("in_", f))
  cat(sprintf("%-10s %-10s %.2f (%.2f-%.2f)\n", s, f, r[1], r[2], r[3]))
}
cat("\ncorrelation (Spearman) of within-group SD with mean coverage among eligible CpGs:", round(cor(e$sd_within, e$mean_cov, method = "spearman"), 3), "\n")
cat("median within-group SD of DML vs other eligible CpGs: ploidy", round(median(e$sd_within[e$dml_ploidy]), 1), "vs", round(median(e$sd_within[!e$dml_ploidy]), 1),
    "; heat", round(median(e$sd_within[e$dml_heat]), 1), "vs", round(median(e$sd_within[!e$dml_heat]), 1), "\n")
cat("expected binomial SD (points) at the median coverage and methylation of eligible CpGs in the 40-60% band:",
    round(100 * sqrt(0.25 / median(e$mean_cov[e$mean_ct >= 40 & e$mean_ct <= 60])), 1), "; observed median within-group SD there:",
    round(median(e$sd_within[e$mean_ct >= 40 & e$mean_ct <= 60]), 1), "\n")

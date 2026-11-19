# Task (c)1: is 1.18 the right comparison value, and by how much is the uncorrected statistic inflated?
# Independent simulation on my own sample of 3,000 eligible CpGs (legacy), own seed.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
d <- readRDS(file.path(S, "sample3000_legacy.rds")); lib <- d$lib
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))[d$picked]
stopifnot(identical(u$start, d$pos$start))
cov <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE]) + 0; Cs <- d$Cs; Ts <- d$Ts
filt <- fread(file.path(R, "04_output/05_differential_methylation/tables/library_filtering_legacy.csv"))
sf <- max(filt$median_coverage) / filt$median_coverage; stopifnot(max(abs(sf - filt$scaling_factor)) < 1e-9)
cat("scaling factors: range", range(sf), "; mean", mean(sf), "; mean heat cohort", mean(sf[1:10]), "; mean control cohort", mean(sf[11:22]), "\n")
raw <- round(sweep(cov, 2, sf, "/")); stopifnot(all(round(sweep(raw, 2, sf, "*")) == cov))
cat("share of sampled CpG-library cells where numCs + numTs differs from the coverage column:", mean(Cs + Ts != cov), "\n")
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat"); X <- cbind(1, trip, heat); X0 <- cbind(1, heat)
fam <- binomial(link = logit)
one <- function(y, w) { o <- suppressWarnings(glm.fit(X, y / w, weights = w, family = fam)); mu <- o$fitted.values
  o0 <- suppressWarnings(glm.fit(X0, y / w, weights = w, family = fam))
  c(phi = sum((y - w * mu)^2 / (w * mu * (1 - mu))) / 19, dev = o0$deviance - o$deviance, mu) }
obs <- t(vapply(1:nrow(Cs), function(i) one(Cs[i, ], Cs[i, ] + Ts[i, ]), numeric(24)))
mu <- obs[, 3:24]
q <- function(x) paste(round(quantile(x, c(.25, .5, .75)), 3), collapse = " / ")
cat("observed dispersion, quartiles:", q(obs[, "phi"]), "; mean", round(mean(obs[, "phi"]), 3), "\n")
set.seed(1212)
sim_rows <- function(scaled) {
  rc <- matrix(rbinom(length(raw), raw, mu), nrow = nrow(raw)); rt <- raw - rc
  if (scaled) { nc <- round(sweep(raw, 2, sf, "*")); sc <- round(nc * (rc / raw)); st <- round(nc * (rt / raw)) } else { sc <- rc; st <- rt }
  t(vapply(1:nrow(sc), function(i) one(sc[i, ], sc[i, ] + st[i, ])[1:2], numeric(2)))
}
for (rep in 1:3) {
  s1 <- sim_rows(TRUE); s0 <- sim_rows(FALSE)
  cat(sprintf("simulation %d: binomial, scaled and rounded: dispersion quartiles %s, mean %.3f; median uncorrected statistic over chi-squared median %.3f; share p < 0.05, 0.01, 0.001 uncorrected: %s\n",
              rep, q(s1[, "phi"]), mean(s1[, "phi"]), median(s1[, "dev"]) / qchisq(.5, 1),
              paste(round(100 * sapply(c(.05, .01, .001), function(a) mean(pchisq(s1[, "dev"], 1, lower.tail = FALSE) < a)), 2), collapse = ", ")))
  cat(sprintf("              binomial, not scaled:        dispersion quartiles %s, mean %.3f; median uncorrected statistic over chi-squared median %.3f; share p < 0.05, 0.01, 0.001 uncorrected: %s\n",
              q(s0[, "phi"]), mean(s0[, "phi"]), median(s0[, "dev"]) / qchisq(.5, 1),
              paste(round(100 * sapply(c(.05, .01, .001), function(a) mean(pchisq(s0[, "dev"], 1, lower.tail = FALSE) < a)), 2), collapse = ", ")))
}
cat("\nObserved counts: median uncorrected statistic over the chi-squared(1) median:", round(median(obs[, "dev"]) / qchisq(.5, 1), 3),
    "; mean statistic (chi-squared mean is 1):", round(mean(obs[, "dev"]), 3), "\n")
# the same under relabelling within cohort (20 relabellings), so that a real ploidy effect cannot be the cause
blocks <- split(1:22, lib$cohort); set.seed(1213)
infl <- vapply(1:20, function(r) { z <- trip; for (b in blocks) z[b] <- trip[b][sample.int(length(b))]; Xz <- cbind(1, z, heat)
  dev <- vapply(1:nrow(Cs), function(i) { y <- Cs[i, ]; w <- y + Ts[i, ]
    suppressWarnings(glm.fit(X0, y / w, weights = w, family = fam))$deviance - suppressWarnings(glm.fit(Xz, y / w, weights = w, family = fam))$deviance }, numeric(1))
  c(median(dev) / qchisq(.5, 1), mean(dev)) }, numeric(2))
cat("20 relabellings: median uncorrected statistic over the chi-squared median: mean", round(mean(infl[1, ]), 3), ", range", round(range(infl[1, ]), 3),
    "; mean statistic: mean", round(mean(infl[2, ]), 3), ", range", round(range(infl[2, ]), 3), "\n")
# does the dispersion rise with coverage (as variation among libraries does, and scaling does not)?
mcov <- rowMeans(Cs + Ts)
cat("Spearman correlation of the dispersion estimate with the mean normalized coverage of the CpG:", round(cor(obs[, "phi"], mcov, method = "spearman"), 3), "\n")
print(round(tapply(obs[, "phi"], cut(mcov, quantile(mcov, 0:4 / 4), include.lowest = TRUE), median), 2))
cat("mean normalized coverage per library at the sampled CpGs: median", median(mcov), "; quartiles", quantile(mcov, c(.25, .75)), "\n")
cat("F(1,19) tail beyond the chi-squared value needed for a single call (31.4176):", pf(31.4176351493246, 1, 19, lower.tail = FALSE),
    "; ratio to 2.08e-8:", pf(31.4176351493246, 1, 19, lower.tail = FALSE) / 2.08082334017924e-08, "\n")
cat("F(1,19) tail beyond the chi-squared quantiles of 0.05, 0.01, 0.001:", sapply(c(.05, .01, .001), function(a) pf(qchisq(a, 1, lower.tail = FALSE), 1, 19, lower.tail = FALSE)), "\n")
cat("t equivalents of the single-call F: 11 v 11:", sqrt(84.1119983932854), "-> standardized difference", sqrt(84.1119983932854) * sqrt(2 / 11),
    "; 6 v 6:", sqrt(250.60665524694), "->", sqrt(250.60665524694) * sqrt(2 / 6), "; 5 v 5:", sqrt(474.470701156092), "->", sqrt(474.470701156092) * sqrt(2 / 5),
    "; 10 v 12:", sqrt(84.1119983932854) * sqrt(1 / 10 + 1 / 12), "\n")

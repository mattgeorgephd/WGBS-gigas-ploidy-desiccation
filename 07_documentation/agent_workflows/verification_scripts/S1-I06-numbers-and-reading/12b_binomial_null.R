# Task (c)1, continued: what the draft's test gives for binomial reads with NO ploidy effect, after the
# scaling and rounding of step 05. Reads are drawn from the proportions fitted by the model WITHOUT ploidy
# (intercept + heat), so the ploidy test is a test of a true null. (In 12_dispersion_reference.R the reads were
# drawn from the full model's fitted proportions, which carry the fitted ploidy differences; its "share p <"
# columns are therefore not null shares and are not used.)
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
d <- readRDS(file.path(S, "sample3000_legacy.rds")); lib <- d$lib; Cs <- d$Cs; Ts <- d$Ts
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))[d$picked]
cov <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE]) + 0
filt <- fread(file.path(R, "04_output/05_differential_methylation/tables/library_filtering_legacy.csv")); sf <- max(filt$median_coverage) / filt$median_coverage
raw <- round(sweep(cov, 2, sf, "/"))
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat"); X <- cbind(1, trip, heat); X0 <- cbind(1, heat); fam <- binomial(link = logit)
mu0 <- t(vapply(1:nrow(Cs), function(i) suppressWarnings(glm.fit(X0, Cs[i, ] / (Cs[i, ] + Ts[i, ]), weights = Cs[i, ] + Ts[i, ], family = fam))$fitted.values, numeric(22)))
tests <- function(sc, st) t(vapply(1:nrow(sc), function(i) { y <- sc[i, ]; w <- y + st[i, ]
  o <- suppressWarnings(glm.fit(X, y / w, weights = w, family = fam)); o0 <- suppressWarnings(glm.fit(X0, y / w, weights = w, family = fam)); mu <- o$fitted.values
  phi <- sum((y - w * mu)^2 / (w * mu * (1 - mu))) / 19; dev <- o0$deviance - o$deviance
  c(phi = phi, dev = dev, p_none = pchisq(dev, 1, lower.tail = FALSE), p_f = if (phi > 1) pf(dev / phi, 1, 19, lower.tail = FALSE) else pchisq(dev, 1, lower.tail = FALSE)) }, numeric(4)))
sh <- function(p) paste(round(100 * c(mean(p < .05), mean(p < .01), mean(p < .001)), 2), collapse = ", ")
set.seed(1214); acc <- list()
for (rep in 1:5) {
  rc <- matrix(rbinom(length(raw), raw, mu0), nrow = nrow(raw)); rt <- raw - rc
  nc <- round(sweep(raw, 2, sf, "*")); s <- tests(round(nc * (rc / raw)), round(nc * (rt / raw))); n <- tests(rc, rt)
  acc[[rep]] <- c(scaled_infl = median(s[, "dev"]) / qchisq(.5, 1), scaled_mean = mean(s[, "dev"]), unscaled_infl = median(n[, "dev"]) / qchisq(.5, 1), unscaled_mean = mean(n[, "dev"]),
                  s05 = mean(s[, "p_none"] < .05), s01 = mean(s[, "p_none"] < .01), s001 = mean(s[, "p_none"] < .001),
                  n05 = mean(n[, "p_none"] < .05), n01 = mean(n[, "p_none"] < .01), n001 = mean(n[, "p_none"] < .001),
                  f05 = mean(s[, "p_f"] < .05), f01 = mean(s[, "p_f"] < .01), f001 = mean(s[, "p_f"] < .001), phi_med = median(s[, "phi"]))
}
a <- do.call(rbind, acc); print(round(a, 4)); cat("means over 5 simulations of 3,000 CpGs:\n"); print(round(colMeans(a), 4))

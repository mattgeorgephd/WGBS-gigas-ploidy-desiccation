# Task (c)2, part 1: sample 3,000 eligible CpGs of the legacy united table, reproduce methylKit's
# corrected p-value (overdispersion = "MN", test = "F", heat as covariate) in base R, and confirm
# against cache/fit_A1_ploidy_mn_f_legacy.rds and against methylKit:::logReg itself.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
pe <- readRDS(file.path(S, "pos_elig_legacy.rds"))
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
stopifnot(identical(u$start, pe$pos$start))
set.seed(606)                                   # my own seed, not the one of the audit script
n_sample <- 3000
elig <- which(pe$elig)
picked <- sort(sample(elig, n_sample))
Cs <- as.matrix(u[picked, paste0("numCs", 1:22), with = FALSE]) + 0
Ts <- as.matrix(u[picked, paste0("numTs", 1:22), with = FALSE]) + 0
saveRDS(list(picked = picked, Cs = Cs, Ts = Ts, pos = pe$pos[picked], lib = lib), file.path(S, "sample3000_legacy.rds"))
# a second, independent sample for a robustness run
set.seed(607); picked2 <- sort(sample(setdiff(elig, picked), n_sample))
saveRDS(list(picked = picked2, Cs = as.matrix(u[picked2, paste0("numCs", 1:22), with = FALSE]) + 0,
             Ts = as.matrix(u[picked2, paste0("numTs", 1:22), with = FALSE]) + 0, pos = pe$pos[picked2], lib = lib), file.path(S, "sample3000b_legacy.rds"))
rm(u); gc(verbose = FALSE)

# methylKit's logReg, rewritten in base R from the source of 1.28.0 (01_methylkit_source.out):
#   full model  ~ treatment + covariates, reduced model ~ covariates (intercept only when there are none)
#   phi = max(1, Pearson chi-squared of the full model / (libraries - parameters))
#   statistic = (deviance reduced - deviance full) / phi
#   test "F" and phi > 1: pf(statistic, ddf, libraries - parameters); otherwise pchisq(statistic, 1)
mk_test <- function(y, w, X_full, X_red) {
  prop <- y / w
  o0 <- suppressWarnings(glm.fit(X_red, prop, weights = w, family = binomial(link = logit)))
  o1 <- suppressWarnings(glm.fit(X_full, prop, weights = w, family = binomial(link = logit)))
  mu <- o1$fitted.values; nprm <- length(o1$coefficients); rdf <- length(w) - nprm
  phi_raw <- sum(((y - w * mu) / sqrt(mu * (w - w * mu)))^2) / rdf
  phi <- if (phi_raw > 1) phi_raw else 1
  dev <- o0$deviance - o1$deviance; ddf <- o0$df.residual - o1$df.residual
  c(dev = dev, phi_raw = phi_raw, ddf = ddf, rdf = rdf,
    p_none = pchisq(dev, 1, lower.tail = FALSE),
    p_mn_chisq = pchisq(dev / phi, 1, lower.tail = FALSE),
    p_mn_f = if (phi > 1) pf(dev / phi, ddf, rdf, lower.tail = FALSE) else pchisq(dev / phi, 1, lower.tail = FALSE))
}
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
X_full <- cbind(1, trip, heat); X_red <- cbind(1, heat)
t0 <- Sys.time()
obs <- t(vapply(seq_len(n_sample), function(i) mk_test(Cs[i, ], Cs[i, ] + Ts[i, ], X_full, X_red), numeric(7)))
cat("seconds for", n_sample, "CpGs (two fits each):", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n")
cache <- function(a) readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_legacy.rds", a)))$fit
relmax <- function(x, y) max(abs(x - y) / pmax(abs(y), .Machine$double.xmin))
c_f <- cache("A1_ploidy_mn_f")$pvalue[picked]; c_chi <- cache("A1_ploidy_mn_chisq")$pvalue[picked]; c_0 <- cache("A0_ploidy")$pvalue[picked]
cat("largest relative difference to the cached p-values: MN F", relmax(obs[, "p_mn_f"], c_f), "; MN chisq", relmax(obs[, "p_mn_chisq"], c_chi),
    "; uncorrected", relmax(obs[, "p_none"], c_0), "\n")
cat("identical to cache: MN F", identical(unname(obs[, "p_mn_f"]), c_f), "; number of CpGs with any difference", sum(obs[, "p_mn_f"] != c_f), "\n")
# against methylKit's own function, 300 CpGs
suppressPackageStartupMessages(library(methylKit))
vars <- data.frame(heat = factor(lib$treatment, levels = c("control", "heat")))
lr <- vapply(1:300, function(i) suppressWarnings(methylKit:::logReg(c(Cs[i, ], Ts[i, ]), vars, treatment = trip, overdispersion = "MN", test = "F"))[["p.value"]], numeric(1))
cat("largest relative difference to methylKit:::logReg at 300 CpGs:", relmax(obs[1:300, "p_mn_f"], lr), "\n")
cat("observed labels, the 3,000 sampled CpGs: share p < 0.05, 0.01, 0.001\n")
sh <- function(p) round(100 * c(mean(p < 0.05), mean(p < 0.01), mean(p < 0.001)), 3)
cat("  uncorrected:", sh(obs[, "p_none"]), "\n  MN chisq:", sh(obs[, "p_mn_chisq"]), "\n  MN F:", sh(obs[, "p_mn_f"]), "\n")
cat("  dispersion (raw, not floored): median", round(median(obs[, "phi_raw"]), 3), "; quartiles", round(quantile(obs[, "phi_raw"], c(.25, .75)), 3),
    "; share above 1.5:", mean(obs[, "phi_raw"] > 1.5), "; share at or below 1:", mean(obs[, "phi_raw"] <= 1), "\n")
# inflation of the uncorrected statistic against its reference distribution
cat("  median uncorrected statistic over the median of chi-squared(1):", round(median(obs[, "dev"]) / qchisq(0.5, 1), 3),
    "; median corrected statistic over the median of F(1,19):", round(median(obs[, "dev"] / pmax(obs[, "phi_raw"], 1)) / qf(0.5, 1, 19), 3), "\n")
saveRDS(obs, file.path(S, "sample3000_observed_tests.rds"))

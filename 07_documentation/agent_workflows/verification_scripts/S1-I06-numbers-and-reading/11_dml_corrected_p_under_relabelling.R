# Are the draft's 1,038 ploidy DML more consistent among animals than the "DML" of relabelled libraries?
# For the observed labels and for the first 150 random labellings of the audit script (seed 20261119):
# the CpGs whose pooled difference exceeds 20 points (all tested CpGs, legacy mode), and at those CpGs
# the uncorrected p-value and the corrected one (MN, F, heat as covariate), computed as methylKit does.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0; ns <- cs + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); rm(u); invisible(gc(verbose = FALSE))
observed <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
blocks <- split(seq_along(observed), lib$cohort)
set.seed(20261119)
labels <- cbind(observed, vapply(seq_len(999L), function(draw) { z <- observed; for (b in blocks) z[b] <- observed[b][sample.int(length(b))]; z }, integer(22)))
labels <- unname(labels[, 1:151])
rs_c <- rowSums(cs); rs_n <- rowSums(ns)
fam <- binomial(link = logit); X_red <- cbind(1, heat)
tests <- function(i, z) {
  y <- cs[i, ]; w <- ns[i, ]; X <- cbind(1, z, heat)
  o0 <- suppressWarnings(glm.fit(X_red, y / w, weights = w, family = fam)); o1 <- suppressWarnings(glm.fit(X, y / w, weights = w, family = fam))
  mu <- o1$fitted.values; phi_raw <- sum(((y - w * mu) / sqrt(mu * (w - w * mu)))^2) / 19; phi <- max(1, phi_raw); dev <- o0$deviance - o1$deviance
  c(p_none = pchisq(dev, 1, lower.tail = FALSE), p_mn_f = if (phi > 1) pf(dev / phi, 1, 19, lower.tail = FALSE) else pchisq(dev, 1, lower.tail = FALSE), phi = phi_raw)
}
t0 <- Sys.time()
res <- rbindlist(lapply(seq_len(ncol(labels)), function(r) {
  z <- labels[, r]; c1 <- drop(cs %*% z); n1 <- drop(ns %*% z)
  eff <- 100 * (c1 / n1 - (rs_c - c1) / (rs_n - n1)); above <- which(abs(eff) > 20)
  tt <- t(vapply(above, tests, numeric(3), z = z))
  data.table(labelling = r - 1L, n_above = length(above), n_unc_p_below_2.3e4 = sum(tt[, "p_none"] < 2.3e-4),
             n_f_05 = sum(tt[, "p_mn_f"] < 0.05), n_f_01 = sum(tt[, "p_mn_f"] < 0.01), n_f_001 = sum(tt[, "p_mn_f"] < 0.001), n_f_0001 = sum(tt[, "p_mn_f"] < 1e-4),
             median_phi = median(tt[, "phi"]), median_abs_diff = median(abs(eff[above])))
}))
cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
fwrite(res, file.path(S, "11_dml_corrected_p_under_relabelling.csv"))
options(width = 200)
o <- res[1]; r <- res[-1]
cat("OBSERVED labels:\n"); print(o)
cat(sprintf("  of the %d CpGs above 20 points: corrected p < 0.05 at %d (%.1f%%), < 0.01 at %d (%.1f%%), < 0.001 at %d (%.1f%%), < 0.0001 at %d (%.2f%%); uncorrected p < 2.3e-4 at %d\n",
            o$n_above, o$n_f_05, 100 * o$n_f_05 / o$n_above, o$n_f_01, 100 * o$n_f_01 / o$n_above, o$n_f_001, 100 * o$n_f_001 / o$n_above, o$n_f_0001, 100 * o$n_f_0001 / o$n_above, o$n_unc_p_below_2.3e4))
r[, `:=`(s05 = n_f_05 / n_above, s01 = n_f_01 / n_above, s001 = n_f_001 / n_above, s0001 = n_f_0001 / n_above, sq = n_unc_p_below_2.3e4 / n_above)]
cat("RELABELLED (150 labellings), share of the CpGs above 20 points with a corrected p below each level:\n")
for (v in c("s05", "s01", "s001", "s0001")) cat(sprintf("  %-6s mean %.1f%% | median %.1f%% | 5th %.1f%% | 95th %.1f%% | labellings with a share at least the observed: %d of %d\n", v, 100 * mean(r[[v]]),
    100 * median(r[[v]]), 100 * quantile(r[[v]], .05), 100 * quantile(r[[v]], .95),
    sum(r[[v]] >= c(s05 = o$n_f_05, s01 = o$n_f_01, s001 = o$n_f_001, s0001 = o$n_f_0001)[[v]] / o$n_above), nrow(r)))
cat(sprintf("  counts: corrected p < 0.01 among CpGs above 20 points: observed %d; relabelled median %g, 5th to 95th %g to %g; labellings with at least the observed count: %d\n",
            o$n_f_01, median(r$n_f_01), quantile(r$n_f_01, .05), quantile(r$n_f_01, .95), sum(r$n_f_01 >= o$n_f_01)))
cat(sprintf("  counts: corrected p < 0.001: observed %d; relabelled median %g, 5th to 95th %g to %g; at least observed: %d\n",
            o$n_f_001, median(r$n_f_001), quantile(r$n_f_001, .05), quantile(r$n_f_001, .95), sum(r$n_f_001 >= o$n_f_001)))
cat(sprintf("  uncorrected p < 2.3e-4 (the largest p with q < 0.01 in the draft's ploidy analysis) among relabelled CpGs above 20 points: %.2f%% overall; smallest share in a labelling %.2f%%\n",
            100 * sum(r$n_unc_p_below_2.3e4) / sum(r$n_above), 100 * min(r$sq)))
cat(sprintf("  median dispersion at the CpGs above 20 points: observed %.2f; relabelled median %.2f (5th to 95th %.2f to %.2f)\n", o$median_phi, median(r$median_phi),
            quantile(r$median_phi, .05), quantile(r$median_phi, .95)))

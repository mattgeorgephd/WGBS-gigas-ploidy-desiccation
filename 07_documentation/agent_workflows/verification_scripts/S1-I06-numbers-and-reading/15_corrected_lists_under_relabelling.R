# How long are the "corrected" lists (overdispersion "MN", test "F", as methylKit computes them) when the ploidy
# labels are reassigned? For each labelling: the CpGs above 20 points (all tested CpGs, legacy mode), and at those
# the corrected p-value. A CpG is counted as "called" when its corrected p-value is at or below the largest
# p-value that has q < 0.01 in the cached fit of the observed labels (a fixed p threshold, so no q-values are
# recomputed), and separately when it is below 0.01 / 480,579 (called whatever the other CpGs do).
suppressPackageStartupMessages(library(data.table)); setDTthreads(1); options(width = 220)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
cache <- function(a) readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_legacy.rds", a)))$fit
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cs_all <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0; ns_all <- cs_all + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); rm(u); invisible(gc(verbose = FALSE))
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat"); fam <- binomial(link = logit)
p_single <- 0.01 / nrow(cs_all)
run <- function(name, analysis, keep, labels, covariate) {
  f <- cache(analysis); p_star <- max(f$pvalue[f$qvalue < 0.01])
  cs <- cs_all[, keep, drop = FALSE]; ns <- ns_all[, keep, drop = FALSE]; rs_c <- rowSums(cs); rs_n <- rowSums(ns)
  rdf <- length(keep) - 2 - !is.null(covariate)
  X0 <- if (is.null(covariate)) NULL else cbind(1, covariate)
  one <- function(z) {
    c1 <- drop(cs %*% z); n1 <- drop(ns %*% z); eff <- 100 * (c1 / n1 - (rs_c - c1) / (rs_n - n1)); above <- which(abs(eff) > 20)
    X <- if (is.null(covariate)) cbind(1, z) else cbind(1, z, covariate)
    tt <- vapply(above, function(i) { y <- cs[i, ]; w <- ns[i, ]
      o <- suppressWarnings(glm.fit(X, y / w, weights = w, family = fam)); mu <- o$fitted.values
      dev <- (if (is.null(X0)) o$null.deviance else suppressWarnings(glm.fit(X0, y / w, weights = w, family = fam))$deviance) - o$deviance
      phi <- sum((y - w * mu)^2 / (w * mu * (1 - mu))) / rdf
      c(if (phi > 1) pf(dev / phi, 1, rdf, lower.tail = FALSE) else pchisq(dev, 1, lower.tail = FALSE), phi) }, numeric(2))
    p <- tt[1, ]; phi <- tt[2, ]
    c(n_above = length(above), called = sum(p <= p_star), called_chisq_route = sum(p <= p_star & phi <= 1), single = sum(p < p_single), single_chisq_route = sum(p < p_single & phi <= 1),
      n_phi_le1 = sum(phi <= 1))
  }
  t0 <- Sys.time(); res <- t(apply(labels, 2, one)); mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  res <- data.table(labelling = seq_len(nrow(res)) - 1L, res); fwrite(res, file.path(S, sprintf("15_corrected_lists_%s.csv", name)))
  o <- res[1]; r <- res[-1]
  cat(sprintf("\n===== %s (%s; %d residual df; %d other labellings; %.1f min) =====\n", name, analysis, rdf, nrow(r), mins))
  cat(sprintf("p threshold of the observed analysis (largest p with q < 0.01): %.3g; cached list: %d CpGs with q < 0.01 and above 20 points\n", p_star, sum(f$qvalue < 0.01 & abs(f$meth.diff) > 20)))
  cat("observed labels:\n"); print(o)
  qs <- function(x) sprintf("median %g, 5th to 95th percentile %g to %g, max %g", median(x), quantile(x, .05), quantile(x, .95), max(x))
  cat("other labellings, called at the observed p threshold:", qs(r$called), "; labellings with at least the observed number:", sum(r$called >= o$called), "of", nrow(r), "\n")
  cat("   of the called CpGs over all other labellings, share called through the chi-squared route (dispersion estimate not above 1):", round(sum(r$called_chisq_route) / max(1, sum(r$called)), 3), "\n")
  cat("other labellings, called on their own (p < 0.01 / 480,579):", qs(r$single), "; observed", o$single, "; labellings with at least the observed number:", sum(r$single >= o$single), "\n")
  cat("   mean number per labelling:", round(mean(r$single), 2), "; expected for 480,579 true null tests at that p: 0.01\n")
  invisible(res)
}
# heat cohort: all 126 distinct groupings, the observed first
hk <- which(lib$cohort == "heat_experiment"); cmb <- combn(10, 5); cmb <- cmb[, cmb[1, ] != 1]
Lh <- apply(cmb, 2, function(s) { z <- integer(10); z[s] <- 1L; z }); io <- colSums(Lh != trip[hk]) == 0; Lh <- cbind(Lh[, io], Lh[, !io])
res_h <- run("heat_cohort", "A4_ploidy_in_heat_mn_f", hk, Lh, NULL)
saveRDS(Lh, file.path(S, "15_labels_heat_cohort.rds"))

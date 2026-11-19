# Direct check that the corrected calls with unchanged p-values are CpGs whose dispersion estimate is not above 1.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1); options(width = 200)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
cache <- function(a) readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_legacy.rds", a)))$fit
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0; ns <- cs + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); pos <- u[, .(chr, start)]; rm(u)
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat"); fam <- binomial(link = logit)
phi_of <- function(i, keep, Xd) { y <- cs[i, keep]; w <- ns[i, keep]; o <- suppressWarnings(glm.fit(Xd, y / w, weights = w, family = fam)); mu <- o$fitted.values
  sum((y - w * mu)^2 / (w * mu * (1 - mu))) / (length(keep) - ncol(Xd)) }
# pooled ploidy and heat (same full model)
for (a in c("A1_ploidy_mn_f", "A1_heat_mn_f")) { f <- cache(a); d <- which(f$qvalue < 0.01 & abs(f$meth.diff) > 20)
  ph <- vapply(d, phi_of, numeric(1), keep = 1:22, Xd = cbind(1, trip, heat)); cat(a, ": dispersion estimates at the", length(d), "calls:", round(sort(ph), 2), "\n") }
hk <- which(heat == 1); f <- cache("A4_ploidy_in_heat_mn_f"); d <- which(f$qvalue < 0.01 & abs(f$meth.diff) > 20)
ph <- vapply(d, phi_of, numeric(1), keep = hk, Xd = cbind(1, trip[hk]))
cat("A4_ploidy_in_heat_mn_f: dispersion estimates at the 125 calls: not above 1 at", sum(ph <= 1), "; quartiles", round(quantile(ph, c(0, .25, .5, .75, 1)), 2), "\n")
pm <- 100 * cs[d, hk] / ns[d, hk]; dip <- pm[, trip[hk] == 0]; tri <- pm[, trip[hk] == 1]
sep <- apply(dip, 1, max) < apply(tri, 1, min) | apply(tri, 1, max) < apply(dip, 1, min)
cat("  of the 125: the five diploids and the five triploids do not overlap at", sum(sep), "; median within-group range (points): diploids", round(median(apply(dip, 1, function(x) diff(range(x)))), 1),
    ", triploids", round(median(apply(tri, 1, function(x) diff(range(x)))), 1), "; median |difference|", round(median(abs(f$meth.diff[d])), 1), "\n")
# all eligible CpGs: how often are the two heat-cohort groups non-overlapping, and how often is the estimate not above 1?
ck <- which(heat == 0); f <- cache("A4_ploidy_in_control_mn_f"); d <- which(f$qvalue < 0.01 & abs(f$meth.diff) > 20)
cat("A4_ploidy_in_control_mn_f: dispersion estimates at the 2 calls:", round(vapply(d, phi_of, numeric(1), keep = ck, Xd = cbind(1, trip[ck])), 2), "\n")
cat("P(chi-squared_8 / 8 * phi <= 1) for a true dispersion phi of 2, 3, 3.67, 5:", round(pchisq(8 / c(2, 3, 3.67, 5), 8), 4), "; with 10 df:", round(pchisq(10 / c(2, 3, 3.67, 5), 10), 4),
    "; with 19 df:", signif(pchisq(19 / c(2, 3, 3.67, 5), 19), 3), "\n")

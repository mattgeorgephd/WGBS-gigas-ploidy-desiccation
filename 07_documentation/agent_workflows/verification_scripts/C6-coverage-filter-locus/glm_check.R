# Independent refit (base R glm, no methylKit) of a sample of CpGs from the corrected united table.
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
stopifnot(identical(lib$legacy_sample_index, 1:22))
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
set.seed(20261119)
for (mode in c("corrected", "legacy")) {
  u <- fread(t5(sprintf("united_cpg_%s.csv.gz", mode)))
  dp <- fread(t5(sprintf("diffmeth_ploidy_%s.csv.gz", mode))); dh <- fread(t5(sprintf("diffmeth_heat_%s.csv.gz", mode)))
  stopifnot(identical(u$chr, dp$chr), identical(u$start, dp$start), identical(u$chr, dh$chr), identical(u$start, dh$start))
  isdml <- which((dp$qvalue < 0.01 & abs(dp$meth.diff) > 20) | (dh$qvalue < 0.01 & abs(dh$meth.diff) > 20))
  border <- which((dp$qvalue < 0.05 & abs(dp$meth.diff) > 18 & abs(dp$meth.diff) <= 20) | (dh$qvalue < 0.05 & abs(dh$meth.diff) > 18 & abs(dh$meth.diff) <= 20))
  idx <- if (mode == "corrected") sort(unique(c(isdml, border, sample(nrow(u), 3000)))) else sort(unique(c(sample(isdml, 300), sample(nrow(u), 700))))
  C <- as.matrix(u[idx, paste0("numCs", 1:22), with = FALSE]); T <- as.matrix(u[idx, paste0("numTs", 1:22), with = FALSE])
  fit1 <- function(cs, ts, g, cv) {
    ok <- TRUE
    f1 <- suppressWarnings(glm(cbind(cs, ts) ~ g + cv, family = binomial)); f0 <- suppressWarnings(glm(cbind(cs, ts) ~ cv, family = binomial))
    p <- pchisq(f0$deviance - f1$deviance, 1, lower.tail = FALSE); if (is.na(p)) p <- 1
    w <- cs + ts
    c(p = p, diff = 100 * (sum(cs[g == 1]) / sum(w[g == 1]) - sum(cs[g == 0]) / sum(w[g == 0])))
  }
  rp <- t(vapply(seq_along(idx), function(i) fit1(C[i, ], T[i, ], trip, heat), numeric(2)))
  rh <- t(vapply(seq_along(idx), function(i) fit1(C[i, ], T[i, ], heat, trip), numeric(2)))
  rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
  cat(sprintf("[%s] %d CpGs refit (%d DML of either list, %d near the 20-point line, rest random)\n", mode, length(idx), sum(idx %in% isdml), sum(idx %in% border)))
  cat(sprintf("   ploidy: max abs diff of meth.diff %.3g; max relative diff of p %.3g; p within 1e-6 relative for %d of %d\n",
              max(abs(rp[, "diff"] - dp$meth.diff[idx])), rel(rp[, "p"], dp$pvalue[idx]), sum(abs(rp[, "p"] - dp$pvalue[idx]) <= 1e-6 * pmax(dp$pvalue[idx], 1e-300)), length(idx)))
  cat(sprintf("   heat:   max abs diff of meth.diff %.3g; max relative diff of p %.3g; p within 1e-6 relative for %d of %d\n",
              max(abs(rh[, "diff"] - dh$meth.diff[idx])), rel(rh[, "p"], dh$pvalue[idx]), sum(abs(rh[, "p"] - dh$pvalue[idx]) <= 1e-6 * pmax(dh$pvalue[idx], 1e-300)), length(idx)))
  # DML status with my p-values but the table's q mapping cannot be recomputed on a subset; check the difference threshold side instead
  cat(sprintf("   |diff| > 20 agreement on refit CpGs: ploidy %s, heat %s\n", all((abs(rp[, "diff"]) > 20) == (abs(dp$meth.diff[idx]) > 20)), all((abs(rh[, "diff"]) > 20) == (abs(dh$meth.diff[idx]) > 20))))
  wp <- which(abs(rp[, "p"] - dp$pvalue[idx]) > 1e-6 * pmax(dp$pvalue[idx], 1e-300)); if (length(wp)) { cat("   ploidy p mismatches (mine, table):\n"); print(head(cbind(mine = rp[wp, "p"], table = dp$pvalue[idx][wp]), 5)) }
  wh <- which(abs(rh[, "p"] - dh$pvalue[idx]) > 1e-6 * pmax(dh$pvalue[idx], 1e-300)); if (length(wh)) { cat("   heat p mismatches (mine, table):\n"); print(head(cbind(mine = rh[wh, "p"], table = dh$pvalue[idx][wh]), 5)) }
  rm(u, dp, dh); gc(verbose = FALSE)
}

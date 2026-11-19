# Task (c)2, part 3: the same level check for the ploidy contrast inside each cohort
# (analyses A4_ploidy_in_heat_mn_f: 5 against 5, 8 residual df; A4_ploidy_in_control_mn_f: 6 against 6, 10 residual df).
# Heat cohort: all 126 distinct groupings (a grouping and its mirror image give the same p-values).
# Control cohort: the observed grouping and 150 random distinct ones of the 462.
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
d <- readRDS(file.path(S, "sample3000_legacy.rds")); lib <- d$lib; n_cpg <- nrow(d$Cs)
fam <- binomial(link = logit)
one_labelling <- function(Cs, W, z) {
  X <- cbind(1, z); rdf <- length(z) - 2
  out <- matrix(NA_real_, nrow(Cs), 4, dimnames = list(NULL, c("p_none", "p_mn_chisq", "p_mn_f", "diff")))
  for (i in seq_len(nrow(Cs))) {
    y <- Cs[i, ]; w <- W[i, ]
    o <- suppressWarnings(glm.fit(X, y / w, weights = w, family = fam)); mu <- o$fitted.values
    phi <- max(1, sum(((y - w * mu) / sqrt(mu * (w - w * mu)))^2) / rdf)
    dev <- o$null.deviance - o$deviance
    out[i, 1] <- pchisq(dev, 1, lower.tail = FALSE); out[i, 2] <- pchisq(dev / phi, 1, lower.tail = FALSE)
    out[i, 3] <- if (phi > 1) pf(dev / phi, 1, rdf, lower.tail = FALSE) else out[i, 2]
    out[i, 4] <- 100 * (sum(y[z == 1]) / sum(w[z == 1]) - sum(y[z == 0]) / sum(w[z == 0]))
  }
  out
}
sh <- function(p) c(mean(p < 0.05), mean(p < 0.01), mean(p < 0.001))
summarise <- function(o) c(none = sh(o[, "p_none"]), mn_chisq = sh(o[, "p_mn_chisq"]), mn_f = sh(o[, "p_mn_f"]), n_above20 = sum(abs(o[, "diff"]) > 20), median_diff = median(o[, "diff"]))
relmax <- function(x, y) max(abs(x - y) / pmax(abs(y), .Machine$double.xmin))
pc <- function(x) round(100 * x, 3)
options(width = 220)
run <- function(cohort, analysis, n_random, seed) {
  keep <- which(lib$cohort == cohort); Cs <- d$Cs[, keep]; W <- Cs + d$Ts[, keep]
  obs <- as.integer(lib$ploidy[keep] == "triploid"); k <- sum(obs)
  all_sets <- combn(length(keep), k)                        # positions of the "triploid" label
  all_sets <- all_sets[, all_sets[1, ] != 1, drop = FALSE]  # one of each mirror pair: library 1 of the cohort is never "triploid"
  labs <- apply(all_sets, 2, function(s) { z <- integer(length(keep)); z[s] <- 1L; z })
  is_obs <- colSums(labs != obs) == 0; stopifnot(sum(is_obs) == 1)
  others <- which(!is_obs)
  if (length(others) > n_random) { set.seed(seed); others <- sort(sample(others, n_random)) }
  labs <- cbind(labs[, is_obs], labs[, others, drop = FALSE])
  cat(sprintf("\n===== %s: %d libraries, %d against %d; %d distinct groupings, %d used besides the observed one =====\n", cohort, length(keep), k, length(keep) - k,
              ncol(all_sets), length(others)))
  t0 <- Sys.time()
  o_obs <- one_labelling(Cs, W, labs[, 1])
  cached <- readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_legacy.rds", analysis)))$fit$pvalue[d$picked]
  cat("largest relative difference of my p-values (observed labels) to the cached fit", analysis, ":", relmax(o_obs[, "p_mn_f"], cached), "\n")
  res <- rbind(summarise(o_obs), t(vapply(seq_len(ncol(labs))[-1], function(r) summarise(one_labelling(Cs, W, labs[, r])), numeric(11))))
  moved <- colSums(labs != obs); moved <- pmin(moved, length(keep) - moved)   # labels moved, up to the mirror image
  res <- cbind(labelling = seq_len(nrow(res)) - 1, moved = moved, res)
  write.csv(res, file.path(S, sprintf("08_relabel_level_%s.csv", cohort)), row.names = FALSE)
  cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
  for (s in c("none", "mn_chisq", "mn_f")) {
    cols <- paste0(s, 1:3); m <- res[-1, cols, drop = FALSE]
    cat(sprintf("%-9s observed %s | regrouped: mean %s | median %s | 5th pct %s | 95th pct %s | max %s | regroupings at or above observed %s of %d\n", s,
                paste(pc(res[1, cols]), collapse = ", "), paste(pc(colMeans(m)), collapse = ", "), paste(pc(apply(m, 2, median)), collapse = ", "),
                paste(pc(apply(m, 2, quantile, 0.05)), collapse = ", "), paste(pc(apply(m, 2, quantile, 0.95)), collapse = ", "), paste(pc(apply(m, 2, max)), collapse = ", "),
                paste(vapply(cols, function(cn) sum(m[, cn] >= res[1, cn]), numeric(1)), collapse = ", "), nrow(m)))
  }
  cat("CpGs above 20 points in the sample: observed", res[1, "n_above20"], "; regrouped median", median(res[-1, "n_above20"]), "; max", max(res[-1, "n_above20"]), "\n")
  cat("MN F share below 0.05 by number of labels moved (up to the mirror image):\n"); print(round(100 * tapply(res[-1, "mn_f1"], res[-1, "moved"], mean), 2)); print(table(res[-1, "moved"]))
  invisible(res)
}
run("heat_experiment", "A4_ploidy_in_heat_mn_f", 1000, 801)
run("control_experiment", "A4_ploidy_in_control_mn_f", 150, 802)

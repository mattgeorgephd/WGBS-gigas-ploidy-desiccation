# Task (c)2, part 2: does the corrected F test hold its level on these data?
# 3,000 eligible CpGs (legacy), ploidy labels reassigned within cohort, group sizes kept.
# For every relabelling: methylKit's p-value for overdispersion "MN" + test "F" with heat as covariate
# (and, from the same fits, the uncorrected test and "MN" + chi-squared), and the tested interaction
# (ploidy x heat against ploidy + heat; "MN" + F with 18 residual degrees of freedom).
args <- commandArgs(trailingOnly = TRUE)
sample_file <- if (length(args) >= 1) args[1] else "sample3000_legacy.rds"
n_relabel <- if (length(args) >= 2) as.integer(args[2]) else 150L
seed <- if (length(args) >= 3) as.integer(args[3]) else 707L
tag <- if (length(args) >= 4) args[4] else "a"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
d <- readRDS(file.path(S, sample_file)); Cs <- d$Cs; Ts <- d$Ts; lib <- d$lib; W <- Cs + Ts; n_cpg <- nrow(Cs)
trip <- as.integer(lib$ploidy == "triploid"); heat <- as.integer(lib$treatment == "heat")
blocks <- split(seq_along(trip), lib$cohort)
fam <- binomial(link = logit)
fit <- function(X, y, w) suppressWarnings(glm.fit(X, y / w, weights = w, family = fam))
pearson <- function(o, y, w) { mu <- o$fitted.values; sum(((y - w * mu) / sqrt(mu * (w - w * mu)))^2) / (length(w) - length(o$coefficients)) }
# reduced model of the ploidy test (intercept + heat): does not depend on the ploidy labels
X_red <- cbind(1, heat)
dev_red <- vapply(seq_len(n_cpg), function(i) fit(X_red, Cs[i, ], W[i, ])$deviance, numeric(1))
one_labelling <- function(z) {
  X_add <- cbind(1, z, heat); X_fac <- cbind(1, z, heat, z * heat)
  out <- matrix(NA_real_, n_cpg, 6, dimnames = list(NULL, c("p_none", "p_mn_chisq", "p_mn_f", "p_int_none", "p_int_mn_f", "diff")))
  for (i in seq_len(n_cpg)) {
    y <- Cs[i, ]; w <- W[i, ]
    oa <- fit(X_add, y, w); of <- fit(X_fac, y, w)
    phi <- max(1, pearson(oa, y, w)); dev <- dev_red[i] - oa$deviance
    out[i, 1] <- pchisq(dev, 1, lower.tail = FALSE)
    out[i, 2] <- pchisq(dev / phi, 1, lower.tail = FALSE)
    out[i, 3] <- if (phi > 1) pf(dev / phi, 1, 19, lower.tail = FALSE) else out[i, 2]
    phi_i <- max(1, pearson(of, y, w)); dev_i <- oa$deviance - of$deviance
    out[i, 4] <- pchisq(dev_i, 1, lower.tail = FALSE)
    out[i, 5] <- if (phi_i > 1) pf(dev_i / phi_i, 1, 18, lower.tail = FALSE) else pchisq(dev_i / phi_i, 1, lower.tail = FALSE)
    out[i, 6] <- 100 * (sum(y[z == 1]) / sum(w[z == 1]) - sum(y[z == 0]) / sum(w[z == 0]))
  }
  out
}
summarise <- function(o) {
  sh <- function(p) c(mean(p < 0.05), mean(p < 0.01), mean(p < 0.001))
  above <- abs(o[, "diff"]) > 20
  c(none = sh(o[, "p_none"]), mn_chisq = sh(o[, "p_mn_chisq"]), mn_f = sh(o[, "p_mn_f"]), int_none = sh(o[, "p_int_none"]), int_mn_f = sh(o[, "p_int_mn_f"]),
    n_above20 = sum(above), n_above20_p_below_2.3e4 = sum(above & o[, "p_none"] < 2.3e-4), median_diff = median(o[, "diff"]))
}
set.seed(seed)
labels <- cbind(trip, vapply(seq_len(n_relabel), function(r) { z <- trip; for (b in blocks) z[b] <- trip[b][sample.int(length(b))]; z }, integer(length(trip))))
t0 <- Sys.time()
res <- t(vapply(seq_len(ncol(labels)), function(r) summarise(one_labelling(labels[, r])), numeric(18)))
cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
# how far each labelling is from the truth: libraries whose label differs, by cohort
moved_heat <- colSums(labels[blocks$heat_experiment, , drop = FALSE] != trip[blocks$heat_experiment])
moved_ctrl <- colSums(labels[blocks$control_experiment, , drop = FALSE] != trip[blocks$control_experiment])
res <- cbind(labelling = c(0, seq_len(n_relabel)), moved_heat = moved_heat, moved_control = moved_ctrl, res)
write.csv(res, file.path(S, sprintf("07_relabel_level_pooled_%s.csv", tag)), row.names = FALSE)
saveRDS(labels, file.path(S, sprintf("07_relabel_labels_%s.rds", tag)))
options(width = 220)
pc <- function(x) round(100 * x, 3)
show <- function(cols, name, rows) {
  m <- res[rows, cols, drop = FALSE]
  cat(sprintf("%-34s mean %s | median %s | 5th pct %s | 95th pct %s | min %s | max %s\n", name,
              paste(pc(colMeans(m)), collapse = ", "), paste(pc(apply(m, 2, median)), collapse = ", "),
              paste(pc(apply(m, 2, quantile, 0.05)), collapse = ", "), paste(pc(apply(m, 2, quantile, 0.95)), collapse = ", "),
              paste(pc(apply(m, 2, min)), collapse = ", "), paste(pc(apply(m, 2, max)), collapse = ", ")))
}
sets <- list(none = c("none1", "none2", "none3"), mn_chisq = c("mn_chisq1", "mn_chisq2", "mn_chisq3"), mn_f = c("mn_f1", "mn_f2", "mn_f3"),
             int_none = c("int_none1", "int_none2", "int_none3"), int_mn_f = c("int_mn_f1", "int_mn_f2", "int_mn_f3"))
cat("\nShares (percent) of the", n_cpg, "sampled eligible CpGs with p < 0.05, 0.01, 0.001\n")
cat("OBSERVED labels:\n"); for (s in names(sets)) cat(sprintf("  %-10s %s\n", s, paste(pc(res[1, sets[[s]]]), collapse = ", ")))
for (k in unique(c(min(30L, n_relabel), n_relabel))) {
  cat(sprintf("\nRELABELLED, first %d relabellings:\n", k))
  for (s in names(sets)) show(sets[[s]], s, 1 + seq_len(k))
}
cat("\nRelabellings (of", n_relabel, ") whose share is at least the observed share, by level 0.05, 0.01, 0.001:\n")
for (s in names(sets)) cat(sprintf("  %-10s %s\n", s, paste(vapply(sets[[s]], function(cn) sum(res[-1, cn] >= res[1, cn]), numeric(1)), collapse = ", ")))
cat("\nCpGs above 20 points per labelling (in the sample): observed", res[1, "n_above20"], "; relabelled median", median(res[-1, "n_above20"]),
    "; of all relabelled CpG-instances above 20 points,", sum(res[-1, "n_above20"]), ", share with uncorrected p < 2.3e-4:",
    round(sum(res[-1, "n_above20_p_below_2.3e4"]) / sum(res[-1, "n_above20"]), 4), "; observed:", res[1, "n_above20_p_below_2.3e4"], "of", res[1, "n_above20"], "\n")
cat("\nCorrelation across relabellings of the MN F share below 0.05 with the number of labels moved (heat cohort, control cohort):",
    round(cor(res[-1, "mn_f1"], res[-1, "moved_heat"]), 3), round(cor(res[-1, "mn_f1"], res[-1, "moved_control"]), 3), "\n")
cat("MN F share below 0.05 by number of labels moved in the heat cohort:\n"); print(round(100 * tapply(res[-1, "mn_f1"], res[-1, "moved_heat"], mean), 2)); print(table(res[-1, "moved_heat"]))

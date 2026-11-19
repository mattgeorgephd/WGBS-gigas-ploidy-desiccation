# 05_code_checks.R  (scratch; claim C4b)
# Checks of my own code:
#   A  the characteristic-function computation of the exact conditional p-value against a direct
#      dynamic-programming computation (no Fourier transform), for a small, a medium and a large term
#   B  the convolution of hypergeometrics (stratified exact test) against brute-force shuffling of gene labels
#      within strata
#   C  sensitivity of the logistic null model to the form of the smooth (knots, degrees of freedom, bins)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(splines)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
P <- readRDS(file.path(scr, "prep.rds"))
source(file.path(scr, "03_stratified.R"))          # make_strata, conv, stratified_exact, build (no side effects when sourced)
set.seed(77)
lg <- fread(file.path(scr, "logistic_results.csv"))

pb_pmf <- function(p, kmax) {                       # Poisson-binomial pmf on 0..kmax by dynamic programming
  f <- c(1, numeric(kmax))
  for (pi in p) f <- f * (1 - pi) + c(0, f[-length(f)]) * pi
  f
}
knots_in <- log(c(2, 4, 8, 20, 50))
null_p <- function(n, y, Z = cbind(1, ns(log(n), knots = knots_in, Boundary.knots = log(c(1, max(n)))))) glm.fit(Z, y, family = binomial())$fitted.values

cat("== A: exact conditional p-value, characteristic functions vs dynamic programming ==\n")
for (case in list(c("heat", "CC", "GO:0043226"), c("heat", "CC", "GO:0030117"), c("heat", "MF", "GO:0001223"), c("ploidy", "BP", "GO:0010467"),
                  c("ploidy", "BP", "GO:0009409"), c("heat", "CC", "GO:0032991"))) {
  B <- build(case[2]); y <- as.numeric(B$u[[paste0("n_dml_", case[1])]] > 0); p0 <- null_p(B$u$n_pool, y)
  x <- as.numeric(B$X[, case[3]]) == 1; M <- sum(y); obs <- sum(y[x])
  ft <- pb_pmf(p0[x], M); fr <- pb_pmf(p0[!x], M)
  joint <- ft * rev(fr)                              # count k inside, M - k outside
  p_dp <- sum(joint[(0:M) >= obs]) / sum(joint)
  cat(sprintf("  %s %s %s: genes %d, DML genes %d | DP %.6g | characteristic function %.6g | ratio %.6f\n", case[1], case[2], case[3], sum(x), obs, p_dp,
              lg[set == case[1] & ontology == case[2] & go_id == case[3], p_exact], p_dp / lg[set == case[1] & ontology == case[2] & go_id == case[3], p_exact]))
}

cat("\n== B: stratified exact test vs brute-force shuffling within strata (200,000 shuffles) ==\n")
for (case in list(c("heat", "CC", "GO:0032991"), c("heat", "CC", "GO:0030117"), c("ploidy", "MF", "GO:0140657"), c("heat", "CC", "GO:0005829"))) {
  B <- build(case[2]); y <- as.numeric(B$u[[paste0("n_dml_", case[1])]] > 0); strata <- make_strata(B$u$n_pool, 50L)
  x <- as.numeric(B$X[, case[3]]); obs <- sum(x * y)
  ex <- stratified_exact(B$X[, case[3], drop = FALSE], y, strata)
  nb <- 200000L; cnt <- numeric(nb)
  idx <- split(seq_along(y), strata)
  for (g in idx) {                                   # shuffle the labels of each stratum: nb x |g| matrix of permuted labels times membership
    if (sum(x[g]) == 0 || sum(y[g]) == 0) next
    mg <- sum(y[g]); ng <- length(g)
    cnt <- cnt + vapply(seq_len(nb), function(i) sum(x[g][sample.int(ng, mg)]), numeric(1))
  }
  cat(sprintf("  %s %s %s: observed %d | exact: mean %.3f sd %.3f p %.4g | shuffles: mean %.3f sd %.3f p %.4g (%d of %d)\n", case[1], case[2], case[3], obs,
              ex$expected, sqrt(ex$variance), ex$p, mean(cnt), sd(cnt), mean(cnt >= obs), sum(cnt >= obs), nb))
}

cat("\n== C: form of the smooth in the logistic null model (per-term refit, signed-root likelihood ratio; exact conditional in brackets) ==\n")
exact_one <- function(p0, x, y) { M <- sum(y); ft <- pb_pmf(p0[x], M); fr <- pb_pmf(p0[!x], M); j <- ft * rev(fr); sum(j[(0:M) >= sum(y[x])]) / sum(j) }
forms <- list(
  "ns log n, knots 2,4,8,20,50 (used)" = function(n) cbind(1, ns(log(n), knots = knots_in, Boundary.knots = log(c(1, max(n))))),
  "linear in log n"                    = function(n) cbind(1, log(n)),
  "ns log n, 3 knots (3,10,40)"        = function(n) cbind(1, ns(log(n), knots = log(c(3, 10, 40)), Boundary.knots = log(c(1, max(n))))),
  "ns log n, 8 knots"                  = function(n) cbind(1, ns(log(n), knots = log(c(2, 3, 5, 8, 13, 22, 40, 80)), Boundary.knots = log(c(1, max(n))))),
  "13 bins (1,2,3,4,5,6-7,...,>100)"   = function(n) model.matrix(~ factor(cut(n, c(0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100, Inf)))),
  "ns of n itself (not log), 5 knots"  = function(n) cbind(1, ns(n, knots = c(2, 4, 8, 20, 50), Boundary.knots = c(1, max(n)))),
  "no adjustment (intercept only)"     = function(n) matrix(1, length(n), 1))
cases <- list(c("heat", "CC", "GO:0043226"), c("heat", "CC", "GO:0043229"), c("heat", "CC", "GO:0032991"), c("heat", "CC", "GO:0030117"),
              c("heat", "MF", "GO:0001223"), c("ploidy", "BP", "GO:0010467"), c("ploidy", "CC", "GO:0005829"), c("ploidy", "MF", "GO:0140657"))
rows <- list()
for (case in cases) {
  B <- build(case[2]); y <- as.numeric(B$u[[paste0("n_dml_", case[1])]] > 0); n <- B$u$n_pool; x <- as.numeric(B$X[, case[3]])
  for (f in names(forms)) {
    Z <- forms[[f]](n); f0 <- glm.fit(Z, y, family = binomial()); f1 <- suppressWarnings(glm.fit(cbind(Z, x), y, family = binomial()))
    b <- f1$coefficients[ncol(Z) + 1]; z <- sign(b) * sqrt(max(f0$deviance - f1$deviance, 0))
    rows[[length(rows) + 1]] <- data.table(set = case[1], ontology = case[2], go_id = case[3], form = f, null_aic = f0$aic,
                                           expected = sum(f0$fitted.values[x == 1]), odds_ratio = exp(b), p_lrt = pnorm(z, lower.tail = FALSE),
                                           p_exact = exact_one(f0$fitted.values, x == 1, y))
  }
}
C <- rbindlist(rows); fwrite(C, file.path(scr, "null_form_sensitivity.csv"))
options(width = 220); print(C[, .(set, ontology, go_id, form, null_aic = round(null_aic, 1), expected = round(expected, 1), OR = round(odds_ratio, 2),
                                  p_lrt = signif(p_lrt, 3), p_exact = signif(p_exact, 3))])

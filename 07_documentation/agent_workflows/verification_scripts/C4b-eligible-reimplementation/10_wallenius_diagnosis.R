# 10_wallenius_diagnosis.R (scratch; claim C4b)
# Why the Wallenius p-value of "organelle" (heat, CC) is larger than the exact one. Four null distributions of the
# number of DML genes in the term, all with the same fitted probabilities p0 (logistic null model) and the same total:
#   A  independent Bernoulli(p0_i) given the total (exact conditional test; Fisher's noncentral hypergeometric
#      with gene-specific odds)
#   B  univariate Wallenius with ONE weight for the term (mean p0 inside / mean p0 outside): what goseq
#      method = "Wallenius" computes
#   C  multivariate Wallenius: genes drawn one after another without replacement with probability proportional to
#      p0_i (gene-specific weights; what goseq method = "Sampling" simulates), by simulation here
#   D  as B but Fisher's noncentral hypergeometric with one odds ratio for the term
suppressPackageStartupMessages({library(data.table); library(Matrix); library(splines); library(BiasedUrn)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
P <- readRDS(file.path(scr, "prep.rds")); source(file.path(scr, "03_stratified.R")); set.seed(5)
pb_pmf <- function(p, kmax) { f <- c(1, numeric(kmax)); for (pi in p) f <- f * (1 - pi) + c(0, f[-length(f)]) * pi; f }
knots_in <- log(c(2, 4, 8, 20, 50))
rows <- list()
for (case in list(c("heat", "CC", "GO:0043226"), c("heat", "CC", "GO:0032991"), c("heat", "CC", "GO:0030117"), c("ploidy", "CC", "GO:0005829"),
                  c("ploidy", "BP", "GO:0010467"), c("heat", "MF", "GO:0005515"))) {
  B <- build(case[2]); y <- as.numeric(B$u[[paste0("n_dml_", case[1])]] > 0); n <- B$u$n_pool
  p0 <- glm.fit(cbind(1, ns(log(n), knots = knots_in, Boundary.knots = log(c(1, max(n))))), y, family = binomial())$fitted.values
  x <- as.numeric(B$X[, case[3]]) == 1; N <- length(y); M <- sum(y); m <- sum(x); obs <- sum(y[x]); kk <- 0:M
  # A
  j <- pb_pmf(p0[x], M) * rev(pb_pmf(p0[!x], M)); j <- j / sum(j)
  A <- c(mean = sum(kk * j), sd = sqrt(sum(kk^2 * j) - sum(kk * j)^2), p = sum(j[kk >= obs]))
  # B
  wgt <- mean(p0[x]) / mean(p0[!x])
  Bm <- meanWNCHypergeo(m, N - m, M, wgt); Bs <- sqrt(varWNCHypergeo(m, N - m, M, wgt))
  Bp <- dWNCHypergeo(obs, m, N - m, M, wgt) + pWNCHypergeo(obs, m, N - m, M, wgt, lower.tail = FALSE)
  # D
  odds <- (mean(p0[x]) / (1 - mean(p0[x]))) / (mean(p0[!x]) / (1 - mean(p0[!x])))
  Dm <- meanFNCHypergeo(m, N - m, M, odds); Ds <- sqrt(varFNCHypergeo(m, N - m, M, odds))
  Dp <- dFNCHypergeo(obs, m, N - m, M, odds) + pFNCHypergeo(obs, m, N - m, M, odds, lower.tail = FALSE)
  # C  (weighted sampling without replacement = the M smallest of Exp(1) / weight)
  nsim <- 20000L; cnt <- integer(nsim); incl <- numeric(N)
  for (b in seq_len(nsim)) { key <- rexp(N) / p0; sel <- key <= sort(key, partial = M)[M]; cnt[b] <- sum(sel[x]); incl <- incl + sel }
  rows[[length(rows) + 1]] <- data.table(set = case[1], ontology = case[2], go_id = case[3], genes = m, observed = obs, total = M,
    A_mean = A["mean"], A_sd = A["sd"], A_p = A["p"], B_weight = wgt, B_mean = Bm, B_sd = Bs, B_p = Bp, D_mean = Dm, D_sd = Ds, D_p = Dp,
    C_mean = mean(cnt), C_sd = sd(cnt), C_z = (obs - mean(cnt)) / sd(cnt), C_p_normal = pnorm((obs - 0.5 - mean(cnt)) / sd(cnt), lower.tail = FALSE),
    C_exceed = sum(cnt >= obs), C_nsim = nsim)
  if (case[3] == "GO:0043226") {
    bins <- cut(n, c(0, 1, 2, 5, 10, 20, 50, 100, Inf), labels = c("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"))
    cat("heat, CC universe: share of genes holding a DML by number of eligible CpGs (observed; fitted p0; inclusion frequency when genes are drawn with weights p0):\n")
    print(data.table(bin = levels(bins), genes = as.integer(table(bins)), observed = round(tapply(y, bins, mean), 3), fitted = round(tapply(p0, bins, mean), 3),
                     weighted_draws = round(tapply(incl / nsim, bins, mean), 3)))
  }
}
R <- rbindlist(rows); fwrite(R, file.path(scr, "wallenius_diagnosis.csv")); options(width = 250)
print(R[, .(set, ontology, go_id, genes, observed, A_mean = round(A_mean, 1), A_sd = round(A_sd, 2), A_p = signif(A_p, 3), B_mean = round(B_mean, 1), B_sd = round(B_sd, 2),
            B_p = signif(B_p, 3), D_mean = round(D_mean, 1), D_sd = round(D_sd, 2), D_p = signif(D_p, 3), C_mean = round(C_mean, 1), C_sd = round(C_sd, 2), C_z = round(C_z, 2),
            C_p_normal = signif(C_p_normal, 3), C_exceed)])

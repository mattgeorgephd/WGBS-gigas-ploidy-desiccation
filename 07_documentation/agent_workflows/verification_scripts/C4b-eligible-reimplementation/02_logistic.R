# 02_logistic.R  (scratch; claim C4b)   -- version 2
# Method (i): per-term logistic regression
#   gene holds a DML ~ natural spline of log(number of eligible CpGs) + term membership
# One-sided for enrichment. Four p-values per term:
#   p_score  null model fitted once; efficient score test (normal approximation)
#   p_lrt    model refitted for EVERY tested term; signed-root likelihood ratio (normal approximation)
#   p_exact  exact conditional test with the fitted null probabilities: the count of DML genes in the
#            term given the total number of DML genes, when gene i holds a DML with probability p0_i
#            independently (Poisson-binomial in the term and outside it; Fisher's noncentral
#            hypergeometric with gene-specific odds). No normal approximation; valid for small terms.
#   p_wrep   a replicate of the goseq calculation done with BiasedUrn directly: ONE weight for the
#            term (mean p0 inside / mean p0 outside) in a univariate Wallenius distribution.
# Benjamini-Hochberg across the terms of an ontology with at least 5 annotated universe genes.
suppressPackageStartupMessages({library(data.table); library(Matrix); library(splines); library(BiasedUrn)})
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scr, "prep.rds"))
suppressPackageStartupMessages(library(GO.db))
term_name <- function(ids) unname(AnnotationDbi::Term(GO.db::GOTERM[ids]))

knots_in <- log(c(2, 4, 8, 20, 50))            # interior knots of the spline of log(n); boundary knots at log(1) and log(max)
null_design <- function(n) cbind(1, ns(log(n), knots = knots_in, Boundary.knots = log(c(1, max(n)))))

# Exact conditional p-values P(S_term >= observed | total), all terms at once, by characteristic functions.
#   p0  null probability per gene; X  genes x terms 0/1; obs observed count per term; M total
exact_conditional <- function(p0, X, obs, M, L = 2048L) {
  stopifnot(M < L / 1.5)
  omega <- 2 * pi * (0:(L - 1)) / L
  lg <- log(outer(p0, exp(1i * omega), function(p, e) 1 - p + p * e))       # genes x L, complex log of each gene's factor
  Lt_re <- as.matrix(crossprod(X, Re(lg))); Lt_im <- as.matrix(crossprod(X, Im(lg)))
  La <- colSums(lg)
  pv <- numeric(ncol(X)); ex <- numeric(ncol(X))
  for (k in seq_len(ncol(X))) {
    lt <- complex(real = Lt_re[k, ], imaginary = Lt_im[k, ])
    f_t <- pmax(Re(fft(exp(lt), inverse = FALSE)) / L, 0)                    # pmf of the count inside the term (index = count + 1)
    f_r <- pmax(Re(fft(exp(La - lt), inverse = FALSE)) / L, 0)               # pmf of the count outside
    kk <- 0:M
    joint <- f_t[kk + 1] * f_r[M - kk + 1]
    pv[k] <- sum(joint[kk >= obs[k]]) / sum(joint)
    ex[k] <- sum(kk * joint) / sum(joint)
  }
  list(p = pmin(pv, 1), expected = ex)
}
# The discrete Fourier transform above uses exp(-i w k); the pmf is recovered from the characteristic function
# sum_k f(k) exp(+i w k) by the forward transform divided by L.

run_one <- function(o, s) {
  long <- P$memb[[o]]
  genes <- sort(unique(long$gene_id))
  tsize <- long[, .N, by = go_id]
  terms <- sort(tsize[N >= 5, go_id])
  X <- sparseMatrix(i = match(long$gene_id, genes), j = match(long$go_id, tsize$go_id), x = 1,
                    dims = c(length(genes), nrow(tsize)), dimnames = list(genes, tsize$go_id))[, terms]
  u <- P$uni[match(genes, gene_id)]
  n <- u$n_pool; y <- as.numeric(u[[paste0("n_dml_", s)]] > 0)
  N <- length(y); M <- sum(y)
  Z <- null_design(n)
  fit0 <- glm.fit(Z, y, family = binomial())
  stopifnot(fit0$converged)
  p0 <- fit0$fitted.values; w <- p0 * (1 - p0)
  m <- as.integer(colSums(X)); obs <- as.integer(crossprod(X, y))

  # ---- score test, all terms at once ----
  U <- as.numeric(crossprod(X, y - p0))
  xwx <- as.numeric(crossprod(X, w))
  A <- as.matrix(crossprod(X, w * Z))                       # terms x q
  V <- xwx - rowSums((A %*% solve(crossprod(Z, w * Z))) * A)
  root <- m == N                                            # every gene is a member: nothing to test
  z_score <- ifelse(root, 0, U / sqrt(pmax(V, 1e-12)))
  p_score <- ifelse(root, 1, pnorm(z_score, lower.tail = FALSE))

  # ---- per-term refit, every term ----
  dev0 <- fit0$deviance; n_fallback <- 0L
  prof_dev <- function(b, x) suppressWarnings(glm.fit(Z, y, family = binomial(), offset = b * x, start = fit0$coefficients)$deviance)
  lrt <- vapply(seq_along(terms), function(k) {
    if (root[k]) return(c(0, 0))
    x <- as.numeric(X[, k])
    f <- suppressWarnings(glm.fit(cbind(Z, x), y, family = binomial(), control = glm.control(maxit = 100)))
    b <- f$coefficients[ncol(Z) + 1]; d <- f$deviance
    if (!is.finite(d) || d > dev0 + 1e-6 || !f$converged) {          # fall back to a one-dimensional profile search
      n_fallback <<- n_fallback + 1L
      op <- optimize(prof_dev, c(-12, 12), x = x)
      b <- op$minimum; d <- op$objective
    }
    c(b, max(dev0 - d, 0))
  }, numeric(2))
  beta <- lrt[1, ]; z_lrt <- sign(beta) * sqrt(lrt[2, ])
  p_lrt <- ifelse(root, 1, pnorm(z_lrt, lower.tail = FALSE))

  # ---- exact conditional test with the fitted null probabilities ----
  ec <- exact_conditional(p0, X, obs, M)
  p_exact <- ifelse(root, 1, ec$p)

  # ---- replicate of the goseq calculation (one mean weight per term, univariate Wallenius) ----
  tot <- sum(p0); sum_in <- as.numeric(crossprod(X, p0))
  p_wrep <- vapply(seq_along(terms), function(k) {
    if (root[k]) return(1)
    wgt <- (sum_in[k] / m[k]) / ((tot - sum_in[k]) / (N - m[k]))
    dWNCHypergeo(obs[k], m[k], N - m[k], M, wgt) + pWNCHypergeo(obs[k], m[k], N - m[k], M, wgt, lower.tail = FALSE)
  }, numeric(1))

  res <- data.table(set = s, ontology = o, go_id = terms, term = term_name(terms),
                    genes_annotated = m, genes_with_dml = obs, expected_null = sum_in, expected_cond = ec$expected, log_or = beta,
                    var_score = V, var_binomial_mean_weight = M * (m / N) * (1 - m / N) * (N - M) / (N - 1),
                    z_score = z_score, p_score = p_score, z_lrt = z_lrt, p_lrt = p_lrt, p_exact = p_exact, p_wrep = p_wrep)
  for (v in c("score", "lrt", "exact", "wrep")) res[, (paste0("q_", v)) := p.adjust(get(paste0("p_", v)), "BH")]
  attr(res, "info") <- data.table(set = s, ontology = o, universe_genes = N, dml_genes = M, terms_tested = length(terms),
                                  null_deviance = fit0$null.deviance, deviance = dev0, lrt_fallbacks = n_fallback)
  res
}

out <- list(); info <- list()
for (s in c("ploidy", "heat")) for (o in c("BP", "CC", "MF")) {
  t0 <- Sys.time()
  r <- run_one(o, s)
  out[[length(out) + 1]] <- r; i <- attr(r, "info"); info[[length(info) + 1]] <- i
  cat(sprintf("%s %s: %d genes, %d with DML, %d terms | q<0.05: score %d, lrt %d, exact %d, wallenius-replicate %d | min q: score %.3g, lrt %.3g, exact %.3g, wrep %.3g | lrt fallbacks %d (%.0f s)\n",
              s, o, i$universe_genes, i$dml_genes, nrow(r), sum(r$q_score < 0.05), sum(r$q_lrt < 0.05), sum(r$q_exact < 0.05), sum(r$q_wrep < 0.05),
              min(r$q_score), min(r$q_lrt), min(r$q_exact), min(r$q_wrep), i$lrt_fallbacks, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- rbindlist(out); info <- rbindlist(info)

# ---- attach step 10's values (reading its output table only) ----
st <- rbindlist(lapply(c("ploidy", "heat"), function(s) {
  fread(file.path(repo, sprintf("04_output/10_go_enrichment_gene_level/tables/go_gene_level_%s_legacy.csv", s)))[universe == "eligible" & tested == TRUE,
    .(set, ontology, go_id, ga10 = genes_annotated, gd10 = genes_with_dml, expected10 = expected, p_fisher10 = p_fisher, q_fisher10 = q_fisher,
      p_wallenius10 = p_wallenius, q_wallenius10 = q_wallenius)]
}))
res <- merge(res, st, by = c("set", "ontology", "go_id"), all = TRUE)
cat(sprintf("\nterms: mine %d, step 10 %d, merged %d; annotated counts identical %s; DML-gene counts identical %s\n", sum(!is.na(res$p_lrt)),
            nrow(st), nrow(res), isTRUE(all(res$genes_annotated == res$ga10)), isTRUE(all(res$genes_with_dml == res$gd10))))
fwrite(res, file.path(scr, "logistic_results.csv"))
fwrite(info, file.path(scr, "logistic_info.csv"))
print(res[, .(spearman_exact_vs_wall10 = cor(p_exact, p_wallenius10, method = "spearman"),
              spearman_wrep_vs_wall10 = cor(p_wrep, p_wallenius10, method = "spearman"),
              med_ratio_exact_to_wall10_top50 = median((p_exact / p_wallenius10)[rank(p_wallenius10) <= 50]),
              med_ratio_wrep_to_wall10_top50 = median((p_wrep / p_wallenius10)[rank(p_wallenius10) <= 50]),
              med_ratio_lrt_to_exact_top50 = median((p_lrt / p_exact)[rank(p_exact) <= 50])), by = .(set, ontology)])

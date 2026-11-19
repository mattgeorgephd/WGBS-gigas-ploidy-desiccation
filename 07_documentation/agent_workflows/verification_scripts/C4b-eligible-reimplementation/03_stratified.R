# 03_stratified.R  (scratch; claim C4b)
# Method (iii), added as an assumption-free check of method (i): an exact stratified test.
# Genes are grouped by their number of eligible CpGs (every value its own stratum while it holds at
# least `min_size` genes, neighbouring values merged above that). Under the null hypothesis that, among
# genes with the same number of eligible CpGs, holding a DML is unrelated to term membership, the number
# of DML genes of a term is a sum over strata of independent hypergeometric counts. Its exact
# distribution is obtained by convolution (the exact form of the Cochran-Mantel-Haenszel test);
# one-sided p = P(count >= observed). This equals the limit of a permutation test that shuffles the
# DML status of genes within strata, with infinitely many permutations.
suppressPackageStartupMessages({library(data.table); library(Matrix)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
P <- readRDS(file.path(scr, "prep.rds"))

make_strata <- function(n, min_size) {
  tab <- table(n); vals <- as.integer(names(tab)); cnt <- as.integer(tab)
  bin <- integer(length(vals)); b <- 1L; acc <- 0L
  for (i in seq_along(vals)) { bin[i] <- b; acc <- acc + cnt[i]; if (acc >= min_size) { b <- b + 1L; acc <- 0L } }
  if (acc > 0L && acc < min_size && b > 1L) bin[bin == b] <- b - 1L
  bin[match(n, vals)]
}
conv <- function(f, h) {
  out <- numeric(length(f) + length(h) - 1L); lf <- length(f)
  for (j in seq_along(h)) if (h[j] > 0) out[j:(j + lf - 1L)] <- out[j:(j + lf - 1L)] + h[j] * f
  out
}
stratified_exact <- function(X, y, strata) {
  G <- max(strata)
  S <- sparseMatrix(i = seq_along(strata), j = strata, x = 1, dims = c(length(strata), G))
  A <- as.matrix(crossprod(X, S))                    # terms x strata: term genes per stratum
  n_g <- as.integer(colSums(S)); m_g <- as.integer(crossprod(S, y))
  obs <- as.integer(crossprod(X, y))
  out <- vapply(seq_len(ncol(X)), function(k) {
    a <- A[k, ]; f <- 1
    for (g in which(a > 0 & m_g > 0)) {
      if (a[g] == n_g[g]) { f <- c(rep(0, m_g[g]), f); next }       # every gene of the stratum is a member: count fixed
      f <- conv(f, dhyper(0:min(a[g], m_g[g]), a[g], n_g[g] - a[g], m_g[g]))
    }
    kk <- seq_along(f) - 1L
    c(sum(f[kk >= obs[k]]), sum(kk * f), sum(kk^2 * f) - sum(kk * f)^2)
  }, numeric(3))
  list(p = pmin(out[1, ], 1), expected = out[2, ], variance = out[3, ], obs = obs, n_strata = G, n_g = n_g, m_g = m_g)
}

build <- function(o) {
  long <- P$memb[[o]]
  genes <- sort(unique(long$gene_id))
  tsize <- long[, .N, by = go_id]; terms <- sort(tsize[N >= 5, go_id])
  X <- sparseMatrix(i = match(long$gene_id, genes), j = match(long$go_id, tsize$go_id), x = 1,
                    dims = c(length(genes), nrow(tsize)), dimnames = list(genes, tsize$go_id))[, terms]
  list(X = X, genes = genes, terms = terms, u = P$uni[match(genes, gene_id)])
}

if (sys.nframe() == 0L) {
  out <- list()
  for (min_size in c(50L, 25L, 100L)) for (o in c("BP", "CC", "MF")) {
    B <- build(o); strata <- make_strata(B$u$n_pool, min_size)
    for (s in c("ploidy", "heat")) {
      t0 <- Sys.time()
      y <- as.numeric(B$u[[paste0("n_dml_", s)]] > 0)
      r <- stratified_exact(B$X, y, strata)
      res <- data.table(set = s, ontology = o, min_stratum = min_size, n_strata = r$n_strata, go_id = B$terms,
                        genes_annotated = as.integer(colSums(B$X)), genes_with_dml = r$obs,
                        expected_strat = r$expected, var_strat = r$variance, p_strat = r$p)
      res[, q_strat := p.adjust(p_strat, "BH")]
      out[[length(out) + 1]] <- res
      cat(sprintf("min stratum %3d (%2d strata) %s %s: %d terms | q<0.05: %d | min q %.3g | smallest p %.3g (%s) (%.0f s)\n", min_size, r$n_strata, s, o,
                  nrow(res), sum(res$q_strat < 0.05), min(res$q_strat), min(res$p_strat), res[which.min(p_strat), go_id],
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  fwrite(rbindlist(out), file.path(scr, "stratified_results.csv"))
}

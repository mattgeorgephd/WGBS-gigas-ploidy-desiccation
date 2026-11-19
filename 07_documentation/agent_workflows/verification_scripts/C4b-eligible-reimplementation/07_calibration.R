# 07_calibration.R  (scratch; claim C4b)
# Error rates of the gene-level tests on null data that keep the real structure of the DML.
# Null data: the DML labels of one set are moved among the eligible CpGs of protein-coding genes
#   P5  rotated within classes of the gene's number of eligible CpGs and within chromosome (clusters kept;
#       the number of DML per class fixed)
#   P1  drawn uniformly among all eligible CpGs (the design of 06_diagnostics/go_random_cpg_control.R)
# and the genes holding a "DML" are recomputed. Such data carry no relation to GO by construction, but they
# keep overlapping genes, shared CpGs and (P5) clustering. Each null data set goes through
#   strat   exact stratified test (03_stratified.R), strata of at least 50 genes
#   exact   exact conditional test with probabilities from the logistic null model (as in 02_logistic.R)
#   wrep    replicate of the goseq calculation: one mean weight per term, univariate Wallenius (BiasedUrn)
# with Benjamini-Hochberg per ontology. Reported: the share of null data sets with at least one term at q < 0.05,
# and the share of all term p-values below 0.05, 0.01, 0.001 and 0.0001.
# Usage: Rscript 07_calibration.R <set> <generator P5|P1> <n_sim> <ontologies,comma> <tests,comma>
suppressPackageStartupMessages({library(data.table); library(Matrix); library(splines); library(BiasedUrn)})
args <- commandArgs(trailingOnly = TRUE)
s <- args[1]; gen <- args[2]; n_sim <- as.integer(args[3]); onts <- strsplit(args[4], ",")[[1]]; tests <- strsplit(args[5], ",")[[1]]
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
P <- readRDS(file.path(scr, "prep.rds"))
source(file.path(scr, "03_stratified.R"))
set.seed(31337 + n_sim + sum(utf8ToInt(paste0(s, gen))))
uni <- P$uni; genes_all <- uni$gene_id
lg_chr <- sprintf("NC_0475%02d.1", 59:68)
n_class <- function(n) cut(n, c(0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100, Inf), labels = FALSE)

# ---- null generator (same construction as in 06_cpg_permutation_adj.R) ----
cp <- P$cpg[eligible == TRUE | get(paste0("dml_", s)) == TRUE]
ovu <- P$ov[row %in% cp$row & gene_id %in% genes_all]; ovu[, n_gene := uni$n_pool[match(gene_id, genes_all)]]
cp[, gene_n := ovu[, .(m = max(n_gene)), by = row][cp, on = "row", m]]
if (gen == "P5") cp <- cp[!is.na(gene_n)]
cp[, chr_blk := ifelse(chr %in% lg_chr, chr, "other")]
cp[, stratum := if (gen == "P5") n_class(gene_n) else 1L]
setorder(cp, stratum, chr, start)
N <- nrow(cp); lab <- cp[[paste0("dml_", s)]]; K <- sum(lab)
ov <- ovu[row %in% cp$row]
C <- sparseMatrix(i = match(ov$row, cp$row), j = match(ov$gene_id, genes_all), x = 1, dims = c(N, length(genes_all)))
blk_key <- if (gen == "P5") paste(cp$stratum, cp$chr_blk) else rep("1", N)
blk_id <- match(blk_key, unique(blk_key)); blk_start <- as.integer(tapply(seq_len(N), blk_id, min)); blk_n <- as.integer(table(blk_id))
dpos <- which(lab); d_blk <- blk_id[dpos]; d_off <- dpos - blk_start[d_blk]
draw_one <- function() {
  if (gen == "P1") return(sample.int(N, K))
  sh <- floor(runif(length(blk_n)) * blk_n)
  blk_start[d_blk] + (d_off + sh[d_blk]) %% blk_n[d_blk]
}

# ---- tests ----
knots_in <- log(c(2, 4, 8, 20, 50))
exact_conditional <- function(p0, X, obs, M) {
  L <- as.integer(2^ceiling(log2(1.6 * M + 200)))
  omega <- 2 * pi * (0:(L - 1)) / L
  lg <- log(outer(p0, exp(1i * omega), function(p, e) 1 - p + p * e))
  Lt_re <- as.matrix(crossprod(X, Re(lg))); Lt_im <- as.matrix(crossprod(X, Im(lg))); La <- colSums(lg)
  vapply(seq_len(ncol(X)), function(k) {
    lt <- complex(real = Lt_re[k, ], imaginary = Lt_im[k, ])
    f_t <- pmax(Re(fft(exp(lt))) / L, 0); f_r <- pmax(Re(fft(exp(La - lt))) / L, 0)
    kk <- 0:M; joint <- f_t[kk + 1] * f_r[M - kk + 1]
    min(1, sum(joint[kk >= obs[k]]) / sum(joint))
  }, numeric(1))
}
B <- lapply(setNames(onts, onts), function(o) { b <- build(o); b$gi <- match(b$genes, genes_all); b$strata <- make_strata(b$u$n_pool, 50L)
  b$m <- as.integer(colSums(b$X)); b$Z <- cbind(1, ns(log(b$u$n_pool), knots = knots_in, Boundary.knots = log(c(1, max(b$u$n_pool))))); b })

rows <- list(); thr <- c(0.05, 0.01, 0.001, 0.0001)
t0 <- Sys.time()
for (i in seq_len(n_sim)) {
  hit <- as.numeric(as.numeric(crossprod(C, sparseVector(x = 1, i = draw_one(), length = N))) > 0)
  for (o in onts) {
    b <- B[[o]]; y <- hit[b$gi]; M <- sum(y); Ng <- length(y); obs <- as.integer(crossprod(b$X, y)); root <- b$m == Ng
    pv <- list()
    if ("strat" %in% tests) pv$strat <- stratified_exact(b$X, y, b$strata)$p
    if (any(c("exact", "wrep") %in% tests)) { p0 <- glm.fit(b$Z, y, family = binomial())$fitted.values }
    if ("exact" %in% tests) pv$exact <- ifelse(root, 1, exact_conditional(p0, b$X, obs, M))
    if ("wrep" %in% tests) {
      tot <- sum(p0); sum_in <- as.numeric(crossprod(b$X, p0))
      pv$wrep <- vapply(seq_along(obs), function(k) {
        if (root[k]) return(1)
        wgt <- (sum_in[k] / b$m[k]) / ((tot - sum_in[k]) / (Ng - b$m[k]))
        dWNCHypergeo(obs[k], b$m[k], Ng - b$m[k], M, wgt) + pWNCHypergeo(obs[k], b$m[k], Ng - b$m[k], M, wgt, lower.tail = FALSE)
      }, numeric(1))
    }
    for (tn in names(pv)) {
      p <- pv[[tn]]; q <- p.adjust(p, "BH")
      rows[[length(rows) + 1]] <- data.table(set = s, generator = gen, sim = i, ontology = o, test = tn, dml_genes = M, terms = length(p),
        n_q05 = sum(q < 0.05), min_p = min(p), min_q = min(q), below_05 = sum(p <= thr[1]), below_01 = sum(p <= thr[2]),
        below_001 = sum(p <= thr[3]), below_0001 = sum(p <= thr[4]),
        big_below_01 = sum(p[b$m >= 200] <= 0.01), big_below_001 = sum(p[b$m >= 200] <= 0.001), big_terms = sum(b$m >= 200))
    }
  }
  if (i %% 20 == 0) { fwrite(rbindlist(rows), file.path(scr, sprintf("calibration_%s_%s_%s.csv", s, gen, paste(tests, collapse = "-"))))
    cat(sprintf("%d of %d (%.1f min)\n", i, n_sim, as.numeric(difftime(Sys.time(), t0, units = "mins")))) }
}
R <- rbindlist(rows); fwrite(R, file.path(scr, sprintf("calibration_%s_%s_%s.csv", s, gen, paste(tests, collapse = "-"))))
options(width = 220)
print(R[, .(sims = .N, share_with_any_q05 = mean(n_q05 > 0), mean_n_q05 = mean(n_q05), median_dml_genes = median(dml_genes),
            frac_p_le_0.05 = sum(below_05) / sum(terms), frac_p_le_0.01 = sum(below_01) / sum(terms),
            frac_p_le_0.001 = sum(below_001) / sum(terms), frac_p_le_0.0001 = sum(below_0001) / sum(terms),
            big_frac_le_0.01 = sum(big_below_01) / sum(big_terms), big_frac_le_0.001 = sum(big_below_001) / sum(big_terms)), by = .(set, generator, ontology, test)])

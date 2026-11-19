# Shared helpers for my checks (scratch only). Re-implements the gene-level test of
# 02_R_functions/go_tests.R so that many simulated gene lists can be tested quickly.
suppressPackageStartupMessages({library(data.table); library(topGO); library(goseq); library(BiasedUrn); library(Matrix)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
sc <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory"
quiet <- function(expr) withCallingHandlers(suppressMessages(expr), warning = function(w) if (grepl("initial point very close", conditionMessage(w))) invokeRestart("muffleWarning"))

# term membership (with ancestors) for an ontology over a set of annotated genes, via topGO as in the pipeline
term_members <- function(ontology, gene_ids, gene2go) {
  start <- factor(c(1L, rep(0L, length(gene_ids) - 1L)), levels = 0:1); names(start) <- gene_ids
  td <- quiet(new("topGOdata", ontology = ontology, allGenes = start, annot = annFUN.gene2GO, gene2GO = gene2go[gene_ids], nodeSize = 1))
  mem <- genesInTerm(td)
  list(genes = genes(td), members = mem)
}
# sparse terms x genes incidence
incidence <- function(members, genes) {
  i <- rep(seq_along(members), lengths(members)); j <- match(unlist(members, use.names = FALSE), genes)
  sparseMatrix(i = i, j = j, x = 1, dims = c(length(members), length(genes)), dimnames = list(names(members), genes))
}
# goseq's Wallenius p-value (over-representation), replicated: weights are the mean PWF in / out of the category
wallenius_p <- function(M, de, pwf) {
  # M terms x genes (annotated genes only); de 0/1 for those genes; pwf for those genes
  N <- ncol(M); n_cat <- as.vector(rowSums(M)); x <- as.vector(M %*% de); num_de <- sum(de); alpha <- sum(pwf)
  avg_w <- as.vector(M %*% pwf) / n_cat
  w <- (avg_w * (N - n_cat)) / (alpha - n_cat * avg_w); w[n_cat == N] <- 1
  p <- vapply(seq_along(x), function(k) {
    dWNCHypergeo(x[k], n_cat[k], N - n_cat[k], num_de, w[k]) + pWNCHypergeo(x[k], n_cat[k], N - n_cat[k], num_de, w[k], lower.tail = FALSE)
  }, numeric(1))
  data.table(go_id = rownames(M), n_cat = n_cat, x = x, weight = w, p_wallenius = p)
}
fisher_p <- function(M, de) {
  N <- ncol(M); n_cat <- as.vector(rowSums(M)); x <- as.vector(M %*% de); num_de <- sum(de)
  phyper(x - 1, n_cat, N - n_cat, num_de, lower.tail = FALSE)
}
# PWF fitted on the whole universe (annotated or not), as go_gene_level_test() does
fit_pwf <- function(de_all, bias_all) { p <- quiet(nullp(de_all, bias.data = bias_all, plot.fit = FALSE)); setNames(p$pwf, rownames(p)) }
bh_count <- function(p, n_cat, min_genes = 5, max_genes = Inf, fdr = 0.05) { k <- n_cat >= min_genes & n_cat <= max_genes; q <- p.adjust(p[k], "BH"); c(tested = sum(k), n_sig = sum(q < fdr), min_q = min(q)) }

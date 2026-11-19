# 04_cpg_permutation.R  (scratch; claim C4b)   -- version 2
# Method (ii): permutation of the DML labels among CpGs, number of DML fixed; for every permutation the genes
# that hold at least one DML are recomputed, and the number of DML genes of each term is compared with the
# observed number.
# Schemes
#   P1  as specified in the task: labels shuffled uniformly among all eligible CpGs (plus the one DML of the
#       set whose mean lies just outside the 10-90% range)
#   P2  labels shuffled uniformly among the eligible CpGs that lie in protein-coding gene bodies; the number of
#       DML in gene bodies is fixed (conditions on the gene-body excess of DML)
#   P3  cluster preserving: the label vector of the genic eligible CpGs, in genome order, is rotated by a random
#       offset within each chromosome (unplaced scaffolds and the mitochondrion form one further block).
#       Spacing of DML along the eligible CpGs is kept, so DML that cluster stay clustered.
#   P4  as P3, among all eligible CpGs (genic and intergenic)
#   P5  stratified and cluster preserving: genic eligible CpGs are grouped by the number of eligible CpGs of their
#       gene (13 classes: 1, 2, 3, 4, 5, 6-7, 8-10, 11-15, 16-20, 21-30, 31-50, 51-100, >100); labels are rotated
#       within class and chromosome, so the number of DML per class is fixed. This is the CpG-level counterpart of
#       "adjusted for the number of eligible CpGs per gene".
#   P6  stratified, uniform shuffle within class (no cluster preservation)
# Statistics per term
#   raw    number of DML genes annotated to the term
#   share  that number divided by the number of DML genes among all annotated genes of the ontology
#          (removes the overall difference in the number of genes hit between the data and the permutations)
# p = (1 + permutations with statistic >= observed) / (1 + B); Benjamini-Hochberg across terms with >= 5 genes.
#
# Usage: Rscript 04_cpg_permutation.R <B> <schemes,comma> <sets,comma> <tag> [file with go_ids to restrict to]
suppressPackageStartupMessages({library(data.table); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
B_total <- if (length(args) >= 1) as.integer(args[1]) else 20000L
schemes <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("P1", "P2", "P3", "P4", "P5", "P6")
sets    <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("ploidy", "heat")
tag     <- if (length(args) >= 4) args[4] else sprintf("B%d", B_total)
restrict <- if (length(args) >= 5) readLines(args[5]) else NULL
chunk   <- 500L
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
P <- readRDS(file.path(scr, "prep.rds"))
set.seed(4101 + B_total)

uni <- P$uni; genes_all <- uni$gene_id
onto <- lapply(setNames(c("BP", "CC", "MF"), c("BP", "CC", "MF")), function(o) {
  long <- P$memb[[o]]; genes <- sort(unique(long$gene_id))
  tsize <- long[, .N, by = go_id]; terms <- sort(tsize[N >= 5, go_id])
  X <- sparseMatrix(i = match(long$gene_id, genes), j = match(long$go_id, tsize$go_id), x = 1,
                    dims = c(length(genes), nrow(tsize)), dimnames = list(genes, tsize$go_id))[, terms]
  keep <- if (is.null(restrict)) seq_along(terms) else which(terms %in% restrict)
  list(X = X[, keep, drop = FALSE], gi = match(genes, genes_all), terms = terms[keep], n_tested = length(terms))
})
lg_chr <- sprintf("NC_0475%02d.1", 59:68)
n_class <- function(n) cut(n, c(0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100, Inf), labels = FALSE)

run_scheme <- function(s, scheme) {
  cp <- P$cpg[eligible == TRUE | get(paste0("dml_", s)) == TRUE]
  ovu <- P$ov[row %in% cp$row & gene_id %in% genes_all]
  ovu[, n_gene := uni$n_pool[match(gene_id, genes_all)]]
  cp[, gene_n := ovu[, .(m = max(n_gene)), by = row][cp, on = "row", m]]            # NA when the CpG is in no protein-coding gene
  if (scheme %in% c("P2", "P3", "P5", "P6")) cp <- cp[!is.na(gene_n)]
  cp[, chr_blk := ifelse(chr %in% lg_chr, chr, "other")]
  cp[, stratum := if (scheme %in% c("P5", "P6")) n_class(gene_n) else 1L]
  setorder(cp, stratum, chr, start)
  N <- nrow(cp); lab <- cp[[paste0("dml_", s)]]; K <- sum(lab)
  ov <- ovu[row %in% cp$row]
  C <- sparseMatrix(i = match(ov$row, cp$row), j = match(ov$gene_id, genes_all), x = 1, dims = c(N, length(genes_all)))
  hit_obs <- as.numeric(crossprod(C, as.numeric(lab)) > 0)
  stopifnot(all(hit_obs == as.numeric(uni[[paste0("n_dml_", s)]] > 0)))
  # blocks: for the rotations (P3, P4, P5) stratum x chromosome; for the uniform shuffles (P1, P2, P6) stratum
  blk_key <- if (scheme %in% c("P3", "P4", "P5")) paste(cp$stratum, cp$chr_blk) else as.character(cp$stratum)
  blk_id <- match(blk_key, unique(blk_key))
  stopifnot(!is.unsorted(blk_id))
  blk_start <- as.integer(tapply(seq_len(N), blk_id, min)); blk_n <- as.integer(table(blk_id))
  dpos <- which(lab); d_blk <- blk_id[dpos]; d_off <- dpos - blk_start[d_blk]            # 0-based position inside the block
  k_blk <- tabulate(d_blk, nbins = length(blk_n))
  draw <- function(b) {
    if (scheme %in% c("P1", "P2", "P6")) {
      return(as.vector(do.call(rbind, lapply(which(k_blk > 0), function(g) {
        matrix(vapply(seq_len(b), function(i) blk_start[g] - 1L + sample.int(blk_n[g], k_blk[g]), integer(k_blk[g])), nrow = k_blk[g])
      }))))
    }
    sh <- matrix(floor(runif(length(blk_n) * b) * rep(blk_n, b)), nrow = length(blk_n))   # blocks x b offsets
    as.vector(blk_start[d_blk] + (d_off + sh[d_blk, , drop = FALSE]) %% blk_n[d_blk])
  }
  acc <- lapply(onto, function(o) list(obs = as.numeric(crossprod(o$X, hit_obs[o$gi])), M_obs = sum(hit_obs[o$gi]),
                                       ex_raw = 0, ex_share = 0, sT = 0, sT2 = 0, M = numeric(0)))
  hit_sum <- numeric(length(genes_all)); done <- 0L
  while (done < B_total) {
    b <- min(chunk, B_total - done)
    D <- sparseMatrix(i = draw(b), j = rep(seq_len(b), each = K), x = 1, dims = c(N, b))
    stopifnot(all(D@x == 1))
    G <- crossprod(C, D); G@x[] <- 1
    hit_sum <- hit_sum + rowSums(G)
    for (o in names(onto)) {
      Go <- G[onto[[o]]$gi, , drop = FALSE]
      T <- as.matrix(crossprod(onto[[o]]$X, Go)); M <- colSums(Go)
      a <- acc[[o]]
      a$ex_raw <- a$ex_raw + rowSums(T >= a$obs)
      a$ex_share <- a$ex_share + rowSums(sweep(T, 2, M, "/") >= a$obs / a$M_obs - 1e-12)
      a$sT <- a$sT + rowSums(T); a$sT2 <- a$sT2 + rowSums(T^2); a$M <- c(a$M, M)
      acc[[o]] <- a
    }
    done <- done + b
  }
  res <- rbindlist(lapply(names(onto), function(o) {
    a <- acc[[o]]
    data.table(set = s, scheme = scheme, ontology = o, go_id = onto[[o]]$terms, genes_annotated = as.integer(colSums(onto[[o]]$X)),
               genes_with_dml = a$obs, perm_mean = a$sT / B_total, perm_sd = sqrt(pmax(a$sT2 / B_total - (a$sT / B_total)^2, 0)),
               exceed_raw = a$ex_raw, exceed_share = a$ex_share, B = B_total,
               p_raw = (1 + a$ex_raw) / (1 + B_total), p_share = (1 + a$ex_share) / (1 + B_total))
  }))
  if (is.null(restrict)) res[, `:=`(q_raw = p.adjust(p_raw, "BH"), q_share = p.adjust(p_share, "BH")), by = ontology]
  tot <- rbindlist(lapply(names(onto), function(o) data.table(set = s, scheme = scheme, ontology = o, pool_cpg = N, dml_labels = K,
    blocks = length(blk_n), dml_genes_observed = acc[[o]]$M_obs, dml_genes_perm_mean = mean(acc[[o]]$M), dml_genes_perm_sd = sd(acc[[o]]$M),
    frac_perm_total_ge_observed = mean(acc[[o]]$M >= acc[[o]]$M_obs), B = B_total)))
  bins <- cut(uni$n_pool, c(0, 1, 2, 5, 10, 20, 50, 100, Inf), labels = c("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"))
  rate <- data.table(set = s, scheme = scheme, cpg_bin = factor(levels(bins), levels = levels(bins)), genes = as.integer(table(bins)),
                     observed_pct = 100 * as.numeric(tapply(hit_obs, bins, mean)), perm_pct = 100 * as.numeric(tapply(hit_sum / B_total, bins, mean)))
  list(res = res, tot = tot, rate = rate)
}

R <- list(); TT <- list(); RT <- list()
for (s in sets) for (scheme in schemes) {
  t0 <- Sys.time(); x <- run_scheme(s, scheme)
  R[[length(R) + 1]] <- x$res; TT[[length(TT) + 1]] <- x$tot; RT[[length(RT) + 1]] <- x$rate
  if (is.null(restrict)) {
    cat(sprintf("%s %s done in %.0f s: %s\n", s, scheme, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                paste(x$res[, sprintf("%s raw q<0.05 %d (min q %.3g), share q<0.05 %d (min q %.3g)", ontology[1], sum(q_raw < 0.05), min(q_raw),
                                      sum(q_share < 0.05), min(q_share)), by = ontology]$V1, collapse = " | ")))
  } else cat(sprintf("%s %s done in %.0f s\n", s, scheme, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  fwrite(rbindlist(R), file.path(scr, sprintf("cpg_perm_results_%s.csv", tag)))
  fwrite(rbindlist(TT), file.path(scr, sprintf("cpg_perm_totals_%s.csv", tag)))
  fwrite(rbindlist(RT), file.path(scr, sprintf("cpg_perm_hit_rate_%s.csv", tag)))
}
options(width = 250)
print(rbindlist(TT)); print(dcast(rbindlist(RT), set + cpg_bin + genes + observed_pct ~ scheme, value.var = "perm_pct"))

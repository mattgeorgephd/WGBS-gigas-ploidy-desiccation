# (d) Robustness of the gene-level result to the eligible range and to the bias variable.
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages({library(GenomicRanges); library(mgcv)})
CP <- readRDS(file.path(sc, "cpg_table_legacy.rds")); st <- readRDS(file.path(sc, "eligible_setup_legacy.rds"))
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
h <- suppressWarnings(findOverlaps(GRanges(CP$chr, IRanges(CP$start + 1L, CP$start + 2L)), GRanges(gt$chr, IRanges(gt$start, gt$end)), ignore.strand = TRUE))
pairs <- data.table(cpg = queryHits(h), gene_id = gt$gene_id[subjectHits(h)])
pairs[, `:=`(mean_pm = CP$mean_pm[cpg], sd_within = CP$sd_within[cpg], mean_cov = CP$mean_cov[cpg], dml_ploidy = CP$dml_ploidy[cpg], dml_heat = CP$dml_heat[cpg])]
stopifnot(all(gt[, .(gene_id, n = n_cpg_tested)][pairs[, .N, by = gene_id], on = "gene_id", n == N]))
# term membership over the widest universe (all protein-coding genes with a tested CpG and GO terms); sub-universes are column subsets
full_file <- file.path(sc, "tested_setup_legacy.rds")
if (!file.exists(full_file)) {
  tmf <- lapply(c(BP = "BP", CC = "CC", MF = "MF"), function(ont) { x <- term_members(ont, gt[n_cpg_tested > 0 & n_go_terms > 0, gene_id], st$gene2go_all); list(genes = x$genes, M = incidence(x$members, x$genes)) })
  saveRDS(tmf, full_file)
} else tmf <- readRDS(full_file)
# check: column subset of the full matrix equals the matrix built on the eligible universe
for (ont in names(tmf)) { g <- st$tm[[ont]]$genes; S <- tmf[[ont]]$M[, g]; S <- S[rowSums(S) > 0, ]; stopifnot(setequal(rownames(S), rownames(st$tm[[ont]]$M)), all(S[rownames(st$tm[[ont]]$M), ] == st$tm[[ont]]$M)) }
cat("column-subset membership equals the membership built by topGO on the eligible universe: TRUE\n")
any_dml <- gt$n_dml_ploidy + gt$n_dml_heat + gt$n_dml_shared > 0
run <- function(bias, label, pwf_method = "goseq") {
  bias_i <- pmax(bias, as.numeric(any_dml) * min(bias[bias > 0])); keep <- bias_i > 0
  uni <- cbind(gt[keep], bias = bias_i[keep])
  rbindlist(lapply(c("ploidy", "heat", "shared"), function(s) {
    de <- as.integer(uni[[paste0("n_dml_", s)]] > 0); names(de) <- uni$gene_id
    pwf <- if (pwf_method == "goseq") fit_pwf(de, uni$bias) else setNames(as.numeric(predict(gam(de ~ s(lb, k = 8), binomial, data.frame(de = de, lb = log2(uni$bias))), type = "response")), uni$gene_id)
    rbindlist(lapply(names(tmf), function(ont) {
      g <- intersect(tmf[[ont]]$genes, uni$gene_id); M <- tmf[[ont]]$M[, g]; M <- M[rowSums(M) > 0, ]
      r <- wallenius_p(M, de[g], pwf[g])
      rbindlist(lapply(c(5, 15, 50), function(k) { b <- bh_count(r$p_wallenius, r$n_cat, k); kk <- r$n_cat >= k
        data.table(analysis = label, set = s, ontology = ont, universe_genes = nrow(uni), annotated = length(g), dml_genes = sum(de[g]), min_genes = k, terms = b[["tested"]], n_q05 = b[["n_sig"]], min_q = round(b[["min_q"]], 4),
                   top = r[kk][which.min(p_wallenius), go_id], top_n = r[kk][which.min(p_wallenius), n_cat]) }))
    }))
  }))
}
count_in <- function(lo, hi) { x <- pairs[mean_pm >= lo & mean_pm <= hi, .N, by = gene_id]; n <- x[gt, on = "gene_id", N]; n[is.na(n)] <- 0L; n }
res <- list()
for (rg in list(c(5, 95), c(20 * 10 / 22, 100 - 20 * 10 / 22), c(10, 90), c(15, 85), c(20, 80), c(25, 75), c(30, 70))) {
  res[[length(res) + 1]] <- run(count_in(rg[1], rg[2]), sprintf("eligible range %.1f to %.1f, bias = number of such CpGs", rg[1], rg[2]))
  cat("done range", rg, "\n")
}
stopifnot(all(count_in(10, 90) == gt$n_cpg_eligible))
# alternative bias variables inside the 10-90 universe
el <- pairs[mean_pm >= 10 & mean_pm <= 90]
for (s in c("ploidy", "heat")) {
  f <- bam(as.formula(sprintf("dml_%s ~ s(mean_pm, k = 10) + s(sd_within, k = 10) + s(log(mean_cov), k = 6)", s)), binomial, unique(el[, .(cpg, mean_pm, sd_within, mean_cov, dml_ploidy, dml_heat)]), discrete = TRUE)
  el[, (paste0("phat_", s)) := as.numeric(predict(f, el, type = "response"))]
  f2 <- bam(as.formula(sprintf("dml_%s ~ s(mean_pm, k = 10)", s)), binomial, unique(el[, .(cpg, mean_pm, sd_within, mean_cov, dml_ploidy, dml_heat)]), discrete = TRUE)
  el[, (paste0("phat_mean_", s)) := as.numeric(predict(f2, el, type = "response"))]
}
gsum <- function(col) { x <- el[, .(v = sum(get(col))), by = gene_id]; v <- x[gt, on = "gene_id", v]; v[is.na(v)] <- 0; v }
# expected number of DML in the gene from the CpG-level landscape model (mean, within-cell SD, coverage); set-specific
alt <- function(col_p, col_h, label) {
  a <- run(gsum(col_p), label, "gam")[set == "ploidy"]; b <- run(gsum(col_h), label, "gam")[set == "heat"]; rbind(a, b)
}
res[[length(res) + 1]] <- run(gt$n_cpg_eligible, "10 to 90, bias = number of eligible CpGs, GAM on log2 scale", "gam")
res[[length(res) + 1]] <- alt("phat_mean_ploidy", "phat_mean_heat", "10 to 90, bias = expected DML from CpG means (GAM weighting)")
res[[length(res) + 1]] <- alt("phat_ploidy", "phat_heat", "10 to 90, bias = expected DML from CpG mean, within-cell SD, coverage (GAM weighting)")
res <- rbindlist(res); res[, term := substr(AnnotationDbi::Term(GO.db::GOTERM[top]), 1, 40)]
fwrite(res, file.path(sc, "range_bias_sensitivity.csv"))
cat("\nTerms with BH q<0.05 (Wallenius), min term size 5 / 15 / 50:\n")
w <- dcast(res, analysis + set + ontology + universe_genes + annotated + dml_genes ~ min_genes, value.var = c("n_q05", "min_q"))
print(w[order(analysis, set, ontology)], nrows = 300)
cat("\nAny analysis with a term at q<0.05 with min 5 genes:\n"); print(res[min_genes == 5 & n_q05 > 0])
cat("\nTop terms where min_q < 0.1 at min 5:\n"); print(res[min_genes == 5 & min_q < 0.1, .(analysis, set, ontology, min_q, top_n, term)])

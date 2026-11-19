# 09_strata_variants.R (scratch; claim C4b)
# Exact stratified test with other stratifications, to see what separates the gene-level result from the
# CpG-level permutation P5: (a) the 13 classes of P5; (b) those classes crossed with chromosome (P5 rotates
# within class and chromosome, so it conditions on the DML density of each chromosome).
suppressPackageStartupMessages({library(data.table); library(Matrix); library(GO.db)})
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scr, "prep.rds")); source(file.path(scr, "03_stratified.R"))
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))
lg_chr <- sprintf("NC_0475%02d.1", 59:68)
n_class <- function(n) cut(n, c(0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100, Inf), labels = FALSE)
out <- list()
for (o in c("BP", "CC", "MF")) {
  B <- build(o); chr <- gt$chr[match(B$genes, gt$gene_id)]; blk <- ifelse(chr %in% lg_chr, chr, "other")
  cls <- n_class(B$u$n_pool)
  variants <- list("27 strata by n (min 50 genes)" = make_strata(B$u$n_pool, 50L), "13 classes of P5" = cls,
                   "13 classes x chromosome" = as.integer(factor(paste(cls, blk))),
                   "6 classes x chromosome" = as.integer(factor(paste(cut(B$u$n_pool, c(0, 1, 2, 5, 10, 20, Inf), labels = FALSE), blk))))
  for (s in c("ploidy", "heat")) {
    y <- as.numeric(B$u[[paste0("n_dml_", s)]] > 0)
    # DML genes per chromosome, for the record
    if (o == "CC") { cat(s, "share of CC-annotated universe genes holding a DML, by chromosome block:\n"); print(round(tapply(y, blk, mean), 3)) }
    for (v in names(variants)) {
      r <- stratified_exact(B$X, y, variants[[v]])
      out[[length(out) + 1]] <- data.table(set = s, ontology = o, variant = v, n_strata = r$n_strata, go_id = B$terms, genes_annotated = as.integer(colSums(B$X)),
                                           genes_with_dml = r$obs, expected = r$expected, sd = sqrt(r$variance), p = r$p)[, q := p.adjust(p, "BH")]
    }
  }
}
R <- rbindlist(out); fwrite(R, file.path(scr, "strata_variants.csv"))
options(width = 220)
print(R[, .(n_strata = n_strata[1], terms = .N, n_q05 = sum(q < 0.05), min_q = signif(min(q), 3), top = go_id[which.min(p)], top_p = signif(min(p), 3)), by = .(set, ontology, variant)])
cat("\norganelle (GO:0043226), heat:\n"); print(R[set == "heat" & go_id == "GO:0043226", .(variant, n_strata, genes_with_dml, expected = round(expected, 1), sd = round(sd, 2), p = signif(p, 3), q = signif(q, 3))])
cat("\nmembrane coat (GO:0030117), heat:\n"); print(R[set == "heat" & go_id == "GO:0030117", .(variant, n_strata, genes_with_dml, expected = round(expected, 1), sd = round(sd, 2), p = signif(p, 3), q = signif(q, 3))])

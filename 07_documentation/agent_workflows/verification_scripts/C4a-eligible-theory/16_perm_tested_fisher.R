# Analysis A (tested universe, Fisher, BH) on permuted-label "DML": does the tested-universe enrichment appear without any ploidy signal?
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages(library(GenomicRanges))
ps <- readRDS(file.path(sc, "perm_sets_ploidy_legacy.rds")); CP <- readRDS(file.path(sc, "cpg_table_legacy.rds")); tmf <- readRDS(file.path(sc, "tested_setup_legacy.rds"))
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
h <- suppressWarnings(findOverlaps(GRanges(CP$chr, IRanges(CP$start + 1L, CP$start + 2L)), GRanges(gt$chr, IRanges(gt$start, gt$end)), ignore.strand = TRUE))
pairs <- data.table(cpg = queryHits(h), gene_id = gt$gene_id[subjectHits(h)]); setkey(pairs, cpg)
prep <- lapply(tmf, function(t) list(genes = t$genes, M = t$M[rowSums(t$M) >= 5, ]))
A <- function(genes) vapply(names(prep), function(ont) { P <- prep[[ont]]; sum(p.adjust(fisher_p(P$M, as.numeric(P$genes %in% genes)), "BH") < 0.05) }, numeric(1))
real <- A(gt[n_dml_ploidy > 0, gene_id]); cat("real ploidy DML, tested universe, Fisher BH q<0.05 (pipeline: 102, 48, 52):", real, "\n")
cat("terms tested:", vapply(prep, function(P) nrow(P$M), 1), " genes:", vapply(prep, function(P) length(P$genes), 1), "\n")
ok <- which(lengths(ps$sets) > 0)
out <- rbindlist(lapply(ok, function(b) { g <- unique(pairs[.(ps$sets[[b]]), nomatch = NULL]$gene_id); a <- A(g); data.table(perm = b, n_dml = length(ps$sets[[b]]), n_genes = length(g), BP = a[["BP"]], CC = a[["CC"]], MF = a[["MF"]]) }))
fwrite(out, file.path(sc, "perm_tested_fisher.csv"))
cat(sprintf("permutations: %d\n", nrow(out)))
for (o in c("BP", "CC", "MF")) cat(sprintf("%s: permuted-label lists give median %g terms (5%%-95%%: %g-%g; max %g); share with >= real (%g): %.3f; share with >= 10 terms: %.3f\n", o, median(out[[o]]), quantile(out[[o]], .05), quantile(out[[o]], .95), max(out[[o]]), real[[o]], mean(out[[o]] >= real[[o]]), mean(out[[o]] >= 10)))
k <- out[n_dml >= 800 & n_dml <= 1300]
for (o in c("BP", "CC", "MF")) cat(sprintf("  lists of 800-1300 loci (n=%d) %s: median %g (5%%-95%%: %g-%g); share >= real: %.3f\n", nrow(k), o, median(k[[o]]), quantile(k[[o]], .05), quantile(k[[o]], .95), mean(k[[o]] >= real[[o]])))

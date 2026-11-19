# Realistic negative control: "DML" from within-cohort permutations of the ploidy labels, passed through
# the location tests and the gene-level GO tests exactly as the real ploidy DML are.
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages(library(GenomicRanges))
ps <- readRDS(file.path(sc, "perm_sets_ploidy_legacy.rds")); cpg <- readRDS(file.path(sc, "cpg_table_legacy.rds"))
st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); uni <- st$uni; tm <- st$tm
stopifnot(nrow(cpg) == length(ps$elig), all(cpg$chr == ps$pos$chr))
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
h <- suppressWarnings(findOverlaps(GRanges(cpg$chr, IRanges(cpg$start + 1L, cpg$start + 2L)), GRanges(gt$chr, IRanges(gt$start, gt$end)), ignore.strand = TRUE))
pairs <- data.table(cpg = queryHits(h), gene_id = gt$gene_id[subjectHits(h)])
setkey(pairs, cpg)
# check: the real ploidy DML reproduce the gene list of the pipeline
real_rows <- which(cpg$dml_ploidy); real_genes <- unique(pairs[.(real_rows), nomatch = NULL]$gene_id)
cat(sprintf("real ploidy DML genes via my mapping: %d (pipeline: %d)\n", length(real_genes), sum(gt$n_dml_ploidy > 0)))
stopifnot(setequal(real_genes, gt[n_dml_ploidy > 0, gene_id]))

loc_or <- function(rows) {
  dml <- seq_len(nrow(cpg)) %in% rows; use <- cpg$eligible | dml
  vapply(c("gene", "CDS", "intron", "TE", "intergenic"), function(t) { tb <- table(factor(dml[use], c(TRUE, FALSE)), factor(cpg[[t]][use], c(TRUE, FALSE))); (tb[1, 1] * tb[2, 2]) / (tb[1, 2] * tb[2, 1]) }, numeric(1))
}
go_run <- function(dml_genes) {
  de_all <- as.integer(uni$gene_id %in% dml_genes); names(de_all) <- uni$gene_id
  pwf <- fit_pwf(de_all, uni$bias)
  rbindlist(lapply(names(tm), function(ont) {
    g <- tm[[ont]]$genes; M <- tm[[ont]]$M; keep <- rowSums(M) >= 5; M <- M[keep, ]
    r <- wallenius_p(M, de_all[g], pwf[g]); pf <- fisher_p(M, de_all[g])
    data.table(ontology = ont, de_genes = sum(de_all[g]), terms = nrow(M), C_wall_q05 = sum(p.adjust(r$p_wallenius, "BH") < 0.05), C_min_q = min(p.adjust(r$p_wallenius, "BH")),
               B_fisher_q05 = sum(p.adjust(pf, "BH") < 0.05), w_p05 = mean(r$p_wallenius < 0.05), w_p01 = mean(r$p_wallenius < 0.01), w_p001 = mean(r$p_wallenius < 0.001),
               # larger terms only (>= 50 genes): where discreteness does not mask miscalibration
               w_p05_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.05), w_p01_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.01), w_p001_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.001),
               C_q05_min15 = sum(p.adjust(r$p_wallenius[r$n_cat >= 15], "BH") < 0.05), C_q05_min100 = sum(p.adjust(r$p_wallenius[r$n_cat >= 100], "BH") < 0.05))
  }))
}
cat("real ploidy DML through my code:\n"); print(go_run(real_genes)); cat("real location ORs (eligible background):\n"); print(round(loc_or(real_rows), 2))
B <- as.integer(commandArgs(TRUE)[1]); if (is.na(B)) B <- 200
ok <- which(lengths(ps$sets) > 0)[seq_len(B)]
res <- list(); loc <- list()
for (b in ok) {
  rows <- ps$sets[[b]]; genes <- unique(pairs[.(rows), nomatch = NULL]$gene_id)
  loc[[length(loc) + 1]] <- data.table(perm = b, n_dml = length(rows), n_genes = length(genes), t(loc_or(rows)))
  res[[length(res) + 1]] <- cbind(perm = b, n_dml = length(rows), go_run(genes))
  if (length(res) %% 25 == 0) cat("done", length(res), "\n")
}
res <- rbindlist(res); loc <- rbindlist(loc)
fwrite(res, file.path(sc, "perm_go_results.csv")); fwrite(loc, file.path(sc, "perm_location_or.csv"))
cat("\nLocation ORs against eligible CpGs for permuted-label 'DML' (real ploidy: gene 2.52, CDS 1.68, intron 0.80, TE 0.54, intergenic 0.30):\n")
print(loc[, lapply(.SD, function(x) round(quantile(x, c(.025, .25, .5, .75, .975)), 2)), .SDcols = c("n_dml", "n_genes", "gene", "CDS", "intron", "TE", "intergenic")])
cat("\nGO, eligible universe, permuted-label 'DML':\n")
print(res[, .(perms = .N, de_genes_median = median(de_genes), C_any_q05 = mean(C_wall_q05 > 0), C_max = max(C_wall_q05), C_median_min_q = round(median(C_min_q), 3), B_fisher_median = median(B_fisher_q05), B_any = mean(B_fisher_q05 > 0), B_max = max(B_fisher_q05),
              w_p05 = round(mean(w_p05), 4), w_p01 = round(mean(w_p01), 4), w_p001 = round(mean(w_p001), 5), w_p05_big = round(mean(w_p05_big), 4), w_p01_big = round(mean(w_p01_big), 4), w_p001_big = round(mean(w_p001_big), 5),
              C_any_min15 = mean(C_q05_min15 > 0), C_any_min100 = mean(C_q05_min100 > 0)), by = ontology])
cat("share of permutations with any Wallenius q<0.05 term in any ontology:", mean(res[, .(any = any(C_wall_q05 > 0)), by = perm]$any), "\n")
cat("quantiles of smallest q (Wallenius, min 5 genes) across permutations, by ontology:\n"); print(res[, as.list(round(quantile(C_min_q, c(.01, .05, .1, .25, .5)), 3)), by = ontology])

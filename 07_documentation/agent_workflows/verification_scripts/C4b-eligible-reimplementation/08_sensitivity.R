# 08_sensitivity.R  (scratch; claim C4b)
# Sensitivity of the exact stratified result (03_stratified.R, strata of at least 50 genes) to choices that
# step 10 fixes:
#   S1  secondary GO IDs remapped to their current primary ID instead of dropped (395 of the 830 IDs that
#       GO.db 3.18.0 does not hold as primary IDs)
#   S2  minimum term size 10 and 20 instead of 5 (fewer tests)
#   S3  a CpG belongs to a gene only when its C lies in the gene body (instead of either base of the CpG)
#   S4  eligibility defined on the pooled methylation level (sum of methylated reads / sum of reads) instead
#       of the mean of the per-library percentages
#   S5  one family per DML set: Benjamini-Hochberg across the three ontologies together
suppressPackageStartupMessages({library(data.table); library(Matrix); library(GO.db); library(AnnotationDbi)})
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scr, "prep.rds"))
source(file.path(scr, "03_stratified.R"))
term_name <- function(ids) unname(Term(GOTERM[ids]))
anc_env <- list(BP = as.list(GOBPANCESTOR), CC = as.list(GOCCANCESTOR), MF = as.list(GOMFANCESTOR))

membership <- function(gg, universe_genes) {        # gg: gene_id, go_id (primary IDs)
  ont_of <- as.data.table(suppressMessages(AnnotationDbi::select(GO.db, keys = unique(gg$go_id), columns = "ONTOLOGY", keytype = "GOID")))
  setnames(ont_of, c("go_id", "ont")); gg <- merge(gg, ont_of, by = "go_id")[!is.na(ont)]
  lapply(setNames(c("BP", "CC", "MF"), c("BP", "CC", "MF")), function(o) {
    d <- unique(gg[ont == o & gene_id %in% universe_genes, .(gene_id, go_id)]); ids <- unique(d$go_id); a <- anc_env[[o]][ids]
    anc <- data.table(go_id = rep(ids, lengths(a)), anc = unlist(a, use.names = FALSE))[anc != "all"]
    unique(rbind(d, merge(d, anc, by = "go_id", allow.cartesian = TRUE)[, .(gene_id, go_id = anc)]))
  })
}
run <- function(memb, uni, label, min_genes = 5L) {
  rbindlist(lapply(c("BP", "CC", "MF"), function(o) {
    long <- memb[[o]]; genes <- sort(unique(long$gene_id)); tsize <- long[, .N, by = go_id]; terms <- sort(tsize[N >= min_genes, go_id])
    X <- sparseMatrix(i = match(long$gene_id, genes), j = match(long$go_id, tsize$go_id), x = 1, dims = c(length(genes), nrow(tsize)),
                      dimnames = list(genes, tsize$go_id))[, terms]
    u <- uni[match(genes, gene_id)]; strata <- make_strata(u$n_pool, 50L)
    rbindlist(lapply(c("ploidy", "heat"), function(s) {
      y <- as.numeric(u[[paste0("n_dml_", s)]] > 0); r <- stratified_exact(X, y, strata)
      data.table(analysis = label, set = s, ontology = o, universe_genes = length(genes), dml_genes = sum(y), go_id = terms,
                 genes_annotated = as.integer(colSums(X)), genes_with_dml = r$obs, expected = r$expected, p = r$p)[, q := p.adjust(p, "BH")]
    }))
  }))
}
summ <- function(R) R[, .(universe_genes = universe_genes[1], dml_genes = dml_genes[1], terms = .N, n_q05 = sum(q < 0.05), min_q = signif(min(q), 3),
                          top = paste0(go_id[which.min(p)], " ", substr(term_name(go_id[which.min(p)]), 1, 30), " (", genes_with_dml[which.min(p)], "/",
                                       genes_annotated[which.min(p)], ", p ", signif(min(p), 2), ")")), by = .(analysis, set, ontology)]
out <- list()

# baseline
out$base <- run(P$memb, P$uni, "baseline (as step 10)")

# S1 secondary IDs remapped
gg <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
syn <- as.list(GOSYNONYM); known <- keys(GO.db)
sec <- unique(gg$go_id[!gg$go_id %in% known & gg$go_id %in% names(syn)])
map <- setNames(vapply(sec, function(i) GOID(syn[[i]]), character(1)), sec)
gg2 <- copy(gg); gg2[go_id %in% sec, go_id := map[go_id]]; gg2 <- unique(gg2)
cat(sprintf("S1: %d secondary IDs remapped to %d primary IDs; gene-GO pairs with a known ID: %d before, %d after\n", length(sec), uniqueN(map),
            gg[go_id %in% known, .N], gg2[go_id %in% known, .N]))
out$s1 <- run(membership(gg2, P$uni$gene_id), P$uni, "S1 secondary IDs remapped")

# S2 minimum term size
out$s2a <- run(P$memb, P$uni, "S2 terms with >= 10 genes", 10L)
out$s2b <- run(P$memb, P$uni, "S2 terms with >= 20 genes", 20L)

# S3 and S4 need the CpG table again
dm <- file.path(repo, "04_output/05_differential_methylation/tables")
united <- fread(file.path(dm, "united_cpg_legacy.csv.gz"), showProgress = FALSE)
nc <- grep("^numCs", names(united), value = TRUE); nt <- sub("numCs", "numTs", nc)
Cm <- as.matrix(united[, ..nc]); Tm <- as.matrix(united[, ..nt])
mean_pm <- rowMeans(100 * Cm / (Cm + Tm)); pooled_pm <- 100 * rowSums(Cm) / (rowSums(Cm) + rowSums(Tm))
cpg <- united[, .(chr, start)]; cpg[, row := .I]; rm(united, Cm, Tm); invisible(gc())
for (s in c("ploidy", "heat")) { d <- fread(file.path(dm, sprintf("dml_%s_legacy.csv", s))); cpg[, (paste0("dml_", s)) := FALSE]; set(cpg, d$united_row, paste0("dml_", s), TRUE) }
g <- fread(file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_gene.gff.gz"), header = FALSE, sep = "\t", quote = "")
g[, `:=`(gene_id = sub("^.*GeneID:([0-9]+).*$", "\\1", V9), biotype = sub("^.*gene_biotype=([^;]+).*$", "\\1", V9))]
pc <- g[biotype == "protein_coding", .(gene_id, chr = V1, start = V4, end = V5)]; setkey(pc, chr, start, end)
gene_counts <- function(elig, c_only) {
  pool <- cpg[elig | dml_ploidy | dml_heat, .(chr, s1 = start + 1L, s2 = if (c_only) start + 1L else start + 2L, row)]
  ov <- foverlaps(pool, pc, by.x = c("chr", "s1", "s2"), by.y = c("chr", "start", "end"), type = "any", nomatch = NULL)[, .(row, gene_id)]
  u <- ov[, .(n_pool = .N, n_dml_ploidy = sum(cpg$dml_ploidy[row]), n_dml_heat = sum(cpg$dml_heat[row])), by = gene_id]
  u
}
gg1 <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
elig_mean <- mean_pm >= 10 & mean_pm <= 90; elig_pool <- pooled_pm >= 10 & pooled_pm <= 90
cat(sprintf("S4: eligible by mean of percentages %d; by pooled level %d; both %d; DML outside the pooled range: ploidy %d, heat %d\n", sum(elig_mean),
            sum(elig_pool), sum(elig_mean & elig_pool), sum(cpg$dml_ploidy & !elig_pool), sum(cpg$dml_heat & !elig_pool)))
u3 <- gene_counts(elig_mean, c_only = TRUE);  out$s3 <- run(membership(gg1, u3$gene_id), u3, "S3 CpG assigned by its C only")
u4 <- gene_counts(elig_pool, c_only = FALSE); out$s4 <- run(membership(gg1, u4$gene_id), u4, "S4 eligibility on pooled methylation")
cat(sprintf("universe genes: baseline %d, S3 %d, S4 %d\n", nrow(P$uni), nrow(u3), nrow(u4)))

R <- rbindlist(out); fwrite(R, file.path(scr, "sensitivity_results.csv"))
options(width = 230); print(summ(R))
# S5 one family per set
b <- out$base[, q_all := p.adjust(p, "BH"), by = set]
cat("\nS5: Benjamini-Hochberg across the three ontologies together (6,378 terms per set):\n")
print(b[, .(terms = .N, n_q05 = sum(q_all < 0.05), min_q = signif(min(q_all), 3), top = go_id[which.min(p)]), by = set])

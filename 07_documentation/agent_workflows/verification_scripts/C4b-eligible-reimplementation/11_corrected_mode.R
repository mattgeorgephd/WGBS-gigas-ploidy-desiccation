# 11_corrected_mode.R (scratch; claim C4b) -- outside the assigned scope (legacy mode), run only to see whether the
# exact stratified result for the heat DML carries over to the DML called with the coverage filter.
suppressPackageStartupMessages({library(data.table); library(Matrix); library(GO.db); library(AnnotationDbi)})
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(scr, "03_stratified.R"))
dm <- file.path(repo, "04_output/05_differential_methylation/tables")
united <- fread(file.path(dm, "united_cpg_corrected.csv.gz"), showProgress = FALSE)
nc <- grep("^numCs", names(united), value = TRUE); nt <- sub("numCs", "numTs", nc)
Cm <- as.matrix(united[, ..nc]); Tm <- as.matrix(united[, ..nt]); mean_pm <- rowMeans(100 * Cm / (Cm + Tm))
cpg <- united[, .(chr, start)]; cpg[, `:=`(row = .I, eligible = mean_pm >= 10 & mean_pm <= 90)]; rm(united, Cm, Tm); invisible(gc())
for (s in c("ploidy", "heat")) { d <- fread(file.path(dm, sprintf("dml_%s_corrected.csv", s))); stopifnot(all(cpg$start[d$united_row] == d$start))
  cpg[, (paste0("dml_", s)) := FALSE]; set(cpg, d$united_row, paste0("dml_", s), TRUE) }
cat("corrected mode: tested CpGs", nrow(cpg), "; eligible", sum(cpg$eligible), "; ploidy DML", sum(cpg$dml_ploidy), "; heat DML", sum(cpg$dml_heat), "\n")
g <- fread(file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_gene.gff.gz"), header = FALSE, sep = "\t", quote = "")
g[, `:=`(gene_id = sub("^.*GeneID:([0-9]+).*$", "\\1", V9), biotype = sub("^.*gene_biotype=([^;]+).*$", "\\1", V9))]
pc <- g[biotype == "protein_coding", .(gene_id, chr = V1, start = V4, end = V5)]; setkey(pc, chr, start, end)
pool <- cpg[eligible | dml_ploidy | dml_heat, .(chr, s1 = start + 1L, s2 = start + 2L, row)]
ov <- foverlaps(pool, pc, by.x = c("chr", "s1", "s2"), by.y = c("chr", "start", "end"), type = "any", nomatch = NULL)[, .(row, gene_id)]
uni <- ov[, .(n_pool = .N, n_dml_ploidy = sum(cpg$dml_ploidy[row]), n_dml_heat = sum(cpg$dml_heat[row])), by = gene_id]
gg <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
ont_of <- as.data.table(suppressMessages(AnnotationDbi::select(GO.db, keys = unique(gg$go_id), columns = "ONTOLOGY", keytype = "GOID"))); setnames(ont_of, c("go_id", "ont"))
gg <- merge(gg, ont_of, by = "go_id")[!is.na(ont)]
anc_env <- list(BP = as.list(GOBPANCESTOR), CC = as.list(GOCCANCESTOR), MF = as.list(GOMFANCESTOR))
st10 <- rbindlist(lapply(c("ploidy", "heat"), function(s) fread(file.path(repo, sprintf("04_output/10_go_enrichment_gene_level/tables/go_gene_level_%s_corrected.csv", s)))[universe == "eligible" & tested == TRUE]))
for (o in c("BP", "CC", "MF")) {
  d <- gg[ont == o & gene_id %in% uni$gene_id, .(gene_id, go_id)]; ids <- unique(d$go_id); a <- anc_env[[o]][ids]
  anc <- data.table(go_id = rep(ids, lengths(a)), anc = unlist(a, use.names = FALSE))[anc != "all"]
  long <- unique(rbind(d, merge(d, anc, by = "go_id", allow.cartesian = TRUE)[, .(gene_id, go_id = anc)]))
  genes <- sort(unique(long$gene_id)); tsize <- long[, .N, by = go_id]; terms <- sort(tsize[N >= 5, go_id])
  X <- sparseMatrix(i = match(long$gene_id, genes), j = match(long$go_id, tsize$go_id), x = 1, dims = c(length(genes), nrow(tsize)), dimnames = list(genes, tsize$go_id))[, terms]
  u <- uni[match(genes, gene_id)]; strata <- make_strata(u$n_pool, 50L)
  for (s in c("ploidy", "heat")) {
    y <- as.numeric(u[[paste0("n_dml_", s)]] > 0); r <- stratified_exact(X, y, strata); q <- p.adjust(r$p, "BH"); k <- which.min(r$p)
    s10 <- st10[set == s & ontology == o]; mm <- match(terms, s10$go_id)
    cat(sprintf("corrected %s %s: %d genes, %d with DML, %d terms (step 10: %d; counts identical %s) | exact stratified q<0.05: %d, min q %.3g, top %s %s (%d of %d, expected %.1f, p %.3g) | step 10 Wallenius min q %.3g\n",
                s, o, length(genes), sum(y), length(terms), nrow(s10), isTRUE(all(s10$genes_with_dml[mm] == r$obs)), sum(q < 0.05), min(q), terms[k],
                Term(GOTERM[[terms[k]]]), r$obs[k], as.integer(sum(X[, k])), r$expected[k], r$p[k], min(s10$q_wallenius)))
  }
}

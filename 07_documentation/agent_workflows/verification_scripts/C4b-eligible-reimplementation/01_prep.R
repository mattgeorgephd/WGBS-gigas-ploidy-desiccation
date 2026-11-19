# 01_prep.R  (scratch; independent re-implementation, claim C4b)
# Builds, without topGO / goseq / the pipeline helpers:
#   - eligible CpGs (mean methylation across the 22 libraries in [10, 90])
#   - CpG-to-gene incidence for protein-coding genes (own parser of the gene GFF, own overlap)
#   - DML labels per CpG (ploidy, heat)
#   - gene x GO-term membership with ancestors from GO.db (own propagation)
# and cross-checks each against the pipeline's tables.
suppressPackageStartupMessages({library(data.table); library(Matrix)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4b-eligible-reimplementation"
dm   <- file.path(repo, "04_output/05_differential_methylation/tables")
ga   <- file.path(repo, "04_output/08_gene_annotation/tables")

# ---- CpGs -------------------------------------------------------------------
united <- fread(file.path(dm, "united_cpg_legacy.csv.gz"), showProgress = FALSE)
cat("united rows:", nrow(united), "\n")
nc <- grep("^numCs", names(united), value = TRUE); nt <- sub("numCs", "numTs", nc)
stopifnot(length(nc) == 22)
pm <- matrix(0, nrow(united), 22)
for (j in 1:22) pm[, j] <- 100 * united[[nc[j]]] / (united[[nc[j]]] + united[[nt[j]]])
stopifnot(!anyNA(pm))
mean_pm <- rowMeans(pm)
pooled_pm <- 100 * rowSums(as.matrix(united[, ..nc])) / (rowSums(as.matrix(united[, ..nc])) + rowSums(as.matrix(united[, ..nt])))
eligible <- mean_pm >= 10 & mean_pm <= 90
cat(sprintf("eligible: %d (%.2f%%); below 10: %.2f%%; above 90: %.2f%%\n", sum(eligible), 100 * mean(eligible),
            100 * mean(mean_pm < 10), 100 * mean(mean_pm > 90)))
cpg <- united[, .(chr, start, end)]
cpg[, `:=`(row = .I, mean_pm = mean_pm, pooled_pm = pooled_pm, eligible = eligible)]
rm(united, pm); invisible(gc())

# ---- DML --------------------------------------------------------------------
for (s in c("ploidy", "heat")) {
  d <- fread(file.path(dm, sprintf("dml_%s_legacy.csv", s)))
  stopifnot(all(cpg$chr[d$united_row] == d$chr), all(cpg$start[d$united_row] == d$start), !anyDuplicated(d$united_row))
  cpg[, (paste0("dml_", s)) := FALSE]
  set(cpg, d$united_row, paste0("dml_", s), TRUE)
  cat(sprintf("%s DML: %d; eligible %d; not eligible %d (mean_pm of those: %s)\n", s, nrow(d), sum(cpg$eligible[d$united_row]),
              sum(!cpg$eligible[d$united_row]), paste(round(cpg$mean_pm[d$united_row][!cpg$eligible[d$united_row]], 2), collapse = ", ")))
}
cat("DML in either list:", sum(cpg$dml_ploidy | cpg$dml_heat), " in both:", sum(cpg$dml_ploidy & cpg$dml_heat),
    " either and not eligible:", sum((cpg$dml_ploidy | cpg$dml_heat) & !cpg$eligible), "\n")

# ---- genes (own parser) -----------------------------------------------------
g <- fread(file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_gene.gff.gz"), header = FALSE, sep = "\t", quote = "")
setnames(g, c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
g[, gene_id := sub("^.*GeneID:([0-9]+).*$", "\\1", attr)]
g[, biotype := sub("^.*gene_biotype=([^;]+).*$", "\\1", attr)]
stopifnot(!anyDuplicated(g$gene_id), all(grepl("^[0-9]+$", g$gene_id)))
cat("gene records:", nrow(g), "\n"); print(g[, .N, by = biotype][order(-N)])
pc <- g[biotype == "protein_coding", .(gene_id, chr, start, end, strand)]
cat("protein-coding genes:", nrow(pc), "\n")

# ---- CpG x gene overlap (own: data.table::foverlaps; the CpG is the two bases C and G) ----
pool <- cpg[eligible | dml_ploidy | dml_heat, .(chr, s1 = start + 1L, s2 = start + 2L, row)]
setkey(pc, chr, start, end)
ov <- foverlaps(pool, pc, by.x = c("chr", "s1", "s2"), by.y = c("chr", "start", "end"), type = "any", nomatch = NULL)
ov <- ov[, .(row, gene_id)]
cat("pool CpGs:", nrow(pool), " CpG-gene pairs:", nrow(ov), " pool CpGs in a pc gene:", uniqueN(ov$row),
    " in two or more:", sum(table(ov$row) > 1), "\n")
# the same for all tested CpGs (only for the cross-check of n_cpg_tested)
allc <- cpg[, .(chr, s1 = start + 1L, s2 = start + 2L, row)]
ov_all <- foverlaps(allc, pc, by.x = c("chr", "s1", "s2"), by.y = c("chr", "start", "end"), type = "any", nomatch = NULL)[, .(row, gene_id)]

per_gene <- merge(pc[, .(gene_id)], ov_all[, .(n_tested = .N), by = gene_id], by = "gene_id", all.x = TRUE)
el <- ov[cpg$eligible[row] == TRUE, .(n_eligible = .N), by = gene_id]
per_gene <- merge(per_gene, el, by = "gene_id", all.x = TRUE)
for (s in c("ploidy", "heat")) {
  x <- ov[cpg[[paste0("dml_", s)]][row] == TRUE, .(n = .N), by = gene_id]; setnames(x, "n", paste0("n_dml_", s))
  per_gene <- merge(per_gene, x, by = "gene_id", all.x = TRUE)
}
for (j in names(per_gene)[-1]) set(per_gene, which(is.na(per_gene[[j]])), j, 0L)

# cross-check against step 08's gene table
gt <- fread(file.path(ga, "gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))
chk <- merge(per_gene, gt[biotype == "protein_coding", .(gene_id, n_cpg_tested, n_cpg_eligible, t_p = n_dml_ploidy, t_h = n_dml_heat, n_go_terms)],
             by = "gene_id", all = TRUE)
cat(sprintf("check vs gene_table_legacy (protein coding, %d genes in mine, %d in step 08, %d merged):\n", nrow(per_gene),
            gt[biotype == "protein_coding", .N], nrow(chk)))
cat("  n_tested identical:", all(chk$n_tested == chk$n_cpg_tested), " n_eligible identical:", all(chk$n_eligible == chk$n_cpg_eligible),
    " n_dml_ploidy identical:", all(chk$n_dml_ploidy == chk$t_p), " n_dml_heat identical:", all(chk$n_dml_heat == chk$t_h), "\n")
cat("  genes with a tested CpG:", sum(per_gene$n_tested > 0), " with an eligible CpG:", sum(per_gene$n_eligible > 0),
    " with ploidy DML:", sum(per_gene$n_dml_ploidy > 0), " with heat DML:", sum(per_gene$n_dml_heat > 0), "\n")

# universe: protein-coding genes with a pool CpG (eligible, or a DML just outside the range)
per_gene[, n_pool := ov[, .N, by = gene_id][per_gene, on = "gene_id", N]]
per_gene[is.na(n_pool), n_pool := 0L]
uni <- per_gene[n_pool > 0]
cat("universe genes (>=1 pool CpG):", nrow(uni), " of them with n_eligible == 0:", sum(uni$n_eligible == 0),
    " step-10 definition pmax(n_eligible, any DML) > 0:", sum(pmax(per_gene$n_eligible, as.integer(per_gene$n_dml_ploidy + per_gene$n_dml_heat > 0)) > 0), "\n")

# ---- GO membership with ancestors (own propagation from GO.db) ------------------
suppressPackageStartupMessages({library(GO.db); library(AnnotationDbi)})
cat("GO.db", as.character(packageVersion("GO.db")), "\n")
gg <- fread(file.path(ga, "gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
cat("gene_go_map rows:", nrow(gg), " genes:", uniqueN(gg$gene_id), " GO IDs:", uniqueN(gg$go_id), "\n")
ont_of <- AnnotationDbi::select(GO.db, keys = unique(gg$go_id), columns = "ONTOLOGY", keytype = "GOID")
ont_of <- as.data.table(ont_of); setnames(ont_of, c("go_id", "ont"))
cat("GO IDs absent from GO.db (dropped):", sum(is.na(ont_of$ont)), "\n")
gg <- merge(gg, ont_of, by = "go_id")[!is.na(ont)]
anc_env <- list(BP = as.list(GOBPANCESTOR), CC = as.list(GOCCANCESTOR), MF = as.list(GOMFANCESTOR))
memb <- list()
for (o in c("BP", "CC", "MF")) {
  d <- gg[ont == o & gene_id %in% uni$gene_id, .(gene_id, go_id)]
  ids <- unique(d$go_id)
  a <- anc_env[[o]][ids]
  anc <- data.table(go_id = rep(ids, lengths(a)), anc = unlist(a, use.names = FALSE))[anc != "all"]
  long <- unique(rbind(d, merge(d, anc, by = "go_id", allow.cartesian = TRUE)[, .(gene_id, go_id = anc)]))
  memb[[o]] <- long
  cat(sprintf("%s: annotated universe genes %d; terms (with ancestors) %d; terms with >= 5 genes %d; gene-term pairs %d\n", o,
              uniqueN(long$gene_id), uniqueN(long$go_id), long[, .N, by = go_id][N >= 5, .N], nrow(long)))
}

saveRDS(list(cpg = cpg[eligible | dml_ploidy | dml_heat], ov = ov, per_gene = per_gene, uni = uni, memb = memb,
             n_tested_cpg = nrow(cpg), n_eligible = sum(cpg$eligible)), file.path(scr, "prep.rds"))
cat("saved prep.rds\n")

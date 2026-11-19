# Independent recomputation of "of the genes that hold a DML and have GO terms, the universe contains ..."
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; out <- commandArgs(trailingOnly = TRUE)[1]
zc <- function(f) sprintf("zcat %s/%s", repo, f)
# gene track parsed directly
g <- fread(cmd = zc("03_input_files/genome_features/cgigas_uk_roslin_v1_gene.gff.gz"), header = FALSE, sep = "\t", quote = "",
           select = c(1, 4, 5, 9), col.names = c("chr", "gs", "ge", "attr"))
g[, gene_id := sub(".*GeneID:([0-9]+).*", "\\1", attr)]
stopifnot(!anyDuplicated(g$gene_id))
setkey(g, chr, gs, ge)
mr <- fread(file.path(out, "all_mrnas.csv"), colClasses = c(gene_id = "character"))     # every mRNA: gene, in_full (has GO), in_uni
mr[, hyphen := grepl("-", parent, fixed = TRUE)]
# full map with true names: which mRNAs have GO terms once NM_ names are repaired?
full <- unique(fread(cmd = zc("03_input_files/annotation/geneid2go.tab.gz"), header = FALSE, sep = "\t", quote = "", col.names = c("id", "go")))
full <- full[grepl("GO:[0-9]{7}", go)]
mr[, num := sub("^[A-Z]{2}_", "", transcript_id)]
bx <- unique(fread(cmd = zc("03_input_files/annotation/20210605-cgigas-roslin-mito-blastx.outfmt6.gz"), header = FALSE, sep = "\t", select = 1)$V1)
mm <- regmatches(bx, regexec("_(mrna|miscrna|ncrna)_([A-Z]{2})_([0-9]+\\.[0-9]+)_[0-9]+$", bx)); mm <- do.call(rbind, mm[lengths(mm) == 4])
true_id <- unique(data.table(map_id = paste0("XM_", mm[, 4]), true = paste0(mm[, 3], "_", mm[, 4])))
full <- merge(full, true_id, by.x = "id", by.y = "map_id", all.x = TRUE)
mr[, has_go_repaired := transcript_id %in% full$true]          # gene route: true IDs, so NM_ mRNAs count
mr[, has_go_legacy_name := in_full]                              # legacy name match (XM_ only)
cat("mRNAs with GO, repaired names:", sum(mr$has_go_repaired), "; by legacy name:", sum(mr$has_go_legacy_name),
    "; NM_ among repaired:", sum(mr$has_go_repaired & grepl("^NM_", mr$transcript_id)), "\n")
gene_flags <- mr[, .(has_go_gene_route = any(has_go_repaired), has_go_legacy_names = any(has_go_legacy_name & !hyphen),
                     in_universe = any(in_uni)), by = gene_id]
legacy_tab <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")
res <- list()
for (s in names(legacy_tab)) {
  d <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/dml_%s_legacy.csv", s)))
  d <- d[, .(chr, cs = start + 1L, ce = start + 2L)]; setkey(d, chr, cs, ce)
  ov <- foverlaps(d, g, by.x = c("chr", "cs", "ce"), by.y = c("chr", "gs", "ge"), type = "any", nomatch = NULL)
  genes <- unique(ov$gene_id)
  gf <- merge(data.table(gene_id = genes), gene_flags, by = "gene_id", all.x = TRUE)
  for (cc in c("has_go_gene_route", "has_go_legacy_names", "in_universe")) set(gf, which(is.na(gf[[cc]])), cc, FALSE)
  # the real legacy annotation table of the draft (from gannet): transcripts and genes that have GOslim-mapped GO terms
  lt <- fread(cmd = zc(sprintf("03_input_files/annotation/legacy/DML-Cov10-20-BOTH-%s.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz", legacy_tab[[s]])),
              header = FALSE, sep = "\t", quote = "", select = 1:2, col.names = c("transcript", "gene_id"), colClasses = list(character = 1:2))
  lt <- unique(lt)
  lt[, tx_in_uni := transcript %in% mr[in_uni == TRUE, transcript_id]]
  lg <- lt[, .(any_listed_tx_in_universe = any(tx_in_uni)), by = gene_id]
  res[[s]] <- data.table(set = s, dml = nrow(d), dml_gene_pairs = nrow(ov), genes_with_dml = length(genes),
    A_genes_go_gene_route = sum(gf$has_go_gene_route), A_in_universe = sum(gf$has_go_gene_route & gf$in_universe),
    B_genes_in_legacy_table = uniqueN(lt$gene_id), B_tx_in_legacy_table = uniqueN(lt$transcript),
    B_genes_with_listed_tx_in_universe = sum(lg$any_listed_tx_in_universe), B_tx_in_universe = sum(lt$tx_in_uni),
    B_genes_all_in_dml_genes = all(lt$gene_id %in% genes),
    C_genes_go_legacy_names = sum(gf$has_go_legacy_names), C_in_universe = sum(gf$has_go_legacy_names & gf$in_universe))
}
res <- rbindlist(res)
res[, `:=`(A_pct = round(100 * A_in_universe / A_genes_go_gene_route, 2),
           B_pct_genes = round(100 * B_genes_with_listed_tx_in_universe / B_genes_in_legacy_table, 2),
           B_pct_tx = round(100 * B_tx_in_universe / B_tx_in_legacy_table, 2),
           C_pct = round(100 * C_in_universe / C_genes_go_legacy_names, 2))]
print(t(res))
fwrite(res, file.path(out, "dml_gene_universe_shares.csv"))

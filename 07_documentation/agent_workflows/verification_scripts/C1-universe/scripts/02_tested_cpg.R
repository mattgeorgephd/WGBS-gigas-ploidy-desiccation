# Alternative explanation test, part 1: did the missing transcripts simply have no covered CpG?
# n_cpg_tested = CpGs with >= 10x in ALL 22 libraries inside the gene body (step 08, legacy mode).
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; out <- commandArgs(trailingOnly = TRUE)[1]
ann <- fread(file.path(out, "annotated_mrnas.csv"), colClasses = c(gene_id = "character"))
gt  <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))
tt  <- fread(file.path(repo, "04_output/08_gene_annotation/tables/transcript_table.csv.gz"), colClasses = c(gene_id = "character"))
# consistency of my own parse with the pipeline's transcript table
m <- merge(ann[, .(transcript_id, gene_id, in_uni)], tt[, .(transcript_id, gene_id2 = gene_id, in_legacy_go_map, in_legacy_universe)], by = "transcript_id")
cat("annotated mRNAs matched in transcript_table:", nrow(m), "; gene ids agree:", all(m$gene_id == m$gene_id2),
    "; in_legacy_go_map all TRUE:", all(m$in_legacy_go_map), "; universe flag agrees:", all(m$in_uni == m$in_legacy_universe), "\n")
cat("transcript_table: in_legacy_go_map TRUE:", sum(tt$in_legacy_go_map), " in_legacy_universe TRUE:", sum(tt$in_legacy_universe), "\n")
stopifnot(!anyDuplicated(gt$gene_id))
ann <- merge(ann, gt[, .(gene_id, n_cpg_tested, n_cpg_eligible, biotype, gene_len = end - start + 1L)], by = "gene_id", all.x = TRUE)
cat("annotated mRNAs without a row in the gene table:", sum(is.na(ann$n_cpg_tested)), "\n")
ann[, status := fifelse(in_uni, "present in universe", fifelse(in_block, "missing: in the block", "missing: 49 others"))]
print(ann[, .(transcripts = .N, genes = uniqueN(gene_id),
              tx_gene_has_tested_cpg = sum(n_cpg_tested > 0), pct_tx = round(100 * mean(n_cpg_tested > 0), 1),
              median_tested_cpg = as.numeric(median(n_cpg_tested)), mean_tested_cpg = round(mean(n_cpg_tested), 1),
              tx_gene_has_eligible_cpg = sum(n_cpg_eligible > 0)), by = status])
g <- unique(ann[, .(gene_id, n_cpg_tested, n_cpg_eligible)])
gs <- ann[, .(any_present = any(in_uni), any_block = any(in_block), all_block = all(in_block), n_tx = .N), by = gene_id]
g <- merge(g, gs, by = "gene_id")
cat("\nGenes with annotated mRNAs:", nrow(g), "\n")
cat("  genes with >= 1 transcript in universe:", sum(g$any_present), "\n")
cat("  genes with no transcript in universe:", sum(!g$any_present), "; of those all transcripts in block:", sum(!g$any_present & g$all_block), "\n")
cat("  genes with transcripts both in universe and in block (split genes):", sum(g$any_present & g$any_block), "\n")
cat("  genes wholly lost (no transcript in universe) that hold a CpG with 10x in all 22 libraries:",
    sum(!g$any_present & g$n_cpg_tested > 0), "of", sum(!g$any_present), sprintf("(%.1f%%)", 100 * mean(g$n_cpg_tested[!g$any_present] > 0)), "\n")
cat("  genes with a transcript in universe that hold such a CpG:",
    sum(g$any_present & g$n_cpg_tested > 0), "of", sum(g$any_present), sprintf("(%.1f%%)", 100 * mean(g$n_cpg_tested[g$any_present] > 0)), "\n")
cat("  tested CpGs inside wholly lost genes (sum over genes; genes can overlap):", sum(g$n_cpg_tested[!g$any_present]), "\n")
cat("  eligible CpGs inside wholly lost genes:", sum(g$n_cpg_eligible[!g$any_present]), "\n")
fwrite(ann, file.path(out, "annotated_mrnas_with_gene_counts.csv"))

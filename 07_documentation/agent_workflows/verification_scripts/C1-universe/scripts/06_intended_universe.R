# What would the intended filter (gene overlaps a CpG with >= 1 read in any of the 22 libraries) have kept?
# c_g = number of union-1x CpGs per gene, taken from the legacy pipeline's own file union_1x-BOTH.GeneIDs (gannet).
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; out <- commandArgs(trailingOnly = TRUE)[1]
cg <- fread(file.path(out, "legacy_union_cpg_per_gene_BOTH.tsv"), header = FALSE, col.names = c("gene_id", "c_union"), colClasses = c("character", "integer"))
cat("legacy union file: genes with a covered CpG:", nrow(cg), "; CpG x gene overlaps:", sum(cg$c_union), "\n")
genes <- fread(file.path(out, "genes.tsv"), header = FALSE, col.names = c("gene_id", "chr", "gs", "ge")); genes[, gene_id := as.character(gene_id)]
cat("genes on track:", nrow(genes), "; without any union CpG:", sum(!genes$gene_id %in% cg$gene_id), "; union genes not on track:", sum(!cg$gene_id %in% genes$gene_id), "\n")
one <- fread(file.path(out, "gene_cov_zr3644_14.tsv"), header = FALSE, col.names = c("gene_id", "chr", "gs", "ge", "c_one", "c_one_5x")); one[, gene_id := as.character(gene_id)]
mr <- fread(file.path(out, "all_mrnas.csv"), colClasses = c(gene_id = "character"))
mr[, hyphen := grepl("-", parent, fixed = TRUE)]
# --- completeness of the table that fed the failing cell: sum over genes of (union CpGs) x (mRNAs with a parsable gene ID)
m_g <- mr[hyphen == FALSE, .(m = .N), by = gene_id]
x <- merge(cg, m_g, by = "gene_id")
cat(sprintf("predicted lines of union_1x-BOTH...transcriptIDs = sum(c_g * m_g) = %d ; notebook cell output: 26098781 ; difference %d\n",
            sum(as.numeric(x$c_union) * x$m), sum(as.numeric(x$c_union) * x$m) - 26098781))
# --- annotated mRNAs and the intended filter
ann <- fread(file.path(out, "annotated_mrnas_go.csv"), colClasses = c(gene_id = "character"))
ann <- merge(ann, cg, by = "gene_id", all.x = TRUE); ann[is.na(c_union), c_union := 0L]
ann <- merge(ann, one[, .(gene_id, c_one, c_one_5x)], by = "gene_id", all.x = TRUE)
ann[, status := fifelse(in_uni, "present in universe", fifelse(in_block, "missing: block", "missing: 49 others"))]
print(ann[, .(transcripts = .N, genes = uniqueN(gene_id),
              gene_has_union_cpg = sum(c_union > 0), pct = round(100 * mean(c_union > 0), 2),
              median_union_cpg = as.numeric(median(c_union)), min_union_cpg = min(c_union),
              gene_has_cpg_in_zr3644_14 = sum(c_one > 0), pct_one = round(100 * mean(c_one > 0), 2),
              median_cpg_zr3644_14 = as.numeric(median(c_one))), by = status])
cat("\nintended filter alone (gene has a union CpG) would keep", sum(ann$c_union > 0), "of", nrow(ann), "annotated mRNAs\n")
cat("...and with a parsable gene id:", sum(ann$c_union > 0 & !ann$hyphen), "\n")
cat("...and with a GO id in the recovered part of the 2017 GOslim table:", sum(ann$c_union > 0 & !ann$hyphen & ann$n_go_rec > 0), "\n")
cat("block transcripts whose gene has zero union CpGs:", ann[in_block == TRUE & c_union == 0, .N], "\n")
print(ann[in_block == TRUE & c_union == 0, .(transcript_id, gene_id, chr, start, end)])
cat("\nThe 49 others, by explanation:\n")
o <- ann[status == "missing: 49 others"]
o[, why := fifelse(hyphen, "gene name with hyphen: gene id lost in legacy parse",
            fifelse(c_union == 0, "gene has no CpG with a read in any library (the intended filter)",
            fifelse(n_go_rec == 0, "no GO id in the recovered 2017 GOslim table", "unexplained")))]
print(o[, .(transcripts = .N, genes = uniqueN(gene_id)), by = why])
print(o[why %in% c("unexplained", "gene has no CpG with a read in any library (the intended filter)"),
        .(transcript_id, gene_id, chr, start, end, c_union, c_one, n_go_rec, why)])
# --- 40 genes drawn at random from the genes that lost every annotated transcript to the block
gs <- ann[, .(all_block = all(in_block), n_tx = .N, c_union = c_union[1], c_one = c_one[1], c_one_5x = c_one_5x[1], chr = chr[1]), by = gene_id]
lost <- gs[all_block == TRUE]
set.seed(20261119); samp <- lost[sample(.N, 40)][order(chr, gene_id)]
cat("\n40 randomly drawn genes whose annotated transcripts are all in the block (seed 20261119):\n")
print(merge(samp, genes[, .(gene_id, gs, ge)], by = "gene_id")[, .(gene_id, chr, gs, ge, n_tx, cpg_zr3644_14 = c_one, cpg_5x_zr3644_14 = c_one_5x, cpg_union22 = c_union)])
cat("of the 40: with a covered CpG in zr3644_14:", sum(samp$c_one > 0), "; with a CpG at >= 5x:", sum(samp$c_one_5x > 0), "; with a union CpG:", sum(samp$c_union > 0), "\n")
cat("all", nrow(lost), "wholly lost genes: with a covered CpG in zr3644_14:", sum(lost$c_one > 0), "; zero:", sum(lost$c_one == 0), "; with union CpG:", sum(lost$c_union > 0), "\n")
fwrite(samp, file.path(out, "sample40_block_genes.csv"))
fwrite(ann, file.path(out, "annotated_mrnas_cov.csv"))

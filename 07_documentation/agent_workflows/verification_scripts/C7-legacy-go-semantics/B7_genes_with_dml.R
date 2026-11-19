# Independent count of genes holding a DML, and of those with GO terms (the reference quantities quoted in I-20).
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
gff <- fread(cmd = paste("zcat", file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_gene.gff.gz")), header = FALSE, sep = "\t", quote = "")
gff[, gene_id := sub(".*GeneID:([0-9]+).*", "\\1", V9)]
genes <- GRanges(gff$V1, IRanges(gff$V4, gff$V5), gene_id = gff$gene_id)
cat("gene records:", length(genes), "distinct gene IDs:", uniqueN(gff$gene_id), "\n")
mr <- fread(cmd = paste("zcat", file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_mRNA.gff.gz")), header = FALSE, sep = "\t", quote = "")
mr[, `:=`(gene_id = sub(".*GeneID:([0-9]+).*", "\\1", V9), tx = sub(".*transcript_id=([^;]+).*", "\\1", V9))]
mr[, num := sub("^[A-Z]{2}_", "", tx)]
gm <- fread(cmd = paste("zcat", file.path(repo, "03_input_files/annotation/geneid2go.tab.gz")), header = FALSE, sep = "\t", col.names = c("id", "go"), quote = "")
gm <- gm[grepl("GO:", go)]
gm[, num := sub("^XM_", "", id)]
# true prefix of each map name from the blastx query IDs: keep only names that are mRNAs (mrna_NM_ or mrna_XM_)
bx <- unique(fread(cmd = paste("zcat", file.path(repo, "03_input_files/annotation/20210605-cgigas-roslin-mito-blastx.outfmt6.gz")), header = FALSE, sep = "\t", select = 1)$V1)
mm <- regmatches(bx, regexec("_(mrna|miscrna|ncrna)_([A-Z]{2})_([0-9]+\\.[0-9]+)_[0-9]+$", bx)); mm <- do.call(rbind, mm[lengths(mm) == 4])
true_tx <- unique(paste0(mm[mm[, 2] == "mrna", 3], "_", mm[mm[, 2] == "mrna", 4]))
tx_with_go <- intersect(true_tx[sub("^[A-Z]{2}_", "", true_tx) %in% gm$num], mr$tx)
genes_with_go <- unique(mr[tx %chin% tx_with_go, gene_id])
legacy_tx_with_go <- intersect(unique(gm$id), mr$tx)                    # the legacy route: name must match as written (XM_ only)
cat("mRNAs", nrow(mr), "; with GO terms (true IDs)", length(tx_with_go), "; genes with GO terms", length(genes_with_go), "\n")
d <- file.path(repo, "08_archive/tables")
for (s in c("PLOIDY", "HEAT", "INTERACTION")) {
  x <- fread(file.path(d, sprintf("DML-getMethylDiff-%s-Cov10-20.csv", s)))
  dml <- GRanges(x$chr, IRanges(x$start + 1L, x$end))           # the lists are BED-like: start is 0-based
  ov <- findOverlaps(dml, genes)
  g <- unique(genes$gene_id[subjectHits(ov)])
  dml0 <- GRanges(x$chr, IRanges(x$start, x$end)); g0 <- unique(genes$gene_id[subjectHits(findOverlaps(dml0, genes))])
  cat(sprintf("%-11s DML %d; in a gene %d; genes holding a DML %d (treating start as 1-based: %d); of these with GO terms %d\n",
              s, nrow(x), length(unique(queryHits(ov))), length(g), length(g0), sum(g %chin% genes_with_go)))
}

suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t05 <- file.path(repo, "04_output/05_differential_methylation/tables"); fdir <- file.path(repo, "03_input_files/genome_features")
g <- fread(file.path(fdir, "cgigas_uk_roslin_v1_gene.gff.gz"), header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5, 7, 9)); setnames(g, c("chr", "s", "e", "strand", "a"))
l <- fread(file.path(fdir, "cgigas_uk_roslin_v1_lncRNA.gff.gz"), header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5, 7, 9)); setnames(l, c("chr", "s", "e", "strand", "a"))
reg <- GRanges("NC_047563.1", IRanges(39633324, 39635522))
gg <- GRanges(g$chr, IRanges(g$s, g$e)); lg <- GRanges(l$chr, IRanges(l$s, l$e))
cat("genes overlapping the region:\n"); print(g[queryHits(findOverlaps(gg, reg)), .(chr, s, e, strand, name = sub(".*Name=([^;]+).*", "\\1", a), biotype = sub(".*gene_biotype=([^;]+).*", "\\1", a))])
cat("lncRNA transcripts overlapping the region:\n"); print(l[queryHits(findOverlaps(lg, reg)), .(chr, s, e, strand, gene = sub(".*;gene=([^;]+).*", "\\1", a))])
for (mode in c("legacy", "corrected")) {
  s <- fread(file.path(t05, sprintf("dml_shared_%s.csv", mode))); u <- fread(file.path(t05, sprintf("united_cpg_%s.csv.gz", mode)), select = c("chr", "start"))
  sg <- GRanges(s$chr, IRanges(s$start + 1L, s$start + 2L)); inl <- overlapsAny(sg, lg); inr <- overlapsAny(sg, reg)
  cat(sprintf("[%s] shared loci %d; in lncRNA %d; in the region %d; in the region and in lncRNA %d; analysed CpGs in the region %d\n", mode, nrow(s), sum(inl), sum(inr), sum(inl & inr), u[chr == "NC_047563.1" & start >= 39633323 & start <= 39635520, .N]))
}

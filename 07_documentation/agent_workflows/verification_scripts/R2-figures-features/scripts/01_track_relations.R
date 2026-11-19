# Independent look at the 13 feature tracks: how they relate to each other.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
dir <- "/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/genome_features"
rd_gff <- function(n) { d <- fread(file.path(dir, sprintf("cgigas_uk_roslin_v1_%s.gff.gz", n)), header = FALSE, sep = "\t", quote = "", select = c(1,3,4,5,7,9), col.names = c("chr","type","start","end","strand","attr")); d }
rd_bed <- function(n) { d <- fread(file.path(dir, sprintf("cgigas_uk_roslin_v1_%s.bed.gz", n)), header = FALSE, sep = "\t", quote = ""); d }
gene <- rd_gff("gene"); mrna <- rd_gff("mRNA"); exon <- rd_gff("exon"); cds <- rd_gff("CDS"); utr <- rd_gff("exonUTR")
lnc <- rd_gff("lncRNA"); up <- rd_gff("upstream"); down <- rd_gff("downstream"); fl <- rd_gff("flanks")
intron <- rd_bed("intron"); noncds <- rd_bed("nonCDS"); inter <- rd_bed("intergenic"); te <- rd_bed("rm.te")
cat("BED column counts: intron", ncol(intron), "nonCDS", ncol(noncds), "intergenic", ncol(inter), "TE", ncol(te), "\n")
print(head(te, 3))
g1 <- function(d) reduce(GRanges(d$chr, IRanges(d$start, d$end)), ignore.strand = TRUE)           # GFF: 1-based closed
b0 <- function(d) { d <- d[V3 > V2]; reduce(GRanges(d$V1, IRanges(d$V2 + 1L, d$V3))) }           # BED read as 0-based half-open
b1 <- function(d) reduce(GRanges(d$V1, IRanges(d$V2, d$V3)))                                      # BED read (wrongly) as 1-based closed
bp <- function(x) sum(as.numeric(width(x)))
ovbp <- function(a, b) { sl <- union(seqlevels(a), seqlevels(b)); seqlevels(a) <- sl; seqlevels(b) <- sl; bp(intersect(a, b, ignore.strand = TRUE)) }
G <- list(gene = g1(gene), mRNA = g1(mrna), exon = g1(exon), CDS = g1(cds), UTR = g1(utr), lnc = g1(lnc),
          up = g1(up[end >= start]), down = g1(down[end >= start]), flanks = g1(fl[end >= start]))
B <- list(intron = b0(intron), nonCDS = b0(noncds), intergenic = b0(inter), TE = b0(te))
B1 <- list(intron = b1(intron), nonCDS = b1(noncds), intergenic = b1(inter), TE = b1(te))
cat("\nbp covered (reduced):\n"); print(sapply(c(G, B), bp))
cat("\n--- BED convention checks (bp of overlap; 0 expected when the reading is right)\n")
cat("intron(0-based) x exon :", ovbp(B$intron, G$exon), "   intron(as 1-based) x exon :", ovbp(B1$intron, G$exon), "\n")
cat("nonCDS(0-based) x exon :", ovbp(B$nonCDS, G$exon), "   nonCDS(as 1-based) x exon :", ovbp(B1$nonCDS, G$exon), "\n")
cat("nonCDS(0-based) x CDS  :", ovbp(B$nonCDS, G$CDS), "\n")
cat("intergenic(0-based) x gene :", ovbp(B$intergenic, G$gene), "   as 1-based:", ovbp(B1$intergenic, G$gene), "\n")
cat("intergenic(0-based) x flanks :", ovbp(B$intergenic, G$flanks), "   as 1-based:", ovbp(B1$intergenic, G$flanks), "\n")
# Is nonCDS the complement of exon, or of CDS?  genome lengths from the fasta index are not here; use per-sequence max end.
allr <- c(G$exon, B$nonCDS); cov_union <- reduce(allr)
cat("\nnonCDS + exon: union bp =", bp(cov_union), " sum of the two =", bp(G$exon) + bp(B$nonCDS), " (equal => disjoint and jointly tiling)\n")
cat("number of union ranges per sequence with exon+nonCDS (1 => gap-free tiling):\n"); print(table(table(as.character(seqnames(cov_union)))))
cu2 <- reduce(c(G$CDS, B$nonCDS)); cat("nonCDS + CDS: union bp =", bp(cu2), " ranges on sequences:", length(cu2), " n seq:", length(unique(as.character(seqnames(cu2)))), "\n")
cat("bp of exon not in CDS (UTR + non-coding exon):", bp(setdiff(G$exon, G$CDS, ignore.strand = TRUE)), "\n")
cat("bp of non-CDS exon covered by the nonCDS track:", ovbp(B$nonCDS, setdiff(G$exon, G$CDS, ignore.strand = TRUE)), "\n")
# exonUTR vs exon minus CDS
utr_expected <- setdiff(G$exon, G$CDS, ignore.strand = TRUE)
cat("\nexonUTR bp:", bp(G$UTR), " exon minus CDS bp:", bp(utr_expected), " overlap:", ovbp(G$UTR, utr_expected), " UTR x CDS bp:", ovbp(G$UTR, G$CDS), "\n")
# intron vs (gene minus exon), and (mRNA+lnc minus exon)
intr_expected <- setdiff(G$gene, G$exon, ignore.strand = TRUE)
cat("intron bp:", bp(B$intron), " gene minus exon bp:", bp(intr_expected), " overlap:", ovbp(B$intron, intr_expected), "\n")
cat("intron x lncRNA transcript bp:", ovbp(B$intron, G$lnc), "   UTR track x lncRNA transcripts bp:", ovbp(G$UTR, G$lnc), "   lncRNA bp:", bp(G$lnc), "\n")
# flanks = upstream + downstream?
cat("flanks bp:", bp(G$flanks), " up U down bp:", bp(reduce(c(G$up, G$down))), " up x down bp:", ovbp(G$up, G$down), "\n")
cat("flanks x gene bp:", ovbp(G$flanks, G$gene), "  up x gene:", ovbp(G$up, G$gene), "  down x gene:", ovbp(G$down, G$gene), "\n")
# intergenic = complement(gene + flanks)?
gf <- reduce(c(G$gene, G$flanks)); tile <- reduce(c(gf, B$intergenic))
cat("intergenic + gene + flanks: union bp =", bp(tile), " sum =", bp(gf) + bp(B$intergenic), " ranges:", length(tile), " sequences:", length(unique(as.character(seqnames(tile)))), "\n")
# TE overlaps
cat("TE bp:", bp(B$TE), " TE x gene bp:", ovbp(B$TE, G$gene), " TE x CDS bp:", ovbp(B$TE, G$CDS), " TE x intergenic bp:", ovbp(B$TE, B$intergenic), "\n")
cat("TE zero-length:", te[V3 <= V2, .N], " intron zero-length:", intron[V3 <= V2, .N], " nonCDS zero-length:", noncds[V3 <= V2, .N], " intergenic zero-length:", inter[V3 <= V2, .N], "\n")
# strand awareness of flanks: match each flank to its gene by ID
gid <- function(a) sub(";.*", "", sub("^ID=", "", a))
gene[, id := gid(attr)]; up[, id := gid(attr)]; down[, id := gid(attr)]
u <- merge(up[end >= start], gene, by = c("chr", "id"), suffixes = c("_f", "_g"))
d <- merge(down[end >= start], gene, by = c("chr", "id"), suffixes = c("_f", "_g"))
cat("\nupstream flank position relative to its gene, by gene strand:\n")
print(u[, .N, by = .(strand_g, where = fifelse(end_f < start_g, "before gene start (lower coord)", fifelse(start_f > end_g, "after gene end (higher coord)", "overlaps gene")))])
cat("downstream flank position relative to its gene, by gene strand:\n")
print(d[, .N, by = .(strand_g, where = fifelse(end_f < start_g, "before gene start (lower coord)", fifelse(start_f > end_g, "after gene end (higher coord)", "overlaps gene")))])
cat("upstream flanks adjacent to gene (+: end_f == start_g - 1; -: start_f == end_g + 1):", u[(strand_g == "+" & end_f == start_g - 1) | (strand_g == "-" & start_f == end_g + 1), .N], "of", nrow(u), "\n")
cat("genes:", nrow(gene), " genes with an upstream record:", uniqueN(up$id), " with a downstream record:", uniqueN(down$id), "\n")
cat("sequences: gene", uniqueN(gene$chr), " TE", uniqueN(te$V1), " intergenic", uniqueN(inter$V1), " nonCDS", uniqueN(noncds$V1), "\n")
cat("mito (NC_001276.1) records: gene", gene[chr == "NC_001276.1", .N], " exon", exon[chr == "NC_001276.1", .N], " CDS", cds[chr == "NC_001276.1", .N], " TE", te[V1 == "NC_001276.1", .N], " intergenic", inter[V1 == "NC_001276.1", .N], " nonCDS", noncds[V1 == "NC_001276.1", .N], "\n")

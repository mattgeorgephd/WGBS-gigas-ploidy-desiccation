# genomic_overlaps.R
# Genome feature tracks and overlap counts, replacing the bedtools steps of the
# legacy notebooks.
#
# Coordinates inside R are 1-based closed intervals (GenomicRanges). GFF files
# already use that convention; BED files and Bismark coverage files use 0-based
# half-open intervals and are converted on read.

# Names of the 13 feature tracks in 03_input_files/genome_features/.
feature_track_files <- function() {
  c(gene = "gene.gff", mRNA = "mRNA.gff", exon = "exon.gff", CDS = "CDS.gff", exon_UTR = "exonUTR.gff",
    intron = "intron.bed", nonCDS = "nonCDS.bed", lncRNA = "lncRNA.gff", upstream_flank = "upstream.gff",
    downstream_flank = "downstream.gff", flanks = "flanks.gff", intergenic = "intergenic.bed", TE = "rm.te.bed")
}

# Read one track as a GRanges object. `keep_attributes = TRUE` keeps the GFF
# attribute column (needed to pull gene and transcript IDs).
read_feature_track <- function(track, keep_attributes = FALSE) {
  file <- feature_track_files()[[track]]
  path <- input_path("genome_features", paste0("cgigas_uk_roslin_v1_", file, ".gz"))
  if (grepl("\\.gff$", file)) {
    d <- data.table::fread(path, header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5, 7, 9),
                           col.names = c("chr", "start", "end", "strand", "attributes"), showProgress = FALSE)
  } else {
    d <- data.table::fread(path, header = FALSE, sep = "\t", quote = "", select = 1:3,
                           col.names = c("chr", "start", "end"), showProgress = FALSE)
    d[, start := start + 1L]                      # BED start is 0-based
    d[, `:=`(strand = "*", attributes = NA_character_)]
  }
  d <- d[end >= start]                            # a few flank records have zero length
  gr <- GenomicRanges::GRanges(d$chr, IRanges::IRanges(d$start, d$end),
                               strand = ifelse(d$strand %in% c("+", "-"), d$strand, "*"))
  if (keep_attributes) S4Vectors::mcols(gr)$attributes <- d$attributes
  gr
}

# CpG positions (0-based start of the C, as in the coverage files and all
# pipeline tables) as two-base GRanges covering the C and the G.
cpg_ranges <- function(chr, start) {
  GenomicRanges::GRanges(chr, IRanges::IRanges(start + 1L, start + 2L))
}

# TRUE for each CpG that overlaps at least one feature by at least one base,
# regardless of strand. Equivalent to `bedtools intersect -u -a cpg -b features`.
overlaps_track <- function(cpg, track_gr) {
  suppressWarnings(IRanges::overlapsAny(cpg, track_gr, ignore.strand = TRUE))
}

# Extract one field from GFF attribute strings, for example
# gff_attribute(x, "GeneID:") or gff_attribute(x, "gene_biotype=").
gff_attribute <- function(attributes, key) {
  m <- regexpr(paste0(key, "[^;,]+"), attributes)
  out <- rep(NA_character_, length(attributes))
  out[m > 0] <- sub(key, "", regmatches(attributes, m), fixed = TRUE)
  out
}

# All (CpG, feature) pairs that overlap by at least one base, regardless of
# strand, as a two-column data.table of row indices. Equivalent to
# `bedtools intersect -wb -a cpg -b features`.
overlap_pairs <- function(cpg, track_gr) {
  h <- suppressWarnings(GenomicRanges::findOverlaps(cpg, track_gr, ignore.strand = TRUE))
  data.table::data.table(cpg = S4Vectors::queryHits(h), feature = S4Vectors::subjectHits(h))
}

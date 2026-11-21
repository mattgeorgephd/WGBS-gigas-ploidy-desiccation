# annotation.R
# Gene, transcript, and GO annotation tables, replacing the cut/tr/join steps of
# the legacy annotation notebook.

# Genes of the gene track: one row per gene record.
read_genes <- function() {
  gr <- read_feature_track("gene", keep_attributes = TRUE)
  a <- S4Vectors::mcols(gr)$attributes
  data.table::data.table(
    gene_id = gff_attribute(a, "GeneID:"), gene = gff_attribute(a, "Name="),
    biotype = gff_attribute(a, "gene_biotype="),
    chr = as.character(GenomicRanges::seqnames(gr)), start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr), strand = as.character(GenomicRanges::strand(gr)))
}

# Transcripts of the mRNA track: one row per transcript, with its gene and product.
#
# legacy_gene_id_lost marks the transcripts that the legacy notebook could not
# join to a gene. It split the attribute string on = ; : , and - and took the
# ninth piece as the gene ID, which is the gene ID only when the gene name holds
# no hyphen (gene-LOC105342885). For names such as gene-Ankyrin-3 the ninth
# piece is the word "GeneID".
read_transcripts <- function() {
  gr <- read_feature_track("mRNA", keep_attributes = TRUE)
  a <- S4Vectors::mcols(gr)$attributes
  product <- gff_attribute(a, "product=")
  decoded <- vapply(unique(product), utils::URLdecode, character(1))
  data.table::data.table(
    transcript_id = gff_attribute(a, "transcript_id="), gene_id = gff_attribute(a, "GeneID:"),
    gene = gff_attribute(a, "gene="), product = unname(decoded[product]),
    chr = as.character(GenomicRanges::seqnames(gr)), start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr), strand = as.character(GenomicRanges::strand(gr)),
    legacy_gene_id_lost = grepl("-", gff_attribute(a, "Parent=gene-"), fixed = TRUE))
}

# Product name of a gene: its transcripts' product names without the
# "transcript variant" suffix, joined when they differ.
gene_products <- function(transcripts) {
  p <- unique(transcripts[, .(gene_id, product = sub(", transcript variant .*$", "", product))])
  p[, .(product = paste(sort(product), collapse = "; "), n_transcripts = .N), by = gene_id][
    transcripts[, .(n_transcripts = .N), by = gene_id], on = "gene_id",
    .(gene_id, product, n_transcripts = i.n_transcripts)]
}

# A frozen transcript-to-GO map in long form: one row per (map_id, go_id).
#   "full"             every transcript with a blastx hit that carries GO terms
#                      (gannet data/geneid2go.tab, June 2021)
#   "legacy_universe"  the subset the legacy topGO analysis used as its gene
#                      universe; a block of transcripts is missing from it
#                      (issue I-39)
# map_id is the transcript name as the legacy map spells it: always XM_<number>,
# whatever the true prefix (see true_transcript_ids()).
read_go_map <- function(which = c("full", "legacy_universe")) {
  which <- match.arg(which)
  f <- if (which == "full") "geneid2go.tab.gz" else "geneid2go-union1x-BOTH.tab.gz"
  d <- unique(data.table::fread(input_path("annotation", f), header = FALSE, sep = "\t",
                                col.names = c("map_id", "go"), showProgress = FALSE))
  long <- d[, .(go_id = trimws(unlist(strsplit(go, ",", fixed = TRUE)))), by = map_id]
  unique(long[grepl("^GO:[0-9]{7}$", go_id)])
}

# True transcript IDs behind the names of the legacy GO map.
# The legacy notebook rebuilt each name from the blastx query ID by keeping the
# number and writing XM_ in front, so curated mRNAs (NM_) and non-coding RNAs
# (XR_) were renamed XM_ and could no longer be matched to the annotation. The
# blastx table still holds the true IDs. Returns map_id, transcript_id, rna_class.
true_transcript_ids <- function() {
  q <- unique(data.table::fread(input_path("annotation", "20210605-cgigas-roslin-mito-blastx.outfmt6.gz"),
                                header = FALSE, sep = "\t", select = 1, col.names = "query", showProgress = FALSE)$query)
  m <- regmatches(q, regexec("_(mrna|miscrna|ncrna)_([A-Z]{2})_([0-9]+\\.[0-9]+)_[0-9]+$", q))
  m <- do.call(rbind, m[lengths(m) == 4])
  out <- unique(data.table::data.table(map_id = paste0("XM_", m[, 4]), transcript_id = paste0(m[, 3], "_", m[, 4]),
                                       rna_class = m[, 2]))
  if (anyDuplicated(out$map_id)) stop("A transcript number occurs under two prefixes")
  out
}

# Legacy GO term table recovered from the legacy annotation tables: GOID, GOterm,
# GOSlim, GOcat (see 03_input_files/annotation/build_legacy_go_terms.R). One row
# per GO ID and GOslim bin: a GO ID filed under several bins has several rows, and
# so has every line of a legacy annotation table that names it.
read_legacy_go_terms <- function() {
  data.table::fread(input_path("annotation", "legacy_go_goslim_terms.tsv.gz"), header = TRUE, sep = "\t",
                    quote = "", na.strings = NULL, showProgress = FALSE)
}

# Read a ten-column DML annotation table (transcript, geneID, chr, start, end,
# meth.diff, GOID, GOterm, GOSlim, GOcat).
#   reader = "exact"   every line becomes a row
#   reader = "legacy"  as the legacy script read it: read.table() with its
#                      default quote characters. A line whose term name holds an
#                      unpaired apostrophe (RNA 3'-end processing) opens or
#                      closes a quoted string, so every line in between
#                      disappears into one field (issue I-40).
dml_annotation_columns <- function() {
  c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
}
read_dml_annotation <- function(path, reader = c("exact", "legacy")) {
  reader <- match.arg(reader)
  # A set without DML in an annotated gene has a table without lines, which
  # neither reader takes: it is returned as a table with the ten columns and no row.
  if (length(readLines(path, n = 1, warn = FALSE)) == 0) {
    return(data.table::data.table(transcript = character(), geneID = integer(), chr = character(), start = integer(),
                                  end = integer(), meth.diff = numeric(), GOID = character(), GOterm = character(),
                                  GOSlim = character(), GOcat = character()))
  }
  if (reader == "legacy") {
    return(data.table::as.data.table(utils::read.table(path, header = FALSE, sep = "\t",
                                                        col.names = dml_annotation_columns())))
  }
  data.table::fread(path, header = FALSE, sep = "\t", quote = "", na.strings = NULL,
                    col.names = dml_annotation_columns(), showProgress = FALSE)
}

# The rows of annotation_summary.csv of step 08 for the analysis `a` (its mode
# column holds the analysis id), one per DML set. A step that reads the
# ten-column tables of step 08 calls this first: without the rows, step 08 has
# not run for the analysis, and the step stops here. `sets` holds the names of
# the DML sets (ploidy, heat, shared).
legacy_annotation_summary <- function(sets, a = config$analysis) {
  rows <- data.table::fread(analysis_input(a, "08_gene_annotation", "tables", "annotation_summary.csv"))
  rows <- rows[rows$mode == a$id, ]
  if (!identical(sort(rows$set), sort(sets))) {
    stop("annotation_summary.csv of step 08 must hold one row of the analysis '", a$id, "' for each DML set (",
         paste(sets, collapse = ", "), ") and holds ",
         if (nrow(rows) == 0) "none" else paste("rows for", paste(rows$set, collapse = ", ")),
         ". Run step 08 for this analysis first.", call. = FALSE)
  }
  rows
}

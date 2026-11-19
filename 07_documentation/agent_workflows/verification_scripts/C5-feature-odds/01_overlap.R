# Independent recomputation of CpG x feature overlaps (legacy mode).
# Does not source anything from the repository.
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
tabdir <- file.path(repo, "04_output/05_differential_methylation/tables")
fdir   <- file.path(repo, "03_input_files/genome_features")

## ---- CpG table -------------------------------------------------------------
u <- fread(file.path(tabdir, "united_cpg_legacy.csv.gz"))
cat("united rows:", nrow(u), " cols:", ncol(u), "\n")
stopifnot(all(u$end - u$start == 2L))
covm <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE])
ncm  <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE])
ntm  <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
cat("coverage == numCs + numTs everywhere:", all(covm == ncm + ntm), "\n")
cat("min coverage:", min(covm), " any NA:", anyNA(covm), "\n")
pm <- 100 * ncm / covm
cpg <- data.table(
  chr = u$chr, start = u$start, end = u$end, row = seq_len(nrow(u)),
  mean_unw    = rowMeans(pm),                       # mean of per-library percentages
  mean_pooled = 100 * rowSums(ncm) / rowSums(covm), # pooled over libraries
  mean_cov    = rowMeans(covm),
  min_cov     = apply(covm, 1, min),
  max_cov     = apply(covm, 1, max)
)
rm(u, covm, ncm, ntm, pm); invisible(gc())

## ---- DML lists -------------------------------------------------------------
dp <- fread(file.path(tabdir, "dml_ploidy_legacy.csv"))
dh <- fread(file.path(tabdir, "dml_heat_legacy.csv"))
ds <- fread(file.path(tabdir, "dml_shared_legacy.csv"))
cat("ploidy DML:", nrow(dp), " heat DML:", nrow(dh), " shared file:", nrow(ds), "\n")
key <- function(d) paste(d$chr, d$start)
cpg[, k := paste(chr, start)]
stopifnot(!anyDuplicated(cpg$k))
cpg[, dml_ploidy := k %in% key(dp)]
cpg[, dml_heat   := k %in% key(dh)]
cpg[, dml_shared := dml_ploidy & dml_heat]
cat("matched ploidy:", sum(cpg$dml_ploidy), " heat:", sum(cpg$dml_heat),
    " in both:", sum(cpg$dml_shared), " shared file all found:", all(key(ds) %in% cpg$k[cpg$dml_shared]), "\n")
# direction of shared loci, from the two lists themselves
mp <- dp$meth.diff[match(cpg$k, key(dp))]
mh <- dh$meth.diff[match(cpg$k, key(dh))]
cpg[, md_ploidy := mp][, md_heat := mh]
cpg[, shared_class := fifelse(!dml_shared, NA_character_,
                       fifelse(md_ploidy > 0 & md_heat > 0, "both hyper",
                        fifelse(md_ploidy < 0 & md_heat < 0, "both hypo", "opposite")))]
print(table(cpg$shared_class, useNA = "ifany"))

## ---- CpG ranges: 0-based start, 2 bases -> 1-based closed [start+1, end] ----
gr <- GRanges(cpg$chr, IRanges(cpg$start + 1L, cpg$end))
stopifnot(all(width(gr) == 2L))

## ---- feature tracks --------------------------------------------------------
read_gff <- function(f) {
  d <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = FALSE, quote = "",
             select = c(1, 3, 4, 5), col.names = c("chr", "type", "s", "e"))
  n0 <- nrow(d)
  bad <- d$e < d$s
  d <- d[!bad]
  list(gr = GRanges(d$chr, IRanges(d$s, d$e)), n = n0, dropped = sum(bad))
}
read_bed <- function(f) {
  d <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = FALSE, quote = "",
             select = 1:3, col.names = c("chr", "s", "e"))
  n0 <- nrow(d)
  bad <- d$e <= d$s            # zero or negative length in half-open coordinates
  d <- d[!bad]
  list(gr = GRanges(d$chr, IRanges(d$s + 1L, d$e)), n = n0, dropped = sum(bad))
}
tracks <- c(
  gene = "cgigas_uk_roslin_v1_gene.gff.gz",
  mRNA = "cgigas_uk_roslin_v1_mRNA.gff.gz",
  exon = "cgigas_uk_roslin_v1_exon.gff.gz",
  CDS  = "cgigas_uk_roslin_v1_CDS.gff.gz",
  exon_UTR = "cgigas_uk_roslin_v1_exonUTR.gff.gz",
  intron = "cgigas_uk_roslin_v1_intron.bed.gz",
  nonCDS = "cgigas_uk_roslin_v1_nonCDS.bed.gz",
  lncRNA = "cgigas_uk_roslin_v1_lncRNA.gff.gz",
  upstream_flank = "cgigas_uk_roslin_v1_upstream.gff.gz",
  downstream_flank = "cgigas_uk_roslin_v1_downstream.gff.gz",
  flanks = "cgigas_uk_roslin_v1_flanks.gff.gz",
  intergenic = "cgigas_uk_roslin_v1_intergenic.bed.gz",
  TE = "cgigas_uk_roslin_v1_rm.te.bed.gz"
)
track_info <- list()
for (tn in names(tracks)) {
  f <- file.path(fdir, tracks[[tn]])
  x <- if (grepl("\\.bed\\.gz$", f)) read_bed(f) else read_gff(f)
  hit <- overlapsAny(gr, x$gr, ignore.strand = TRUE)
  set(cpg, j = paste0("in_", tn), value = hit)
  track_info[[tn]] <- data.table(track = tn, records = x$n, dropped_zero_len = x$dropped,
                                 bp_covered = sum(as.numeric(width(reduce(x$gr)))),
                                 cpg_hits = sum(hit))
  cat(sprintf("%-18s records %7d dropped %4d  CpG hits %7d\n", tn, x$n, x$dropped, sum(hit)))
}
track_info <- rbindlist(track_info)
fwrite(track_info, file.path(scr, "track_info.csv"))
saveRDS(cpg, file.path(scr, "cpg_features.rds"))
cat("saved\n")

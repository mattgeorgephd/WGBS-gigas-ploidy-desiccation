# Same question as script 07, this time with the pipeline's own helpers (run from the scratch copy of the repository).
suppressPackageStartupMessages({library(here); library(data.table); library(GenomicRanges)})
source(here("run_config.R")); invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
te <- read_feature_track("TE")                                   # pipeline reading: BED start + 1
te_true <- GenomicRanges::shift(te, 0); start(te_true) <- start(te) - 1L   # RepeatMasker begin is already 1-based
cat("first TE record as read by the pipeline:", as.character(seqnames(te)[1]), start(te)[1], end(te)[1], " width", width(te)[1], "(RepeatMasker: 4109-4128, 20 bases)\n")
for (mode in config$modes) {
  u <- fread(output_path("05_differential_methylation", "tables", sprintf("united_cpg_%s.csv.gz", mode)), select = c("chr", "start"))
  gr <- cpg_ranges(u$chr, u$start)
  a <- overlaps_track(gr, te); b <- overlaps_track(gr, te_true)
  dp <- fread(output_path("05_differential_methylation", "tables", sprintf("dml_ploidy_%s.csv", mode))); dh <- fread(output_path("05_differential_methylation", "tables", sprintf("dml_heat_%s.csv", mode)))
  cat(sprintf("[%s] analysed CpGs in TE: pipeline %d, with 1-based starts %d (gained %d, lost %d); ploidy DML %d -> %d; heat DML %d -> %d\n", mode, sum(a), sum(b), sum(b & !a), sum(a & !b),
              sum(a[dp$united_row]), sum(b[dp$united_row]), sum(a[dh$united_row]), sum(b[dh$united_row])))
  w <- which(b & !a)[1:3]; print(data.table(u[w], G_position_1based = u$start[w] + 2L))
  hit <- findOverlaps(gr[w], te_true); print(data.table(cpg = queryHits(hit), te_start = start(te_true)[subjectHits(hit)], te_end = end(te_true)[subjectHits(hit)]))
}

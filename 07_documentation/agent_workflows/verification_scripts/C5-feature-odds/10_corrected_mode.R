# Same independent computation on the corrected-mode tables (side check; the claim is about legacy mode).
suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tabdir <- file.path(repo, "04_output/05_differential_methylation/tables"); fdir <- file.path(repo, "03_input_files/genome_features")
u <- fread(file.path(tabdir, "united_cpg_corrected.csv.gz"))
nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
cpg <- data.table(chr = u$chr, start = u$start, mean_ct = rowMeans(100 * nc / (nc + nt)))
rm(u, nc, nt); invisible(gc())
cpg[, k := paste(chr, start)]
for (s in c("ploidy", "heat")) { d <- fread(file.path(tabdir, sprintf("dml_%s_corrected.csv", s))); set(cpg, j = paste0("dml_", s), value = cpg$k %in% paste(d$chr, d$start)) }
cpg[, elig := mean_ct >= 10 & mean_ct <= 90]
cat("corrected mode: tested", nrow(cpg), " eligible", sum(cpg$elig), " ploidy DML", sum(cpg$dml_ploidy), " heat DML", sum(cpg$dml_heat),
    " DML outside eligible:", sum(cpg$dml_ploidy & !cpg$elig), sum(cpg$dml_heat & !cpg$elig), "\n")
gr <- GRanges(cpg$chr, IRanges(cpg$start + 1L, cpg$start + 2L))
rd <- function(file, bed) {
  d <- fread(cmd = paste("zcat", shQuote(file.path(fdir, file))), sep = "\t", header = FALSE, quote = "", select = if (bed) 1:3 else c(1, 4, 5))
  setnames(d, c("chr", "s", "e")); if (bed) d[, s := s + 1L]; d <- d[e >= s]
  suppressWarnings(overlapsAny(gr, GRanges(d$chr, IRanges(d$s, d$e))))
}
tr <- list(gene = c("cgigas_uk_roslin_v1_gene.gff.gz", FALSE), intron = c("cgigas_uk_roslin_v1_intron.bed.gz", TRUE), CDS = c("cgigas_uk_roslin_v1_CDS.gff.gz", FALSE),
           TE = c("cgigas_uk_roslin_v1_rm.te.bed.gz", TRUE), intergenic = c("cgigas_uk_roslin_v1_intergenic.bed.gz", TRUE))
pa <- fread(file.path(repo, "04_output/07_genome_features/tables/feature_association_tests_corrected.csv"))
for (t in names(tr)) {
  f <- rd(tr[[t]][1], as.logical(tr[[t]][2]))
  for (s in c("ploidy", "heat")) for (b in c("all analysed CpGs", "eligible CpGs")) {
    dml <- cpg[[paste0("dml_", s)]]; bgv <- if (b == "eligible CpGs") cpg$elig else rep(TRUE, nrow(cpg))
    ft <- fisher.test(matrix(c(sum(dml & f), sum(bgv & !dml & f), sum(dml & !f), sum(bgv & !dml & !f)), 2))
    p <- pa[background == b & set == paste(s, "DML") & track == t]
    cat(sprintf("%-10s %-6s %-18s mine %.4f (%.2f-%.2f)  pipeline %.4f  equal %s\n", t, s, b, ft$estimate, ft$conf.int[1], ft$conf.int[2], p$odds_ratio,
                isTRUE(all.equal(unname(ft$estimate), p$odds_ratio))))
  }
}

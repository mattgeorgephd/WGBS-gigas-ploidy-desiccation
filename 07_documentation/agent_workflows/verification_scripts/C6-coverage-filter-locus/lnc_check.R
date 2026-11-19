suppressPackageStartupMessages(library(data.table))
options(width = 220)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
gf <- function(f) file.path(repo, "03_input_files/genome_features", f)
gene <- fread(cmd = paste("zcat", gf("cgigas_uk_roslin_v1_gene.gff.gz")), header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5, 9), col.names = c("chr", "s", "e", "attr"))
gene[, biotype := sub(".*gene_biotype=([^;]+).*", "\\1", attr)][, name := sub(".*Name=([^;]+).*", "\\1", attr)]
cat("gene records:", nrow(gene), "; biotypes:\n"); print(table(gene$biotype))
lnc <- gene[biotype == "lncRNA", .(chr, s, e, name)]; setkey(lnc, chr, s, e)
lnctrack <- fread(cmd = paste("zcat", gf("cgigas_uk_roslin_v1_lncRNA.gff.gz")), header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5), col.names = c("chr", "s", "e")); setkey(lnctrack, chr, s, e)
te <- fread(cmd = paste("zcat", gf("cgigas_uk_roslin_v1_rm.te.bed.gz")), header = FALSE, sep = "\t", select = 1:3, col.names = c("chr", "s", "e")); te[, s := s + 1L]; setkey(te, chr, s, e)
hit <- function(d, track) { q <- d[, .(chr, s = start + 1L, e = start + 2L, i = .I)]; setkey(q, chr, s, e)
  o <- foverlaps(q, track, type = "any", nomatch = NULL, which = TRUE); seq_len(nrow(d)) %in% q$i[o$xid] }
L0 <- 39633323L; L1 <- 39635520L
for (mode in c("legacy", "corrected")) {
  for (set in c("ploidy", "heat", "shared")) {
    d <- fread(t5(sprintf("dml_%s_%s.csv", set, mode)))
    h <- hit(d, lnc); h2 <- hit(d, lnctrack)
    atloc <- d$chr == "NC_047563.1" & d$start >= L0 & d$start <= L1
    cat(sprintf("[%s %s] DML %d; in lncRNA genes (gene.gff, biotype lncRNA): %d; in lncRNA transcript track: %d; at the locus: %d; lncRNA DML at the locus: %d; lncRNA DML elsewhere: %d\n",
                mode, set, nrow(d), sum(h), sum(h2), sum(atloc), sum(h & atloc), sum(h & !atloc)))
  }
}
# genes overlapping the locus span
cat("\ngenes overlapping NC_047563.1:", L0 + 1, "-", L1 + 2, ":\n"); print(gene[chr == "NC_047563.1" & s <= L1 + 2 & e >= L0 + 1, .(name, biotype, s, e)])
# ---- characterise the CpGs the filter removes ----
z <- readRDS("cand_matrix.rds"); cand <- z$cand; k2 <- z$k2; nlo <- z$nlib_over
cand[, `:=`(removed = !k2, nlib_over = nlo, te = hit(cand, te))]
cat("\nshare of CpGs inside an annotated transposable element (RepeatMasker track): kept", sprintf("%.1f%%", 100 * cand[removed == FALSE, mean(te)]),
    "; removed", sprintf("%.1f%%", 100 * cand[removed == TRUE, mean(te)]), "\n")
print(cand[removed == TRUE, .(n = .N, pct_in_TE = round(100 * mean(te), 1), mito = sum(chr == "NC_001276.1"), unplaced = sum(grepl("^NW_", chr))),
           by = .(libraries_over = cut(nlib_over, c(0, 1, 2, 5, 10, 21, 22), labels = c("1", "2", "3-5", "6-10", "11-21", "22")))][order(libraries_over)])
cat("sequence class of removed vs kept:\n")
cls <- function(chr) fifelse(chr == "NC_001276.1", "mitochondrial", fifelse(grepl("^NW_", chr), "unplaced scaffold", "linkage group"))
print(dcast(cand[, .N, by = .(class = cls(chr), removed)], class ~ removed, value.var = "N", fill = 0)[, pct_removed := round(100 * `TRUE` / (`TRUE` + `FALSE`), 1)][])
# clustering of removed CpGs: runs within 1 kb
r <- cand[removed == TRUE][order(chr, start)]
r[, g := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(start) > 1000))]
rs <- r[, .(chr = chr[1], from = min(start), to = max(start), n = .N, med_nlib = median(nlib_over)), by = g]
cat("\nremoved CpGs:", nrow(r), "in", nrow(rs), "runs (neighbours within 1 kb); isolated singletons:", sum(rs$n == 1), "; in runs of >= 10:", rs[n >= 10, sum(n)], sprintf("(%.1f%%)", 100 * rs[n >= 10, sum(n)] / nrow(r)), "\n")
cat("largest runs of removed CpGs:\n"); print(rs[order(-n)][1:12, .(chr, from, to, span = to - from, n, med_nlib)])
cat("rank of the I-09 locus among runs by number of removed CpGs:", which(rs[order(-n)]$chr == "NC_047563.1" & rs[order(-n)]$from <= L1 & rs[order(-n)]$to >= L0), "\n")
print(rs[chr == "NC_047563.1" & from <= L1 + 2000 & to >= L0 - 2000])

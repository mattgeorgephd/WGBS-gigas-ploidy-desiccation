# Independent recomputation of the I-39 numbers from raw inputs (not from step 08/09 tables).
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
out  <- commandArgs(trailingOnly = TRUE)[1]

# ---- mRNA track, parsed directly from the GFF --------------------------------
gff <- fread(cmd = sprintf("zcat %s/03_input_files/genome_features/cgigas_uk_roslin_v1_mRNA.gff.gz", repo),
             header = FALSE, sep = "\t", quote = "", col.names = c("chr","src","type","start","end","score","strand","phase","attr"))
cat("mRNA records:", nrow(gff), "\n")
gff[, transcript_id := sub(".*transcript_id=([^;]+).*", "\\1", attr)]
gff[, gene_id := sub(".*GeneID:([0-9]+).*", "\\1", attr)]
gff[, parent := sub(".*Parent=gene-([^;]+).*", "\\1", attr)]
stopifnot(!anyDuplicated(gff$transcript_id))
cat("prefix table of mRNA transcript ids:\n"); print(table(substr(gff$transcript_id, 1, 3)))
cat("digits in accession number:\n"); print(table(nchar(sub("^[A-Z]{2}_([0-9]+)\\..*", "\\1", gff$transcript_id))))
cat("version digits:\n"); print(table(sub(".*\\.", "", gff$transcript_id)))

# ---- full map and universe map ------------------------------------------------
full <- fread(cmd = sprintf("zcat %s/03_input_files/annotation/geneid2go.tab.gz", repo), header = FALSE, sep = "\t", quote = "",
              col.names = c("id", "go"))
uni  <- fread(cmd = sprintf("zcat %s/03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz", repo), header = FALSE, sep = "\t", quote = "",
              col.names = c("id", "go"))
cat("full map lines:", nrow(full), " names:", uniqueN(full$id), " names with a GO id:", uniqueN(full$id[grepl("GO:[0-9]{7}", full$go)]), "\n")
cat("universe lines:", nrow(uni), " names:", uniqueN(uni$id), " names with a GO id:", uniqueN(uni$id[grepl("GO:[0-9]{7}", uni$go)]), "\n")
cat("universe lines without any GO id:", sum(!grepl("GO:[0-9]{7}", uni$go)), "\n")
cat("universe names not in full map:", sum(!unique(uni$id) %in% full$id), "\n")
# is the universe map exactly the full-map lines of its ids (the awk filter of the notebook)?
chk <- full[id %in% uni$id]
cat("full-map lines whose id is in the universe:", nrow(chk), " identical to universe file (same order):",
    identical(chk$id, uni$id) && identical(chk$go, uni$go), "\n")

with_go <- unique(full$id[grepl("GO:[0-9]{7}", full$go)])
gff[, in_full := transcript_id %in% with_go]
gff[, in_full_any := transcript_id %in% full$id]
gff[, in_uni := transcript_id %in% uni$id]
cat("mRNAs on track with GO terms in full map (annotated):", sum(gff$in_full), "\n")
cat("mRNAs on track named in full map at all (with or without GO):", sum(gff$in_full_any), "\n")
cat("universe names on mRNA track:", sum(unique(uni$id) %in% gff$transcript_id), "of", uniqueN(uni$id), "\n")
cat("universe names that are annotated mRNAs:", sum(gff$in_uni & gff$in_full), "\n")
ann <- gff[in_full == TRUE]
cat("annotated mRNAs missing from the universe:", sum(!ann$in_uni), "of", nrow(ann), sprintf("(%.2f%%)", 100*mean(!ann$in_uni)), "\n")
cat("annotated mRNAs by prefix:\n"); print(table(substr(ann$transcript_id,1,3)))

runs_of <- function(d, label) {
  r <- rle(d$in_uni)
  miss <- r$lengths[!r$values]
  cat(sprintf("%-38s runs of missing: %5d; longest %5d; missing outside longest %5d; runs of present: %d\n", label,
              length(miss), max(miss), sum(miss) - max(miss), sum(r$values)))
  invisible(r)
}
# ---- orderings ---------------------------------------------------------------
ann[, acc_num := as.numeric(sub("^[A-Z]{2}_([0-9]+)\\..*", "\\1", transcript_id))]
ann[, ver := as.integer(sub(".*\\.", "", transcript_id))]
o_byte <- ann[base::order(ann$transcript_id, method = "radix")]            # C-locale byte order
r <- runs_of(o_byte, "byte order of the ID")
i <- which.max(ifelse(r$values, 0L, r$lengths)); e <- sum(r$lengths[seq_len(i)]); s <- e - r$lengths[i] + 1L
cat("  longest block:", o_byte$transcript_id[s], "to", o_byte$transcript_id[e], " n =", r$lengths[i], "\n")
cat("  last present before:", o_byte$transcript_id[s-1], "  first present after:", o_byte$transcript_id[e+1], "\n")
cat("  position of block in ordered list:", s, "to", e, "of", nrow(o_byte), "\n")
blk <- o_byte$transcript_id[s:e]
fwrite(data.table(transcript_id = blk), file.path(out, "block_ids.txt"), col.names = FALSE)
o_num <- ann[order(acc_num, ver)]
runs_of(o_num, "numeric accession order")
cat("  numeric order identical to byte order:", identical(o_num$transcript_id, o_byte$transcript_id), "\n")
# all annotated transcripts with accession number in [first, last] of the block
lo <- ann[transcript_id == blk[1], acc_num]; hi <- ann[transcript_id == blk[length(blk)], acc_num]
cat("  annotated mRNAs with accession number in block range:", ann[acc_num >= lo & acc_num <= hi, .N],
    "  of them in universe:", ann[acc_num >= lo & acc_num <= hi, sum(in_uni)], "\n")
# the same block relative to ALL names of the full map with GO terms (incl. XR_/NM_ written as XM_)
allnames <- data.table(id = with_go); allnames[, acc := as.numeric(sub("^XM_([0-9]+)\\..*", "\\1", id))]
cat("  full-map names (any RNA class) with GO terms in block range:", allnames[acc >= lo & acc <= hi, .N],
    "  of them in universe:", allnames[acc >= lo & acc <= hi, sum(id %in% uni$id)], "\n")
cat("  universe ids below block:", sum(ann$in_uni & ann$acc_num < lo), "  above block:", sum(ann$in_uni & ann$acc_num > hi), "\n")
cat("  annotated mRNAs below block:", sum(ann$acc_num < lo), "  above block:", sum(ann$acc_num > hi), "\n")
cat("  missing below block:", sum(!ann$in_uni & ann$acc_num < lo), "  missing above block:", sum(!ann$in_uni & ann$acc_num > hi), "\n")

# other "natural" orderings
chr_levels <- unique(gff$chr)                                     # GFF order of sequences
o_gen <- ann[order(match(chr, chr_levels), start, end)]
runs_of(o_gen, "genomic order (GFF sequence, start)")
o_gff <- ann                                                       # the file order of the GFF itself
runs_of(o_gff, "order of the mRNA GFF file")
o_gene <- ann[base::order(ann$gene_id, ann$transcript_id, method = "radix")]
runs_of(o_gene, "gene ID (string), then ID")
o_genen <- ann[order(as.numeric(gene_id), transcript_id)]
runs_of(o_genen, "gene ID (numeric), then ID")
o_ver <- ann[order(ver, acc_num)]
runs_of(o_ver, "version, then accession")
o_file <- ann[order(match(transcript_id, unique(full$id)))]
runs_of(o_file, "order of first appearance in full map")

# series composition of the block
ann[, series := fifelse(acc_num < 12000000, "XM_011 (2014-15)", fifelse(acc_num < 30000000, "XM_0199-XM_020 (2016-17)", "XM_034 (2020)"))]
ann[, in_block := transcript_id %in% blk]
cat("\nseries x status (annotated mRNAs):\n")
print(ann[, .(n = .N, in_universe = sum(in_uni), in_block = sum(in_block), missing_outside_block = sum(!in_uni & !in_block),
              min_acc = min(acc_num), max_acc = max(acc_num)), by = series][order(min_acc)])

# ---- per-sequence shares ------------------------------------------------------
lg <- sprintf("NC_0475%02d.1", 59:68)
ann[, sequence := fifelse(chr %in% lg, chr, "unplaced scaffolds / other")]
by_chr <- ann[, .(annotated = .N, in_universe = sum(in_uni)), by = sequence][order(sequence)]
by_chr <- rbind(by_chr, data.table(sequence = "all", annotated = nrow(ann), in_universe = sum(ann$in_uni)))
by_chr[, pct := round(100 * in_universe / annotated, 2)]
print(by_chr)
cat("other sequence names:\n"); print(ann[!chr %in% lg, .N, by = .(prefix = substr(chr, 1, 3))])
ref <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables/legacy_universe_by_chromosome.csv"))
cmp <- merge(by_chr[, .(sequence = sub(" / other", "", sequence), annotated, in_universe)], ref[, .(sequence, mRNAs_with_GO_terms, in_legacy_universe)], by = "sequence")
cat("per-sequence counts equal to the pipeline table:", all(cmp$annotated == cmp$mRNAs_with_GO_terms & cmp$in_universe == cmp$in_legacy_universe), "\n")

# the 49 others
others <- ann[in_uni == FALSE & in_block == FALSE]
others[, hyphen_gene := grepl("-", parent, fixed = TRUE)]
cat("\nmissing outside the block:", nrow(others), "; of them in a gene whose name holds a hyphen:", sum(others$hyphen_gene), "\n")
print(others[, .(transcript_id, gene_id, parent, chr, start, end)])
fwrite(ann[, .(transcript_id, gene_id, parent, chr, start, end, strand, acc_num, ver, in_uni, in_block, series)], file.path(out, "annotated_mrnas.csv"))
fwrite(gff[, .(transcript_id, gene_id, parent, chr, start, end, strand, in_full, in_full_any, in_uni)], file.path(out, "all_mrnas.csv"))

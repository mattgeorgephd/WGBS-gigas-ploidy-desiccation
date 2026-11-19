suppressPackageStartupMessages({library(data.table); library(stringi)})
out <- commandArgs(trailingOnly = TRUE)[1]; dl <- commandArgs(trailingOnly = TRUE)[2]
ann <- fread(file.path(out, "annotated_mrnas.csv"), colClasses = c(gene_id = "character"))
runs <- function(ids, label) {
  x <- ann$in_uni[match(ids, ann$transcript_id)]; r <- rle(x); miss <- r$lengths[!r$values]
  i <- which.max(ifelse(r$values, 0L, r$lengths)); e <- sum(r$lengths[seq_len(i)]); s <- e - r$lengths[i] + 1L
  cat(sprintf("%-46s runs of missing %5d; longest %5d (%s to %s); missing outside it %5d\n", label, length(miss), max(miss), ids[s], ids[e], sum(miss) - max(miss)))
}
ids <- ann$transcript_id
runs(sort(ids, method = "radix"), "C / POSIX byte order")
runs(stri_sort(ids, locale = "en_US"), "ICU en_US collation")
runs(stri_sort(ids, opts_collator = stri_opts_collator(locale = "en_US", alternate_shifted = TRUE)), "ICU en_US, punctuation ignored")
runs(stri_sort(ids, opts_collator = stri_opts_collator(locale = "en_US", numeric = TRUE)), "ICU en_US, numeric (natural) collation")
runs(ids[order(gsub("[^A-Za-z0-9]", "", ids), method = "radix")], "punctuation stripped, then bytes (old glibc)")
runs(ids[order(as.numeric(sub("^XM_", "", ids)))], "numeric value of the accession.version")
runs(ids[order(stri_reverse(ids), method = "radix")], "reversed string (control: unrelated order)")
set.seed(1); runs(sample(ids), "random order (control)")
# whole lines as the legacy sort saw them: ID <tab> gene id ... ; key ties broken by later fields, so ID order rules
# the legacy unique-ID list (order of first appearance in the big table)
u <- readLines(file.path(dl, "data_gene_annotations_BOTH_union_1x-BOTH.transcriptIDs.unique"))
cat("unique-ID list from gannet: n =", length(u), "; duplicates:", anyDuplicated(u) > 0, "; all 14 characters XM_:", all(grepl("^XM_[0-9]{9}\\.[0-9]$", u)),
    "; strictly increasing in byte order:", !is.unsorted(u, strictly = TRUE), "; same set as universe:", setequal(u, ann$transcript_id[ann$in_uni]), "\n")
i <- which(u == "XM_011455140.3"); cat("neighbours in the list around the jump:", u[(i - 1):(i + 2)], "\n")
# accession number versus genome position: is accession order collinear with position inside the XM_034 series?
a <- ann[acc_num > 3e7]
lg <- sprintf("NC_0475%02d.1", 59:68)
cat("\nXM_034 series (2020 annotation): accession range per sequence\n")
print(a[, .(n = .N, min_acc = min(acc_num), median_acc = median(acc_num), max_acc = max(acc_num), in_block = sum(in_block), in_universe = sum(in_uni)),
        by = .(sequence = fifelse(chr %in% lg, chr, "unplaced"))][order(median_acc)])
b <- ann[acc_num < 3e7]
cat("\nolder series (XM_011, XM_0199, XM_020): by sequence\n")
print(b[, .(n = .N, in_block = sum(in_block), in_universe = sum(in_uni), pct_in_universe = round(100 * mean(in_uni), 1)),
        by = .(sequence = fifelse(chr %in% lg, chr, "unplaced"))][order(sequence)])

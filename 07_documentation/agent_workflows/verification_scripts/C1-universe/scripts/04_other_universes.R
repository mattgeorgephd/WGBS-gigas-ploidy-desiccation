# The legacy notebooks built the same kind of universe for the heat-only (RONIT) and the
# control-only (YAAMINI) data. The 22-library union contains the CpGs of both subsets, so under
# the intended filter the BOTH universe must contain (nearly) every transcript of either subset.
suppressPackageStartupMessages(library(data.table))
out <- commandArgs(trailingOnly = TRUE)[1]; dl <- commandArgs(trailingOnly = TRUE)[2]
rd <- function(f) fread(file.path(dl, f), header = FALSE, sep = "\t", quote = "", col.names = c("id", "go"))
both <- rd("data_gene_annotations_BOTH_geneid2go-union1x-BOTH.tab")
ron  <- rd("data_gene_annotations_RONIT_geneid2go-union1x-RONIT.tab")
yaa  <- rd("data_gene_annotations_YAAMINI_geneid2go-union1x-YAAMINI.tab")
ann <- fread(file.path(out, "annotated_mrnas_go.csv"), colClasses = c(gene_id = "character"))
ids <- list(BOTH = unique(both$id), RONIT = unique(ron$id), YAAMINI = unique(yaa$id))
cat("unique transcript names: ", paste(names(ids), lengths(ids), collapse = "; "), "\n")
for (n in names(ids)) cat(sprintf("%-8s names that are annotated mRNAs: %d; names not among the annotated mRNAs: %d; annotated mRNAs absent: %d\n", n,
    sum(ids[[n]] %in% ann$transcript_id), sum(!ids[[n]] %in% ann$transcript_id), sum(!ann$transcript_id %in% ids[[n]])))
cat("BOTH names absent from RONIT:", sum(!ids$BOTH %in% ids$RONIT), " absent from YAAMINI:", sum(!ids$BOTH %in% ids$YAAMINI), "\n")
cat("RONIT names absent from BOTH:", sum(!ids$RONIT %in% ids$BOTH), "; YAAMINI names absent from BOTH:", sum(!ids$YAAMINI %in% ids$BOTH),
    "; in RONIT or YAAMINI but not BOTH:", sum(!union(ids$RONIT, ids$YAAMINI) %in% ids$BOTH), "\n")
ann[, `:=`(in_ron = transcript_id %in% ids$RONIT, in_yaa = transcript_id %in% ids$YAAMINI)]
ann[, status := fifelse(in_uni, "present in BOTH", fifelse(in_block, "missing: block", "missing: 49 others"))]
print(ann[, .(n = .N, in_RONIT = sum(in_ron), in_YAAMINI = sum(in_yaa), in_either = sum(in_ron | in_yaa), in_neither = sum(!in_ron & !in_yaa)), by = status])
cat("\nBlock transcripts by presence in the subset universes:\n"); print(ann[in_block == TRUE, .N, by = .(in_ron, in_yaa)])
cat("\nThe 49 others:\n"); print(ann[status == "missing: 49 others", .N, by = .(hyphen, in_ron, in_yaa)])
# run structure of the YAAMINI and RONIT universes in accession order
setorder(ann, acc_num)
for (n in c("in_ron", "in_yaa")) {
  r <- rle(ann[[n]]); miss <- r$lengths[!r$values]
  cat(sprintf("%s: missing annotated mRNAs %d in %d runs; longest run %d\n", n, sum(miss), length(miss), max(miss)))
  i <- which.max(ifelse(r$values, 0L, r$lengths)); e <- sum(r$lengths[seq_len(i)]); s <- e - r$lengths[i] + 1L
  cat("   longest run:", ann$transcript_id[s], "to", ann$transcript_id[e], "\n")
}
# names in the subset universes that are not annotated mRNAs: what are they?
for (n in c("RONIT", "YAAMINI")) { x <- setdiff(ids[[n]], ann$transcript_id); cat(n, "names not annotated mRNAs:", length(x), head(x, 8), "\n") }
# are the lines of the subset universe maps the same lines as in BOTH (same GO strings)?
m <- merge(unique(both), unique(ron), by = "id", suffixes = c("_both", "_ron"))
cat("ids in BOTH and RONIT maps with identical GO string:", sum(m$go_both == m$go_ron), "of", nrow(m), "\n")
fwrite(ann[, .(transcript_id, gene_id, chr, in_uni, in_block, in_ron, in_yaa)], file.path(out, "universe_membership_three_sets.csv"))

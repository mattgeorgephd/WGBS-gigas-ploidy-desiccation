# Logical test that needs no coverage data: the legacy route assigned CpGs to GENES and then
# expanded each gene to all of its mRNAs. Two mRNAs of one gene with the same GO list are
# therefore both in the intended universe or both out of it.
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; out <- commandArgs(trailingOnly = TRUE)[1]
ann <- fread(file.path(out, "annotated_mrnas.csv"), colClasses = c(gene_id = "character"))
full <- unique(fread(cmd = sprintf("zcat %s/03_input_files/annotation/geneid2go.tab.gz", repo), header = FALSE, sep = "\t", quote = "",
              col.names = c("id", "go")))
cat("names with more than one distinct GO string in full map:", full[, .N, by = id][N > 1, .N], "\n")
gos <- full[, .(go = paste(sort(unique(trimws(unlist(strsplit(go, ","))))), collapse = ",")), by = id]
ann <- merge(ann, gos, by.x = "transcript_id", by.y = "id", all.x = TRUE)
stopifnot(!anyNA(ann$go))
ann[, hyphen := grepl("-", parent, fixed = TRUE)]
# sibling present in universe, same gene
ann[, gene_has_present := any(in_uni), by = gene_id]
pres <- unique(ann[in_uni == TRUE, .(gene_id, go)])[, key := paste(gene_id, go, sep = "|")]
ann[, same_go_sibling_present := paste(gene_id, go, sep = "|") %in% pres$key]
# any shared GO id with a present sibling
pres_long <- unique(ann[in_uni == TRUE, .(go_id = unlist(strsplit(go, ","))), by = gene_id])[, k := paste(gene_id, go_id)]
miss_long <- ann[in_uni == FALSE, .(go_id = unlist(strsplit(go, ","))), by = .(transcript_id, gene_id)]
miss_long[, shared := paste(gene_id, go_id) %in% pres_long$k]
sh <- miss_long[, .(shares_go_with_present_sibling = any(shared)), by = transcript_id]
ann <- merge(ann, sh, by = "transcript_id", all.x = TRUE)
ann[, status := fifelse(in_uni, "present", fifelse(in_block, "block", "others49"))]
print(ann[in_uni == FALSE, .(transcripts = .N, genes = uniqueN(gene_id),
            in_gene_with_present_sibling = sum(gene_has_present),
            identical_go_list_as_present_sibling = sum(same_go_sibling_present),
            shares_a_go_id_with_present_sibling = sum(shares_go_with_present_sibling)), by = status])
# the 49 others in detail
legacy_terms <- fread(cmd = sprintf("zcat %s/03_input_files/annotation/legacy_go_goslim_terms.tsv.gz", repo), header = TRUE, sep = "\t", quote = "")
cat("recovered legacy GOslim table: GO ids:", uniqueN(legacy_terms$GOID), "\n")
o <- ann[status == "others49"]
o[, n_go := lengths(strsplit(go, ","))]
o[, n_go_in_recovered_goslim := vapply(strsplit(go, ","), function(v) sum(v %in% legacy_terms$GOID), integer(1))]
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))
o <- merge(o, gt[, .(gene_id, n_cpg_tested)], by = "gene_id")
setorder(o, acc_num)
print(o[, .(transcript_id, gene_id, parent, hyphen, gene_has_present, n_go, n_go_in_recovered_goslim, n_cpg_tested, go = substr(go, 1, 40))])
cat("\n49 others: hyphen gene:", sum(o$hyphen), "; not hyphen:", sum(!o$hyphen),
    "; not hyphen and no GO id in recovered GOslim table:", sum(!o$hyphen & o$n_go_in_recovered_goslim == 0),
    "; not hyphen, with a GO id known to be in the GOslim table:", sum(!o$hyphen & o$n_go_in_recovered_goslim > 0), "\n")
# how often do annotated transcripts in general have no GO id in the recovered GOslim table?
ann[, n_go_rec := vapply(strsplit(go, ","), function(v) sum(v %in% legacy_terms$GOID), integer(1))]
print(ann[, .(n = .N, no_go_in_recovered_table = sum(n_go_rec == 0)), by = status])
fwrite(o, file.path(out, "others49.csv"))
fwrite(ann, file.path(out, "annotated_mrnas_go.csv"))

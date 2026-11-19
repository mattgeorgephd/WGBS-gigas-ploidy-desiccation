suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
ms <- readLines("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/manuscript/manuscript_gdoc.txt")
rows <- ms[grepl("^\\| (ploidy|heat|interaction)-DML \\|", ms)]
p <- lapply(strsplit(rows, "|", fixed = TRUE), function(v) trimws(gsub("\\\\", "", v[-1])))
m <- as.data.table(do.call(rbind, lapply(p, function(v) v[1:8])))
setnames(m, c("dml_set", "gene_id", "gene_name", "chromosome", "start", "end", "ms_diff", "category"))
m[, `:=`(start = as.integer(start), end = as.integer(end), ms_diff = as.numeric(ms_diff), row = .I)]
t1 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables/table_1_legacy.csv"), colClasses = c(gene_id = "character"))
t1[, row9 := .I]
cat("manuscript rows:", nrow(m), "; step 09 rows:", nrow(t1), "\n")
x <- merge(m, t1, by = c("dml_set", "gene_id", "chromosome", "start"), all = TRUE, suffixes = c("_ms", "_s9"))
setorder(x, row)
x[, `:=`(name_same = gene_name_ms == gene_name_s9, end_same = end_ms == end_s9, rounds_to_ms = round(methylation_difference) == ms_diff, same_order = row == row9)]
options(width = 230)
print(x[, .(row, dml_set, gene_id, start, ms_diff, regenerated = methylation_difference, rounds_to_ms, name_same, end_same, same_order, category)])
cat("unmatched rows:", sum(is.na(x$row) | is.na(x$row9)), "; names identical:", sum(x$name_same), "; ends identical:", sum(x$end_same),
    "; rounded value equals manuscript:", sum(x$rounds_to_ms), "; same row order:", sum(x$same_order), "\n")
cat("rows where |regenerated| < 20.5 and manuscript prints 20 or 21:\n"); print(x[abs(methylation_difference) < 21, .(dml_set, gene_name_ms, start, ms_diff, methylation_difference)])

cat("\n==== step 09 Fisher tables against the archive, every row\n")
arch_dir <- c(ploidy = "BOTH-PLOIDY", heat = "BOTH-HEAT", shared = "BOTH-COMBO"); lab <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")
fmt <- function(p) ifelse(p < 1e-30, "< 1e-30", ifelse(p < 1e-4, format(p, digits = 2, scientific = TRUE), NA))
for (s in names(arch_dir)) for (o in c("BP", "CC", "MF")) {
  a <- readLines(file.path(repo, "08_archive/GOterms", arch_dir[[s]], sprintf("%s-%s-FisherTestResults.csv", lab[[s]], o)))[-1]
  parts <- strsplit(a, ",", fixed = TRUE); a_go <- vapply(parts, `[`, "", 1); tl <- t(vapply(parts, function(v) tail(v, 4), character(4)))
  f9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_fisher_%s_%s.csv", s, o)))
  ap <- as.numeric(sub("< ", "", tl[, 4]))
  cat(sprintf("%-6s %s: rows %d/%d; GO IDs in the same order %s; Annotated equal %s; Significant equal %s; Expected equal %s; max relative difference of p against the printed archive value (p >= 1e-30): %.3g\n",
              s, o, length(a_go), nrow(f9), identical(a_go, f9$go_id), all(as.integer(tl[, 1]) == f9$annotated), all(as.integer(tl[, 2]) == f9$significant),
              all(as.numeric(tl[, 3]) == f9$expected), max(abs(f9$p_fisher - ap)[tl[, 4] != "< 1e-30"] / ap[tl[, 4] != "< 1e-30"])))
}
cat("\n==== regenerated annotation tables (step 08) against the gannet originals\n")
cols <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
for (s in names(arch_dir)) {
  orig <- readLines(file.path("/home/user/gannet_cache/legacy_goannot", sprintf("DML-Cov10-20-BOTH-%s.GOAnnot", lab[[s]])))
  regen <- readLines(gzfile(file.path(repo, "04_output/08_gene_annotation/tables", sprintf("dml_go_annotation_%s_legacy.tsv.gz", s))))
  o2 <- fread(text = orig, header = FALSE, sep = "\t", quote = "", col.names = cols, na.strings = NULL)
  r2 <- fread(text = regen, header = FALSE, sep = "\t", quote = "", col.names = cols, na.strings = NULL)
  same_struct <- nrow(o2) == nrow(r2) && all(o2$transcript == r2$transcript & o2$geneID == r2$geneID & o2$chr == r2$chr & o2$start == r2$start & o2$GOID == r2$GOID &
                                               o2$GOterm == r2$GOterm & o2$GOSlim == r2$GOSlim & o2$GOcat == r2$GOcat)
  cat(sprintf("%-6s lines original %d regenerated %d; byte-identical lines %d; same rows in the same order apart from meth.diff digits: %s; max |meth.diff difference| %.3g\n",
              s, length(orig), length(regen), sum(orig == regen[seq_along(orig)]), same_struct, if (same_struct) max(abs(o2$meth.diff - r2$meth.diff)) else NA))
}

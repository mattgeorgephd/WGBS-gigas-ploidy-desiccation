suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
for (o in c("CC", "MF")) {
  g <- fread(file.path(repo, "08_archive/GOterms/BOTH-PLOIDY", sprintf("GENELIST-%s-EnrichedGO-DMLOnly.csv", o)))
  s9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_dml_in_top10_genes_ploidy_%s.csv", o)))
  cat("\n", o, ": archive cols:", names(g), "| step09 cols:", names(s9), "\n")
  m <- merge(g[, .(geneID, chr, start, end, md_arch = meth.diff)], s9[, .(geneID, chr, start, end, md_s9 = meth.diff, gene, product)],
             by = c("geneID", "chr", "start", "end"), all = TRUE)
  cat("rows archive", nrow(g), "step09", nrow(s9), "merged", nrow(m), "; unmatched:", sum(is.na(m$md_arch) | is.na(m$md_s9)),
      "; max abs diff", max(abs(m$md_arch - m$md_s9), na.rm = TRUE), "\n")
  print(m[is.na(md_arch) | is.na(md_s9) | abs(md_arch - md_s9) > 1e-7])
  print(sapply(g, class)); print(sapply(s9, class))
  print(head(g[, sprintf("%.10f", meth.diff)])); print(head(s9[, sprintf("%.10f", meth.diff)]))
}

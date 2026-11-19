S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
x <- readRDS(file.path(S, "tables.rds")); universe <- x$universe
for (s in c("ploidy", "heat")) {
  p <- readRDS(file.path(S, sprintf("pattern_%s.rds", s))); g <- p$g; tx <- p$tx
  g$cls <- cut(g$n_tx, c(0, 1, 2, 5, Inf), labels = c("1", "2", "3 to 5", "more than 5"))
  a <- aggregate(cbind(genes = 1, kept = as.integer(tx_kept > 0), tested = as.integer(tx_tested > 0), in_univ = as.integer(tx_univ > 0)) ~ cls, data = g, FUN = sum)
  a$pct_kept <- round(100 * a$kept / a$genes, 1); a$pct_tested <- round(100 * a$tested / a$genes, 1); a$pct_in_universe_exact <- round(100 * a$in_univ / a$genes, 1)
  cat("\n", s, ": genes by number of transcripts in the table (isoforms holding the DML)\n", sep = ""); print(a, row.names = FALSE)
  # transcripts per tested gene: legacy vs complete read
  tested_tx <- tx[tx$kept > 0 & tx$in_universe, ]; all_tx <- tx[tx$in_universe, ]
  cat(sprintf("transcripts of interest per gene in the test: legacy reader %.2f (%d transcripts, %d genes); complete read %.2f (%d transcripts, %d genes)\n",
              nrow(tested_tx) / length(unique(tested_tx$geneID)), nrow(tested_tx), length(unique(tested_tx$geneID)),
              nrow(all_tx) / length(unique(all_tx$geneID)), nrow(all_tx), length(unique(all_tx$geneID))))
  # single-row accounting
  cat(sprintf("lines of transcripts that are lost entirely: %d; lines lost from transcripts that are kept: %d\n",
              sum(tx$lines[tx$kept == 0]), sum(tx$lines[tx$kept > 0] - tx$kept[tx$kept > 0])))
}

# Would the comments "#662 / #178 / #214" (written BEFORE the script's Freq != 0 filter) arise with factor columns (R < 4.0 default)?
f <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader/gannet_dl/data_BOTH_PLOIDY.GOAnnot"
cn <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
have_dplyr <- requireNamespace("dplyr", quietly = TRUE)
for (saf in c(FALSE, TRUE)) {
  d <- suppressWarnings(read.table(f, header = FALSE, sep = "\t", col.names = cn, stringsAsFactors = saf))
  n <- sapply(c("P", "C", "F"), function(k) {
    sub <- if (have_dplyr) dplyr::filter(d, d$GOcat == k) else d[which(d$GOcat == k), ]
    nrow(as.data.frame(table(sub$GOID)))
  })
  cat(sprintf("stringsAsFactors = %-5s rows %d; nrow(as.data.frame(table(GOID))) for P/C/F before any Freq filter: %s (dplyr used: %s)\n", saf, nrow(d), paste(n, collapse = "/"), have_dplyr))
}

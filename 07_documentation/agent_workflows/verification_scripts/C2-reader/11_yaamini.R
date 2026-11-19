f <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader/gannet_dl/data_YAAMINI_PLOIDY.GOAnnot"
cn <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
w <- character(0)
a <- withCallingHandlers(read.table(f, header = FALSE, sep = "\t", col.names = cn), warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") })
b <- read.table(f, header = FALSE, sep = "\t", col.names = cn, quote = "", comment.char = "")
cat(sprintf("control-only (YAAMINI) ploidy table: lines %d; read.table defaults %d rows (%.1f%%), %d of %d transcripts, %d of %d genes; warnings: %s\n",
            nrow(b), nrow(a), 100 * nrow(a) / nrow(b), length(unique(a$transcript)), length(unique(b$transcript)),
            length(unique(a$geneID)), length(unique(b$geneID)), if (length(w)) paste(unique(w), collapse = " | ") else "none"))

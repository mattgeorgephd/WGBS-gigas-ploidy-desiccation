# 04_pattern.R (scratch) Which rows, transcripts and genes survive the legacy reader?
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
x <- readRDS(file.path(S, "tables.rds")); E <- x$E; L <- x$L; universe <- x$universe
options(width = 230)
lgs <- sprintf("NC_0475%02d.1", 59:68)

for (s in names(E)) {
  e <- E[[s]]; l <- L[[s]]
  cat("\n==================== ", s, " ====================\n", sep = "")
  # records of the legacy reader in terms of file lines (from the independent state machine)
  rec <- read.delim(file.path(S, sprintf("records_%s.tsv", s)))
  stopifnot(nrow(rec) == nrow(l))
  # the transcript column of legacy row i must be the transcript on file line rec$first_line[i]
  stopifnot(identical(e$transcript[rec$first_line], l$transcript), identical(e$GOID[rec$first_line], l$GOID))
  # ... and its GOcat is that of the LAST line of the record (chimeric when the record spans lines)
  multi <- rec$last_line > rec$first_line
  last_cat <- e$GOcat[rec$last_line]
  eof_open <- s == "ploidy"
  if (eof_open) { ok <- identical(last_cat[-nrow(rec)], l$GOcat[-nrow(rec)]) } else ok <- identical(last_cat, l$GOcat)
  cat(sprintf("legacy rows = state-machine records: %d; GOcat of each legacy row equals GOcat of the record's last file line: %s\n", nrow(rec), ok))
  chim <- sum(multi & e$GOID[rec$first_line] != e$GOID[rec$last_line])
  cat(sprintf("multi-line (merged) rows: %d; of these, GOID (first line) and GOSlim/GOcat (last line) come from different GO terms: %d; transcript differs between first and last line: %d; gene differs: %d\n",
              sum(multi), chim, sum(e$transcript[rec$first_line] != e$transcript[rec$last_line]), sum(e$geneID[rec$first_line] != e$geneID[rec$last_line])))
  cat(sprintf("merged rows whose GOcat differs from the true GOcat of their GOID line: %d\n", sum(e$GOcat[rec$first_line] != l$GOcat)))

  # line status: start of a record (kept as a row) or absorbed
  e$kept <- FALSE; e$kept[rec$first_line] <- TRUE
  cat(sprintf("lines kept as row starts: %d of %d; absorbed: %d\n", sum(e$kept), nrow(e), sum(!e$kept)))
  blk <- rec[multi, ]; blk$absorbed <- blk$last_line - blk$first_line
  blk <- blk[order(-blk$absorbed), ]
  cat(sprintf("absorbed lines by merged record, largest first: %s\n", paste(blk$absorbed, collapse = " ")))
  cat(sprintf("share of all absorbed lines in the 3 largest merged records: %.1f%%; in the 7 largest: %.1f%%\n",
              100 * sum(head(blk$absorbed, 3)) / sum(blk$absorbed), 100 * sum(head(blk$absorbed, 7)) / sum(blk$absorbed)))

  # is the file sorted by transcript?  (C-locale order)
  cat(sprintf("file sorted by transcript accession (C order): %s; transcripts contiguous: %s\n",
              !is.unsorted(e$transcript, strictly = FALSE) || identical(e$transcript, sort(e$transcript, method = "radix")),
              length(rle(e$transcript)$values) == length(unique(e$transcript))))

  # ---- per transcript ----
  tx <- aggregate(cbind(lines = 1, kept = kept) ~ transcript + geneID + chr, data = transform(e, kept = as.integer(kept)), FUN = sum)
  tx$status <- ifelse(tx$kept == 0, "lost", ifelse(tx$kept == tx$lines, "all rows kept", "some rows kept"))
  tx$in_universe <- tx$transcript %in% universe
  cat("transcripts by status:\n"); print(table(tx$status))
  stopifnot(sum(tx$kept > 0) == length(unique(l$transcript)))
  # ---- per gene ----
  g <- aggregate(cbind(n_tx = 1, tx_kept = as.integer(kept > 0), tx_univ = as.integer(in_universe), tx_tested = as.integer(kept > 0 & in_universe)) ~ geneID, data = tx, FUN = sum)
  g$status <- ifelse(g$tx_kept == 0, "no transcript kept (gene lost)", ifelse(g$tx_kept == g$n_tx, "every transcript kept", "some transcripts kept"))
  cat("genes by status:\n"); print(table(g$status))
  cat(sprintf("genes: %d; with >=1 transcript kept: %d; with >=1 kept transcript in the universe (tested): %d; with >=1 transcript in universe (exact reader): %d\n",
              nrow(g), sum(g$tx_kept > 0), sum(g$tx_tested > 0), sum(g$tx_univ > 0)))
  cat(sprintf("mean transcripts per gene (in this table): genes kept %.2f, genes lost %.2f\n", mean(g$n_tx[g$tx_kept > 0]), mean(g$n_tx[g$tx_kept == 0])))

  # ---- reader x universe, transcripts ----
  tb <- table(kept_by_reader = tx$kept > 0, in_universe = tx$in_universe)
  cat("transcripts: kept by reader x in legacy universe\n"); print(addmargins(tb))
  if (all(dim(tb) == 2)) { ft <- fisher.test(tb); cat(sprintf("odds ratio %.2f (95%% CI %.2f to %.2f), p = %.3g\n", ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)) }
  cat(sprintf("share in universe: kept transcripts %.1f%%, lost transcripts %.1f%%\n", 100 * mean(tx$in_universe[tx$kept > 0]),
              if (any(tx$kept == 0)) 100 * mean(tx$in_universe[tx$kept == 0]) else NA))

  # ---- accession structure ----
  tx <- tx[order(tx$transcript), ]
  tx$series <- sub("^(XM_[0-9]{3}).*$", "\\1", tx$transcript)
  cat("transcripts kept by accession series (first three digits):\n")
  print(cbind(all = table(tx$series), kept = table(factor(tx$series[tx$kept > 0], levels = sort(unique(tx$series)))),
              in_universe = table(factor(tx$series[tx$in_universe], levels = sort(unique(tx$series))))))
  r <- rle(tx$kept > 0); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  runs <- data.frame(kept = r$values, n_tx = r$lengths, first = tx$transcript[starts], last = tx$transcript[ends])
  cat(sprintf("runs of consecutive transcripts (accession order): %d runs; kept runs %d, lost runs %d\n", nrow(runs), sum(runs$kept), sum(!runs$kept)))
  cat("largest lost runs:\n"); print(head(runs[!runs$kept, ][order(-runs$n_tx[!runs$kept]), ], 8), row.names = FALSE)
  cat("largest kept runs:\n"); print(head(runs[runs$kept, ][order(-runs$n_tx[runs$kept]), ], 8), row.names = FALSE)

  # ---- by chromosome: genes ----
  gchr <- unique(tx[, c("geneID", "chr")]); gchr$seq <- ifelse(gchr$chr %in% lgs, gchr$chr, "scaffolds")
  gchr <- merge(gchr, g[, c("geneID", "tx_kept", "tx_tested", "tx_univ")], by = "geneID")
  bychr <- aggregate(cbind(genes = 1, kept = as.integer(tx_kept > 0), tested = as.integer(tx_tested > 0), in_univ_exact = as.integer(tx_univ > 0)) ~ seq, data = gchr, FUN = sum)
  bychr$pct_kept <- round(100 * bychr$kept / bychr$genes, 1); bychr$pct_tested <- round(100 * bychr$tested / bychr$genes, 1)
  bychr$pct_univ_exact <- round(100 * bychr$in_univ_exact / bychr$genes, 1)
  cat("genes with GO terms holding a DML, by sequence: kept by reader, tested (kept and in universe), in universe had the reader been exact\n"); print(bychr, row.names = FALSE)
  if (any(bychr$kept < bychr$genes)) {
    ct <- suppressWarnings(chisq.test(rbind(bychr$kept, bychr$genes - bychr$kept)))
    cat(sprintf("heterogeneity of reader retention across sequences: chi-squared = %.1f, df = %d, p = %.3g\n", ct$statistic, ct$parameter, ct$p.value))
  }
  ct2 <- suppressWarnings(chisq.test(rbind(bychr$tested, bychr$genes - bychr$tested)))
  cat(sprintf("heterogeneity of 'tested' across sequences: chi-squared = %.1f, df = %d, p = %.3g\n", ct2$statistic, ct2$parameter, ct2$p.value))
  write.csv(bychr, file.path(S, sprintf("retention_by_chromosome_%s.csv", s)), row.names = FALSE)

  # ---- DML and direction ----
  e$dml <- paste(e$chr, e$start); l$dml <- paste(l$chr, l$start)
  d_all <- unique(e[, c("dml", "meth.diff")]); d_all$meth.diff <- as.numeric(d_all$meth.diff)
  d_all$kept <- d_all$dml %in% l$dml
  cat(sprintf("DML in the table: %d; in rows kept: %d (%.1f%%). Hypermethylated share: all %.1f%%, kept %.1f%%, lost %.1f%%\n", nrow(d_all), sum(d_all$kept), 100 * mean(d_all$kept),
              100 * mean(d_all$meth.diff > 0), 100 * mean(d_all$meth.diff[d_all$kept] > 0), if (any(!d_all$kept)) 100 * mean(d_all$meth.diff[!d_all$kept] > 0) else NA))

  # ---- transcripts whose rows hold a quote-toggling line ----
  napos <- lengths(regmatches(e$GOterm, gregexpr("'", e$GOterm, fixed = TRUE)))
  dq <- grepl('^".*"$', e$GOterm)
  e$toggle <- napos %% 2 == 1 & !dq                         # odd number of apostrophes outside a double-quoted name
  e$toggle_dq <- napos %% 2 == 1 & dq
  cat(sprintf("lines with an apostrophe: %d; with an odd number outside double quotes (togglers): %d; odd number inside a double-quoted name: %d; even number: %d\n",
              sum(napos > 0), sum(e$toggle), sum(e$toggle_dq), sum(napos > 0 & napos %% 2 == 0)))
  tt <- unique(e$transcript[e$toggle]); tg <- unique(e$geneID[e$toggle])
  cat(sprintf("transcripts carrying a toggling term: %d (genes %d); kept by the reader: %d (%.0f%%) against %.0f%% of all other transcripts\n",
              length(tt), length(tg), sum(tt %in% l$transcript), 100 * mean(tt %in% l$transcript),
              100 * mean(setdiff(unique(e$transcript), tt) %in% l$transcript)))
  cat("toggling terms (lines):\n"); print(sort(table(paste(e$GOID[e$toggle | e$toggle_dq], e$GOterm[e$toggle | e$toggle_dq])), decreasing = TRUE))
  saveRDS(list(tx = tx, g = g, e = e), file.path(S, sprintf("pattern_%s.rds", s)))
}

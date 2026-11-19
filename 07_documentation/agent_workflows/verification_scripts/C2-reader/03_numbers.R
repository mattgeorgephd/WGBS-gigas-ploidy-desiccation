# 03_numbers.R  (scratch) Recompute every number of I-40 from the legacy tables, independently
# of the pipeline code (base R only; nothing from 02_R_functions is sourced).
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
ann <- file.path(repo, "03_input_files/annotation")
files <- c(
  ploidy = file.path(ann, "legacy/DML-Cov10-20-BOTH-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz"),
  heat   = file.path(ann, "legacy/DML-Cov10-20-BOTH-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz"),
  shared = file.path(ann, "legacy/DML-Cov10-20-BOTH-COMBO.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz"))
cn <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
options(width = 220)

# ---- exact read without read.table at all: readLines + strsplit ----
exact_read <- function(f) {
  l <- readLines(f)
  p <- strsplit(paste0(l, "\t<END>"), "\t", fixed = TRUE)     # sentinel keeps trailing empty fields
  stopifnot(all(lengths(p) == 11), all(vapply(p, `[`, "", 11) == "<END>"))
  d <- as.data.frame(do.call(rbind, p)[, 1:10], stringsAsFactors = FALSE); names(d) <- cn
  d$line <- seq_along(l)
  d
}
# ---- the legacy call, verbatim arguments ----
legacy_read <- function(f) suppressWarnings(read.table(f, header = FALSE, sep = "\t", col.names = cn))

# ---- universe: first column of the map, read with readLines (not topGO::readMappings) ----
ul <- readLines(file.path(ann, "geneid2go-union1x-BOTH.tab.gz"))
u_names <- gsub(" ", "", sub("\t.*$", "", ul))          # readMappings strips blanks from names
universe <- unique(u_names)
cat(sprintf("universe map: %d lines, %d distinct transcript names\n", length(ul), length(universe)))

# ---- independent transcript -> gene map (legacy intermediate: geneID <tab> transcriptID) ----
tg <- read.table(file.path(ann, "legacy/cgigas_uk_roslin_v1_mRNA.transcriptID.geneID.gz"), sep = "\t", quote = "",
                 comment.char = "", col.names = c("gene", "transcript"), colClasses = "character")
tx2gene <- setNames(tg$gene, tg$transcript)
cat(sprintf("transcript-gene map: %d rows, %d transcripts, %d genes\n", nrow(tg), length(unique(tg$transcript)), length(unique(tg$gene))))

res <- list(); E <- list(); L <- list()
for (s in names(files)) {
  e <- exact_read(files[[s]]); l <- legacy_read(files[[s]])
  E[[s]] <- e; L[[s]] <- l
  # gene IDs in the table agree with the independent map?
  chk <- unique(e[, c("transcript", "geneID")])
  agree <- sum(tx2gene[chk$transcript] == chk$geneID, na.rm = TRUE); miss <- sum(is.na(tx2gene[chk$transcript]))
  kept_tx <- unique(l$transcript); kept_gene <- unique(as.character(l$geneID))
  tested_tx <- intersect(kept_tx, universe)
  tested_gene_tablecol <- unique(as.character(l$geneID[l$transcript %in% universe]))
  tested_gene_map <- unique(unname(tx2gene[tested_tx]))
  all_tx <- unique(e$transcript); all_gene <- unique(e$geneID)
  exact_tested_tx <- intersect(all_tx, universe)
  exact_tested_gene <- unique(e$geneID[e$transcript %in% universe])
  res[[s]] <- data.frame(
    set = s, lines = nrow(e), rows_legacy = nrow(l), pct_rows_read = round(100 * nrow(l) / nrow(e), 2),
    pct_rows_skipped = round(100 * (1 - nrow(l) / nrow(e)), 2),
    tx_all = length(all_tx), tx_kept = length(kept_tx), genes_all = length(all_gene), genes_kept = length(kept_gene),
    pct_genes_kept = round(100 * length(kept_gene) / length(all_gene), 1),
    tx_tested = length(tested_tx), genes_tested = length(tested_gene_tablecol), genes_tested_via_map = length(tested_gene_map),
    pct_genes_tested = round(100 * length(tested_gene_tablecol) / length(all_gene), 1),
    genes_in_universe_exact_reader = length(exact_tested_gene), tx_in_universe_exact_reader = length(exact_tested_tx),
    pct_genes_in_universe_exact = round(100 * length(exact_tested_gene) / length(all_gene), 1),
    tx_gene_pairs_agree_with_map = agree, tx_missing_from_map = miss, tx_gene_pairs = nrow(chk))
}
res <- do.call(rbind, res)
print(t(res))
write.csv(res, file.path(S, "numbers_recomputed.csv"), row.names = FALSE)
saveRDS(list(E = E, L = L, universe = universe, tx2gene = tx2gene), file.path(S, "tables.rds"))

# ---- per-ontology rows and distinct GO IDs, legacy reader, against the script's comments ----
cat("\nRows by GOcat and distinct GO IDs (legacy reader) / (exact reader)\n")
for (s in names(files)) {
  l <- L[[s]]; e <- E[[s]]
  cat(sprintf("%-6s legacy rows P/C/F = %d/%d/%d  other GOcat rows = %d (values: %s);  distinct GOIDs P/C/F = %d/%d/%d\n", s,
              sum(l$GOcat == "P", na.rm = TRUE), sum(l$GOcat == "C", na.rm = TRUE), sum(l$GOcat == "F", na.rm = TRUE),
              sum(!(l$GOcat %in% c("P", "C", "F"))), paste(sprintf("'%s'x%d", names(table(l$GOcat[!(l$GOcat %in% c("P", "C", "F"))], useNA = "ifany")),
                                                                    as.integer(table(l$GOcat[!(l$GOcat %in% c("P", "C", "F"))], useNA = "ifany"))), collapse = ","),
              length(unique(l$GOID[which(l$GOcat == "P")])), length(unique(l$GOID[which(l$GOcat == "C")])), length(unique(l$GOID[which(l$GOcat == "F")]))))
  cat(sprintf("%-6s exact  rows P/C/F = %d/%d/%d  other GOcat rows = %d (values: %s);  distinct GOIDs P/C/F = %d/%d/%d\n", s,
              sum(e$GOcat == "P"), sum(e$GOcat == "C"), sum(e$GOcat == "F"), sum(!(e$GOcat %in% c("P", "C", "F"))),
              paste(sprintf("'%s'x%d", names(table(e$GOcat[!(e$GOcat %in% c("P", "C", "F"))])), as.integer(table(e$GOcat[!(e$GOcat %in% c("P", "C", "F"))]))), collapse = ","),
              length(unique(e$GOID[e$GOcat == "P"])), length(unique(e$GOID[e$GOcat == "C"])), length(unique(e$GOID[e$GOcat == "F"]))))
}

# ---- last row of the ploidy table as the legacy reader returned it ----
cat("\nLast two rows of the ploidy table, legacy reader (GOterm truncated to 60 characters):\n")
lp <- L$ploidy; tl <- tail(lp, 2); tl$GOterm <- substr(gsub("\n", "\\\\n", tl$GOterm), 1, 60); print(tl)
cat("is.na / empty GOSlim, GOcat in the last row: ", is.na(tail(lp$GOSlim, 1)), is.na(tail(lp$GOcat, 1)),
    identical(tail(lp$GOSlim, 1), ""), identical(tail(lp$GOcat, 1), ""), "\n")
cat("rows with NA anywhere (legacy reader): ", sapply(L, function(d) sum(!complete.cases(d))), "\n")

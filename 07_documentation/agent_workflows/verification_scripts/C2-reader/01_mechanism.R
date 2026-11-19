# 01_mechanism.R  (scratch; independent of the pipeline helpers)
# Read the three legacy DML annotation tables with different read.table settings
# and record rows, transcripts, genes, and warnings.
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
files <- c(
  ploidy = file.path(repo, "03_input_files/annotation/legacy/DML-Cov10-20-BOTH-PLOIDY.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz"),
  heat   = file.path(repo, "03_input_files/annotation/legacy/DML-Cov10-20-BOTH-HEAT.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz"),
  shared = file.path(repo, "03_input_files/annotation/legacy/DML-Cov10-20-BOTH-COMBO.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz"))
cn <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")

try_read <- function(label, expr) {
  w <- character(0)
  val <- tryCatch(withCallingHandlers(expr, warning = function(m) { w <<- c(w, conditionMessage(m)); invokeRestart("muffleWarning") }),
                  error = function(e) structure(conditionMessage(e), class = "read_error"))
  if (inherits(val, "read_error")) {
    return(list(df = NULL, row = data.frame(reader = label, rows = NA, transcripts = NA, genes = NA, dml = NA,
                                            multi_line_rows = NA, warnings = paste("ERROR:", val))))
  }
  list(df = val, row = data.frame(
    reader = label, rows = nrow(val), transcripts = length(unique(val$transcript)), genes = length(unique(val$geneID)),
    dml = length(unique(paste(val$chr, val$start))),
    multi_line_rows = sum(grepl("\n", val$GOterm, fixed = TRUE)),
    warnings = if (length(w)) paste(unique(w), collapse = " | ") else ""))
}

out <- list(); keep <- list()
for (s in names(files)) {
  f <- files[[s]]
  n_lines <- length(readLines(f))
  readers <- list(
    "A defaults (as legacy script)"              = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn)),
    "B quote = '' (quoting off)"                 = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, quote = "")),
    "C quote = '\"' only"                         = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, quote = "\"")),
    "D quote = \"'\" only"                        = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, quote = "'")),
    "E defaults + comment.char = ''"             = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, comment.char = "")),
    "F defaults + fill = TRUE"                   = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, fill = TRUE)),
    "G defaults + fill = TRUE + comment.char=''" = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, fill = TRUE, comment.char = "")),
    "H quote = '' + comment.char = '' + fill"    = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, quote = "", comment.char = "", fill = TRUE)),
    "I defaults + stringsAsFactors = TRUE"       = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, stringsAsFactors = TRUE)),
    "J read.delim(header = FALSE) defaults"      = quote(read.delim(f, header = FALSE, col.names = cn)),
    "K defaults + allowEscapes = TRUE"           = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, allowEscapes = TRUE)),
    "L defaults + strip.white = TRUE"            = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, strip.white = TRUE)),
    "M defaults, colClasses = 'character'"       = quote(read.table(f, header = FALSE, sep = "\t", col.names = cn, colClasses = "character")))
  for (lab in names(readers)) {
    r <- try_read(lab, eval(readers[[lab]]))
    out[[length(out) + 1]] <- cbind(set = s, lines_in_file = n_lines, r$row)
    if (substr(lab, 1, 1) %in% c("A", "B", "C", "D")) keep[[paste(s, substr(lab, 1, 1))]] <- r$df
  }
}
res <- do.call(rbind, out)
res$pct_rows <- round(100 * res$rows / res$lines_in_file, 2)
options(width = 250)
print(res, row.names = FALSE)
write.csv(res, file.path(S, "mechanism_reader_comparison.csv"), row.names = FALSE)
saveRDS(keep, file.path(S, "reads_ABCD.rds"))

# Column classes after the default read (what type.convert made of them)
for (s in names(files)) {
  cat("\n", s, " default read column classes: ", paste(sapply(keep[[paste(s, "A")]], class), collapse = ", "), "\n", sep = "")
}
cat("\nR version: ", R.version.string, "\n")

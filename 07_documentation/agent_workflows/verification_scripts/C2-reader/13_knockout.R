# Knock-out experiment: which lines cause the loss?
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
cn <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
rd <- function(lines) { tf <- tempfile(tmpdir = S, fileext = ".tsv"); on.exit(unlink(tf)); writeLines(lines, tf)
  suppressWarnings(read.table(tf, header = FALSE, sep = "\t", col.names = cn)) }
for (s in c("ploidy", "heat", "shared")) {
  l <- readLines(file.path(S, sprintf("legacy_%s.tsv", s)))
  term <- vapply(strsplit(l, "\t", fixed = TRUE), `[`, "", 8)
  napos <- lengths(regmatches(term, gregexpr("'", term, fixed = TRUE)))
  dq <- grepl('^".*"$', term)
  odd_out <- napos %% 2 == 1 & !dq; odd_in <- napos %% 2 == 1 & dq; even <- napos > 0 & napos %% 2 == 0
  base <- rd(l)
  a <- l; a[even] <- gsub("'", "^", a[even], fixed = TRUE)                 # neutralise paired apostrophes only
  b <- l; b[odd_out] <- gsub("'", "^", b[odd_out], fixed = TRUE)           # neutralise unpaired apostrophes outside double quotes
  c2 <- l; c2[odd_in] <- gsub("'", "^", c2[odd_in], fixed = TRUE)          # neutralise unpaired apostrophes inside double-quoted names
  ra <- rd(a); rb <- rd(b); rc <- rd(c2)
  cat(sprintf("%-6s lines %d | apostrophe lines: paired %d, unpaired outside double quotes %d, unpaired inside a double-quoted name %d\n", s, length(l), sum(even), sum(odd_out), sum(odd_in)))
  cat(sprintf("       rows: as is %d | paired apostrophes removed %d (same rows start on same transcript+GOID: %s) | unpaired-outside removed %d | unpaired-inside removed %d\n",
              nrow(base), nrow(ra), identical(paste(base$transcript, base$GOID), paste(ra$transcript, ra$GOID)), nrow(rb), nrow(rc)))
}

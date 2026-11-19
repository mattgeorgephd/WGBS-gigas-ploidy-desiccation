# Shared helper: parse an archived, unquoted topGO GenTable CSV from both ends.
# Line layout: GO.ID,Term,Annotated,Significant,Expected,classic
# Term may contain commas; "classic" may be "< 1e-30" or "<1e-30".
parse_archive <- function(f) {
  l <- readLines(f)
  hdr <- l[1]; l <- l[-1]
  parts <- strsplit(l, ",", fixed = TRUE)
  n <- lengths(parts)
  stopifnot(all(n >= 6))
  d <- data.frame(
    go_id       = vapply(parts, function(p) p[1], ""),
    term        = vapply(parts, function(p) paste(p[2:(length(p) - 4)], collapse = ","), ""),
    annotated   = vapply(parts, function(p) p[length(p) - 3], ""),
    significant = vapply(parts, function(p) p[length(p) - 2], ""),
    expected    = vapply(parts, function(p) p[length(p) - 1], ""),
    classic     = vapply(parts, function(p) p[length(p)], ""),
    n_fields    = n,
    stringsAsFactors = FALSE)
  # sanity of the both-ends parse
  stopifnot(all(grepl("^GO:[0-9]{7}$", d$go_id)))
  stopifnot(all(grepl("^[0-9]+$", d$annotated)), all(grepl("^[0-9]+$", d$significant)))
  stopifnot(all(grepl("^[0-9.eE+-]+$", d$expected)))
  d$annotated <- as.integer(d$annotated); d$significant <- as.integer(d$significant)
  d$expected_num <- as.numeric(d$expected)
  d$rank <- seq_len(nrow(d))
  attr(d, "header") <- hdr
  d
}

# Numeric value of the printed p; NA for "< 1e-30"-style entries.
classic_numeric <- function(x) suppressWarnings(as.numeric(x))

# What GenTable prints for a vector of p-values (topGO: format.pval, digits 2, eps 1e-30).
gentable_p <- function(p) format.pval(p, digits = 2, eps = 1e-30)

# What GenTable prints for the term name (numChar = 40).
short_term <- function(x, numChar = 40) paste0(substr(x, 1, numChar), ifelse(nchar(x) > numChar, "...", ""))

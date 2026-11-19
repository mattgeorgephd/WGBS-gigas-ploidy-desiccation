# sample_sheet.R
# The one place where sample groups are defined.

# Read the WGBS library sheet. Rows are in legacy order (legacy_sample_index 1 to
# 22), which is the column order of every methylation matrix in the pipeline.
#
# The sheet is checked before it is used, because a label that is not one of the
# expected ones would become NA and drop its library from a sample set without
# any error: the steps would then write counts for a smaller or an empty set.
# The function stops unless the sheet has 22 rows, a different library name and
# a different coverage file in every row, the sample indices 1 to 22, and no
# other label than those listed in `allowed`, each of them for at least one
# library. If a label of the sheet is ever renamed, change `allowed` here,
# sample_sets() below and the two comparisons with "triploid" and "heat" in
# fit_mode() of step 05 together.
read_libraries <- function() {
  allowed <- list(ploidy = c("diploid", "triploid"), treatment = c("control", "heat"))
  lib <- utils::read.csv(input_path("metadata", "wgbs_libraries.csv"), stringsAsFactors = FALSE)
  bad_sheet <- function(...) stop("wgbs_libraries.csv: ", sprintf(...), call. = FALSE)
  absent <- setdiff(c("library_id", "legacy_sample_index", "cov_file", names(allowed)), names(lib))
  if (length(absent) > 0) bad_sheet("column(s) not found: %s", paste(absent, collapse = ", "))
  if (nrow(lib) != 22) bad_sheet("expected 22 libraries, found %d row(s)", nrow(lib))
  for (column in c("library_id", "cov_file")) {
    if (anyNA(lib[[column]]) || anyDuplicated(lib[[column]])) bad_sheet("column '%s' must hold a different value in every row", column)
  }
  lib <- lib[order(lib$legacy_sample_index), ]
  if (!isTRUE(all.equal(lib$legacy_sample_index, 1:22))) {
    bad_sheet("column 'legacy_sample_index' must hold the whole numbers 1 to 22, each once")
  }
  for (column in names(allowed)) {
    quoted <- function(labels) paste(sprintf("'%s'", labels), collapse = ", ")
    unexpected <- setdiff(unique(lib[[column]]), allowed[[column]])
    if (length(unexpected) > 0) {
      bad_sheet("column '%s' may hold only %s; found %s in %d row(s)", column, quoted(allowed[[column]]),
                quoted(unexpected), sum(!lib[[column]] %in% allowed[[column]]))
    }
    unused <- setdiff(allowed[[column]], lib[[column]])
    if (length(unused) > 0) bad_sheet("column '%s' has no library labelled %s", column, quoted(unused))
    lib[[column]] <- factor(lib[[column]], levels = allowed[[column]])
  }
  lib$group     <- factor(paste(lib$ploidy, lib$treatment),
                          levels = c("diploid control", "diploid heat", "triploid control", "triploid heat"))
  rownames(lib) <- NULL
  lib
}

# Named list of row indices into the library sheet for the sample sets used in
# the coverage summaries.
sample_sets <- function(lib) {
  list(
    all      = seq_len(nrow(lib)),
    diploid  = which(lib$ploidy == "diploid"),
    triploid = which(lib$ploidy == "triploid"),
    control  = which(lib$treatment == "control"),
    heat     = which(lib$treatment == "heat")
  )
}

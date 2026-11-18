# sample_sheet.R
# The one place where sample groups are defined.

# Read the WGBS library sheet. Rows are in legacy order (legacy_sample_index 1 to
# 22), which is the column order of every methylation matrix in the pipeline.
read_libraries <- function() {
  lib <- utils::read.csv(input_path("metadata", "wgbs_libraries.csv"), stringsAsFactors = FALSE)
  lib <- lib[order(lib$legacy_sample_index), ]
  stopifnot(nrow(lib) == 22, !anyDuplicated(lib$library_id))
  lib$ploidy    <- factor(lib$ploidy, levels = c("diploid", "triploid"))
  lib$treatment <- factor(lib$treatment, levels = c("control", "heat"))
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

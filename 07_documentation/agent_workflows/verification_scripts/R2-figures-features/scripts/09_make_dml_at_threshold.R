# Emulate what step 05 writes for another value of config$dml_difference, from its own per-CpG tables
# (step 05 itself is not rerun). Same filters and merge as write_mode() in 05_differential_methylation.R.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE); root <- args[1]; th <- as.numeric(args[2])
tab <- function(f) file.path(root, "04_output/05_differential_methylation/tables", f)
for (mode in c("legacy", "corrected")) {
  dml <- list()
  for (contrast in c("ploidy", "heat")) {
    a <- fread(tab(sprintf("diffmeth_%s_%s.csv.gz", contrast, mode)))
    d <- a[qvalue < 0.01 & abs(meth.diff) > th]; d[, direction := fifelse(meth.diff > 0, "hyper", "hypo")]
    fwrite(d, tab(sprintf("dml_%s_%s.csv", contrast, mode))); dml[[contrast]] <- d
  }
  shared <- merge(dml$ploidy, dml$heat, by = c("chr", "start", "end", "strand", "united_row"), suffixes = c("_ploidy", "_heat"))
  shared[, class := fifelse(meth.diff_ploidy > 0 & meth.diff_heat > 0, "both_hyper", fifelse(meth.diff_ploidy < 0 & meth.diff_heat < 0, "both_hypo", "opposite"))]
  setorder(shared, united_row); fwrite(shared, tab(sprintf("dml_shared_%s.csv", mode)))
  cat(sprintf("threshold %g, %s: ploidy %d, heat %d, shared %d\n", th, mode, nrow(dml$ploidy), nrow(dml$heat), nrow(shared)))
}

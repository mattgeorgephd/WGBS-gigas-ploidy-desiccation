suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
for (mode in c("legacy", "corrected")) for (contrast in c("ploidy", "heat")) {
  d <- fread(t5(sprintf("diffmeth_%s_%s.csv.gz", contrast, mode)))
  cat(sprintf("[%s %s] |diff| > 20: %d; of these q < 0.01: %d; q >= 0.01: %d; largest q among |diff| > 20: %.3g; CpGs with q < 0.01 at any difference: %d (%.1f%%)\n", mode, contrast,
              sum(abs(d$meth.diff) > 20), sum(abs(d$meth.diff) > 20 & d$qvalue < 0.01), sum(abs(d$meth.diff) > 20 & d$qvalue >= 0.01),
              max(d$qvalue[abs(d$meth.diff) > 20]), sum(d$qvalue < 0.01), 100 * mean(d$qvalue < 0.01)))
}

suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t5 <- file.path(repo, "04_output/05_differential_methylation/tables")
source(file.path(repo, "02_R_functions/methylkit_helpers.R")); config <- list(dml_difference = 20)
for (md in c("legacy", "corrected")) {
  u <- fread(file.path(t5, sprintf("united_cpg_%s.csv.gz", md)))
  mp <- mean_percent_methylation(u); el <- dml_eligible(mp, config)
  cat(sprintf("[%s] tested %d ; eligible (mean 10-90%%): %d (%.1f%%) ; mean < 10%%: %.1f%%\n", md, nrow(u), sum(el), 100 * mean(el), 100 * mean(mp < 10)))
  for (ct in c("ploidy", "heat")) {
    d <- fread(file.path(t5, sprintf("dml_%s_%s.csv", ct, md)))
    out <- d[!el[united_row]]
    cat(sprintf("   %s DML: %d ; not 'eligible': %d", ct, nrow(d), nrow(out)))
    if (nrow(out)) cat(" -> ", paste(sprintf("%s:%d mean %.2f%% diff %.1f", out$chr, out$start, mp[out$united_row], out$meth.diff), collapse = " ; "))
    cat("\n")
  }
}

suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- file.path(repo, "04_output/05_differential_methylation/tables"); ar <- file.path(repo, "08_archive/tables")
for (ct in c("PLOIDY", "HEAT")) {
  old <- fread(file.path(ar, sprintf("DML-getMethylDiff-%s-Cov10-20.csv", ct))); setnames(old, 1, "row")
  new <- fread(file.path(t5, sprintf("dml_%s_legacy.csv", tolower(ct))))
  cat(sprintf("%s: legacy rows %d, pipeline rows %d ; same order of (row, chr, start, end): %s ; max |meth.diff diff| %.3g ; max rel p diff %.3g ; max rel q diff %.3g\n", ct, nrow(old), nrow(new),
              identical(old$row, new$united_row) && identical(old$chr, new$chr) && identical(old$start, new$start) && identical(old$end, new$end),
              max(abs(old$meth.diff - new$meth.diff)), max(abs(old$pvalue - new$pvalue) / old$pvalue), max(abs(old$qvalue - new$qvalue) / old$qvalue)))
}
old <- fread(file.path(ar, "DML-getMethylDiff-INTERACTION-Cov10-20.csv")); setnames(old, 1, "row")
new <- fread(file.path(t5, "dml_shared_legacy.csv"))
cat(sprintf("INTERACTION: legacy rows %d, pipeline %d ; legacy 'row' range %d..%d\n", nrow(old), nrow(new), min(old$row), max(old$row)))
m <- match(paste(old$chr, old$start), paste(new$chr, new$start))
cat("  all legacy loci in pipeline shared set:", !anyNA(m), "; legacy meth.diff equals pipeline heat difference: max abs diff", max(abs(old$meth.diff - new$meth.diff_heat[m])), "\n")
cat("  legacy 'row' equals united_row:", sum(old$row == new$united_row[m]), "of", nrow(old), "\n")
hd <- fread(file.path(t5, "dml_heat_legacy.csv"))
cat("  legacy 'row' equals the locus's position within the heat DML table (the legacy select() index):", sum(old$row == match(paste(old$chr, old$start), paste(hd$chr, hd$start))), "of", nrow(old), "\n")
print(head(old, 3)); print(head(new[m][, .(chr, start, united_row, meth.diff_heat)], 3))

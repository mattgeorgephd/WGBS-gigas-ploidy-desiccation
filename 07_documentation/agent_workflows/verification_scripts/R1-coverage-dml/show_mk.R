suppressPackageStartupMessages(library(methylKit))
show_method <- function(f, sig) { cat("\n=====", f, sig, "=====\n"); print(getMethod(f, sig)) }
show_method("filterByCoverage", "methylRaw")
show_method("normalizeCoverage", "methylRawList")
show_method("unite", "methylRawList")
cat("\n===== .procBismarkCoverage =====\n"); print(methylKit:::.procBismarkCoverage)
cat("\n===== .readTableFast =====\n"); print(methylKit:::.readTableFast)

suppressPackageStartupMessages(library(methylKit))
ns <- asNamespace("methylKit")
for (nm in c(".readTableFast")) { cat("\n-- ", nm, "\n"); try(print(get(nm, envir = ns))) }
cat("==== args(unite) ====\n"); print(args(unite))
cat("==== unite methylRawList ====\n"); print(getMethod("unite", "methylRawList")@.Data)
cat("==== normalizeCoverage methylRawList ====\n"); print(getMethod("normalizeCoverage", "methylRawList")@.Data)
cat("==== quantile default type ====\n"); print(formals(stats:::quantile.default)$type)

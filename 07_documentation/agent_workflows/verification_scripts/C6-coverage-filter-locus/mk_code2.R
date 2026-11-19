suppressPackageStartupMessages(library(methylKit))
cat("==== args(methRead) ====\n"); print(args(methRead))
cat("==== showMethods(methRead) ====\n"); showMethods("methRead")
ns <- asNamespace("methylKit")
cat("==== methRead character-signature method ====\n")
m <- findMethods("methRead")
for (i in seq_along(m)) { cat("\n#### signature: ", paste(m[[i]]@defined, collapse=", "), "\n"); print(m[[i]]@.Data) }
for (nm in c(".procBismarkCoverage", ".check.pipeline.list", ".structureAMPoutput", ".structureGeneric", "fread.gzipped")) {
  cat("\n-- ", nm, "\n"); try(print(get(nm, envir = ns)))
}

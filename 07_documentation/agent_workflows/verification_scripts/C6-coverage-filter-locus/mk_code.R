suppressPackageStartupMessages(library(methylKit))
cat("methylKit version: ", as.character(packageVersion("methylKit")), "\n")
cat("R: ", R.version.string, "\n\n")
cat("==== args(filterByCoverage) ====\n"); print(args(filterByCoverage))
cat("==== showMethods(filterByCoverage) ====\n"); showMethods("filterByCoverage")
cat("==== getMethod filterByCoverage methylRaw ====\n"); print(getMethod("filterByCoverage", "methylRaw"))
cat("==== getMethod filterByCoverage methylRawList ====\n"); print(getMethod("filterByCoverage", "methylRawList"))
cat("==== existsMethod methylRawDB ====\n"); print(existsMethod("filterByCoverage", "methylRawDB"))
cat("==== getMethod filterByCoverage methylRawDB ====\n"); try(print(getMethod("filterByCoverage", "methylRawDB")))
cat("==== getMethod filterByCoverage methylRawListDB ====\n"); try(print(getMethod("filterByCoverage", "methylRawListDB")))
cat("==== internal .filterByCoverage? ====\n")
for (nm in c(".filterByCoverage", ".CpG.dinuc.unify", ".setMethylDBNames")) {
  cat("--", nm, "\n"); try(print(get(nm, envir = asNamespace("methylKit"))))
}

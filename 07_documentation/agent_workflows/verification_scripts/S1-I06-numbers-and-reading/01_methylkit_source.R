# Print the source of the functions that calculateDiffMeth dispatches to (methylKit 1.28.0).
suppressPackageStartupMessages(library(methylKit))
cat("methylKit version:", as.character(packageVersion("methylKit")), "\n")
cat("\n==== logReg ====\n"); print(methylKit:::logReg)
cat("\n==== calculateDiffMeth methylBase method ====\n")
print(selectMethod("calculateDiffMeth", "methylBase"))
cat("\n==== .calculateDiffMeth? ====\n")
for (f in c(".calculateDiffMeth", "fast.fisher", "SLIMfunc", "QValuesfun", "mc.fish", "estimateShrinkageMN", "glm.fit.fast")) {
  if (exists(f, envir = asNamespace("methylKit"), inherits = FALSE)) { cat("\n==== ", f, " ====\n"); print(get(f, envir = asNamespace("methylKit"))) }
}

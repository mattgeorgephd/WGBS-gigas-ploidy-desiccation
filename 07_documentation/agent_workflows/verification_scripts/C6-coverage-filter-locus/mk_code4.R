suppressPackageStartupMessages(library(methylKit))
ns <- asNamespace("methylKit")
cat("==== args(calculateDiffMeth) ====\n"); print(args(calculateDiffMeth))
cat("==== calculateDiffMeth methylBase ====\n"); print(getMethod("calculateDiffMeth", "methylBase")@.Data)
for (nm in c("logReg", "p.adjusted", "SLIMfunc", "QValuesfun")) { cat("\n-- ", nm, "\n"); try(print(get(nm, envir = ns))) }
cat("==== getMethylDiff methylDiff ====\n"); print(getMethod("getMethylDiff", "methylDiff")@.Data)

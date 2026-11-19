# S4 check: structure of a cached fit (read-only)
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
root <- "/home/user/WGBS-gigas-ploidy-desiccation"
f <- readRDS(file.path(root, "06_diagnostics/dml_sensitivity/cache/fit_A2_interaction_mn_f_legacy.rds"))
print(names(f)); str(f, max.level = 1)
print(head(f$fit))
print(f[setdiff(names(f), "fit")])

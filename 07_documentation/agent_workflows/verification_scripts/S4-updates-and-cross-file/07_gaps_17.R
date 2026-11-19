# S4 check: margin by which the "apart" group is separated at the 17 CpGs of the heat-cohort corrected list,
# from the united table (unrounded percentages), and at NC_047565.1:35030701.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
root <- "/home/user/WGBS-gigas-ploidy-desiccation"
out <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S4-updates-and-cross-file"
load(file.path(out, "united_legacy_matrices.RData"))
h <- fread(file.path(root, "06_diagnostics/dml_sensitivity/dml_A4_ploidy_in_heat_mn_f_legacy.csv"))
k <- match(paste(h$chr, h$start), paste(chr, pos)); stopifnot(!anyNA(k), nrow(h) == 125)
M <- P[k, ]
gap <- function(M, S) { a <- M[, S, drop = FALSE]; b <- M[, -S, drop = FALSE]
  hi <- apply(a, 1, min) - apply(b, 1, max); lo <- apply(b, 1, min) - apply(a, 1, max); data.table(higher = hi, lower = lo) }
gt <- gap(M, 6:10); gd <- gap(M, 1:5)
tA <- gt$higher > 0 | gt$lower > 0; dA <- gd$higher > 0 | gd$lower > 0
cat("triploid heat apart:", sum(tA), "(higher", sum(gt$higher > 0), ", lower", sum(gt$lower > 0), ") ; diploid heat apart:", sum(dA), "(higher", sum(gd$higher > 0), ", lower", sum(gd$lower > 0), ")\n")
cat("gaps, triploid apart:", paste(round(sort(pmax(gt$higher, gt$lower)[tA]), 2), collapse = " "), "\n")
cat("gaps, diploid apart:", paste(round(sort(pmax(gd$higher, gd$lower)[dA]), 2), collapse = " "), "\n")
g <- c(pmax(gt$higher, gt$lower)[tA], pmax(gd$higher, gd$lower)[dA])
cat("of the 17: gap below 5 points:", sum(g < 5), "(largest of these", round(max(g[g < 5]), 2), "); gap above 10 points:", sum(g > 10), "(range", round(min(g[g > 10]), 1), "to", round(max(g[g > 10]), 1), ")\n")
i <- which(chr == "NC_047565.1" & pos == 35030701); x <- P[i, ]
cat(sprintf("NC_047565.1:35030701: triploid heat %.2f to %.2f ; largest other %.2f ; gap %.1f\n", min(x[6:10]), max(x[6:10]), max(x[-(6:10)]), min(x[6:10]) - max(x[-(6:10)])))
# dispersion-free description: how many of the 125 have the 12 control libraries all inside the range spanned by the ten heat libraries' two groups, etc. (not needed)
# corrected-mode pooled lists: is NC_047565.1:35030701 in both?
a <- fread(file.path(root, "06_diagnostics/dml_sensitivity/dml_A1_ploidy_mn_f_corrected.csv")); b <- fread(file.path(root, "06_diagnostics/dml_sensitivity/dml_A1_heat_mn_f_corrected.csv"))
cat("corrected mode: CpGs in both corrected lists:", paste(intersect(paste(a$chr, a$start), paste(b$chr, b$start)), collapse = "; "), "\n")

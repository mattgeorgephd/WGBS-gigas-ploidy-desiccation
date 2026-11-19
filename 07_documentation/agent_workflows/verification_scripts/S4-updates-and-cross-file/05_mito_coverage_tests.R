# S4 check: relative coverage of the mitochondrial CpGs by cohort and ploidy (exploratory by-product of 04).
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
out <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S4-updates-and-cross-file"
load(file.path(out, "united_legacy_matrices.RData"))
mt <- chr == "NC_001276.1"
nuc <- !mt
rel <- colMeans(N[mt, ]) / apply(N[nuc, ], 2, median)
rel2 <- colSums(N[mt, ]) / colSums(N[nuc, ]) * sum(nuc) / sum(mt)     # mean over mean
names(rel) <- names(rel2) <- lib$library_id
print(round(rbind(mean_over_median = rel, mean_over_mean = rel2), 2))
w <- function(a, b) wilcox.test(a, b, exact = TRUE)$p.value
cat(sprintf("cohort (heat against control): p = %.2g ; ranges %.1f to %.1f and %.1f to %.1f\n", w(rel[1:10], rel[11:22]), min(rel[1:10]), max(rel[1:10]), min(rel[11:22]), max(rel[11:22])))
cat(sprintf("ploidy inside the heat cohort: diploid %.1f to %.1f, triploid %.1f to %.1f, exact Wilcoxon p = %.3g\n", min(rel[1:5]), max(rel[1:5]), min(rel[6:10]), max(rel[6:10]), w(rel[1:5], rel[6:10])))
cat(sprintf("ploidy inside the control cohort: diploid %.1f to %.1f, triploid %.1f to %.1f, exact Wilcoxon p = %.3g\n", min(rel[11:16]), max(rel[11:16]), min(rel[17:22]), max(rel[17:22]), w(rel[11:16], rel[17:22])))
cat(sprintf("same with mean over mean: heat cohort p = %.3g ; control cohort p = %.3g\n", w(rel2[1:5], rel2[6:10]), w(rel2[11:16], rel2[17:22])))
mp <- 100 * colSums(C[mt, ]) / colSums(N[mt, ])
cat(sprintf("mitochondrial percent methylation by ploidy: heat cohort p = %.3g ; control cohort p = %.3g\n", w(mp[1:5], mp[6:10]), w(mp[11:16], mp[17:22])))

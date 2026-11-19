# Would the permuted-label "DML" also pass the q < 0.01 filter? LRT p-values (methylKit's test, no overdispersion) for them.
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; sc <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory"
set.seed(4101)   # same seed and stream as 05_perm_counts.R, so the first permutations are the same label sets
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"; stratum <- paste(lib$treatment, lib$acute_heat_shock)
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); W <- C + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); storage.mode(C) <- "double"; storage.mode(W) <- "double"; rm(u); gc(verbose = FALSE)
wdiff <- function(g) 100 * (rowSums(C[, g, drop = FALSE]) / rowSums(W[, g, drop = FALSE]) - rowSums(C[, !g, drop = FALSE]) / rowSums(W[, !g, drop = FALSE]))
lrt_p <- function(rows, g) vapply(rows, function(r) { y <- cbind(C[r, ], W[r, ] - C[r, ]); f1 <- suppressWarnings(glm(y ~ g + heat, binomial)); f0 <- suppressWarnings(glm(y ~ heat, binomial)); pchisq(f0$deviance - f1$deviance, 1, lower.tail = FALSE) }, numeric(1))
dq <- fread(file.path(repo, "04_output/05_differential_methylation/tables/diffmeth_ploidy_legacy.csv.gz"))
real <- which(abs(dq$meth.diff) > 20); chk <- sample(real, 200); pp <- lrt_p(chk, trip)
cat(sprintf("my LRT against the pipeline p-values on 200 real ploidy DML: max |log10 ratio| = %.2e\n", max(abs(log10(pp) - log10(dq$pvalue[chk])))))
cat(sprintf("real ploidy DML p-values: median %.2g; 90th pct %.2g; max %.2g; p cutoff at q<0.01: %.2g; CpGs with q<0.01: %d\n", median(dq$pvalue[real]), quantile(dq$pvalue[real], .9), max(dq$pvalue[real]), max(dq[qvalue < 0.01, pvalue]), sum(dq$qvalue < 0.01)))
set.seed(4101)
for (b in 1:6) {
  g <- logical(22); for (s in unique(stratum)) { i <- which(stratum == s); g[sample(i, sum(trip[i]))] <- TRUE }
  hit <- which(abs(wdiff(g)) > 20); p <- lrt_p(hit, g)
  cat(sprintf("permutation %d: %d CpGs pass 20 points; LRT p: median %.2g; 90th pct %.2g; max %.2g; share above 2.3e-4: %.4f; above 1e-5: %.4f\n", b, length(hit), median(p), quantile(p, .9), max(p), mean(p > 2.3e-4), mean(p > 1e-5)))
}

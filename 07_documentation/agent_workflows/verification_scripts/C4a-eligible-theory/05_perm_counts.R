# Label-permutation null for the ploidy contrast: how many CpGs pass the 20-point rule when ploidy labels
# are shuffled within cohort (and within the acute-heat-shock stratum)? And where do they fall?
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
sc <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory"
mode <- "legacy"; set.seed(4101)
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
stratum <- paste(lib$treatment, lib$acute_heat_shock)
print(table(stratum, lib$ploidy))
u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); W <- C + as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
storage.mode(C) <- "double"; storage.mode(W) <- "double"
pos <- u[, .(chr, start)]; rm(u); gc(verbose = FALSE)
m <- rowMeans(100 * C / W); elig <- m >= 10 & m <= 90
wdiff <- function(g) 100 * (rowSums(C[, g, drop = FALSE]) / rowSums(W[, g, drop = FALSE]) - rowSums(C[, !g, drop = FALSE]) / rowSums(W[, !g, drop = FALSE]))
obs <- wdiff(trip); cat(sprintf("observed ploidy: |weighted diff| > 20 at %d CpGs (%d eligible)\n", sum(abs(obs) > 20), sum(abs(obs) > 20 & elig)))
dq <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/diffmeth_ploidy_%s.csv.gz", mode)))
cat(sprintf("real data: largest p-value with q < 0.01: %.3g; smallest p-value with q >= 0.01: %.3g; CpGs with q<0.01: %d (%.1f%%); with |diff|>20: %d; with both: %d\n", max(dq[qvalue < 0.01, pvalue]), min(dq[qvalue >= 0.01, pvalue]),
            sum(dq$qvalue < 0.01), 100 * mean(dq$qvalue < 0.01), sum(abs(dq$meth.diff) > 20), sum(dq$qvalue < 0.01 & abs(dq$meth.diff) > 20)))
pcut <- max(dq[qvalue < 0.01, pvalue]); rm(dq)
perm_labels <- function(stratified = TRUE) {
  g <- logical(22)
  if (stratified) for (s in unique(stratum)) { i <- which(stratum == s); g[sample(i, sum(trip[i]))] <- TRUE }
  else for (s in unique(lib$treatment)) { i <- which(lib$treatment == s); g[sample(i, sum(trip[i]))] <- TRUE }
  g
}
lrt_p <- function(rows, g) vapply(rows, function(r) {
  y <- cbind(C[r, ], W[r, ] - C[r, ])
  f1 <- suppressWarnings(glm(y ~ g + heat, binomial)); f0 <- suppressWarnings(glm(y ~ heat, binomial))
  pchisq(f0$deviance - f1$deviance, 1, lower.tail = FALSE)
}, numeric(1))
B <- 300; out <- vector("list", B); sets <- vector("list", B)
for (b in seq_len(B)) {
  g <- perm_labels(TRUE)
  if (all(g == trip) || all(g == !trip)) next
  d <- wdiff(g); hit <- which(abs(d) > 20)
  pmaxp <- NA_real_
  if (b <= 10) { pp <- lrt_p(hit, g); pmaxp <- max(pp) }
  out[[b]] <- data.table(perm = b, n_pass = length(hit), n_pass_eligible = sum(elig[hit]), n_hyper = sum(d[hit] > 0), overlap_with_observed = sum(abs(obs[hit]) > 20),
                         same_as_true = sum(g == trip), max_lrt_p = pmaxp)
  sets[[b]] <- hit
  if (b %% 50 == 0) cat("perm", b, "\n")
}
out <- rbindlist(out)
cat(sprintf("\nStratified permutations: %d. CpGs passing |diff|>20: median %d, mean %.0f, range %d to %d; 5%%/95%% quantiles %d / %d\n", nrow(out), as.integer(median(out$n_pass)), mean(out$n_pass), min(out$n_pass), max(out$n_pass),
            as.integer(quantile(out$n_pass, .05)), as.integer(quantile(out$n_pass, .95))))
cat(sprintf("share of permutations with at least as many as observed (1038): %.3f\n", mean(out$n_pass >= sum(abs(obs) > 20))))
cat(sprintf("largest LRT p-value among permuted 'DML' in the first 10 permutations: %.3g (real-data p cutoff at q<0.01: %.3g)\n", max(out$max_lrt_p, na.rm = TRUE), pcut))
cat(sprintf("share of passing CpGs outside the eligible range: %.4f\n", 1 - sum(out$n_pass_eligible) / sum(out$n_pass)))
cat("by number of libraries keeping their true label (22 = identity, 11 = fully mixed):\n"); print(out[, .(perms = .N, median_pass = median(n_pass)), by = same_as_true][order(same_as_true)])
# unstratified version (within cohort only)
out2 <- rbindlist(lapply(seq_len(200), function(b) { g <- perm_labels(FALSE); d <- wdiff(g); data.table(n_pass = sum(abs(d) > 20), same_as_true = sum(g == trip)) }))
cat(sprintf("\nWithin-cohort permutations (acute shock not held balanced): median %d, mean %.0f, range %d to %d; share >= observed %.3f\n", as.integer(median(out2$n_pass)), mean(out2$n_pass), min(out2$n_pass), max(out2$n_pass), mean(out2$n_pass >= sum(abs(obs) > 20))))
saveRDS(list(summary = out, sets = sets, pos = pos, elig = elig), file.path(sc, "perm_sets_ploidy_legacy.rds"))

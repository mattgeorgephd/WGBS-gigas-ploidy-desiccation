# Task (a): recompute the numbers of I-06 (Rating to end) from the cached fits and the united tables.
# Read-only on the repository. One thread.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
stopifnot(identical(lib$legacy_sample_index, 1:22))
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
cache <- function(a, m) readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_%s.rds", a, m)))$fit
ana <- fread(file.path(R, "06_diagnostics/dml_sensitivity/analyses.csv"))
cnt <- fread(file.path(R, "06_diagnostics/dml_sensitivity/dml_counts.csv"))
out <- list()
for (mode in c("legacy", "corrected")) {
  u <- fread(file.path(R, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
  Cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); Ts <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
  pm <- 100 * Cs / (Cs + Ts)
  mean_pm <- rowMeans(pm)
  elig <- mean_pm >= 10 & mean_pm <= 90
  cat(sprintf("\n[%s] CpGs %d; eligible %d\n", mode, nrow(u), sum(elig)))
  pos <- u[, .(chr, start)]
  saveRDS(list(pos = pos, elig = elig, mean_pm = mean_pm), file.path(S, sprintf("pos_elig_%s.rds", mode)))
  pooled <- function(cols) 100 * rowSums(Cs[, cols, drop = FALSE]) / (rowSums(Cs[, cols, drop = FALSE]) + rowSums(Ts[, cols, drop = FALSE]))
  # interaction contrast as the script defines it
  inter <- (pooled(which(trip & heat)) - pooled(which(trip & !heat))) - (pooled(which(!trip & heat)) - pooled(which(!trip & !heat)))
  ids <- ana[grepl(mode, modes), analysis_id]
  for (a in ids) {
    f <- cache(a, mode); stopifnot(nrow(f) == nrow(u))
    eff <- if (grepl("interaction", a)) inter else f$meth.diff
    sig <- f$qvalue < 0.01; dml <- sig & abs(eff) > 20
    row <- data.table(analysis_id = a, mode = mode, n = nrow(f), n_above = sum(abs(eff) > 20), n_sig = sum(sig), n_dml = sum(dml),
                      n_hyper = sum(dml & eff > 0), n_hypo = sum(dml & eff < 0),
                      p05 = mean(f$pvalue[elig] < 0.05), p01 = mean(f$pvalue[elig] < 0.01), p001 = mean(f$pvalue[elig] < 0.001),
                      med_eff_elig = median(eff[elig]),
                      # counts behind the shares, and the same shares over ALL tested CpGs
                      n05 = sum(f$pvalue[elig] < 0.05), n01 = sum(f$pvalue[elig] < 0.01), n001 = sum(f$pvalue[elig] < 0.001),
                      all05 = mean(f$pvalue < 0.05), all01 = mean(f$pvalue < 0.01), all001 = mean(f$pvalue < 0.001))
    out[[length(out) + 1]] <- row
    # the DML list file
    lf <- fread(file.path(R, "06_diagnostics/dml_sensitivity", sprintf("dml_%s_%s.csv", a, mode)))
    same <- nrow(lf) == sum(dml) && identical(lf$chr, pos$chr[dml]) && identical(lf$start, pos$start[dml])
    if (!same) cat("  LIST FILE DIFFERS:", a, mode, "\n")
  }
  rm(u, Cs, Ts, pm); gc(verbose = FALSE)
}
res <- rbindlist(out)
m <- merge(res, cnt, by = c("analysis_id", "mode"), sort = FALSE)
cat("\nAgreement of my recomputation with dml_counts.csv (all must be TRUE):\n")
print(m[, .(n = all(n == n_cpg_tested), above = all(n_above == n_effect_above_cutoff), sig = all(n_sig == n_significant), dml = all(n_dml.x == n_dml.y),
            hyper = all(n_hyper.x == n_hyper.y), hypo = all(n_hypo.x == n_hypo.y),
            p05 = max(abs(p05 - share_eligible_p_below_0_05)), p01 = max(abs(p01 - share_eligible_p_below_0_01)), p001 = max(abs(p001 - share_eligible_p_below_0_001)),
            med = max(abs(med_eff_elig - median_effect_size_eligible)))])
options(width = 250)
print(res[, .(analysis_id, mode, n_above, n_sig, n_dml, n_hyper, n_hypo, p05 = round(100 * p05, 2), p01 = round(100 * p01, 2), p001 = round(100 * p001, 3),
              n05, n01, n001, med_eff_elig = round(med_eff_elig, 2), all05 = round(100 * all05, 2), all01 = round(100 * all01, 2), all001 = round(100 * all001, 3))])
fwrite(res, file.path(S, "03_numbers_from_cache.csv"))
# statistic needed for a single call
cat("\nStatistic needed for a single call (q cutoff / number tested):\n")
for (mode in c("legacy", "corrected")) {
  mm <- c(legacy = 480579, corrected = 454555)[[mode]]; p1 <- 0.01 / mm
  cat(sprintf("[%s] p = %.4g; chisq(1) %.3f; F(1,19) %.3f; F(1,18) %.3f; F(1,10) %.3f; F(1,8) %.3f\n", mode, p1,
              qchisq(p1, 1, lower.tail = FALSE), qf(p1, 1, 19, lower.tail = FALSE), qf(p1, 1, 18, lower.tail = FALSE),
              qf(p1, 1, 10, lower.tail = FALSE), qf(p1, 1, 8, lower.tail = FALSE)))
}
cat("choose(10,5)*choose(12,6) =", choose(10, 5) * choose(12, 6), "; choose(12,6) =", choose(12, 6), "; choose(10,5) =", choose(10, 5),
    "; without acute: choose(6,3)*choose(12,6) =", choose(6, 3) * choose(12, 6), "\n")

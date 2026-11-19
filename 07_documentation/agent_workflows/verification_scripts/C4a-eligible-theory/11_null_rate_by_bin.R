suppressPackageStartupMessages(library(data.table)); options(width = 200)
sc <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory"
ps <- readRDS(file.path(sc, "perm_sets_ploidy_legacy.rds")); CP <- readRDS(file.path(sc, "cpg_table_legacy.rds"))
sets <- ps$sets[lengths(ps$sets) > 0]; B <- length(sets)
hits <- tabulate(unlist(sets), nbins = nrow(CP))          # times each CpG passed the 20-point rule over B permutations
CP[, null_rate := hits / B]
cat(sprintf("permutations: %d; CpGs passing in at least one permutation: %d; mean number passing per permutation: %.0f\n", B, sum(hits > 0), sum(hits) / B))
CP[, mean_bin := cut(mean_pm, c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100), right = FALSE, include.lowest = TRUE)]
cat("\nNull (label-permutation) probability of passing the 20-point rule, and the real ploidy DML rate, by mean methylation:\n")
print(CP[, .(cpgs = .N, null_pct = round(100 * mean(null_rate), 3), real_ploidy_pct = round(100 * mean(dml_ploidy), 3), real_heat_pct = round(100 * mean(dml_heat), 3), expected_null_dml = round(sum(null_rate), 1), real_ploidy = sum(dml_ploidy)), by = mean_bin][order(mean_bin)])
el <- CP[eligible == TRUE]
el[, sd_bin := cut(sd_within, c(0, 5, 10, 15, 20, 25, 30, 100), right = FALSE)]
cat("\nSame, by within-cell SD (eligible CpGs):\n")
print(el[, .(cpgs = .N, null_pct = round(100 * mean(null_rate), 3), real_ploidy_pct = round(100 * mean(dml_ploidy), 3), real_heat_pct = round(100 * mean(dml_heat), 3), expected_null = round(sum(null_rate), 1), real_ploidy = sum(dml_ploidy), real_heat = sum(dml_heat)), by = sd_bin][order(sd_bin)])
cat("\nNull pass probability by genomic feature within the eligible set (this is what a 'fair' background would have to match):\n")
for (t in c("gene", "CDS", "intron", "TE", "intergenic")) cat(sprintf("  %-10s in feature: %.3f%%   outside: %.3f%%   ratio %.2f\n", t, 100 * mean(el$null_rate[el[[t]]]), 100 * mean(el$null_rate[!el[[t]]]), mean(el$null_rate[el[[t]]]) / mean(el$null_rate[!el[[t]]])))
# real ploidy DML: observed share in feature vs share expected from null rates (weighted background)
cat("\nShare of loci in each feature: real ploidy DML, eligible CpGs (unweighted), eligible CpGs weighted by null pass probability:\n")
for (t in c("gene", "CDS", "intron", "TE", "intergenic")) cat(sprintf("  %-10s real ploidy DML %.1f%%; real heat DML %.1f%%; eligible %.1f%%; null-weighted eligible %.1f%%\n", t, 100 * mean(CP[dml_ploidy == TRUE][[t]]), 100 * mean(CP[dml_heat == TRUE][[t]]), 100 * mean(el[[t]]), 100 * weighted.mean(el[[t]], el$null_rate)))
cat(sprintf("\ncorrelation of null pass rate with real DML status among eligible: ploidy AUC-like: mean null rate at ploidy DML %.3f vs at other eligible %.4f; at heat DML %.3f\n", mean(el$null_rate[el$dml_ploidy]), mean(el$null_rate[!el$dml_ploidy]), mean(el$null_rate[el$dml_heat])))
cat(sprintf("real ploidy DML that passed in >= 10%% of permutations: %d of %d; that never passed in %d permutations: %d\n", sum(CP$null_rate[CP$dml_ploidy] >= 0.10), sum(CP$dml_ploidy), B, sum(CP$null_rate[CP$dml_ploidy] == 0)))
# heat DML direction near the lower bound
h <- fread("/home/user/WGBS-gigas-ploidy-desiccation/04_output/05_differential_methylation/tables/dml_heat_legacy.csv"); h[, mean_pm := CP$mean_pm[united_row]]
h[, mb := cut(mean_pm, c(0, 10, 11, 12, 15, 20, 50, 80, 85, 88, 89, 90, 100), right = FALSE)]
cat("\nHeat DML direction by mean methylation (10 heat v 12 control: a hyper DML needs mean >= 9.09, a hypo DML >= 10.91; mirrored at the top):\n"); print(dcast(h[, .N, by = .(mb, direction)], mb ~ direction, value.var = "N", fill = 0))

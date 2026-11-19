# Where methylKit's "MN" estimate of the dispersion is not above 1, logReg sets it to 1 and uses the chi-squared
# test: the "corrected" p-value is then the uncorrected p-value. How many of the corrected DML are of that kind?
suppressPackageStartupMessages(library(data.table)); setDTthreads(1); options(width = 220)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
cache <- function(a, m = "legacy") readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_%s.rds", a, m)))$fit
pe <- readRDS(file.path(S, "pos_elig_legacy.rds"))
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0; ns <- cs + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); pos <- u[, .(chr, start)]; rm(u)
pm <- 100 * cs / ns
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
grp <- list(dc = which(!trip & !heat), dh = which(!trip & heat), tc = which(trip & !heat), th = which(trip & heat))
pairs <- list(c("A1_ploidy_mn_f", "A0_ploidy"), c("A1_heat_mn_f", "A0_heat"), c("A2_interaction_mn_f", "A2_interaction"),
              c("A4_ploidy_in_control_mn_f", "A4_ploidy_in_control"), c("A4_ploidy_in_heat_mn_f", "A4_ploidy_in_heat"))
pool <- function(cols) 100 * rowSums(cs[, cols, drop = FALSE]) / rowSums(ns[, cols, drop = FALSE])
inter <- (pool(grp$th) - pool(grp$tc)) - (pool(grp$dh) - pool(grp$dc))
for (p in pairs) {
  f <- cache(p[1]); f0 <- cache(p[2])
  eff <- if (grepl("interaction", p[1])) inter else f$meth.diff
  same <- f$pvalue == f0$pvalue                        # the correction changed nothing: dispersion estimate not above 1
  sig <- f$qvalue < 0.01; dml <- sig & abs(eff) > 20
  libs <- if (grepl("in_control", p[1])) which(!heat) else if (grepl("in_heat", p[1])) which(heat) else 1:22
  # a group of the analysis whose libraries are all at 0% or all at 100%
  gs <- if (grepl("in_control", p[1])) grp[c("dc", "tc")] else if (grepl("in_heat", p[1])) grp[c("dh", "th")] else if (grepl("ploidy", p[1])) list(which(!trip), which(trip)) else
        if (grepl("heat", p[1])) list(which(!heat), which(heat)) else grp
  flat <- Reduce(`|`, lapply(gs, function(g) rowSums(cs[, g, drop = FALSE]) == 0 | rowSums(cs[, g, drop = FALSE]) == rowSums(ns[, g, drop = FALSE])))
  cat(sprintf("%-28s all CpGs with p unchanged by the correction: %d (%.2f%%); eligible: %.2f%% | q < 0.01: %d, of them unchanged %d | DML: %d, of them unchanged %d; DML with a group entirely at 0%% or 100%%: %d\n",
              p[1], sum(same), 100 * mean(same), 100 * mean(same[pe$elig]), sum(sig), sum(sig & same), sum(dml), sum(dml & same), sum(dml & flat)))
  if (sum(dml) <= 20) print(cbind(pos[dml], p_corrected = f$pvalue[dml], p_uncorrected = f0$pvalue[dml], unchanged = same[dml], effect = round(eff[dml], 1), group_all_0_or_100 = flat[dml]))
}
# the heat-cohort corrected list in more detail: how many of the 125 have a group of five entirely at 0% or 100%, or are called through the chi-squared route
f <- cache("A4_ploidy_in_heat_mn_f"); f0 <- cache("A4_ploidy_in_heat"); dml <- f$qvalue < 0.01 & abs(f$meth.diff) > 20; same <- f$pvalue == f0$pvalue
th0 <- rowSums(cs[, grp$th]) == 0; dh0 <- rowSums(cs[, grp$dh]) == 0; th1 <- rowSums(cs[, grp$th]) == rowSums(ns[, grp$th]); dh1 <- rowSums(cs[, grp$dh]) == rowSums(ns[, grp$dh])
cat("\nheat-cohort corrected list (125): triploids all 0%:", sum(dml & th0), "; diploids all 0%:", sum(dml & dh0), "; triploids all 100%:", sum(dml & th1), "; diploids all 100%:", sum(dml & dh1),
    "; called with the chi-squared reference (dispersion estimate not above 1):", sum(dml & same), "; called with the F reference:", sum(dml & !same), "\n")
cat("  F statistic needed with q < 0.01 when 220 CpGs are significant together: p about", signif(max(f$pvalue[f$qvalue < 0.01]), 3), "-> F(1,8) of", round(qf(max(f$pvalue[f$qvalue < 0.01]), 1, 8, lower.tail = FALSE), 1),
    "or chi-squared of", round(qchisq(max(f$pvalue[f$qvalue < 0.01]), 1, lower.tail = FALSE), 1), "\n")
cat("  of the 220 significant: unchanged", sum(f$qvalue < 0.01 & same), "\n")
# the same question for all eligible CpGs: how often is a 5-library group entirely 0 or 100
cat("  eligible CpGs with a heat-cohort ploidy group entirely at 0% or 100%:", sum(pe$elig & (th0 | dh0 | th1 | dh1)), "of", sum(pe$elig), "\n")

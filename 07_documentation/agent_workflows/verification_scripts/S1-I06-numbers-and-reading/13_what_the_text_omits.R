# Task (d) and (c)5: things in the caches and tables that the text of I-06 does not report.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1); options(width = 220)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
cache <- function(a, m = "legacy") readRDS(file.path(R, "06_diagnostics/dml_sensitivity/cache", sprintf("fit_%s_%s.rds", a, m)))$fit
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0; ns <- cs + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); pos <- u[, .(chr, start)]; rm(u)
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
pool <- function(cols) 100 * rowSums(cs[, cols, drop = FALSE]) / rowSums(ns[, cols, drop = FALSE])
d_heat <- pool(which(trip & heat)) - pool(which(!trip & heat)); d_ctrl <- pool(which(trip & !heat)) - pool(which(!trip & !heat))
a0 <- cache("A0_ploidy"); dml <- a0$qvalue < 0.01 & abs(a0$meth.diff) > 20
cat("draft's ploidy DML:", sum(dml), "\n")
x <- data.table(pooled = a0$meth.diff[dml], heat = d_heat[dml], ctrl = d_ctrl[dml])
cat("  within the heat cohort: same sign as the pooled difference", sum(sign(x$heat) == sign(x$pooled)), "; above 20 points with that sign", sum(sign(x$heat) == sign(x$pooled) & abs(x$heat) > 20),
    "; median |difference|", round(median(abs(x$heat)), 1), "\n")
cat("  within the control cohort: same sign as the pooled difference", sum(sign(x$ctrl) == sign(x$pooled)), "; above 20 points with that sign", sum(sign(x$ctrl) == sign(x$pooled) & abs(x$ctrl) > 20),
    "; opposite sign", sum(sign(x$ctrl) == -sign(x$pooled)), "; median |difference|", round(median(abs(x$ctrl)), 1), "\n")
cat("  above 20 points with the pooled sign in both cohorts:", sum(sign(x$heat) == sign(x$pooled) & abs(x$heat) > 20 & sign(x$ctrl) == sign(x$pooled) & abs(x$ctrl) > 20),
    "; in the heat cohort only:", sum(sign(x$heat) == sign(x$pooled) & abs(x$heat) > 20 & !(sign(x$ctrl) == sign(x$pooled) & abs(x$ctrl) > 20)),
    "; in the control cohort only:", sum(!(sign(x$heat) == sign(x$pooled) & abs(x$heat) > 20) & sign(x$ctrl) == sign(x$pooled) & abs(x$ctrl) > 20),
    "; in neither:", sum(!(sign(x$heat) == sign(x$pooled) & abs(x$heat) > 20) & !(sign(x$ctrl) == sign(x$pooled) & abs(x$ctrl) > 20)), "\n")
cat("  larger |difference| in the heat cohort than in the control cohort:", sum(abs(x$heat) > abs(x$ctrl)), "\n")
# DML of the within-cohort uncorrected fits among the draft's list
h <- cache("A4_ploidy_in_heat"); cc <- cache("A4_ploidy_in_control")
hd <- h$qvalue < 0.01 & abs(h$meth.diff) > 20; cd <- cc$qvalue < 0.01 & abs(cc$meth.diff) > 20
cat("  of the 1,038: DML in the heat-only fit", sum(dml & hd), "; in the control-only fit", sum(dml & cd), "; in both", sum(dml & hd & cd), "; in neither", sum(dml & !hd & !cd),
    "; in both with the same sign", sum(dml & hd & cd & sign(h$meth.diff) == sign(cc$meth.diff)), "\n")
cat("  heat-only DML:", sum(hd), "; control-only DML:", sum(cd), "; in both:", sum(hd & cd), "; of these with the same sign:", sum(hd & cd & sign(h$meth.diff) == sign(cc$meth.diff)), "\n")
# corrected lists
hf <- cache("A4_ploidy_in_heat_mn_f"); cf <- cache("A4_ploidy_in_control_mn_f"); pf_ <- cache("A1_ploidy_mn_f"); tf <- cache("A1_heat_mn_f"); intf <- cache("A2_interaction_mn_f")
hfd <- hf$qvalue < 0.01 & abs(hf$meth.diff) > 20; cfd <- cf$qvalue < 0.01 & abs(cf$meth.diff) > 20
pfd <- pf_$qvalue < 0.01 & abs(pf_$meth.diff) > 20; tfd <- tf$qvalue < 0.01 & abs(tf$meth.diff) > 20
cat("\ncorrected heat-cohort ploidy list:", sum(hfd), "; hyper", sum(hfd & hf$meth.diff > 0), "; hypo", sum(hfd & hf$meth.diff < 0), "; among the draft's 1,038:", sum(hfd & dml),
    "; among the uncorrected heat-only list:", sum(hfd & hd), "; contains the two corrected pooled ploidy CpGs:", sum(hfd & pfd), "; contains corrected heat (cohort) CpGs:", sum(hfd & tfd), "\n")
bh <- p.adjust(hf$pvalue, "BH")
cat("  with Benjamini-Hochberg:", sum(bh < 0.01 & abs(hf$meth.diff) > 20), "; q-values of the 125: range", signif(range(hf$qvalue[hfd]), 3), "; p-values: largest", signif(max(hf$pvalue[hfd]), 3),
    "; SLIM share of true nulls", round(median((hf$qvalue / bh)[bh < 1]), 4), "\n")
cat("  with q < 0.05 instead of 0.01 (same 20 points):", sum(hf$qvalue < 0.05 & abs(hf$meth.diff) > 20), "; pooled ploidy with q < 0.05:", sum(pf_$qvalue < 0.05 & abs(pf_$meth.diff) > 20),
    "; with q < 0.1:", sum(pf_$qvalue < 0.1 & abs(pf_$meth.diff) > 20), "; cohort contrast q < 0.05:", sum(tf$qvalue < 0.05 & abs(tf$meth.diff) > 20), "; q < 0.1:", sum(tf$qvalue < 0.1 & abs(tf$meth.diff) > 20), "\n")
cat("  control-cohort corrected CpGs:\n"); print(cbind(pos[cfd], cf[cfd], in_draft_list = dml[cfd], pooled_corrected_q = pf_$qvalue[cfd]))
cat("  pooled corrected ploidy CpGs:\n"); print(cbind(pos[pfd], pf_[pfd], heat_cohort_diff = d_heat[pfd], control_cohort_diff = d_ctrl[pfd], heat_only_q = hf$qvalue[pfd], control_only_q = cf$qvalue[pfd]))
cat("  the 24 CpGs with corrected q < 0.01 for ploidy: |difference| quartiles", round(quantile(abs(pf_$meth.diff[pf_$qvalue < 0.01]), c(0, .25, .5, .75, 1)), 1), "\n")
cat("  the 287 CpGs with corrected q < 0.01 for heat: |difference| quartiles", round(quantile(abs(tf$meth.diff[tf$qvalue < 0.01]), c(0, .25, .5, .75, 1)), 1), "\n")
# corrected p at the draft's DML, within each cohort
cat("\nOf the draft's 1,038 ploidy DML: corrected p < 0.05 in the heat-only fit", sum(hf$pvalue[dml] < 0.05), "; in the control-only fit", sum(cf$pvalue[dml] < 0.05),
    "; in both", sum(hf$pvalue[dml] < 0.05 & cf$pvalue[dml] < 0.05), "; both and same sign", sum(hf$pvalue[dml] < 0.05 & cf$pvalue[dml] < 0.05 & sign(d_heat[dml]) == sign(d_ctrl[dml])), "\n")
# heat list: corrected p at the draft's 2,351
a0h <- cache("A0_heat"); dmlh <- a0h$qvalue < 0.01 & abs(a0h$meth.diff) > 20
cat("Of the draft's", sum(dmlh), "heat DML: corrected p < 0.05 at", sum(tf$pvalue[dmlh] < 0.05), sprintf("(%.1f%%)", 100 * mean(tf$pvalue[dmlh] < 0.05)), "; < 0.01 at", sum(tf$pvalue[dmlh] < 0.01),
    sprintf("(%.1f%%)", 100 * mean(tf$pvalue[dmlh] < 0.01)), "; < 0.001 at", sum(tf$pvalue[dmlh] < 0.001), sprintf("(%.1f%%)", 100 * mean(tf$pvalue[dmlh] < 0.001)),
    "; < 0.0001 at", sum(tf$pvalue[dmlh] < 1e-4), sprintf("(%.2f%%)", 100 * mean(tf$pvalue[dmlh] < 1e-4)), "\n")
cat("Of the draft's 1,038 ploidy DML (from the cache): corrected p < 0.05 at", sum(pf_$pvalue[dml] < 0.05), "; < 0.01 at", sum(pf_$pvalue[dml] < 0.01), "; < 0.001 at", sum(pf_$pvalue[dml] < 0.001),
    "; < 0.0001 at", sum(pf_$pvalue[dml] < 1e-4), "\n")
# regroupings of the control cohort against the full-data shares
w <- read.csv(file.path(S, "08_relabel_level_control_experiment.csv")); full <- c(0.0670586438607208, 0.0121105681124007, 0.00167990226023213)
cat("\ncontrol cohort, 150 regroupings (3,000 CpGs): with a corrected-F share at least the full-data share (6.7%, 1.2%, 0.17%):", sapply(1:3, function(k) sum(w[-1, paste0("mn_f", k)] >= full[k])), "\n")
w2 <- read.csv(file.path(S, "08_relabel_level_heat_experiment.csv")); fullh <- c(0.178313989004276, 0.0532223579718998, 0.00978924862553452)
cat("heat cohort, 125 regroupings: with a corrected-F share at least the full-data share (18%, 5.3%, 0.98%):", sapply(1:3, function(k) sum(w2[-1, paste0("mn_f", k)] >= fullh[k])),
    "; farthest regroupings (two exchanges, n =", sum(w2$moved[-1] == 4), "): mean share", round(100 * colMeans(w2[-1, ][w2$moved[-1] == 4, paste0("mn_f", 1:3)]), 3), "\n")
cat("control cohort regroupings: mean share", round(100 * colMeans(w[-1, paste0("mn_f", 1:3)]), 3), "\n")

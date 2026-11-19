# Combine the two relabelling runs (two independent samples of 3,000 eligible CpGs, 150 relabellings each).
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
a <- read.csv(file.path(S, "07_relabel_level_pooled_a.csv")); b <- read.csv(file.path(S, "07_relabel_level_pooled_b.csv"))
r <- rbind(cbind(run = "a", a[-1, ]), cbind(run = "b", b[-1, ])); o <- rbind(a[1, ], b[1, ])
pc <- function(x) round(100 * x, 2)
desc <- function(cols, rows = rep(TRUE, nrow(r))) {
  m <- r[rows, cols]
  cat(sprintf("  n = %d | mean %s | median %s | 5th %s | 95th %s | max %s\n", nrow(m), paste(pc(colMeans(m)), collapse = ", "), paste(pc(apply(m, 2, median)), collapse = ", "),
              paste(pc(apply(m, 2, quantile, .05)), collapse = ", "), paste(pc(apply(m, 2, quantile, .95)), collapse = ", "), paste(pc(apply(m, 2, max)), collapse = ", ")))
}
for (s in c("none", "mn_chisq", "mn_f", "int_none", "int_mn_f")) { cat(s, "(all 300 relabellings)\n"); desc(paste0(s, 1:3)) }
# pooled over CpG-instances: standard error of the mean share treating relabellings as the units
for (s in c("none", "mn_chisq", "mn_f", "int_mn_f")) {
  m <- r[, paste0(s, 1:3)]
  cat(sprintf("%s: mean share and its standard error over 300 relabellings: %s\n", s,
              paste(sprintf("%.3f%% (se %.3f)", 100 * colMeans(m), 100 * apply(m, 2, sd) / sqrt(nrow(m))), collapse = "; ")))
}
far <- r$moved_heat %in% c(4, 6); near <- !far
cat("\nMN F, relabellings with two exchanges in the heat cohort (farthest from its real grouping):\n"); desc(paste0("mn_f", 1:3), far)
cat("MN F, relabellings that keep the real heat-cohort grouping or change it by one exchange:\n"); desc(paste0("mn_f", 1:3), near)
cat("uncorrected, far:\n"); desc(paste0("none", 1:3), far); cat("uncorrected, near:\n"); desc(paste0("none", 1:3), near)
cat("MN chisq, far:\n"); desc(paste0("mn_chisq", 1:3), far)
cat("interaction MN F, far:\n"); desc(paste0("int_mn_f", 1:3), far); cat("interaction MN F, near:\n"); desc(paste0("int_mn_f", 1:3), near)
# against the shares of the full data (dml_counts.csv)
full <- list("ploidy, pooled (A1_ploidy_mn_f)" = c(0.0573915699450214, 0.0123549175320709, 0.00135919364691509),
             "heat (A1_heat_mn_f)" = c(0.110094685400122, 0.0355222968845449, 0.00673488087965791),
             "ploidy among controls (A4_ploidy_in_control_mn_f)" = c(0.0670586438607208, 0.0121105681124007, 0.00167990226023213))
for (n in names(full)) cat(sprintf("%s: relabellings (of 300, pooled design, MN F) with a share at least that of the full data: %s; among the %d far ones: %s\n", n,
                                   paste(sapply(1:3, function(k) sum(r[[paste0("mn_f", k)]] >= full[[n]][k])), collapse = ", "), sum(far),
                                   paste(sapply(1:3, function(k) sum(r[[paste0("mn_f", k)]][far] >= full[[n]][k])), collapse = ", ")))
cat("\nuncorrected test: relabellings (of 300) with a share at least that of the full data for ploidy (33.4%, 21.1%, 11.7%):",
    sum(r$none1 >= 0.333720219914478), sum(r$none2 >= 0.210949908368968), sum(r$none3 >= 0.117058643860721), "\n")
cat("CpGs above 20 points: observed", o$n_above20, "; relabelled median by run:", tapply(r$n_above20, r$run, median), "; share of relabelled instances with uncorrected p < 2.3e-4:",
    round(sum(r$n_above20_p_below_2.3e4) / sum(r$n_above20), 4), "\n")

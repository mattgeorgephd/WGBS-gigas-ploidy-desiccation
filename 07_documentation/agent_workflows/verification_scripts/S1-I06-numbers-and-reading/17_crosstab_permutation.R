suppressPackageStartupMessages(library(data.table)); options(width = 200)
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
for (mode in c("legacy", "corrected")) {
  r <- fread(file.path(S, sprintf("09_permutation_reference_%s.csv", mode))); obs <- r$n_above[1]; r <- r[labelling > 0]
  r[, heat_dist := pmin(moved_heat, 10 - moved_heat)]
  cat("\n", mode, ": observed", obs, "\n")
  print(r[, .(labellings = .N, median = as.numeric(median(n_above)), reach_observed = sum(n_above >= obs)), keyby = .(exchanges_in_heat_cohort = heat_dist / 2, high_split_imbalance = high_imbalance)])
  cat("reach observed: total", sum(r$n_above >= obs), "; with at most one exchange in the heat cohort", sum(r$n_above >= obs & r$heat_dist <= 2), "of", sum(r$heat_dist <= 2),
      "; with two exchanges", sum(r$n_above >= obs & r$heat_dist == 4), "of", sum(r$heat_dist == 4),
      "; with four or five high libraries in one group", sum(r$n_above >= obs & r$high_imbalance >= 3), "of", sum(r$high_imbalance >= 3),
      "; two exchanges AND 3:2 split", sum(r$n_above >= obs & r$heat_dist == 4 & r$high_imbalance == 1), "of", sum(r$heat_dist == 4 & r$high_imbalance == 1),
      "; median of the latter", median(r[heat_dist == 4 & high_imbalance == 1, n_above]), "\n")
}

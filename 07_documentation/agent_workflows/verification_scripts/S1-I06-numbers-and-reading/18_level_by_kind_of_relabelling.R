# The corrected F shares of the two relabelling runs, by kind of relabelling:
# exchanges in the heat cohort (0, 1, 2) and split of the five high-level libraries of I-49 (3:2, 4:1, 5:0).
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- readRDS(file.path(S, "sample3000_legacy.rds"))$lib
high <- lib$library_id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15"); trip <- as.integer(lib$ploidy == "triploid"); hk <- lib$cohort == "heat_experiment"
out <- list()
for (tag in c("a", "b")) {
  L <- readRDS(file.path(S, sprintf("07_relabel_labels_%s.rds", tag))); r <- read.csv(file.path(S, sprintf("07_relabel_level_pooled_%s.csv", tag)))
  mv <- colSums(L[hk, ] != trip[hk]); r$exch <- pmin(mv, 10 - mv) / 2; r$imb <- abs(colSums(L[high, ]) - colSums(1L - L[high, ])); r$run <- tag
  r$obs05 <- r$mn_f1[1]; r$obs01 <- r$mn_f2[1]; r$obs001 <- r$mn_f3[1]; r$obs_u05 <- r$none1[1]; r$obs_above <- r$n_above20[1]
  out[[tag]] <- r[-1, ]
}
r <- do.call(rbind, out); pc <- function(x) round(100 * x, 2)
cat("all relabellings:", nrow(r), "\n")
agg <- function(sel, name) cat(sprintf("%-46s n = %3d | corrected F: mean %s, median %s, 95th pct %s | at least the observed share of its run: %d, %d, %d | uncorrected mean %s | CpGs above 20 points: median %g, at least observed %d\n",
  name, sum(sel), paste(pc(colMeans(r[sel, c("mn_f1", "mn_f2", "mn_f3")])), collapse = ", "), paste(pc(apply(r[sel, c("mn_f1", "mn_f2", "mn_f3")], 2, median)), collapse = ", "),
  paste(pc(apply(r[sel, c("mn_f1", "mn_f2", "mn_f3")], 2, quantile, .95)), collapse = ", "),
  sum(r$mn_f1[sel] >= r$obs05[sel]), sum(r$mn_f2[sel] >= r$obs01[sel]), sum(r$mn_f3[sel] >= r$obs001[sel]), paste(pc(colMeans(r[sel, c("none1", "none2", "none3")])), collapse = ", "),
  median(r$n_above20[sel]), sum(r$n_above20[sel] >= r$obs_above[sel])))
agg(rep(TRUE, nrow(r)), "all")
agg(r$exch <= 1, "real heat-cohort grouping, or one exchange")
agg(r$exch == 2, "two exchanges in the heat cohort")
agg(r$imb == 1, "high-level libraries split 3:2")
agg(r$imb >= 3, "high-level libraries split 4:1 or 5:0")
agg(r$exch == 2 & r$imb == 1, "two exchanges and 3:2 split")
agg(r$exch == 2 & r$imb >= 3, "two exchanges and 4:1 or 5:0 split")
agg(r$exch <= 1 & r$imb == 1, "at most one exchange and 3:2 split")
cat("observed shares (run a, run b): 0.05:", pc(unique(r$obs05)), "; 0.01:", pc(unique(r$obs01)), "; 0.001:", pc(unique(r$obs001)), "\n")

# Control cohort: the 150 regroupings of 08_relabel_level_within_cohort.R rebuilt with the same seed, classified by
# whether the two high-level control libraries (zr3644_14, zr3644_15; both triploid) fall into the same group.
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- readRDS(file.path(S, "sample3000_legacy.rds"))$lib
keep <- which(lib$cohort == "control_experiment"); obs <- as.integer(lib$ploidy[keep] == "triploid")
all_sets <- combn(12, 6); all_sets <- all_sets[, all_sets[1, ] != 1, drop = FALSE]
labs <- apply(all_sets, 2, function(s) { z <- integer(12); z[s] <- 1L; z }); is_obs <- colSums(labs != obs) == 0
others <- which(!is_obs); set.seed(802); others <- sort(sample(others, 150)); labs <- cbind(labs[, is_obs], labs[, others])
r <- read.csv(file.path(S, "08_relabel_level_control_experiment.csv"))
moved <- colSums(labs != obs); moved <- pmin(moved, 12 - moved); stopifnot(all(moved == r$moved))     # the rebuilt labellings are those of the run
hi <- which(lib$library_id[keep] %in% c("zr3644_14", "zr3644_15")); together <- labs[hi[1], ] == labs[hi[2], ]
pc <- function(x) round(100 * x, 2); cols <- c("mn_f1", "mn_f2", "mn_f3")
for (k in c(TRUE, FALSE)) { sel <- together[-1] == k; m <- r[-1, ][sel, ]
  cat(sprintf("high-level control libraries %s: n = %d | corrected F mean %s, median %s, 95th pct %s | at least observed (7.2, 1.6, 0.167 on the sample): %d, %d, %d | CpGs above 20: median %g (observed %g), at least observed %d\n",
              if (k) "in the same group" else "in different groups", nrow(m), paste(pc(colMeans(m[, cols])), collapse = ", "), paste(pc(apply(m[, cols], 2, median)), collapse = ", "),
              paste(pc(apply(m[, cols], 2, quantile, .95)), collapse = ", "), sum(m$mn_f1 >= r$mn_f1[1]), sum(m$mn_f2 >= r$mn_f2[1]), sum(m$mn_f3 >= r$mn_f3[1]),
              median(m$n_above20), r$n_above20[1], sum(m$n_above20 >= r$n_above20[1]))) }
cat("observed grouping has them together:", together[1], "\n")
# heat cohort: by the split of the three high-level heat libraries (zr3534_3, 4, 5; all diploid)
keep <- which(lib$cohort == "heat_experiment"); obs <- as.integer(lib$ploidy[keep] == "triploid")
s <- combn(10, 5); s <- s[, s[1, ] != 1, drop = FALSE]; labs <- apply(s, 2, function(x) { z <- integer(10); z[x] <- 1L; z }); io <- colSums(labs != obs) == 0; labs <- cbind(labs[, io], labs[, !io])
r <- read.csv(file.path(S, "08_relabel_level_heat_experiment.csv")); moved <- colSums(labs != obs); moved <- pmin(moved, 10 - moved); stopifnot(all(moved == r$moved))
hi <- which(lib$library_id[keep] %in% c("zr3534_3", "zr3534_4", "zr3534_5")); k3 <- colSums(labs[hi, ]); split3 <- pmax(k3, 3 - k3)     # 3 = all three in one group, 2 = split 2:1
for (k in c(3, 2)) for (mv in c(2, 4)) { sel <- split3[-1] == k & r$moved[-1] == mv; if (!any(sel)) next; m <- r[-1, ][sel, ]
  cat(sprintf("heat cohort, three high libraries %s, %d exchange(s): n = %d | corrected F mean %s, max %s | CpGs above 20 (sample): median %g, max %g (observed %g)\n",
              if (k == 3) "together" else "split 2:1", mv / 2, nrow(m), paste(pc(colMeans(m[, cols])), collapse = ", "), paste(pc(apply(m[, cols], 2, max)), collapse = ", "),
              median(m$n_above20), max(m$n_above20), r$n_above20[1])) }

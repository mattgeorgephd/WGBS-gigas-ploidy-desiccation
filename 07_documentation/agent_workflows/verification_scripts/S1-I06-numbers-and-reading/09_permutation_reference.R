# Task (c)3: the relabelling reference of item 5, recomputed with the labellings of the audit script
# (same seed and sampling code), then looked at for what the audit text does not report.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "legacy"
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
u <- fread(file.path(R, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0
ns <- cs + as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
rm(u); invisible(gc(verbose = FALSE))
observed <- as.integer(lib$ploidy == "triploid")
blocks <- split(seq_along(observed), lib$cohort)          # control_experiment first, then heat_experiment, as in the audit script
cat("block order:", names(blocks), "\n")
set.seed(20261119)
labels <- cbind(observed, vapply(seq_len(999L), function(draw) { z <- observed; for (b in blocks) z[b] <- observed[b][sample.int(length(b))]; z }, integer(22)))
labels <- unname(labels)
rs_c <- rowSums(cs); rs_n <- rowSums(ns)
count_above <- function(L) {
  n <- integer(ncol(L))
  for (from in seq(1L, ncol(L), by = 40L)) {
    cols <- from:min(from + 39L, ncol(L))
    c1 <- cs %*% L[, cols, drop = FALSE]; n1 <- ns %*% L[, cols, drop = FALSE]
    eff <- 100 * (c1 / n1 - (rs_c - c1) / (rs_n - n1))
    n[cols] <- as.integer(colSums(abs(eff) > 20))
  }
  n
}
t0 <- Sys.time(); n_above <- count_above(labels); cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
q <- quantile(n_above, c(0, .05, .25, .5, .75, .95, 1), names = FALSE)
cat(sprintf("[%s] observed %d; min %g, q05 %g, q25 %g, median %g, q75 %g, q95 %g, max %g\n", mode, n_above[1], q[1], q[2], q[3], q[4], q[5], q[6], q[7]))
cat("labellings (of 1,000, the observed one included) with at least the observed count:", sum(n_above >= n_above[1]),
    "; random labellings (of 999) with at least:", sum(n_above[-1] >= n_above[1]), "; strictly more:", sum(n_above[-1] > n_above[1]),
    "; exactly equal:", sum(n_above[-1] == n_above[1]), "\n")
cat("median and quantiles over the 999 random labellings only: median", median(n_above[-1]), "; 5th", quantile(n_above[-1], .05), "; 95th", quantile(n_above[-1], .95), "\n")
# identity of the labellings
keyf <- function(z) paste(z, collapse = ""); keys <- apply(labels, 2, keyf); mirror <- apply(1L - labels, 2, keyf)
cat("distinct labellings among the 1,000:", length(unique(keys)), "; distinct up to the mirror image:", length(unique(pmin(keys, mirror))),
    "; random labellings equal to the observed one:", sum(keys[-1] == keys[1]), "; equal to its mirror image:", sum(keys[-1] == mirror[1]), "\n")
# labels moved (up to the mirror image), and the split of the five high-level libraries (I-49)
moved_heat <- colSums(labels[blocks$heat_experiment, ] != observed[blocks$heat_experiment]); moved_ctrl <- colSums(labels[blocks$control_experiment, ] != observed[blocks$control_experiment])
high <- lib$library_id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15")
imbalance <- abs(colSums(labels[high, ]) - colSums(1L - labels[high, ]))      # |high libraries in group 1 minus high libraries in group 0|: 1, 3 or 5
obs_n <- n_above[1]   # the observed count, under a name that no column carries (inside a grouped call n_above[1] would be the first count of the group)
res <- data.table(labelling = 0:999, n_above = n_above, moved_heat, moved_ctrl, high_imbalance = imbalance)
fwrite(res, file.path(S, sprintf("09_permutation_reference_%s.csv", mode)))
cat("\nBy the split of the five high-level libraries between the two groups (observed labelling: imbalance", imbalance[1], "):\n")
print(res[labelling > 0, .(labellings = .N, median = as.numeric(median(n_above)), q05 = quantile(n_above, .05), q95 = quantile(n_above, .95), max = as.numeric(max(n_above)),
                           at_least_observed = sum(n_above >= obs_n)), keyby = high_imbalance])
same <- res[labelling > 0 & high_imbalance == imbalance[1]]
cat(sprintf("Among the %d random labellings with the same imbalance as the observed one: median %g; at least observed %d; p = (%d + 1)/(%d + 1) = %.3f\n",
            nrow(same), median(same$n_above), sum(same$n_above >= n_above[1]), sum(same$n_above >= n_above[1]), nrow(same), (sum(same$n_above >= n_above[1]) + 1) / (nrow(same) + 1)))
cat("\nBy labels moved in the heat cohort (0 = heat cohort keeps its true labels or their mirror; 10 = all exchanged):\n")
res[, heat_dist := pmin(moved_heat, 10 - moved_heat)]; res[, ctrl_dist := pmin(moved_ctrl, 12 - moved_ctrl)]
print(res[labelling > 0, .(labellings = .N, median = as.numeric(median(n_above)), max = as.numeric(max(n_above)), at_least_observed = sum(n_above >= obs_n)), keyby = heat_dist])
cat("By labels moved in the control cohort:\n")
print(res[labelling > 0, .(labellings = .N, median = as.numeric(median(n_above)), max = as.numeric(max(n_above)), at_least_observed = sum(n_above >= obs_n)), keyby = ctrl_dist])
cat("Spearman correlation of the count with heat_dist:", round(cor(res$n_above[-1], res$heat_dist[-1], method = "spearman"), 3),
    "; with ctrl_dist:", round(cor(res$n_above[-1], res$ctrl_dist[-1], method = "spearman"), 3),
    "; with high_imbalance:", round(cor(res$n_above[-1], res$high_imbalance[-1], method = "spearman"), 3), "\n")
# Hybrid labellings: one cohort keeps its true labels, the other is relabelled (all 924 or all 252 regroupings).
all_labellings <- function(idx, k) { cmb <- combn(length(idx), k); apply(cmb, 2, function(s) { z <- integer(22); z[idx[s]] <- 1L; z }) }
Lc <- all_labellings(blocks$control_experiment, 6); Lc <- Lc + observed * (lib$cohort == "heat_experiment")      # heat true, control regrouped
Lh <- all_labellings(blocks$heat_experiment, 5);    Lh <- Lh + observed * (lib$cohort == "control_experiment")   # control true, heat regrouped
nc <- count_above(Lc); nh <- count_above(Lh)
cat(sprintf("\nHeat cohort keeps its true labels, control cohort regrouped in all %d ways: median %g, 5th to 95th %g to %g, min %g, max %g; regroupings with at least the observed %d: %d\n",
            ncol(Lc), median(nc), quantile(nc, .05), quantile(nc, .95), min(nc), max(nc), n_above[1], sum(nc >= n_above[1])))
cat(sprintf("Control cohort keeps its true labels, heat cohort regrouped in all %d ways: median %g, 5th to 95th %g to %g, min %g, max %g; regroupings with at least the observed %d: %d\n",
            ncol(Lh), median(nh), quantile(nh, .05), quantile(nh, .95), min(nh), max(nh), n_above[1], sum(nh >= n_above[1])))
fwrite(data.table(which = rep(c("control_regrouped", "heat_regrouped"), c(length(nc), length(nh))), n_above = c(nc, nh)), file.path(S, sprintf("09_hybrid_labellings_%s.csv", mode)))

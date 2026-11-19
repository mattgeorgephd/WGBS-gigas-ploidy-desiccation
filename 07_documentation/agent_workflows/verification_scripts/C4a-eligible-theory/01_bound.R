# Check (a): is the 10%-90% bound exact? How many DML fall outside, and how far?
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
sc <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory"
mode <- commandArgs(TRUE)[1]; if (is.na(mode)) mode <- "legacy"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))
stopifnot(all(lib$legacy_sample_index == 1:22))
trip <- lib$ploidy == "triploid"; heat <- lib$treatment == "heat"
cat("group sizes: ploidy", table(lib$ploidy), " heat", table(lib$treatment), "\n")
print(table(lib$ploidy, lib$treatment))
u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
cat("rows:", nrow(u), "\n")
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); Tt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE])
cov <- C + Tt
covcol <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE]); cat(sprintf("cells where coverage != numCs+numTs: %d of %d (max abs diff %d)\n", sum(cov != covcol), length(cov), max(abs(cov - covcol)))); rm(covcol)
pm <- 100 * C / cov
m <- rowMeans(pm)                                   # unweighted mean of per-library percentages (mean_percent_methylation)
mw <- 100 * rowSums(C) / rowSums(cov)               # coverage-weighted pooled over all libraries
cat(sprintf("mean<10: %d (%.2f%%); mean>90: %d (%.2f%%); eligible [10,90]: %d (%.2f%%)\n", sum(m < 10), 100*mean(m < 10), sum(m > 90), 100*mean(m > 90), sum(m >= 10 & m <= 90), 100*mean(m >= 10 & m <= 90)))
# group pooled (weighted) and unweighted means
wp <- function(g) 100 * rowSums(C[, g, drop = FALSE]) / rowSums(cov[, g, drop = FALSE])
up <- function(g) rowMeans(pm[, g, drop = FALSE])
d_ploidy_w <- wp(trip) - wp(!trip); d_heat_w <- wp(heat) - wp(!heat)
d_ploidy_u <- up(trip) - up(!trip); d_heat_u <- up(heat) - up(!heat)
dm <- list(ploidy = fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/diffmeth_ploidy_%s.csv.gz", mode))),
           heat   = fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/diffmeth_heat_%s.csv.gz", mode))))
cat("max |meth.diff - weighted pooled diff|: ploidy", max(abs(dm$ploidy$meth.diff - d_ploidy_w)), " heat", max(abs(dm$heat$meth.diff - d_heat_w)), "\n")
cat("max |meth.diff - unweighted diff|: ploidy", max(abs(dm$ploidy$meth.diff - d_ploidy_u)), " heat", max(abs(dm$heat$meth.diff - d_heat_u)), "\n")
cat("cor(weighted, unweighted diff): ploidy", cor(d_ploidy_w, d_ploidy_u), " heat", cor(d_heat_w, d_heat_u), "\n")
dml <- list(ploidy = fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/dml_ploidy_%s.csv", mode))),
            heat   = fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/dml_heat_%s.csv", mode))))
for (s in names(dml)) {
  r <- dml[[s]]$united_row
  stopifnot(all(u$chr[r] == dml[[s]]$chr), all(u$start[r] == dml[[s]]$start))
  mm <- m[r]
  cat(sprintf("\n[%s] %d DML; mean range %.3f to %.3f; outside [10,90]: %d (below %d, above %d)\n", s, length(r), min(mm), max(mm), sum(mm < 10 | mm > 90), sum(mm < 10), sum(mm > 90)))
  out <- which(mm < 10 | mm > 90)
  if (length(out)) print(data.table(set = s, chr = dml[[s]]$chr[out], start = dml[[s]]$start[out], meth.diff = dml[[s]]$meth.diff[out], unweighted_mean = mm[out], weighted_mean = mw[r][out],
        unweighted_diff = (if (s == "ploidy") d_ploidy_u else d_heat_u)[r][out], g1_w = (if (s == "ploidy") wp(trip) else wp(heat))[r][out], g0_w = (if (s == "ploidy") wp(!trip) else wp(!heat))[r][out]))
  # how close to the bound
  cat(sprintf("  DML with mean in [10,12): %d; [12,15): %d; [15,20): %d; [20,80]: %d; (80,85]: %d; (85,88]: %d; (88,90]: %d\n", sum(mm >= 10 & mm < 12), sum(mm >= 12 & mm < 15), sum(mm >= 15 & mm < 20), sum(mm >= 20 & mm <= 80), sum(mm > 80 & mm <= 85), sum(mm > 85 & mm <= 88), sum(mm > 88 & mm <= 90)))
  # unweighted difference of DML
  du <- (if (s == "ploidy") d_ploidy_u else d_heat_u)[r]
  cat(sprintf("  DML whose UNWEIGHTED group difference is <= 20: %d (min |unweighted diff| %.2f)\n", sum(abs(du) <= 20), min(abs(du))))
  dml[[s]][, `:=`(mean_pm = mm, d_unw = du)]
}
# all DML loci, distinct
allr <- unique(c(dml$ploidy$united_row, dml$heat$united_row))
cat(sprintf("\nDML rows summed over two lists: %d; distinct loci: %d; distinct loci outside [10,90]: %d\n", nrow(dml$ploidy) + nrow(dml$heat), length(allr), sum(m[allr] < 10 | m[allr] > 90)))
# Theoretical bounds
cat("\nTheoretical lower bound on the unweighted mean given |unweighted diff| > 20:\n")
cat("  ploidy 11 v 11: ", 20 * 11 / 22, "  heat 10 v 12: smaller group high:", 20 * 10 / 22, " larger group high:", 20 * 12 / 22, "\n")
# CpGs in the grey zone
cat(sprintf("tested CpGs with mean in [9.0909,10): %d; in (90, 90.909]: %d\n", sum(m >= 20*10/22 & m < 10), sum(m > 90 & m <= 100 - 20*10/22)))
# CpGs that reach |meth.diff|>20 (regardless of q) by eligibility
for (s in names(dm)) {
  big <- abs(dm[[s]]$meth.diff) > 20
  cat(sprintf("[%s] CpGs with |meth.diff|>20: %d; of them outside [10,90]: %d; min mean %.2f max mean %.2f\n", s, sum(big), sum(big & (m < 10 | m > 90)), min(m[big]), max(m[big])))
  print(data.table(chr = u$chr, start = u$start, meth.diff = dm[[s]]$meth.diff, q = dm[[s]]$qvalue, mean = m)[big & (m < 10 | m > 90)])
}
# Coverage by cohort at united CpGs
cat("\nMedian per-CpG coverage (after normalization) by library:\n"); print(round(apply(cov, 2, median)))
cat("mean coverage heat cohort / control cohort (per CpG ratio, quantiles):\n"); rr <- rowMeans(cov[, heat]) / rowMeans(cov[, !heat]); print(quantile(rr, c(.01,.05,.25,.5,.75,.95,.99)))
cat("max share of a group's coverage held by one library (ploidy groups), quantiles:\n")
mx <- pmax(apply(cov[, trip], 1, max) / rowSums(cov[, trip]), apply(cov[, !trip], 1, max) / rowSums(cov[, !trip])); print(quantile(mx, c(.5,.9,.99,.999,1)))
saveRDS(list(m = m, mw = mw, d_ploidy_w = d_ploidy_w, d_heat_w = d_heat_w, d_ploidy_u = d_ploidy_u, d_heat_u = d_heat_u,
             q_ploidy = dm$ploidy$qvalue, q_heat = dm$heat$qvalue, chr = u$chr, start = u$start,
             mean_cov = rowMeans(cov), sd_pm = apply(pm, 1, sd),
             m_heat = up(heat), m_ctrl = up(!heat), m_trip = up(trip), m_dip = up(!trip)), file.path(sc, sprintf("cpg_stats_%s.rds", mode)))

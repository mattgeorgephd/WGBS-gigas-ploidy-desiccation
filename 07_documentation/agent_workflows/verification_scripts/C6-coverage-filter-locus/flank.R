suppressPackageStartupMessages(library(data.table))
options(width = 250)
z <- readRDS("locus.rds"); reg <- z$reg; tested <- z$tested
pl <- fread("locus_per_library.csv")
cut <- fread("cutoffs_mine.csv")
reg <- merge(reg, cut[, .(library_id, cutoff = cutoff_mine, median_all)], by = "library_id")
reg <- merge(reg, pl[, .(library_id, ploidy, treatment, state = fifelse(n_over > 0, "high", "low"))], by = "library_id")
L0 <- 39633323L; L1 <- 39635520L
reg[, zone := fifelse(start >= L0 & start <= L1, "locus", fifelse(start >= L0 - 5000 & start < L0, "up5kb", fifelse(start > L1 & start <= L1 + 5000, "down5kb",
              fifelse(start < L0, "up_5to100kb", "down_5to100kb"))))]
reg[, ratio := cov / median_all]
cat("CpG sites with reads in at least one library: locus", reg[zone == "locus", uniqueN(start)], "; up 5 kb", reg[zone == "up5kb", uniqueN(start)],
    "; down 5 kb", reg[zone == "down5kb", uniqueN(start)], "\n")
s <- reg[, .(n_cpg = .N, median_cov = as.numeric(median(cov)), ratio = round(median(cov) / median_all[1], 2), n_over = sum(cov >= cutoff),
             meth = round(100 * sum(nC) / sum(cov), 1)), by = .(library_id, state, ploidy, treatment, zone)]
w <- dcast(s, library_id + state + ploidy + treatment ~ zone, value.var = c("n_cpg", "median_cov", "ratio", "n_over", "meth"))
w <- w[match(pl$library_id, library_id)]
cat("\nPer library: number of covered CpGs, median coverage, ratio to library median, CpGs at or above cutoff, pooled methylation\n")
print(w[, .(library_id, state, ploidy, treatment, n_up = n_cpg_up5kb, n_loc = n_cpg_locus, n_dn = n_cpg_down5kb,
            cov_up = median_cov_up5kb, cov_loc = median_cov_locus, cov_dn = median_cov_down5kb,
            r_up = ratio_up5kb, r_loc = ratio_locus, r_dn = ratio_down5kb, r_far_up = ratio_up_5to100kb, r_far_dn = ratio_down_5to100kb,
            over_up = n_over_up5kb, over_loc = n_over_locus, over_dn = n_over_down5kb,
            m_up = meth_up5kb, m_loc = meth_locus, m_dn = meth_down5kb)])
cat("\nBy state (high = at least one of the 54 tested CpGs at or above cutoff): median across libraries of the ratio to library median\n")
print(s[, .(n_lib = .N, ratio_median = median(ratio), ratio_min = min(ratio), ratio_max = max(ratio), n_over_total = sum(n_over)), by = .(state, zone)][order(state, zone)])
# all CpGs in the locus interval (not only the 54 tested), per library
cat("\nAll covered CpGs inside the locus interval (not only the 54 tested):\n")
print(reg[zone == "locus", .(n_cpg = .N, n_ge10 = sum(cov >= 10), median_cov = as.numeric(median(cov)), n_over = sum(cov >= cutoff)), by = .(library_id, state)][match(pl$library_id, library_id)])
# fine profile: 500 bp bins, +-12 kb
reg[, bin := (start - L0) %/% 500L]
prof <- reg[start >= L0 - 12000 & start <= L1 + 12000, .(r = median(ratio)), by = .(library_id, state, bin)][
  , .(n_lib = .N, ratio_median = round(median(r), 2), ratio_min = round(min(r), 2), ratio_max = round(max(r), 2)), by = .(state, bin)]
pw <- dcast(prof, bin ~ state, value.var = c("n_lib", "ratio_median", "ratio_min", "ratio_max"))
pw[, from := L0 + bin * 500L]
cat("\n500 bp bins (bin 0 starts at the first tested CpG; the locus is bins 0 to 4): median over libraries of the bin's median coverage / library median\n")
print(pw[, .(bin, from, n_high = n_lib_high, high_median = ratio_median_high, high_min = ratio_min_high, high_max = ratio_max_high,
             n_low = n_lib_low, low_median = ratio_median_low, low_min = ratio_min_low, low_max = ratio_max_low)], nrows = 200)
# coarser: 5 kb bins +-100 kb
reg[, bin5 := (start - L0) %/% 5000L]
p5 <- reg[, .(r = median(ratio), n = .N), by = .(library_id, state, bin5)][, .(ratio_median = round(median(r), 2), ratio_max = round(max(r), 2), cpg_median = median(n)), by = .(state, bin5)]
p5w <- dcast(p5, bin5 ~ state, value.var = c("ratio_median", "ratio_max", "cpg_median"))
cat("\n5 kb bins across +-100 kb:\n"); print(p5w, nrows = 100)
# boundaries of the high block: in high libraries, first and last CpG at or above cutoff in a contiguous sense
hb <- reg[state == "high" & cov >= cutoff & start >= L0 - 20000 & start <= L1 + 20000, .(first_over = min(start), last_over = max(start), n_over = .N), by = library_id]
cat("\nIn the 13 high libraries, CpGs at or above cutoff within +-20 kb:\n"); print(hb)
cat("positions at or above cutoff in at least 5 libraries, +-100 kb:\n")
po <- reg[cov >= cutoff, .N, by = start][N >= 5][order(start)]
cat("n =", nrow(po), "; range", range(po$start), "; outside the locus interval:", po[start < L0 | start > L1, .N], "\n")
print(po[start < L0 | start > L1])

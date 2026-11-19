suppressPackageStartupMessages(library(data.table))
z <- readRDS("locus.rds"); reg <- z$reg
pl <- fread("locus_per_library.csv"); cut <- fread("cutoffs_mine.csv")
reg <- merge(reg, cut[, .(library_id, cutoff = cutoff_mine, median_all)], by = "library_id")
reg <- merge(reg, pl[, .(library_id, state = fifelse(n_over > 0, "high", "low"))], by = "library_id")
B0 <- 39633065L; B1 <- 39636662L      # first and last CpG at or above cutoff in the high libraries
reg[, zone := fifelse(start >= B0 & start <= B1, "block", fifelse(start >= B0 - 5000 & start < B0, "up5kb_outside_block",
              fifelse(start > B1 & start <= B1 + 5000, "down5kb_outside_block", "far")))]
s <- reg[zone != "far", .(n_cpg = .N, median_cov = as.numeric(median(cov)), ratio = median(cov) / median_all[1], n_over = sum(cov >= cutoff), n_ge10 = sum(cov >= 10)), by = .(library_id, state, zone)]
print(s[, .(n_lib = .N, cpg_sites_median = median(n_cpg), ratio_median = round(median(ratio), 2), ratio_min = round(min(ratio), 2), ratio_max = round(max(ratio), 2),
            pairs_over = sum(n_over), pairs = sum(n_cpg)), by = .(zone, state)][order(zone, state)])
cat("block length:", B1 + 2 - B0, "bp; CpG sites with reads in the block:", reg[zone == "block", uniqueN(start)], "\n")
cat("all 22 libraries, block: ratio range", round(range(s[zone == "block", ratio]), 2), "\n")
cat("all 22 libraries, flanks outside block: ratio range up", round(range(s[zone == "up5kb_outside_block", ratio]), 2), " down", round(range(s[zone == "down5kb_outside_block", ratio]), 2), "\n")
# local contrast: block median / flank median, per library
w <- dcast(s, library_id + state ~ zone, value.var = "median_cov")
w[, contrast := block / ((up5kb_outside_block + down5kb_outside_block) / 2)]
print(w[order(state, contrast)])
cat("block-to-flank coverage contrast: low libraries", round(range(w[state == "low", contrast]), 1), "; high libraries", round(range(w[state == "high", contrast]), 1), "\n")

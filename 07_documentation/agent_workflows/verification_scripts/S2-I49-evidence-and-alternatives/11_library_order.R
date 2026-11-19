# 11_library_order.R
# Does the level follow the order of the library numbers?  I-49 notes that library numbers run diploids first,
# triploids second.  Here: are libraries with neighbouring numbers alike in level (both projects), and where do
# the libraries of unusual profile sit in the numbering?
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 200)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
sm <- fread(file.path(repo, "03_input_files/bismark_reports/061021-big_bismark_summary_report.txt"), sep = "\t")
setnames(sm, gsub(" ", "_", names(sm)))
sm[, lib := sub("_R1_bismark_bt2_pe.bam", "", File)]
sm[, level := 100 * Methylated_CpGs / (Methylated_CpGs + Unmethylated_CpGs)]
get <- function(prefix, n) { x <- sm[grepl(paste0("^", prefix, "_"), lib)]; x[, num := as.integer(sub(paste0(prefix, "_"), "", lib))]; setorder(x, num); stopifnot(identical(x$num, 1:n)); x }
z <- get("zr3644", 24); h <- get("zr3534", 10)
cat("zr3534, level by library number:\n"); print(round(setNames(h$level, h$num), 2))
cat("zr3644, level by library number:\n"); print(round(setNames(z$level, z$num), 2))

# statistic: mean absolute difference in level between libraries with neighbouring numbers (small = neighbours alike)
mad_adj <- function(x, pairs) mean(abs(x[pairs[, 1]] - x[pairs[, 2]]))
perm_test <- function(x, blocks, pairs, B = 100000, seed = 1) {
  set.seed(seed); obs <- mad_adj(x, pairs)
  ref <- replicate(B, { y <- x; for (b in blocks) y[b] <- sample(x[b]); mad_adj(y, pairs) })
  c(observed = obs, reference_median = median(ref), p_neighbours_more_alike = mean(ref <= obs + 1e-12))
}
adj <- function(idx) cbind(idx[-length(idx)], idx[-1])
cat("\n-- zr3644: neighbours within the four blocks of six (1-6, 7-12, 13-18, 19-24), levels permuted within block --\n")
blocks <- list(1:6, 7:12, 13:18, 19:24); pairs <- do.call(rbind, lapply(blocks, adj))
print(round(perm_test(z$level, blocks, pairs), 4))
cat("-- zr3644: the two blocks used in the study only (1-6, 13-18) --\n")
blocks2 <- list(1:6, 13:18); pairs2 <- do.call(rbind, lapply(blocks2, adj))
print(round(perm_test(z$level, blocks2, pairs2), 4))
cat("-- zr3534: neighbours within the two blocks of five (1-5 diploid, 6-10 triploid), permuted within block --\n")
blocksh <- list(1:5, 6:10); pairsh <- do.call(rbind, lapply(blocksh, adj))
# exact enumeration: 5! x 5! = 14,400 arrangements
library(utils)
perms <- function(v) { if (length(v) == 1) return(matrix(v, 1)); do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], perms(v[-i])))) }
p5 <- perms(1:5)
obsh <- mad_adj(h$level, pairsh)
refh <- as.vector(outer(seq_len(nrow(p5)), seq_len(nrow(p5)), Vectorize(function(a, b) { y <- h$level; y[1:5] <- h$level[p5[a, ]]; y[6:10] <- h$level[5 + p5[b, ]]; mad_adj(y, pairsh) })))
cat("observed", round(obsh, 4), "; median of the 14,400 arrangements", round(median(refh), 4), "; arrangements with neighbours at least as alike:", sum(refh <= obsh + 1e-12), "of", length(refh), "=", round(mean(refh <= obsh + 1e-12), 4), "\n")
cat("-- zr3534 diploid block alone (1-5): arrangements of the five values in which the three highest are consecutive --\n")
cat("  3 x 3! x 2! = 36 of 120 =", 36 / 120, "\n")
cat("-- chance that the k highest libraries of a block carry consecutive numbers --\n")
cat("  three highest of the ten zr3534 libraries consecutive: 8 of", choose(10, 3), "=", round(8 / choose(10, 3), 4), "\n")
cat("  two highest of the 24 zr3644 libraries adjacent: 23 of", choose(24, 2), "=", round(23 / choose(24, 2), 4), "; of the 12 used in the study (1-6, 13-18): 10 of", choose(12, 2), "=", round(10 / choose(12, 2), 4), "\n")
cat("  both:", round(8 / choose(10, 3) * 23 / choose(24, 2), 4), "(24 libraries) or", round(8 / choose(10, 3) * 10 / choose(12, 2), 4), "(12 libraries)\n")

cat("\n==== position of the libraries of unusual profile ====\n")
Rd <- readRDS(file.path(scratch, "03_correlations.rds"))$Rd
Mt <- as.matrix(fread(file.path(scratch, "04_thinned_correlation_22.csv"))[, -1])
id <- rownames(Rd)
high <- which(id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15")); nh <- setdiff(1:22, high)
typ <- sapply(nh, function(j) mean(Mt[j, setdiff(nh, j)])); names(typ) <- id[nh]
typd <- sapply(nh, function(j) mean(Rd[j, setdiff(nh, j)])); names(typd) <- id[nh]
print(data.table(library = id[nh], mean_r_with_other_16_equal_coverage = round(typ, 4), rank = rank(typ), noise_corrected = round(typd, 4), rank_nc = rank(typd))[order(rank)])
first2 <- c("zr3534_1", "zr3534_2", "zr3644_1", "zr3644_2")
cat("ranks (1 = least like the others) of the libraries numbered 1 and 2 of each project among the 17: equal coverage", rank(typ)[first2], "; noise-corrected", rank(typd)[first2], "\n")
cat("chance that four named libraries are the four lowest of 17: 1 of", choose(17, 4), "; that they are all among the five lowest: ", choose(13, 1), "of", choose(17, 5), "=", round(13 / choose(17, 5), 4), "\n")
# Wilcoxon: libraries numbered 1-2 against the other 13 not high
cat("exact Wilcoxon, typicality of the four against the other 13: p =", signif(wilcox.test(typ[first2], typ[setdiff(names(typ), first2)], exact = TRUE)$p.value, 3), "\n")

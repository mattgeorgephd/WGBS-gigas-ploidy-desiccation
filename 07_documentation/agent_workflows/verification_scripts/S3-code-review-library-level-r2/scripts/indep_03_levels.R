# indep_03_levels.R
# Independent recomputation of parts 1, 2 and 7 of library_methylation_level.R:
# level_group_tests.csv (all measures; the exact Wilcoxon p-value also by my own enumeration of rank sums),
# level_correlations.csv (with the exact null distribution of Spearman's statistic computed by dynamic programming for n = 6, 10, 12),
# pair_similarity.csv (subsets from bit patterns, mean correlation from row and column sums of the full matrix).
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
out <- file.path(S, "out"); committed <- file.path(repo, "06_diagnostics/library_methylation_level")
lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE)
bis <- read.csv(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"), stringsAsFactors = FALSE)
stopifnot(identical(bis$library_id, lib$library_id), identical(lib$legacy_sample_index, 1:22))
meas <- list(pct_mCpG_genome = 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG))
chh <- 100 * bis$meth_CHH / (bis$meth_CHH + bis$unmeth_CHH); chg <- 100 * bis$meth_CHG / (bis$meth_CHG + bis$unmeth_CHG)
for (mode in c("legacy", "corrected")) {
  u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
  pm <- lapply(1:22, function(j) 100 * u[[paste0("numCs", j)]] / (u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]]))
  mean_pm <- Reduce(`+`, pm) / 22; el <- mean_pm >= 10 & mean_pm <= 90
  meas[[paste0("mean_pct_tested_", mode)]] <- sapply(pm, mean)
  meas[[paste0("mean_pct_eligible_", mode)]] <- sapply(pm, function(x) mean(x[el]))
  rm(u, pm); invisible(gc())
}
meas[["pct_mCHH_genome"]] <- chh
lev <- fread(file.path(committed, "library_level.csv"))
cat("library_level.csv against my values (6 significant digits):\n")
for (nm in names(meas)) cat(sprintf("  %-30s %s\n", nm, isTRUE(all.equal(signif(meas[[nm]], 6), lev[[nm]], tolerance = 1e-12))))
cat(sprintf("  %-30s %s\n", "pct_mCHG_genome", isTRUE(all.equal(signif(chg, 6), lev$pct_mCHG_genome, tolerance = 1e-12))))
el_file <- read.csv(file.path(repo, "04_output/03_global_methylation/tables/elisa_percent_5mC.csv"), stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(el_file$animal_id))
elisa <- el_file$percent_5mC[match(lib$animal_id, el_file$animal_id)]; elisa[lib$animal_id == "" | is.na(lib$animal_id)] <- NA
cat(sprintf("  %-30s %s ; animals with both: %s ; of these above the top standard of the ELISA curve: %s\n", "elisa_percent_5mC",
            isTRUE(all.equal(signif(elisa, 6), lev$elisa_percent_5mC, tolerance = 1e-12)), paste(lib$animal_id[!is.na(elisa)], collapse = " "),
            paste(el_file$animal_id[el_file$animal_id %in% lib$animal_id & el_file$above_top_standard], collapse = " ")))
cat("  level_rank equals rank of -level:", all(lev$level_rank == rank(-meas$pct_mCpG_genome)), "\n")
for (mode in c("legacy", "corrected")) {
  pcs <- read.csv(file.path(repo, sprintf("04_output/04_cpg_coverage/tables/pca_scores_10x_%s.csv", mode)), stringsAsFactors = FALSE)
  stopifnot(identical(pcs$library_id, lib$library_id))
  for (k in 1:4) { meas[[sprintf("pc%d_10x_%s", k, mode)]] <- pcs[[paste0("PC", k)]]
    stopifnot(isTRUE(all.equal(signif(pcs[[paste0("PC", k)]], 6), lev[[sprintf("pc%d_10x_%s", k, mode)]], tolerance = 1e-12))) }
}
cat("  principal components 1 to 4 of both modes equal the step 04 tables at 6 significant digits: TRUE\n")

# ---- part 2: group tests ----------------------------------------------------------------------------
w_cache <- list()
wilcox_exact_own <- function(x, y) {          # two-sided exact p from every assignment of the ranks; no ties allowed
  r <- rank(c(x, y)); stopifnot(!anyDuplicated(r)); n <- length(r); k <- length(x)
  key <- paste(n, k); if (is.null(w_cache[[key]])) w_cache[[key]] <<- colSums(combn(n, k))
  w <- sum(r[seq_len(k)]); all_w <- w_cache[[key]]; e <- k * (n + 1) / 2
  c(p = mean(abs(all_w - e) >= abs(w - e) - 1e-9), assignments = length(all_w), as_extreme = sum(abs(all_w - e) >= abs(w - e) - 1e-9))
}
sets <- list(
  "ploidy inside the heat cohort" = list(a = lib$ploidy == "diploid" & lib$treatment == "heat", b = lib$ploidy == "triploid" & lib$treatment == "heat"),
  "ploidy inside the control cohort" = list(a = lib$ploidy == "diploid" & lib$treatment == "control", b = lib$ploidy == "triploid" & lib$treatment == "control"),
  "cohort" = list(a = lib$treatment == "control", b = lib$treatment == "heat"),
  "ploidy, all 22 libraries" = list(a = lib$ploidy == "diploid", b = lib$ploidy == "triploid"),
  "acute heat shock inside the heat cohort" = list(a = lib$treatment == "heat" & lib$acute_heat_shock != "yes", b = lib$acute_heat_shock == "yes"))
tests_committed <- fread(file.path(committed, "level_group_tests.csv"))
mine <- rbindlist(lapply(c("pct_mCpG_genome", "mean_pct_tested_legacy", "mean_pct_eligible_legacy", "mean_pct_tested_corrected", "mean_pct_eligible_corrected", "pct_mCHH_genome"), function(ms)
  rbindlist(lapply(names(sets), function(nm) {
    x <- meas[[ms]][sets[[nm]]$a]; y <- meas[[ms]][sets[[nm]]$b]; own <- wilcox_exact_own(x, y)
    data.table(measure = ms, comparison = nm, n_a = length(x), n_b = length(y), min_a = min(x), median_a = median(x), max_a = max(x),
               min_b = min(y), median_b = median(y), max_b = max(y), difference_of_means_b_minus_a = mean(y) - mean(x),
               a_libraries_above_every_b_library = sum(x > max(y)), b_libraries_above_every_a_library = sum(y > max(x)),
               complete_separation = max(x) < min(y) | max(y) < min(x),
               wilcoxon_p = wilcox.test(x, y, exact = TRUE)$p.value, wilcoxon_p_own_enumeration = own[["p"]], assignments = own[["assignments"]], as_extreme = own[["as_extreme"]],
               ties = anyDuplicated(c(x, y)) > 0)
  }))))
stopifnot(identical(mine$measure, tests_committed$measure), identical(mine$comparison, tests_committed$comparison))
cat("\nlevel_group_tests.csv (30 rows):\n")
for (col in c("n_a", "n_b", "min_a", "median_a", "max_a", "min_b", "median_b", "max_b", "difference_of_means_b_minus_a", "a_libraries_above_every_b_library",
              "b_libraries_above_every_a_library", "wilcoxon_p"))
  cat(sprintf("  %-36s %s\n", col, isTRUE(all.equal(signif(as.numeric(mine[[col]]), 6), as.numeric(tests_committed[[col]]), tolerance = 1e-12))))
cat(sprintf("  %-36s %s\n", "complete_separation", identical(mine$complete_separation, tests_committed$complete_separation)))
cat(sprintf("  %-36s %s\n", "wilcoxon_p_is_exact (no ties)", identical(!mine$ties, tests_committed$wilcoxon_p_is_exact)))
cat(sprintf("  %-36s %s (largest difference %.2e)\n", "R's exact p = my enumeration", isTRUE(all.equal(mine$wilcoxon_p, mine$wilcoxon_p_own_enumeration, tolerance = 1e-10)),
            max(abs(mine$wilcoxon_p - mine$wilcoxon_p_own_enumeration))))
options(width = 220)
print(mine[measure == "pct_mCpG_genome", .(comparison, n_a, n_b, median_a = signif(median_a, 6), median_b = signif(median_b, 6), diff = signif(difference_of_means_b_minus_a, 6),
                                             a_above = a_libraries_above_every_b_library, b_above = b_libraries_above_every_a_library, complete_separation,
                                             wilcoxon_p = signif(wilcoxon_p, 6), own = signif(wilcoxon_p_own_enumeration, 6), as_extreme, assignments)])
fwrite(mine, file.path(out, "level_group_tests_recomputed.csv"))

# ---- part 2: correlations ------------------------------------------------------------------------------
# Exact null distribution of sum(d^2) for n untied pairs: dynamic programming over the set of ranks already used.
spearman_exact_dist <- function(n) {
  maxd <- n * (n^2 - 1) / 3
  cur <- list(); cur[["0"]] <- c(1, rep(0, maxd))           # key: bit mask of used ranks; value: counts by partial sum of d^2
  for (pos in 1:n) {
    nxt <- new.env()
    for (key in names(cur)) {
      mask <- as.integer(key); v <- cur[[key]]
      for (r in 1:n) if (bitwAnd(mask, bitwShiftL(1L, r - 1L)) == 0L) {
        nk <- as.character(bitwOr(mask, bitwShiftL(1L, r - 1L))); sh <- (pos - r)^2
        add <- c(rep(0, sh), v[seq_len(maxd + 1 - sh)])
        if (is.null(nxt[[nk]])) nxt[[nk]] <- add else nxt[[nk]] <- nxt[[nk]] + add
      }
    }
    cur <- as.list(nxt)
  }
  stopifnot(length(cur) == 1); cnt <- cur[[1]]; stopifnot(abs(sum(cnt) - factorial(n)) < 0.5)
  cnt                                                         # cnt[q + 1] = number of permutations with sum(d^2) = q
}
dist_cache <- list()
spearman_p_exact_own <- function(x, y) {
  n <- length(x); stopifnot(!anyDuplicated(x), !anyDuplicated(y), n <= 12)
  if (is.null(dist_cache[[as.character(n)]])) dist_cache[[as.character(n)]] <<- spearman_exact_dist(n)
  cnt <- dist_cache[[as.character(n)]]; q <- round(sum((rank(x) - rank(y))^2)); maxd <- n * (n^2 - 1) / 3
  lower <- sum(cnt[1:(q + 1)]) / factorial(n); upper <- sum(cnt[(q + 1):(maxd + 1)]) / factorial(n)
  min(1, 2 * min(lower, upper))                                # two-sided, as cor.test defines it
}
L <- meas$pct_mCpG_genome; heat <- lib$treatment == "heat"; ctrl <- lib$treatment == "control"
pairs <- c(lapply(c("legacy", "corrected"), function(mode) c(sprintf("pc%d_10x_%s", 1:4, mode), paste0("mean_pct_tested_", mode), paste0("mean_pct_eligible_", mode))))
rows <- list()
addrow <- function(yname, x, y, set) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rows[[length(rows) + 1]] <<- data.table(y = yname, libraries = set, n = n, pearson = cor(x, y), spearman = cor(rank(x), rank(y)),
    p_from_cor.test_exact_TRUE = suppressWarnings(cor.test(x, y, method = "spearman", exact = TRUE)$p.value),
    p_exact_own = if (n <= 12) spearman_p_exact_own(x, y) else NA_real_,
    p_t_approximation = suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value), sum_d2 = round(sum((rank(x) - rank(y))^2)))
}
for (yn in unlist(pairs)) addrow(yn, L, meas[[yn]], "all 22")
addrow("pct_mCHH_genome", L, chh, "all 22"); addrow("pct_mCHH_genome", L[heat], chh[heat], "heat cohort"); addrow("pct_mCHH_genome", L[ctrl], chh[ctrl], "control cohort")
addrow("elisa_percent_5mC", L, elisa, "animals with both a WGBS library and an ELISA value")
cors <- rbindlist(rows); cor_committed <- fread(file.path(committed, "level_correlations.csv"))
stopifnot(identical(cors$y, cor_committed$y), identical(cors$libraries, cor_committed$libraries))
cat("\nlevel_correlations.csv (16 rows):\n")
for (col in c("n", "pearson", "spearman")) cat(sprintf("  %-28s %s\n", col, isTRUE(all.equal(signif(as.numeric(cors[[col]]), 6), as.numeric(cor_committed[[col]]), tolerance = 1e-12))))
cat(sprintf("  %-28s %s\n", "spearman_p_exact = cor.test", isTRUE(all.equal(signif(cors$p_from_cor.test_exact_TRUE, 6), cor_committed$spearman_p_exact, tolerance = 1e-12))))
cors[, relative_error_of_committed_p := p_from_cor.test_exact_TRUE / p_exact_own - 1]
print(cors[, .(y, libraries = substr(libraries, 1, 14), n, spearman = signif(spearman, 6), sum_d2, p_committed = signif(p_from_cor.test_exact_TRUE, 6), p_exact_own = signif(p_exact_own, 6),
               relative_error = signif(relative_error_of_committed_p, 3), p_t_approximation = signif(p_t_approximation, 6))])
fwrite(cors, file.path(out, "level_correlations_recomputed.csv"))

# ---- part 7: similarity of the libraries of one group --------------------------------------------------
sim_committed <- fread(file.path(committed, "pair_similarity.csv"))
popcount <- function(x) { n <- 0L; while (any(x > 0L)) { n <- n + bitwAnd(x, 1L); x <- bitwShiftR(x, 1L) }; n }
sim <- list()
for (mode in c("legacy", "corrected")) {
  cm <- read.csv(file.path(repo, sprintf("04_output/04_cpg_coverage/tables/library_correlation_10x_%s.csv", mode)), check.names = FALSE, stringsAsFactors = FALSE)
  R <- as.matrix(cm[, -1]); stopifnot(identical(cm$library_id, lib$library_id), identical(colnames(R), lib$library_id), isTRUE(all.equal(R, t(R), check.attributes = FALSE, tolerance = 0)),
                                      all(diag(R) == 1))
  for (coh in c("heat", "control")) {
    members <- which(lib$treatment == coh); n <- length(members)
    for (pl in c("diploid", "triploid")) {
      own <- which(lib$treatment == coh & lib$ploidy == pl); k <- length(own); other <- setdiff(members, own)
      masks <- 0:(2^n - 1); masks <- masks[popcount(masks) == k]
      ind <- sapply(masks, function(mk) as.numeric(bitwAnd(mk, bitwShiftL(1L, 0:(n - 1))) != 0L))     # members x subsets, 0/1
      Rm <- R[members, members]
      ref <- (colSums(ind * (Rm %*% ind)) - k) / (k * (k - 1))                                         # (sum of the block - diagonal) / ordered pairs
      own_col <- which(colSums(ind[match(own, members), , drop = FALSE]) == k); stopifnot(length(own_col) == 1)
      obs <- ref[own_col]
      pr <- which(upper.tri(Rm), arr.ind = TRUE); pr_r <- Rm[pr]; pr_own <- members[pr[, 1]] %in% own & members[pr[, 2]] %in% own
      best <- which.max(replace(ref, own_col, -Inf)); best_members <- members[ind[, best] == 1]
      sim[[length(sim) + 1]] <- data.table(mode = mode, cohort = paste(coh, "cohort"), group = pl, n_libraries = k, n_pairs = sum(pr_own), mean_r_within_group = obs,
        min_r_within_group = min(pr_r[pr_own]), max_r_within_group = max(pr_r[pr_own]),
        mean_r_within_other_group = (sum(R[other, other]) - length(other)) / (length(other) * (length(other) - 1)),
        mean_r_between_groups = mean(R[own, other]), largest_mean_r_of_another_subset = max(ref[-own_col]),
        pairs_of_group_among_the_n_pairs_most_similar_in_cohort = sum(pr_own[order(pr_r, decreasing = TRUE)[seq_len(sum(pr_own))]]),
        subsets_of_this_size = length(ref), subsets_at_or_above_observed = sum(ref >= obs), p_enumeration = sum(ref >= obs) / length(ref),
        # extra: is the count sensitive to the tolerance, and which subset is the largest other one
        other_subsets_within_1e_9_of_observed = sum(abs(ref[-own_col] - obs) < 1e-9), subsets_at_or_below_observed = sum(ref <= obs),
        largest_other_subset = paste(lib$library_id[best_members], collapse = " "), own_libraries_in_it = sum(best_members %in% own),
        smallest_gap_between_two_pair_correlations = min(diff(sort(pr_r))))
    }
  }
}
sim <- rbindlist(sim)
stopifnot(identical(sim$mode, sim_committed$mode), identical(sim$cohort, sim_committed$cohort), identical(sim$group, sim_committed$group))
cat("\npair_similarity.csv (8 rows):\n")
for (col in names(sim_committed)[4:15]) cat(sprintf("  %-58s %s\n", col, isTRUE(all.equal(signif(as.numeric(sim[[col]]), 6), as.numeric(sim_committed[[col]]), tolerance = 1e-12))))
print(sim[, .(mode, cohort, group, mean_r_within_group = signif(mean_r_within_group, 6), largest_other = signif(largest_mean_r_of_another_subset, 6), largest_other_subset, own_libraries_in_it,
              at_or_above = subsets_at_or_above_observed, at_or_below = subsets_at_or_below_observed, subsets = subsets_of_this_size, p = signif(p_enumeration, 6),
              near_ties = other_subsets_within_1e_9_of_observed, min_gap = signif(smallest_gap_between_two_pair_correlations, 2))])
fwrite(sim, file.path(out, "pair_similarity_recomputed.csv"))

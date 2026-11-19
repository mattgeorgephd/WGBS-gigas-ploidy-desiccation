# feature_odds_matched.R
#
# Purpose   Are DML more likely than comparable CpGs to lie in gene bodies,
#           coding sequence, introns, transposable elements, or intergenic
#           sequence? "Comparable" is the point. Step 07 compares DML with all
#           tested CpGs and with "eligible" CpGs (mean methylation across the
#           22 libraries between 10% and 90%; dml_eligible()). Inside the
#           eligible range the chance that a CpG exceeds the 20-point
#           difference without any treatment effect still depends on its mean
#           methylation and on how much animals of one group differ, and both
#           differ between genome features. This script makes the comparison
#           fairer step by step and shows what is left (issue I-41).
#           For each mode, for the ploidy and the heat DML list, and for the
#           tracks gene body, CDS, intron, exon UTR, lncRNA, transposable
#           element and intergenic (tracks, coordinates and overlap rule of
#           step 07), it computes the odds ratio of a DML lying in the track:
#             1. against all tested CpGs and against eligible CpGs, with Fisher
#                95% intervals. These are the values of step 07; the script
#                stops if they differ from feature_association_tests_<mode>.csv.
#             2. against eligible CpGs in the same 5-point bin of mean
#                methylation (Mantel-Haenszel odds ratio).
#             3. against eligible CpGs in the same methylation bin and the same
#                decile of within-group standard deviation (SD).
#             4. against CpGs with a mean of 5% to 95%, and of 15% to 85%, to
#                show how much comparison 1 depends on the range.
#             5. transposable elements only: against eligible CpGs in the top
#                quartile of mean coverage, and in the lower three quartiles.
#             6. within gene bodies: CDS against intron.
#             7. comparisons 1 (eligible), 2 and 3 again with an interval from a
#                bootstrap over 100 kb blocks, because neighbouring CpGs are
#                not independent observations.
#             8. ploidy only: the odds ratios of comparison 1 (eligible) and of
#                comparison 3 for CpG lists made with shuffled ploidy labels.
#                This reference needs no model: it is what a list of this kind
#                shows when the labels carry no information about ploidy.
#           It also writes the within-group SD of DML and of the other eligible
#           CpGs (medians, overall and inside each SD decile).
#
# Run       Rscript 06_diagnostics/feature_odds_matched.R     (about 4 minutes; about 1.5 GB of memory)
#
# Reads     04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           04_output/05_differential_methylation/tables/dml_{ploidy,heat}_<mode>.csv
#           04_output/05_differential_methylation/tables/diffmeth_{ploidy,heat}_<mode>.csv.gz
#           04_output/07_genome_features/tables/feature_association_tests_<mode>.csv
#           03_input_files/genome_features/      seven feature tracks
#           03_input_files/metadata/wgbs_libraries.csv
# Writes    06_diagnostics/feature_odds_matched.csv
#               one row per mode, DML list, track, comparison and kind of interval
#           06_diagnostics/feature_odds_matched_within_group_sd.csv
#               within-group SD of DML and of other eligible CpGs, overall and by decile
#           06_diagnostics/feature_odds_label_permutation.csv
#               per track: the odds ratios of the real list (comparisons 1 and 3) and their
#               distribution over the label permutations; share of eligible CpGs called
#               inside and outside the track
#           06_diagnostics/feature_odds_label_permutation_counts.csv
#               difference rule and q-value rule in the real data; CpGs called per
#               permutation; share of eligible CpGs called, by bin of mean methylation
#
# Notes     - Convention of step 07, kept throughout: every DML of a list is
#             compared with the CpGs of the background that are not in that
#             list. The few DML outside the background stay in the DML row.
#           - Within-group SD of a CpG: the pooled standard deviation of the
#             per-library percent methylation around the means of the four
#             groups (diploid control, diploid heat, triploid control, triploid
#             heat), 18 degrees of freedom. Differences between the groups do
#             not enter it. Deciles are those of the eligible CpGs. Mean and SD
#             come from the same 22 libraries as the contrasts.
#           - Methylation bins are closed on the left: 10 to 15 means at least
#             10 and below 15. A stratum is used when it holds a DML and
#             another CpG, and CpGs inside and outside the track; the other
#             strata add nothing to a Mantel-Haenszel odds ratio. The output
#             gives the number of strata used and of DML in them.
#           - Unstratified odds ratios are the conditional estimates of
#             fisher.test(), as in step 07; stratified ones are the common odds
#             ratio of mantelhaen.test() with its interval. Both intervals
#             treat CpGs as independent. Deciles are a coarse match; the SD
#             table shows what difference in SD is left inside each decile, and
#             n_sd_groups below sets a finer grouping.
#           - Block interval (comparison 7): the sequences are cut into 100 kb
#             windows; each of 2,000 resamples draws, with replacement, as many
#             windows as hold a compared CpG; the interval runs from the 2.5th
#             to the 97.5th percentile of the odds ratio recomputed on each
#             resample. A bootstrap is used and not a delete-one-block
#             jackknife, because a jackknife interval rests on a normal
#             distribution of the log odds ratio, and that fails where one
#             window holds most of the DML of a track (in legacy mode 53 of
#             the 65 ploidy DML in lncRNA lie in the locus of issue I-09).
#           - Label permutation (comparison 8): the ploidy labels are shuffled
#             within cohort (5 and 5 in the heat cohort, 6 and 6 in the control
#             cohort), the read-pooled difference triploid minus diploid is
#             recomputed at every tested CpG from the normalized counts
#             (methylKit's meth.diff: pooled methylated reads over pooled reads
#             of each group), and the CpGs whose absolute difference exceeds
#             config$dml_difference are called. 1,000 permutations. A second
#             scheme, reported beside the first, also keeps the acute heat
#             shock balanced (3 and 3, and 2 and 2, in the heat cohort; issue
#             I-04), as the real labels do. With 1,000 permutations a share near
#             0.13 has a standard error of 0.01.
#           - Comparison 3 under shuffled labels: for every permutation the
#             within-group SD and its deciles are computed again, around the
#             means of the four groups that the shuffled labels define (label
#             by cohort). A CpG called under shuffled labels differs between
#             the shuffled groups, and that difference has to stay out of its
#             SD exactly as the real group difference stays out of the SD of a
#             real DML; with the SD of the real groups the shuffled lists would
#             be matched to CpGs that are too variable.
#           - THE PERMUTATION APPLIES THE DIFFERENCE RULE ONLY, NOT THE Q-VALUE
#             RULE: refitting the logistic regression at every CpG for every
#             permutation is not feasible here. Evidence that the q-value rule
#             does not bind in the real data, counted from
#             diffmeth_<contrast>_<mode>.csv.gz and written to the counts
#             table: every CpG above the difference threshold for ploidy also
#             has q below config$dml_qvalue (1,038 of 1,038 in legacy mode, 952
#             of 952 in corrected mode; the largest q among them is about half
#             the cutoff), and the script stops if that is not so. For heat one
#             CpG of 2,352 fails the q-value rule in legacy mode and none in
#             corrected mode. What the q-values would be under shuffled labels
#             is not computed, so a count per permutation is an upper bound on
#             what the full rule would call.
#           - For heat no label permutation exists: heat coincides with cohort
#             (issue I-01), so no shuffle of the heat labels keeps the cohorts
#             apart. The output says so; nothing is put in its place.
#           - A difference of exactly the threshold is not called. The rule is
#             evaluated on whole numbers (see difference_above_threshold()), so
#             that rounding cannot call such a tie.
#           - A CpG whose per-library percentages span no more than the
#             threshold cannot be called under any labelling, because a pooled
#             percentage is a weighted mean of per-library percentages. Such
#             CpGs are skipped in the permutation loop; the result is the same
#             as evaluating every tested CpG.
#           - The counts table also gives the number of CpGs called when the
#             ploidy labels are exchanged in one cohort only. That labelling is
#             the ploidy by cohort contrast and is one of the possible
#             permutations. It calls more than three times as many CpGs as the
#             real labels, about as many as the largest counts among the
#             permutations (see also issues I-02 and I-07).
#           - Random numbers: set.seed(20261119) before the label permutations
#             (1,000 per scheme, the same labellings for both modes) and
#             set.seed(20261120) before the block resamples of each mode
#             (2,000). Two runs give identical files.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges"))
setDTthreads(config$cores)
t_start <- Sys.time()

# ---- Settings of this diagnostic -------------------------------------------------------
tracks                <- c("gene", "CDS", "intron", "exon_UTR", "lncRNA", "TE", "intergenic")  # names of step 07
contrasts             <- c("ploidy", "heat")
methylation_bin_width <- 5          # percentage points
n_sd_groups           <- 10         # deciles of within-group SD among eligible CpGs
range_lower_bounds    <- c(5, 15)   # comparison 4: mean methylation from the bound to 100 minus the bound
coverage_quantile     <- 0.75       # comparison 5: CpGs above this quantile of mean coverage are the top quartile
coverage_track        <- "TE"
block_size            <- 1e5        # base pairs
n_boot                <- 2000       # block resamples
n_perm                <- 1000       # label permutations per scheme
perm_chunk            <- 25         # permutations evaluated at a time (memory only)
seed_permutation      <- 20261119
seed_bootstrap        <- 20261120
out_file <- function(f) here("06_diagnostics", f)

lib <- read_libraries()
n_lib <- nrow(lib)
is_triploid <- lib$ploidy == "triploid"
is_heat <- lib$treatment == "heat"
real_labels <- list(ploidy = as.numeric(is_triploid), heat = as.numeric(is_heat))
in_first_cohort <- as.numeric(lib$cohort == lib$cohort[1])
# Heat coincides with cohort: each cohort holds one treatment. This is why the
# heat labels cannot be permuted (issue I-01).
stopifnot(length(unique(lib$cohort)) == 2,
          all(tapply(as.character(lib$treatment), lib$cohort, function(x) length(unique(x))) == 1))
no_heat_permutation <- paste("No label permutation exists for heat: heat coincides with cohort (issue I-01),",
                             "so the heat labels cannot be shuffled within cohort.")

# ---- Label permutations (drawn once, used for both modes) -----------------------------
perm_strata <- list(`within cohort` = lib$cohort,
                    `within cohort and acute heat shock` = paste(lib$cohort, lib$acute_heat_shock))
stopifnot(all(vapply(perm_strata, function(s) all(tapply(is_triploid, s, function(x) any(x) && !all(x))), logical(1))))
draw_labels <- function(strata) {
  labels <- numeric(n_lib)
  for (s in unique(strata)) {
    members <- which(strata == s)
    labels[members[sample.int(length(members), sum(is_triploid[members]))]] <- 1
  }
  labels
}
set.seed(seed_permutation)
label_draws <- lapply(perm_strata, function(s) vapply(seq_len(n_perm), function(i) draw_labels(s), numeric(n_lib)))
for (scheme in names(perm_strata)) {      # every draw keeps the number of triploid labels in every stratum
  stopifnot(all(rowsum(label_draws[[scheme]], perm_strata[[scheme]]) ==
                as.vector(rowsum(real_labels$ploidy, perm_strata[[scheme]]))))
}

# ---- Functions ---------------------------------------------------------------------------
# The four groups that a labelling defines: label by cohort. Libraries by
# groups, entries 0 or 1. For the real ploidy labels these are the groups of the
# sample sheet (diploid control, diploid heat, triploid control, triploid heat).
group_membership <- function(labels) {
  cbind(labels * in_first_cohort, labels * (1 - in_first_cohort),
        (1 - labels) * in_first_cohort, (1 - labels) * (1 - in_first_cohort))
}
real_groups <- rowsum(group_membership(real_labels$ploidy), as.integer(lib$group))   # sheet groups by columns
stopifnot(all(rowSums(real_groups > 0) == 1), all(colSums(real_groups > 0) == 1), all(colSums(real_groups) >= 2))
# Pooled within-group SD of each row of x (CpGs by libraries) around the means
# of the groups in `membership`; sum_sq = rowSums(x^2).
pooled_within_sd <- function(x, sum_sq, membership) {
  group_sums <- x %*% membership
  within_ss <- sum_sq - drop(group_sums^2 %*% (1 / colSums(membership)))
  sqrt(pmax(within_ss, 0) / (nrow(membership) - ncol(membership)))
}
# Counts of the 2 x 2 table of step 07: the DML of a list against the other
# CpGs of the background, by overlap with a track. Order: DML in the track, DML
# outside it, other CpGs in the track, other CpGs outside it.
table_counts <- function(dml, background, feature) {
  other <- background & !dml
  c(sum(dml & feature), sum(dml & !feature), sum(other & feature), sum(other & !feature))
}
# Conditional odds ratio of fisher.test() for counts in that order, with or
# without its 95% interval.
fisher_odds <- function(k, interval = TRUE) {
  ft <- fisher.test(matrix(c(k[1], k[3], k[2], k[4]), 2), conf.int = interval)
  if (interval) c(unname(ft$estimate), ft$conf.int[1], ft$conf.int[2]) else unname(ft$estimate)
}
# Mantel-Haenszel odds ratio from stratum counts. Each argument is a matrix
# with one row per table set (here: per resample) and one column per stratum.
mantel_haenszel <- function(a, b, c, d) {
  n <- a + b + c + d
  rowSums(a * d / n, na.rm = TRUE) / rowSums(b * c / n, na.rm = TRUE)   # an emptied stratum (0 / 0) adds nothing
}
percentile_interval <- function(x) {
  if (!all(is.finite(x))) return(c(NA_real_, NA_real_))
  quantile(x, c(0.025, 0.975), names = FALSE)
}
# TRUE where the read-pooled percent methylation of the libraries labelled 1
# differs from that of the libraries labelled 0 by more than `threshold` points.
# counts_c, counts_w: methylated and total reads, CpGs by libraries. labels:
# libraries by labellings, entries 0 or 1. With p1 = c1 / w1 and p0 = c0 / w0,
# 100 |p1 - p0| > threshold is evaluated as 100 |c1 w0 - c0 w1| > threshold w1 w0.
# All terms are whole numbers that a double holds exactly, so a difference of
# exactly the threshold is never called through rounding.
difference_above_threshold <- function(counts_c, counts_w, labels, threshold) {
  c1 <- counts_c %*% labels; w1 <- counts_w %*% labels
  c0 <- rowSums(counts_c) - c1; w0 <- rowSums(counts_w) - w1
  100 * abs(c1 * w0 - c0 * w1) > threshold * w1 * w0
}
odds_row <- function(mode_name, set_name, track_name, comparison, estimate, estimate_method, interval_method,
                     k, n_strata, notes) {
  data.table(mode = mode_name, dml_set = set_name, track = track_name, comparison = comparison,
             interval_method = interval_method, odds_ratio = estimate[1], ci_low = estimate[2], ci_high = estimate[3],
             n_strata_used = n_strata, estimate_method = estimate_method,
             n_dml = k[1] + k[2], n_dml_in_track = k[1], n_other_cpg = k[3] + k[4], n_other_cpg_in_track = k[3],
             notes = notes)
}
count_row <- function(mode_name, set_name, scheme, quantity, value, bin = "all", bin_width = NA_real_, notes = "") {
  data.table(mode = mode_name, dml_set = set_name, permutation_scheme = scheme, quantity = quantity,
             mean_methylation_bin = bin, bin_width = bin_width, value = as.numeric(value), notes = notes)
}

method_fisher <- "conditional maximum likelihood (fisher.test)"
method_mh     <- "Mantel-Haenszel (mantelhaen.test)"
interval_fisher <- "exact conditional (fisher.test); CpGs treated as independent"
interval_mh     <- "Mantel-Haenszel (mantelhaen.test); CpGs treated as independent"
interval_block  <- sprintf("percentile bootstrap over %g kb blocks, %d resamples", block_size / 1e3, n_boot)
label_all      <- "all tested CpGs"
label_eligible <- "eligible CpGs"
step07_background <- setNames(c("all analysed CpGs", "eligible CpGs"), c(label_all, label_eligible))  # names in step 07
label_bin      <- sprintf("eligible CpGs in the same %g-point bin of mean methylation", methylation_bin_width)
label_bin_sd   <- sprintf("%s and the same %s of within-group SD", label_bin,
                          if (n_sd_groups == 10) "decile" else sprintf("quantile group (%d groups)", n_sd_groups))
stopifnot(coverage_quantile == 0.75)       # the two labels below say "quartile"
label_cov_high <- "eligible CpGs in the top quartile of mean coverage"
label_cov_low  <- "eligible CpGs in the lower three quartiles of mean coverage"
label_cds      <- "eligible CpGs in gene bodies that lie in CDS or in an intron, not both: CDS against intron"

track_gr <- lapply(setNames(tracks, tracks), function(tr) { log_msg("reading track %s", tr); read_feature_track(tr) })
invisible(gc(verbose = FALSE))

# ---- One mode ----------------------------------------------------------------------------
analyse_mode <- function(mode_name) {
  tab05 <- function(f) output_path("05_differential_methylation", "tables", f)
  odds <- list(); counts <- list(); perm_rows <- list(); sd_rows <- list()
  matched_real <- list()                         # comparison 3 for the real lists, by DML list and track

  log_msg("[%s] reading the tested CpGs", mode_name)
  united <- fread(tab05(sprintf("united_cpg_%s.csv.gz", mode_name)))
  stopifnot(sum(grepl("^numCs", names(united))) == n_lib, !anyNA(united))
  n_cpg <- nrow(united)
  mean_pm <- mean_percent_methylation(united)
  eligible <- dml_eligible(mean_pm, config)
  num_c <- as.matrix(united[, paste0("numCs", seq_len(n_lib)), with = FALSE])
  storage.mode(num_c) <- "double"
  num_w <- num_c + as.matrix(united[, paste0("numTs", seq_len(n_lib)), with = FALSE])
  mean_coverage <- rowMeans(as.matrix(united[, paste0("coverage", seq_len(n_lib)), with = FALSE]))
  position <- data.table(chr = united$chr, start = united$start)
  rm(united); gc(verbose = FALSE)
  # The cross-multiplied difference rule is exact as long as these products fit a double.
  stopifnot(min(num_w) > 0, 100 * max(rowSums(num_w))^2 < 2^53)

  # Per-library percent methylation: within-group SD and span.
  pm <- 100 * num_c / num_w
  sd_within <- pooled_within_sd(pm, rowSums(pm^2), group_membership(real_labels$ploidy))
  pm_high <- pm[, 1]; pm_low <- pm[, 1]
  for (j in seq_len(n_lib)[-1]) { pm_high <- pmax(pm_high, pm[, j]); pm_low <- pmin(pm_low, pm[, j]) }
  pm_span <- pm_high - pm_low
  pm_eligible <- pm[eligible, ]                  # kept for the SD under shuffled labels
  rm(pm, pm_high, pm_low); gc(verbose = FALSE)

  # Strata of the matched comparisons.
  bin_breaks <- seq(0, 100, by = methylation_bin_width)
  meth_bin <- findInterval(mean_pm, bin_breaks, rightmost.closed = TRUE)
  sd_probs <- seq(0, 1, length.out = n_sd_groups + 1)
  sd_cuts <- quantile(sd_within[eligible], sd_probs, names = FALSE)
  sd_group <- findInterval(sd_within, sd_cuts[-c(1, n_sd_groups + 1)]) + 1L
  strata <- setNames(list(meth_bin, (meth_bin - 1L) * n_sd_groups + sd_group), c(label_bin, label_bin_sd))
  high_coverage <- mean_coverage > quantile(mean_coverage[eligible], coverage_quantile, names = FALSE)

  # ---- DML lists; the difference rule and the q-value rule in the real data ----------
  is_dml <- list()
  for (contrast in contrasts) {
    dml_list <- fread(tab05(sprintf("dml_%s_%s.csv", contrast, mode_name)))
    stopifnot(!anyDuplicated(dml_list$united_row), identical(dml_list$chr, position$chr[dml_list$united_row]),
              identical(dml_list$start, position$start[dml_list$united_row]))
    is_dml[[contrast]] <- seq_len(n_cpg) %in% dml_list$united_row

    tested <- fread(tab05(sprintf("diffmeth_%s_%s.csv.gz", contrast, mode_name)),
                    select = c("chr", "start", "qvalue", "meth.diff"))
    stopifnot(identical(tested$chr, position$chr), identical(tested$start, position$start))
    lab <- real_labels[[contrast]]
    pooled <- 100 * (drop(num_c %*% lab) / drop(num_w %*% lab) - drop(num_c %*% (1 - lab)) / drop(num_w %*% (1 - lab)))
    # methylKit's meth.diff is the read-pooled difference, and the exact rule of
    # this script calls the same CpGs as a comparison of meth.diff with the threshold.
    above <- abs(tested$meth.diff) > config$dml_difference
    stopifnot(max(abs(pooled - tested$meth.diff)) < 1e-8,
              identical(drop(difference_above_threshold(num_c, num_w, matrix(lab), config$dml_difference)), above),
              identical(above & tested$qvalue < config$dml_qvalue, is_dml[[contrast]]))
    if (contrast == "ploidy" && !identical(above, is_dml[[contrast]])) {
      stop("[", mode_name, "] ", sum(above & !is_dml[[contrast]]), " CpGs pass the difference rule for ploidy and fail ",
           "the q-value rule. The label permutation of this script applies the difference rule only and ",
           "would not be comparable with the DML list.", call. = FALSE)
    }
    real <- function(quantity, value) count_row(mode_name, contrast, "real labels", quantity, value)
    counts[[length(counts) + 1]] <- rbind(
      real("tested CpGs", n_cpg),
      real(sprintf("CpGs with absolute difference above %g points", config$dml_difference), sum(above)),
      real(sprintf("of these, CpGs with q-value below %g", config$dml_qvalue), sum(above & tested$qvalue < config$dml_qvalue)),
      real(sprintf("of these, CpGs with q-value at or above %g (removed by the q-value rule)", config$dml_qvalue),
           sum(above & tested$qvalue >= config$dml_qvalue)),
      real("largest q-value among CpGs with absolute difference above the threshold", signif(max(tested$qvalue[above]), 4)),
      real(sprintf("CpGs with q-value below %g, any difference", config$dml_qvalue), sum(tested$qvalue < config$dml_qvalue)),
      real("DML in the list of step 05", sum(is_dml[[contrast]])),
      real("DML outside the eligible range", sum(is_dml[[contrast]] & !eligible)))
    rm(tested, dml_list, pooled, above)
  }
  log_msg("[%s] %d tested CpGs, %d eligible, %d ploidy DML, %d heat DML", mode_name, n_cpg, sum(eligible),
          sum(is_dml$ploidy), sum(is_dml$heat))

  # ---- Overlap with the tracks, as in step 07 ---------------------------------------------
  cpg_gr <- cpg_ranges(position$chr, position$start)
  in_track <- vapply(tracks, function(tr) overlaps_track(cpg_gr, track_gr[[tr]]), logical(n_cpg))
  rm(cpg_gr)
  step07 <- as.data.frame(fread(output_path("07_genome_features", "tables",
                                           sprintf("feature_association_tests_%s.csv", mode_name))))

  # ---- Blocks and block resamples -----------------------------------------------------------
  # Blocks are numbered over the CpGs that enter a comparison with eligible CpGs.
  compared <- eligible | is_dml$ploidy | is_dml$heat
  block_key <- paste(position$chr[compared], position$start[compared] %/% block_size)
  block_of <- integer(n_cpg)
  block_of[compared] <- match(block_key, unique(block_key))
  n_blocks <- max(block_of)
  set.seed(seed_bootstrap)
  boot_weights <- t(vapply(seq_len(n_boot), function(i) tabulate(sample.int(n_blocks, n_blocks, replace = TRUE), n_blocks),
                           numeric(n_blocks)))             # resamples by blocks: times each block is drawn
  # Cell counts of every resample (resamples by cells), for the CpGs `rows` whose
  # cells are `cell`. Only the blocks that hold a CpG of a cell are multiplied.
  resample_counts <- function(rows, cell, n_cell) {
    by_block <- matrix(tabulate((cell - 1L) * n_blocks + block_of[rows], n_blocks * n_cell), n_blocks, n_cell)
    resampled <- matrix(0, n_boot, n_cell)
    for (j in seq_len(n_cell)) {
      filled <- which(by_block[, j] > 0L)
      if (length(filled) > 0) resampled[, j] <- boot_weights[, filled, drop = FALSE] %*% by_block[filled, j]
    }
    resampled
  }
  block_note <- sprintf("%d blocks", n_blocks)
  log_msg("[%s] %d blocks of %g kb hold the eligible CpGs and DML", mode_name, n_blocks, block_size / 1e3)

  # ---- Comparisons 1 to 7 -----------------------------------------------------------------------
  for (set_name in contrasts) {
    dml <- is_dml[[set_name]]
    outside_note <- function(background) {
      n_out <- sum(dml & !background)
      if (n_out == 0) "" else sprintf("DML of the list outside the background, counted as DML: %d", n_out)
    }
    for (track_name in tracks) {
      feature <- in_track[, track_name]
      add <- function(...) odds[[length(odds) + 1]] <<- odds_row(mode_name, set_name, track_name, ...)

      # 1. All tested CpGs and eligible CpGs; must equal step 07.
      backgrounds <- setNames(list(rep(TRUE, n_cpg), eligible), c(label_all, label_eligible))
      for (bg in names(backgrounds)) {
        k <- table_counts(dml, backgrounds[[bg]], feature)
        est <- fisher_odds(k)
        ref <- step07[step07$background == step07_background[[bg]] & step07$set == paste(set_name, "DML") &
                        step07$track == track_name, ]
        if (nrow(ref) != 1 || ref$dml_in_track != k[1] || ref$dml_total != k[1] + k[2] ||
            !isTRUE(all.equal(c(ref$odds_ratio, ref$ci_low, ref$ci_high), est, tolerance = 1e-9))) {
          stop(sprintf("[%s] %s DML, track %s, %s: the odds ratio differs from step 07", mode_name, set_name, track_name, bg),
               call. = FALSE)
        }
        add(bg, est, method_fisher, interval_fisher, k, 1L, outside_note(backgrounds[[bg]]))
      }

      # 7 for comparison 1 (eligible): the same estimate, interval from block resamples.
      use <- which(eligible | dml)
      cell4 <- 1L + 2L * (!dml[use]) + (!feature[use])       # 1 DML in, 2 DML out, 3 other in, 4 other out
      k <- tabulate(cell4, 4L)
      resampled <- resample_counts(use, cell4, 4L)
      boot_est <- apply(resampled, 1, fisher_odds, interval = FALSE)
      add(label_eligible, c(fisher_odds(k, interval = FALSE), percentile_interval(boot_est)), method_fisher, interval_block,
          k, 1L, paste(c(block_note, if (!all(is.finite(boot_est))) sprintf("%d resamples without a finite odds ratio, no interval",
                                                                             sum(!is.finite(boot_est)))), collapse = "; "))

      # 2 and 3. Eligible CpGs in the same stratum; 7 for both.
      for (strat in names(strata)) {
        stratum <- strata[[strat]][use]
        levels_present <- sort(unique(stratum))
        stratum_index <- match(stratum, levels_present)
        by_stratum <- matrix(tabulate((stratum_index - 1L) * 4L + cell4, 4L * length(levels_present)), ncol = 4, byrow = TRUE)
        informative <- by_stratum[, 1] + by_stratum[, 2] > 0 & by_stratum[, 3] + by_stratum[, 4] > 0 &
          by_stratum[, 1] + by_stratum[, 3] > 0 & by_stratum[, 2] + by_stratum[, 4] > 0
        used <- by_stratum[informative, , drop = FALSE]
        stopifnot(nrow(used) >= 2)
        mh <- mantelhaen.test(array(t(used[, c(1, 3, 2, 4)]), c(2, 2, nrow(used))))
        stopifnot(isTRUE(all.equal(unname(mh$estimate),
                                   mantel_haenszel(t(used[, 1]), t(used[, 2]), t(used[, 3]), t(used[, 4])))))
        k_used <- colSums(used)
        strata_note <- sprintf("%d of %d strata that hold a CpG are used; %d of %d DML lie in them", nrow(used),
                               length(levels_present), k_used[1] + k_used[2], sum(dml))
        add(strat, c(unname(mh$estimate), mh$conf.int[1], mh$conf.int[2]), method_mh, interval_mh, k_used, nrow(used), strata_note)
        if (strat == label_bin_sd) matched_real[[set_name]][track_name] <- unname(mh$estimate)

        in_used <- informative[stratum_index]
        cell <- (match(stratum_index[in_used], which(informative)) - 1L) * 4L + cell4[in_used]
        resampled <- resample_counts(use[in_used], cell, 4L * nrow(used))
        pick <- function(j) resampled[, seq(j, ncol(resampled), by = 4L), drop = FALSE]
        boot_est <- mantel_haenszel(pick(1), pick(2), pick(3), pick(4))
        add(strat, c(unname(mh$estimate), percentile_interval(boot_est)), method_mh, interval_block, k_used, nrow(used),
            paste(c(strata_note, block_note, if (!all(is.finite(boot_est))) sprintf("%d resamples without a finite odds ratio, no interval",
                                                                                 sum(!is.finite(boot_est)))), collapse = "; "))
      }

      # 4. Other ranges of mean methylation. dml_eligible() gives the range of a
      #    difference threshold of twice the lower bound.
      for (lower in range_lower_bounds) {
        in_range <- dml_eligible(mean_pm, list(dml_difference = 2 * lower))
        k <- table_counts(dml, in_range, feature)
        add(sprintf("CpGs with mean methylation %g%% to %g%%", lower, 100 - lower), fisher_odds(k), method_fisher,
            interval_fisher, k, 1L, outside_note(in_range))
      }

      # 5. By mean coverage, for one track.
      if (track_name == coverage_track) {
        for (top in c(TRUE, FALSE)) {
          k <- table_counts(dml & high_coverage == top, eligible & high_coverage == top, feature)
          add(if (top) label_cov_high else label_cov_low, fisher_odds(k), method_fisher, interval_fisher, k, 1L,
              "quartiles of mean normalized coverage among eligible CpGs; DML are split by the same cut")
        }
      }

      # 6. Within gene bodies: CDS against intron.
      if (track_name == "CDS") {
        one_of_two <- in_track[, "gene"] & xor(in_track[, "CDS"], in_track[, "intron"])
        k <- table_counts(dml & one_of_two, eligible & one_of_two, feature)
        add(label_cds, fisher_odds(k), method_fisher, interval_fisher, k, 1L, "")
      }
    }

    # Within-group SD of DML and of the other eligible CpGs, overall and by decile.
    other <- eligible & !dml
    sd_summary <- function(label, keep, from, to) {
      data.table(mode = mode_name, dml_set = set_name, sd_group = label, sd_from = from, sd_to = to,
                 n_dml = sum(dml & keep), n_other_eligible = sum(other & keep),
                 median_sd_dml = median(sd_within[dml & keep]), median_sd_other_eligible = median(sd_within[other & keep]),
                 pct_dml_in_gene_body = if (any(dml & keep)) 100 * mean(in_track[dml & keep, "gene"]) else NA_real_,
                 pct_other_eligible_in_gene_body = 100 * mean(in_track[other & keep, "gene"]))
    }
    sd_rows[[length(sd_rows) + 1]] <- rbind(
      sd_summary("all", rep(TRUE, n_cpg), sd_cuts[1], sd_cuts[n_sd_groups + 1]),
      rbindlist(lapply(seq_len(n_sd_groups), function(g) sd_summary(sprintf("%d of %d", g, n_sd_groups), sd_group == g,
                                                                    sd_cuts[g], sd_cuts[g + 1]))))
  }
  rm(boot_weights); gc(verbose = FALSE)

  # ---- Comparison 8: ploidy labels shuffled ----------------------------------------------------
  log_msg("[%s] label permutations", mode_name)
  candidate <- which(pm_span > config$dml_difference)        # the only CpGs that any labelling can call
  cand_c <- num_c[candidate, ]; cand_w <- num_w[candidate, ]
  called_real <- candidate[drop(difference_above_threshold(cand_c, cand_w, matrix(real_labels$ploidy), config$dml_difference))]
  stopifnot(identical(called_real, which(is_dml$ploidy)))
  # One labelling of special interest: the ploidy labels exchanged in one cohort
  # only, which is the ploidy by cohort contrast. Both schemes can draw it.
  exchanged <- ifelse(in_first_cohort == 1, 1 - real_labels$ploidy, real_labels$ploidy)
  counts[[length(counts) + 1]] <- count_row(
    mode_name, "ploidy", "ploidy labels exchanged in one cohort only", "CpGs called",
    sum(difference_above_threshold(cand_c, cand_w, matrix(exchanged), config$dml_difference)),
    notes = "the ploidy by cohort contrast; one of the possible labellings of both permutation schemes")
  cand_eligible <- as.numeric(eligible[candidate])
  n_eligible <- sum(eligible)
  cand_track <- in_track[candidate, ] + 0
  eligible_in_track <- colSums(in_track[eligible, ])
  real_odds <- vapply(tracks, function(tr) fisher_odds(table_counts(is_dml$ploidy, eligible, in_track[, tr]), interval = FALSE),
                      numeric(1))
  # Comparison 3 for the list that one labelling calls. The within-group SD and
  # its deciles are those of the four groups that this labelling defines (label
  # by cohort), as they are for the real labels: a CpG called under shuffled
  # labels differs between the shuffled groups, and that difference must stay
  # out of its SD as the real difference stays out of the SD of a real DML.
  # Called CpGs outside the eligible range have no comparison CpGs and add nothing.
  eligible_sum_sq <- rowSums(pm_eligible^2)
  stratum_offset <- (meth_bin[eligible] - 1L) * n_sd_groups
  n_strata_all <- max(meth_bin[eligible]) * n_sd_groups
  track_cell <- 2L - in_track[eligible, ]                    # 1 inside the track, 2 outside
  cand_row_in_eligible <- match(candidate, which(eligible))  # NA for a candidate outside the eligible range
  matched_odds_for <- function(called_candidates, labels) {
    is_called <- logical(n_eligible)
    is_called[cand_row_in_eligible[called_candidates & !is.na(cand_row_in_eligible)]] <- TRUE
    sd_now <- pooled_within_sd(pm_eligible, eligible_sum_sq, group_membership(labels))
    cuts <- quantile(sd_now, sd_probs[-c(1, n_sd_groups + 1)], names = FALSE)
    stratum <- stratum_offset + findInterval(sd_now, cuts) + 1L
    vapply(seq_along(tracks), function(j) {
      cell <- (stratum - 1L) * 4L + 2L * (!is_called) + track_cell[, j]      # cells 1 to 4 as in comparison 3
      m <- matrix(tabulate(cell, 4L * n_strata_all), ncol = 4, byrow = TRUE)
      mantel_haenszel(t(m[, 1]), t(m[, 2]), t(m[, 3]), t(m[, 4]))
    }, numeric(1))
  }
  # With the real labels this is comparison 3 of the ploidy list. The tolerance
  # allows for a CpG at a decile boundary that falls on the other side when the
  # SD is computed in a separate matrix product; an error in the method would
  # show as a far larger difference.
  stopifnot(isTRUE(all.equal(matched_odds_for(candidate %in% called_real, real_labels$ploidy),
                             unname(matched_real$ploidy[tracks]), tolerance = 1e-3)))
  in_bin <- function(width) findInterval(mean_pm, seq(0, 100, by = width), rightmost.closed = TRUE)
  bin_label <- function(index, width) sprintf("%g to %g", (index - 1) * width, index * width)
  bin_widths <- unique(c(methylation_bin_width, 2 * methylation_bin_width))
  for (width in bin_widths) {
    bins <- in_bin(width)
    for (b in sort(unique(bins[eligible]))) {
      here_bin <- eligible & bins == b
      counts[[length(counts) + 1]] <- rbind(
        count_row(mode_name, "ploidy", "real labels", "eligible CpGs in the bin", sum(here_bin), bin_label(b, width), width),
        count_row(mode_name, "ploidy", "real labels", "eligible CpGs that are ploidy DML, percent of the eligible CpGs in the bin",
                  100 * mean(is_dml$ploidy[here_bin]), bin_label(b, width), width))
    }
  }

  for (scheme in names(label_draws)) {
    draws <- label_draws[[scheme]]
    n_called <- numeric(n_perm); n_called_eligible <- numeric(n_perm)
    called_in_track <- matrix(0, n_perm, length(tracks)); called_eligible_in_track <- matrix(0, n_perm, length(tracks))
    times_called <- numeric(length(candidate))
    matched_perm <- matrix(0, n_perm, length(tracks))
    for (first in seq(1, n_perm, by = perm_chunk)) {
      cols <- first:min(n_perm, first + perm_chunk - 1)
      called <- difference_above_threshold(cand_c, cand_w, draws[, cols, drop = FALSE], config$dml_difference) + 0
      for (i in seq_along(cols)) matched_perm[cols[i], ] <- matched_odds_for(called[, i] == 1, draws[, cols[i]])
      n_called[cols] <- colSums(called)
      n_called_eligible[cols] <- drop(crossprod(cand_eligible, called))
      called_in_track[cols, ] <- t(crossprod(cand_track, called))
      called_eligible_in_track[cols, ] <- t(crossprod(cand_track * cand_eligible, called))
      times_called <- times_called + rowSums(called)
    }
    rm(called)
    # Odds ratio of comparison 1 (eligible) for every permuted list.
    perm_odds <- vapply(seq_along(tracks), function(j) vapply(seq_len(n_perm), function(p) {
      other_in <- eligible_in_track[j] - called_eligible_in_track[p, j]
      fisher_odds(c(called_in_track[p, j], n_called[p] - called_in_track[p, j], other_in,
                    n_eligible - n_called_eligible[p] - other_in), interval = FALSE)
    }, numeric(1)), numeric(n_perm))
    times_all <- numeric(n_cpg); times_all[candidate] <- times_called
    same_as_real <- sum(colSums(draws == real_labels$ploidy) %in% c(0, n_lib))
    scheme_note <- sprintf("%d permutations; %d of them equal the real labels or their mirror image", n_perm, same_as_real)
    for (j in seq_along(tracks)) {
      inside <- eligible & in_track[, j]; outside <- eligible & !in_track[, j]
      perm_rows[[length(perm_rows) + 1]] <- data.table(
        mode = mode_name, dml_set = "ploidy", permutation_scheme = scheme, track = tracks[j], n_permutations = n_perm,
        odds_ratio_real = real_odds[[j]], perm_median = median(perm_odds[, j]),
        perm_percentile_2.5 = quantile(perm_odds[, j], 0.025, names = FALSE),
        perm_percentile_97.5 = quantile(perm_odds[, j], 0.975, names = FALSE),
        share_perm_at_least_real = mean(perm_odds[, j] >= real_odds[[j]]),
        matched_odds_ratio_real = matched_real$ploidy[[tracks[j]]], matched_perm_median = median(matched_perm[, j]),
        matched_perm_percentile_2.5 = quantile(matched_perm[, j], 0.025, names = FALSE),
        matched_perm_percentile_97.5 = quantile(matched_perm[, j], 0.975, names = FALSE),
        matched_share_perm_at_least_real = mean(matched_perm[, j] >= matched_real$ploidy[[tracks[j]]]),
        pct_eligible_called_in_track = 100 * sum(times_all[inside]) / (sum(inside) * n_perm),
        pct_eligible_called_outside_track = 100 * sum(times_all[outside]) / (sum(outside) * n_perm),
        notes = scheme_note)
    }
    perm_count <- function(quantity, value) count_row(mode_name, "ploidy", scheme, quantity, value, notes = scheme_note)
    counts[[length(counts) + 1]] <- rbind(
      perm_count("permutations", n_perm),
      perm_count("CpGs called per permutation, median", median(n_called)),
      perm_count("CpGs called per permutation, 5th percentile", quantile(n_called, 0.05, names = FALSE)),
      perm_count("CpGs called per permutation, 95th percentile", quantile(n_called, 0.95, names = FALSE)),
      perm_count("CpGs called per permutation, smallest", min(n_called)),
      perm_count("CpGs called per permutation, largest", max(n_called)),
      perm_count("share of permutations that call at least as many CpGs as the real labels", mean(n_called >= sum(is_dml$ploidy))),
      perm_count("called CpGs outside the eligible range, percent of all called CpGs", 100 * (1 - sum(n_called_eligible) / sum(n_called))))
    for (width in bin_widths) {
      bins <- in_bin(width)
      for (b in sort(unique(bins[eligible]))) {
        here_bin <- eligible & bins == b
        counts[[length(counts) + 1]] <- count_row(
          mode_name, "ploidy", scheme, "eligible CpGs called per permutation, percent of the eligible CpGs in the bin",
          100 * sum(times_all[here_bin]) / (sum(here_bin) * n_perm), bin_label(b, width), width, scheme_note)
      }
    }
    log_msg("[%s] %s: median %g CpGs called per permutation (real labels %d)", mode_name, scheme, median(n_called),
            sum(is_dml$ploidy))
  }
  # Heat: no permutation.
  heat_odds <- vapply(tracks, function(tr) fisher_odds(table_counts(is_dml$heat, eligible, in_track[, tr]), interval = FALSE),
                      numeric(1))
  perm_rows[[length(perm_rows) + 1]] <- data.table(
    mode = mode_name, dml_set = "heat", permutation_scheme = "none", track = tracks, n_permutations = 0L,
    odds_ratio_real = heat_odds,
    perm_median = NA_real_, perm_percentile_2.5 = NA_real_, perm_percentile_97.5 = NA_real_, share_perm_at_least_real = NA_real_,
    matched_odds_ratio_real = unname(matched_real$heat[tracks]), matched_perm_median = NA_real_,
    matched_perm_percentile_2.5 = NA_real_, matched_perm_percentile_97.5 = NA_real_, matched_share_perm_at_least_real = NA_real_,
    pct_eligible_called_in_track = NA_real_, pct_eligible_called_outside_track = NA_real_,
    notes = no_heat_permutation)
  counts[[length(counts) + 1]] <- count_row(mode_name, "heat", "none", "CpGs called per permutation", NA_real_,
                                            notes = no_heat_permutation)

  list(odds = rbindlist(odds), sd = rbindlist(sd_rows), perm = rbindlist(perm_rows), counts = rbindlist(counts))
}

results <- lapply(config$modes, analyse_mode)
odds_tab  <- rbindlist(lapply(results, `[[`, "odds"))
sd_tab    <- rbindlist(lapply(results, `[[`, "sd"))
perm_tab  <- rbindlist(lapply(results, `[[`, "perm"))
count_tab <- rbindlist(lapply(results, `[[`, "counts"))

# ---- Write ---------------------------------------------------------------------------------
round_columns <- function(tab, columns, digits) {
  tab[, (columns) := lapply(.SD, round, digits), .SDcols = columns]
}
round_columns(odds_tab, c("odds_ratio", "ci_low", "ci_high"), 4)
round_columns(sd_tab, c("sd_from", "sd_to", "median_sd_dml", "median_sd_other_eligible", "pct_dml_in_gene_body",
                        "pct_other_eligible_in_gene_body"), 2)
perm_tab[, ratio_in_to_outside := pct_eligible_called_in_track / pct_eligible_called_outside_track]
setcolorder(perm_tab, c(setdiff(names(perm_tab), "notes"), "notes"))
round_columns(perm_tab, setdiff(names(perm_tab), c("mode", "dml_set", "permutation_scheme", "track", "n_permutations", "notes")), 4)
count_tab[, value := signif(value, 6)]
fwrite(odds_tab, out_file("feature_odds_matched.csv"))
fwrite(sd_tab, out_file("feature_odds_matched_within_group_sd.csv"))
fwrite(perm_tab, out_file("feature_odds_label_permutation.csv"))
fwrite(count_tab, out_file("feature_odds_label_permutation_counts.csv"))

# ---- Summary for the log ---------------------------------------------------------------------
options(width = 250)
show <- odds_tab[, .(mode, dml_set, track, comparison = substr(comparison, 1, 62),
                     interval = fifelse(interval_method == interval_block, "block bootstrap", "CpGs independent"),
                     odds_ratio, ci_low, ci_high, n_strata_used, n_dml)]
cat("\nOdds ratio of a DML lying in a gene body, by comparison:\n"); print(show[track == "gene"], nrows = 200)
cat("\nOdds ratios against eligible CpGs in the same methylation bin and SD decile:\n")
print(show[odds_tab$comparison == label_bin_sd, !"comparison"], nrows = 200)
cat("\nLabel permutation, ploidy:\n")
print(perm_tab[dml_set == "ploidy", !c("dml_set", "n_permutations", "notes", "pct_eligible_called_in_track",
                                       "pct_eligible_called_outside_track")], nrows = 200)
cat("\nMedian within-group SD (percentage points):\n"); print(sd_tab[sd_group == "all"])
cat("\nDifference rule, q-value rule, and CpGs called per permutation:\n")
print(count_tab[mean_methylation_bin == "all", .(mode, dml_set, permutation_scheme, quantity = substr(quantity, 1, 80),
                                                 value = sprintf("%.6g", value))], nrows = 300)
log_msg("feature_odds_matched finished in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))

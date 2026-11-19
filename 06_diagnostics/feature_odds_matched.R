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
#                show how much comparison 1 depends on the range. Each range is
#                given twice: with every DML of the list, as in comparison 1,
#                and with the DML inside the range only.
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
#                Comparison 1 is computed for every labelling within cohort,
#                comparison 3 for 1,000 labellings drawn at random.
#           It also writes the within-group SD of DML and of the other eligible
#           CpGs (medians, overall and inside each SD decile).
#
# Run       Rscript 06_diagnostics/feature_odds_matched.R     (about 10 minutes; about 2.1 GB of memory)
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
#               per track and permutation scheme: the odds ratios of the real list; their
#               distribution over every labelling (comparison 1) and over the random
#               labellings (comparison 3, with Monte Carlo ranges); share of eligible CpGs
#               called inside and outside the track
#           06_diagnostics/feature_odds_label_permutation_counts.csv
#               difference rule, q-value rule and exact ties in the real data; CpGs called
#               per labelling; share of eligible CpGs called, by bin of mean methylation
#
# Notes     - Convention of step 07: every DML of a list is compared with the
#             CpGs of the background that are not in that list. The few DML
#             outside the background stay in the DML row, and the notes column
#             gives their number. For the range 15% to 85% they are several
#             dozen, compared with a background that holds no CpG at their
#             methylation level; this is why comparison 4 is also given with
#             the DML inside the range only. Everything else keeps the
#             convention.
#           - Within-group SD of a CpG: the pooled standard deviation of the
#             per-library percent methylation around the means of the four
#             groups (diploid control, diploid heat, triploid control, triploid
#             heat), 18 degrees of freedom. Differences between the groups do
#             not enter it. Mean and SD come from the same 22 libraries as the
#             contrasts. Deciles are ten groups of equal size by rank of the SD
#             among the eligible CpGs; equal SDs keep the order of the rows.
#             Ranks are used and not cut points, because a cut point can
#             coincide with the SD of a CpG, and its decile would then depend
#             on the last bit of a number. A DML outside the eligible range
#             gets the decile of the eligible CpG next below it in SD.
#           - Methylation bins are closed on the left: 10 to 15 means at least
#             10 and below 15. A stratum is used when it holds a DML and
#             another CpG, and CpGs inside and outside the track; the other
#             strata add nothing to a Mantel-Haenszel odds ratio. The output
#             gives the number of strata used and of DML in them.
#           - Unstratified odds ratios are the conditional estimates of
#             fisher.test(), as in step 07; stratified ones are the common odds
#             ratio of mantelhaen.test() with its interval. Both intervals
#             treat CpGs as independent. fisher.test() finds its estimate and
#             its limits by a root search with the default tolerance of
#             uniroot(), about 0.0001 on the odds ratio where it is below 1 and
#             on its inverse where it is above 1. They are therefore accurate
#             to about four significant digits, and the last digit written can
#             differ from a recomputation with other software. The method is
#             kept because equality with step 07 is the purpose of comparison
#             1. Deciles are a coarse match; the SD table shows what difference
#             in SD is left inside each decile, and n_sd_groups below sets a
#             finer grouping.
#           - Block interval (comparison 7): the sequences are cut into 100 kb
#             windows; each of 10,000 resamples draws, with replacement, as many
#             windows as hold a compared CpG; the interval runs from the 2.5th
#             to the 97.5th percentile of the odds ratio recomputed on each
#             resample. A bootstrap is used and not a delete-one-block
#             jackknife, because a jackknife interval rests on a normal
#             distribution of the log odds ratio, and that fails where one
#             window holds most of the DML of a track (in legacy mode 53 of
#             the 65 ploidy DML in lncRNA lie in the locus of issue I-09). The
#             ends of a block interval carry resampling error. The notes column
#             of every block row gives a 95% range for each end, taken from
#             the order statistics of the resamples (quantile_mc_range());
#             digits of an end beyond that range are not determined.
#           - Label permutation (comparison 8): the ploidy labels are shuffled
#             within cohort (5 and 5 in the heat cohort, 6 and 6 in the control
#             cohort), the read-pooled difference triploid minus diploid is
#             recomputed from the normalized counts (methylKit's meth.diff:
#             pooled methylated reads over pooled reads of each group), and the
#             CpGs whose absolute difference exceeds config$dml_difference are
#             called. A second scheme, reported beside the first, also keeps
#             the acute heat shock balanced (3 and 3, and 2 and 2, in the heat
#             cohort; issue I-04), as the real labels do.
#           - Every labelling is evaluated, not a random sample of them:
#             252 x 924 = 232,848 labellings within cohort, of which 110,880
#             are balanced on the acute heat shock. A labelling and its mirror
#             image (all labels exchanged) call the same CpGs, so one of each
#             pair is computed. The number of CpGs called per labelling, its
#             percentiles, the share of labellings that reach the real number,
#             the share of eligible CpGs called by methylation bin and by
#             track, and the distribution of the odds ratio of comparison 1
#             (eligible) are therefore exact for each scheme and carry no
#             Monte Carlo error. The real labelling and its mirror image are
#             among the labellings and are counted in every share.
#           - The odds ratio of a labelling is the cross-product ratio
#             (a d / b c) of its 2 x 2 table. The conditional estimate of
#             fisher.test() for every labelling and track would add several
#             hundred thousand root searches per mode, and for the small odds
#             ratios of some labellings its tolerance is coarse. The
#             permutation table gives both estimates for the real list, and
#             the share of labellings at or above the real value compares
#             cross-product ratios on both sides.
#           - Comparison 3 under shuffled labels is computed for 1,000 random
#             labellings per scheme, because it needs the within-group SD and
#             its deciles again for every labelling: they are computed around
#             the means of the four groups that the shuffled labels define
#             (label by cohort). A CpG called under shuffled labels differs
#             between the shuffled groups, and that difference has to stay out
#             of its SD exactly as the real group difference stays out of the
#             SD of a real DML; with the SD of the real groups the shuffled
#             lists would be matched to CpGs that are too variable. These
#             columns carry Monte Carlo error: every percentile comes with a
#             95% range from order statistics and every share with its
#             standard error. A difference between the two schemes that lies
#             inside those ranges is sampling error. Each random labelling is
#             also looked up in the enumeration, and the script stops unless
#             both give the same counts.
#           - THE PERMUTATION APPLIES THE DIFFERENCE RULE ONLY, NOT THE Q-VALUE
#             RULE: refitting the logistic regression at every CpG for every
#             labelling is not feasible here. Evidence that the q-value rule
#             does not bind in the real data, counted from
#             diffmeth_<contrast>_<mode>.csv.gz and written to the counts
#             table: every CpG above the difference threshold for ploidy also
#             has q below config$dml_qvalue (1,038 of 1,038 in legacy mode, 952
#             of 952 in corrected mode; the largest q among them is about half
#             the cutoff), and the script stops if that is not so. For heat one
#             CpG of 2,352 fails the q-value rule in legacy mode and none in
#             corrected mode. What the q-values would be under shuffled labels
#             is not computed, so a count per labelling is an upper bound on
#             what the full rule would call.
#           - For heat no label permutation exists: heat coincides with cohort
#             (issue I-01), so no shuffle of the heat labels keeps the cohorts
#             apart. The output says so; nothing is put in its place.
#           - A difference of exactly the threshold is not called. The rule is
#             evaluated on whole numbers (difference_sides()), so that rounding
#             cannot call such a tie. Step 05 compares the stored difference, a
#             double, with the threshold (methylKit::getMethylDiff), so there a
#             tie is called only if its stored value was rounded upwards. The
#             script counts the exact ties of each contrast and those among
#             them that step 05 counted as above the threshold, writes both to
#             the counts table, and requires its own rule and that of step 05
#             to agree at every other CpG. The labellings are compared with
#             the list that the exact rule calls under the real labels.
#           - Only CpGs that some labelling can call are evaluated: those whose
#             per-library percentages span more than the threshold, because a
#             pooled percentage is a weighted mean of per-library percentages.
#             In the enumeration a CpG is evaluated for one choice of labels in
#             the first cohort only if some choice in the second cohort could
#             call it. That test compares the largest value that the left side
#             of the cross-multiplied rule can take with the smallest value of
#             its right side. It is made on whole numbers, so it cannot drop a
#             CpG through rounding. The random labellings are evaluated without
#             it, and their lookup in the enumeration checks it.
#           - The counts table also gives the number of CpGs called when the
#             ploidy labels are exchanged in one cohort only. That labelling is
#             the ploidy by cohort contrast and is one of the labellings of
#             both schemes. It calls more than three times as many CpGs as the
#             real labels (see also issues I-02 and I-07).
#           - Random numbers: set.seed(20261119) before the random labellings
#             (1,000 per scheme, the same labellings for both modes) and
#             set.seed(20261120) before the block resamples of each mode
#             (10,000). Two runs give identical files.

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
n_boot                <- 10000      # block resamples
n_perm                <- 1000       # random labellings per scheme, for comparison 3 under shuffled labels
perm_chunk            <- 25         # random labellings evaluated at a time (memory only)
enum_rows             <- 8000       # enumeration of the labellings: CpGs evaluated at a time (memory only)
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

# ---- Label permutations: the random labellings (drawn once, used for both modes) ------
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

# ---- Label permutations: every labelling within cohort ---------------------------------
# A labelling within cohort is a choice of the libraries labelled triploid in
# the first cohort of the sample sheet together with a choice in the second.
# labels_a and labels_b hold the choices (libraries of one cohort by choices,
# entries 0 or 1); the enumeration numbers the labellings first by the column of
# labels_a, then by the column of labels_b.
rows_a <- which(in_first_cohort == 1); rows_b <- which(in_first_cohort == 0)
subsets_as_labels <- function(n, k) {
  chosen <- combn(n, k)
  labels <- matrix(0, n, ncol(chosen))
  labels[cbind(as.vector(chosen), rep(seq_len(ncol(chosen)), each = k))] <- 1
  labels
}
labels_a <- subsets_as_labels(length(rows_a), sum(is_triploid[rows_a]))
labels_b <- subsets_as_labels(length(rows_b), sum(is_triploid[rows_b]))
# With both cohorts split in halves, the mirror image of a labelling (all labels
# exchanged) is a labelling too and calls the same CpGs. One of each pair is
# evaluated: the one that labels the first library of the first cohort diploid.
mirror_pairs <- 2 * sum(is_triploid[rows_a]) == length(rows_a) && 2 * sum(is_triploid[rows_b]) == length(rows_b)
if (mirror_pairs) labels_a <- labels_a[, labels_a[1, ] == 0, drop = FALSE]
labellings_per_row <- if (mirror_pairs) 2 else 1         # labellings that one evaluated labelling stands for
# The schemes: TRUE for the evaluated labellings that keep the number of
# triploid labels of every stratum of the scheme.
keeps_strata <- function(labels, rows, strata) {
  colSums(rowsum(labels, strata[rows]) != as.vector(rowsum(real_labels$ploidy[rows], strata[rows]))) == 0
}
in_scheme <- lapply(perm_strata, function(s) rep(keeps_strata(labels_a, rows_a, s), each = ncol(labels_b)) &
                                              rep(keeps_strata(labels_b, rows_b, s), times = ncol(labels_a)))
stopifnot(all(labellings_per_row * vapply(in_scheme, sum, numeric(1)) ==     # as many as the binomial coefficients say
                vapply(perm_strata, function(s) prod(choose(table(s), tapply(is_triploid, s, sum))), numeric(1))))
# Row of each labelling (libraries by labellings, entries 0 or 1) in the enumeration.
subset_code <- function(labels) drop(2^(seq_len(nrow(labels)) - 1) %*% labels)
labelling_row <- function(labels) {
  if (mirror_pairs) { mirrored <- labels[rows_a[1], ] == 1; labels[, mirrored] <- 1 - labels[, mirrored] }
  index_a <- match(subset_code(labels[rows_a, , drop = FALSE]), subset_code(labels_a))
  index_b <- match(subset_code(labels[rows_b, , drop = FALSE]), subset_code(labels_b))
  stopifnot(!anyNA(index_a), !anyNA(index_b))
  (index_a - 1L) * ncol(labels_b) + index_b
}
real_row <- labelling_row(matrix(real_labels$ploidy))
draw_rows <- lapply(label_draws, labelling_row)
stopifnot(all(vapply(names(perm_strata), function(s) in_scheme[[s]][real_row] && all(in_scheme[[s]][draw_rows[[s]]]), logical(1))))

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
# Group of each value among n_groups groups of equal size (to within one), by
# rank. Equal values keep the order in which they stand.
rank_group <- function(values, n_groups) {
  as.integer(ceiling(n_groups * rank(values, ties.method = "first") / length(values)))
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
# Monte Carlo error of a percentile of random draws: 95% range for the q-th
# quantile of the distribution the draws come from, from order statistics. The
# number of draws at or below that quantile is binomial, so the quantile lies
# between the two order statistics below with a probability of at least 0.95,
# whatever the distribution.
quantile_mc_range <- function(draws, q) {
  n <- length(draws)
  lower <- qbinom(0.025, n, q); upper <- qbinom(0.975, n, q) + 1
  if (!all(is.finite(draws)) || lower < 1 || upper > n) return(c(NA_real_, NA_real_))
  sort(draws)[c(lower, upper)]
}
# Percentile interval of the block resamples of an odds ratio, with a note on
# the resampling error of its two ends.
block_interval <- function(boot_est) {
  if (!all(is.finite(boot_est))) {
    return(list(ends = c(NA_real_, NA_real_),
                note = sprintf("%d resamples without a finite odds ratio, no interval", sum(!is.finite(boot_est)))))
  }
  low <- quantile_mc_range(boot_est, 0.025); high <- quantile_mc_range(boot_est, 0.975)
  list(ends = quantile(boot_est, c(0.025, 0.975), names = FALSE),
       note = sprintf("resampling error, 95%% range of each end from order statistics: ci_low %.4f to %.4f, ci_high %.4f to %.4f",
                      low[1], low[2], high[1], high[2]))
}
# The two sides of the difference rule, for CpGs (rows) and labellings (columns
# of `labels`, libraries by labellings, entries 0 or 1). counts_c, counts_w:
# methylated and total reads, CpGs by libraries. With p1 = c1 / w1 the
# read-pooled share of the libraries labelled 1 and p0 = c0 / w0 that of the
# others, 100 |p1 - p0| against the threshold is 100 |c1 w0 - c0 w1| against
# threshold w1 w0. For a whole-number threshold all terms are whole numbers that
# a double holds exactly, so "above" and "equal" are decided without rounding.
difference_sides <- function(counts_c, counts_w, labels, threshold) {
  c1 <- counts_c %*% labels; w1 <- counts_w %*% labels
  c0 <- rowSums(counts_c) - c1; w0 <- rowSums(counts_w) - w1
  list(difference = 100 * abs(c1 * w0 - c0 * w1), limit = threshold * w1 * w0)
}
difference_above_threshold <- function(counts_c, counts_w, labels, threshold) {
  sides <- difference_sides(counts_c, counts_w, labels, threshold)
  sides$difference > sides$limit
}
difference_equals_threshold <- function(counts_c, counts_w, labels, threshold) {
  sides <- difference_sides(counts_c, counts_w, labels, threshold)
  sides$difference == sides$limit
}
# Smallest and largest value in each row of a matrix.
row_extremes <- function(x) {
  low <- x[, 1]; high <- x[, 1]
  for (j in seq_len(ncol(x))[-1]) { low <- pmin(low, x[, j]); high <- pmax(high, x[, j]) }
  list(low = low, high = high)
}
# Counts of called CpGs that describe one list: all, eligible, in each track,
# eligible and in each track. One column per CpG property, in this order.
list_features <- function(is_eligible, track_overlap) {
  features <- cbind(1, is_eligible, track_overlap, track_overlap & is_eligible) + 0
  colnames(features) <- c("called", "called_eligible", paste0("called_in_", tracks), paste0("called_eligible_in_", tracks))
  features
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

method_fisher <- sprintf("conditional maximum likelihood (fisher.test; root search with tolerance %.1g)", .Machine$double.eps^0.25)
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
  threshold <- config$dml_difference

  log_msg("[%s] reading the tested CpGs", mode_name)
  united <- fread(tab05(sprintf("united_cpg_%s.csv.gz", mode_name)))
  stopifnot(sum(grepl("^numCs", names(united))) == n_lib, !anyNA(united))
  n_cpg <- nrow(united)
  mean_pm <- mean_percent_methylation(united)
  eligible <- dml_eligible(mean_pm, config)
  n_eligible <- sum(eligible)
  num_c <- as.matrix(united[, paste0("numCs", seq_len(n_lib)), with = FALSE])
  storage.mode(num_c) <- "double"
  num_w <- num_c + as.matrix(united[, paste0("numTs", seq_len(n_lib)), with = FALSE])
  mean_coverage <- rowMeans(as.matrix(united[, paste0("coverage", seq_len(n_lib)), with = FALSE]))
  position <- data.table(chr = united$chr, start = united$start)
  rm(united); gc(verbose = FALSE)
  # The cross-multiplied difference rule is exact for a whole-number threshold,
  # as long as these products fit a double.
  stopifnot(threshold == round(threshold), min(num_w) > 0, 100 * max(rowSums(num_w))^2 < 2^53)

  # Per-library percent methylation: within-group SD and span.
  pm <- 100 * num_c / num_w
  sd_within <- pooled_within_sd(pm, rowSums(pm^2), group_membership(real_labels$ploidy))
  pm_range <- row_extremes(pm)
  pm_span <- pm_range$high - pm_range$low
  pm_eligible <- pm[eligible, ]                  # kept for the SD under shuffled labels
  rm(pm, pm_range); gc(verbose = FALSE)

  # Strata of the matched comparisons. SD groups are groups of equal size by
  # rank among the eligible CpGs; a CpG outside the eligible range gets the
  # group of the eligible CpG next below it in SD (the first group if there is none).
  bin_breaks <- seq(0, 100, by = methylation_bin_width)
  meth_bin <- findInterval(mean_pm, bin_breaks, rightmost.closed = TRUE)
  sd_group <- integer(n_cpg)
  sd_group[eligible] <- rank_group(sd_within[eligible], n_sd_groups)
  n_at_or_below <- findInterval(sd_within[!eligible], sort(sd_within[eligible]))
  sd_group[!eligible] <- as.integer(ceiling(n_sd_groups * pmax(n_at_or_below, 1L) / n_eligible))
  strata <- setNames(list(meth_bin, (meth_bin - 1L) * n_sd_groups + sd_group), c(label_bin, label_bin_sd))
  high_coverage <- mean_coverage > quantile(mean_coverage[eligible], coverage_quantile, names = FALSE)

  # ---- DML lists; the difference rule and the q-value rule in the real data ----------
  is_dml <- list(); called_exact <- list()
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
    # methylKit's meth.diff is the read-pooled difference. Step 05 compares the
    # stored value with the threshold; this script decides on whole numbers.
    # The two rules can differ only at a CpG whose difference is exactly the
    # threshold, where step 05 calls the CpG if the stored value was rounded upwards.
    above <- abs(tested$meth.diff) > threshold
    called_exact[[contrast]] <- drop(difference_above_threshold(num_c, num_w, matrix(lab), threshold))
    exact_tie <- drop(difference_equals_threshold(num_c, num_w, matrix(lab), threshold))
    stopifnot(max(abs(pooled - tested$meth.diff)) < 1e-8,
              identical(called_exact[[contrast]][!exact_tie], above[!exact_tie]),
              identical(above & tested$qvalue < config$dml_qvalue, is_dml[[contrast]]))
    if (contrast == "ploidy" && !identical(above, is_dml[[contrast]])) {
      stop("[", mode_name, "] ", sum(above & !is_dml[[contrast]]), " CpGs pass the difference rule for ploidy and fail ",
           "the q-value rule. The label permutation of this script applies the difference rule only and ",
           "would not be comparable with the DML list.", call. = FALSE)
    }
    real <- function(quantity, value) count_row(mode_name, contrast, "real labels", quantity, value)
    counts[[length(counts) + 1]] <- rbind(
      real("tested CpGs", n_cpg),
      real(sprintf("CpGs with absolute difference above %g points", threshold), sum(above)),
      real(sprintf("of these, CpGs with q-value below %g", config$dml_qvalue), sum(above & tested$qvalue < config$dml_qvalue)),
      real(sprintf("of these, CpGs with q-value at or above %g (removed by the q-value rule)", config$dml_qvalue),
           sum(above & tested$qvalue >= config$dml_qvalue)),
      real("largest q-value among CpGs with absolute difference above the threshold", signif(max(tested$qvalue[above]), 4)),
      real(sprintf("CpGs with q-value below %g, any difference", config$dml_qvalue), sum(tested$qvalue < config$dml_qvalue)),
      real(sprintf("CpGs with absolute difference of exactly %g points (exact ties, not called by this script)", threshold),
           sum(exact_tie)),
      real("of these, CpGs that step 05 counted as above the threshold (stored difference rounded upwards)",
           sum(exact_tie & above)),
      real("CpGs called by the whole-number rule of this script (for ploidy, the number the label permutations are compared with)",
           sum(called_exact[[contrast]])),
      real("DML in the list of step 05", sum(is_dml[[contrast]])),
      real("DML outside the eligible range", sum(is_dml[[contrast]] & !eligible)))
    rm(tested, dml_list, pooled, above, exact_tie)
  }
  log_msg("[%s] %d tested CpGs, %d eligible, %d ploidy DML, %d heat DML", mode_name, n_cpg, n_eligible,
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
  boot_weights <- matrix(0, n_boot, n_blocks)              # resamples by blocks: times each block is drawn
  for (i in seq_len(n_boot)) boot_weights[i, ] <- tabulate(sample.int(n_blocks, n_blocks, replace = TRUE), n_blocks)
  # Cell counts of every resample (resamples by cells), for the CpGs `rows` whose
  # cells are `cell`. For a cell that fewer than a tenth of the blocks hold, only
  # those blocks are multiplied; this is faster, and the copy of their weights is
  # small. All terms are whole numbers, so both forms give the same counts.
  resample_counts <- function(rows, cell, n_cell) {
    by_block <- matrix(tabulate((cell - 1L) * n_blocks + block_of[rows], n_blocks * n_cell), n_blocks, n_cell)
    resampled <- matrix(0, n_boot, n_cell)
    for (j in seq_len(n_cell)) {
      filled <- which(by_block[, j] > 0L)
      if (10 * length(filled) > n_blocks) {
        resampled[, j] <- boot_weights %*% by_block[, j]
      } else if (length(filled) > 0) {
        resampled[, j] <- boot_weights[, filled, drop = FALSE] %*% by_block[filled, j]
      }
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
      boot <- block_interval(apply(resampled, 1, fisher_odds, interval = FALSE))
      add(label_eligible, c(fisher_odds(k, interval = FALSE), boot$ends), method_fisher, interval_block, k, 1L,
          paste(c(block_note, boot$note), collapse = "; "))

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
        boot <- block_interval(mantel_haenszel(pick(1), pick(2), pick(3), pick(4)))
        add(strat, c(unname(mh$estimate), boot$ends), method_mh, interval_block, k_used, nrow(used),
            paste(c(strata_note, block_note, boot$note), collapse = "; "))
      }

      # 4. Other ranges of mean methylation. dml_eligible() gives the range of a
      #    difference threshold of twice the lower bound. First with every DML of
      #    the list (the convention of step 07), then with the DML inside the range only.
      for (lower in range_lower_bounds) {
        in_range <- dml_eligible(mean_pm, list(dml_difference = 2 * lower))
        range_label <- sprintf("CpGs with mean methylation %g%% to %g%%", lower, 100 - lower)
        k <- table_counts(dml, in_range, feature)
        add(paste0(range_label, ", every DML of the list"), fisher_odds(k), method_fisher, interval_fisher, k, 1L,
            outside_note(in_range))
        k <- table_counts(dml & in_range, in_range, feature)
        add(paste0(range_label, ", DML inside the range only"), fisher_odds(k), method_fisher, interval_fisher, k, 1L,
            sprintf("DML of the list outside the range, left out: %d", sum(dml & !in_range)))
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
    # sd_from and sd_to: smallest and largest SD among the eligible CpGs of the group.
    other <- eligible & !dml
    sd_summary <- function(label, keep) {
      data.table(mode = mode_name, dml_set = set_name, sd_group = label,
                 sd_from = min(sd_within[eligible & keep]), sd_to = max(sd_within[eligible & keep]),
                 n_dml = sum(dml & keep), n_other_eligible = sum(other & keep),
                 median_sd_dml = median(sd_within[dml & keep]), median_sd_other_eligible = median(sd_within[other & keep]),
                 pct_dml_in_gene_body = if (any(dml & keep)) 100 * mean(in_track[dml & keep, "gene"]) else NA_real_,
                 pct_other_eligible_in_gene_body = 100 * mean(in_track[other & keep, "gene"]))
    }
    sd_rows[[length(sd_rows) + 1]] <- rbind(
      sd_summary("all", rep(TRUE, n_cpg)),
      rbindlist(lapply(seq_len(n_sd_groups), function(g) sd_summary(sprintf("%d of %d", g, n_sd_groups), sd_group == g))))
  }
  rm(boot_weights); gc(verbose = FALSE)

  # ---- Comparison 8: ploidy labels shuffled ----------------------------------------------------
  # The list that the exact rule calls under the real labels. It is the list of
  # step 05 unless step 05 called an exact tie (counted above).
  called_real <- called_exact$ploidy
  n_called_real <- sum(called_real)
  eligible_in_track <- colSums(in_track[eligible, ])
  real_odds <- vapply(tracks, function(tr) fisher_odds(table_counts(is_dml$ploidy, eligible, in_track[, tr]), interval = FALSE),
                      numeric(1))
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
  # One labelling of special interest: the ploidy labels exchanged in one cohort
  # only, which is the ploidy by cohort contrast.
  exchanged <- ifelse(in_first_cohort == 1, 1 - real_labels$ploidy, real_labels$ploidy)
  n_called_exchanged <- sum(difference_above_threshold(num_c, num_w, matrix(exchanged), threshold))
  counts[[length(counts) + 1]] <- count_row(
    mode_name, "ploidy", "ploidy labels exchanged in one cohort only", "CpGs called", n_called_exchanged,
    notes = if (mirror_pairs) "the ploidy by cohort contrast; one of the labellings of both permutation schemes"
            else "the ploidy by cohort contrast")

  # ---- 8a. Every labelling within cohort ---------------------------------------------------------
  log_msg("[%s] every labelling of the ploidy labels within cohort", mode_name)
  # A pooled percentage is a weighted mean of per-library percentages, so only a
  # CpG whose per-library percentages span more than the threshold can be called.
  candidate <- which(pm_span > threshold)
  total_c <- rowSums(num_c); total_w <- rowSums(num_w)
  # Called CpGs are counted by pattern of eligibility and track overlap, and the
  # pattern counts are added up to the features of list_features().
  cand_features <- list_features(eligible[candidate], in_track[candidate, ])
  pattern_key <- eligible[candidate] + 2 * drop(in_track[candidate, ] %*% 2^(seq_along(tracks) - 1))
  pattern <- match(pattern_key, sort(unique(pattern_key)))
  pattern_features <- cand_features[match(seq_len(max(pattern)), pattern), , drop = FALSE]
  n_first <- ncol(labels_a); n_second <- ncol(labels_b)
  per_labelling <- matrix(0, n_first * n_second, ncol(cand_features), dimnames = list(NULL, colnames(cand_features)))
  times_called <- matrix(0, length(candidate), length(in_scheme), dimnames = list(NULL, names(in_scheme)))
  for (slice in split(seq_along(candidate), ceiling(seq_along(candidate) / enum_rows))) {
    # The difference rule, cross-multiplied as in difference_sides() and split by
    # cohort. With C, W the methylated and total reads of all libraries and c1, w1
    # those of the libraries labelled 1, c1 w0 - c0 w1 = c1 W - C w1, and both c1
    # and w1 are the sum of a first-cohort part and a second-cohort part. The rule
    # is 100 |c1 W - C w1| > threshold w1 (W - w1).
    cpgs <- candidate[slice]
    slice_total_c <- total_c[cpgs]; slice_total_w <- total_w[cpgs]
    part_a_w <- num_w[cpgs, rows_a] %*% labels_a
    part_a_x <- 100 * ((num_c[cpgs, rows_a] %*% labels_a) * slice_total_w - slice_total_c * part_a_w)
    part_b_w <- num_w[cpgs, rows_b] %*% labels_b
    part_b_x <- 100 * ((num_c[cpgs, rows_b] %*% labels_b) * slice_total_w - slice_total_c * part_b_w)
    range_w <- row_extremes(part_b_w); range_x <- row_extremes(part_b_x)
    slice_pattern <- pattern[slice]
    for (first in seq_len(n_first)) {
      # With this choice in the first cohort, can any choice in the second cohort
      # call the CpG? The left side is at most its value at the smallest or the
      # largest second-cohort part; the right side is at least its value at the
      # smallest or the largest number of reads labelled 1, because w1 (W - w1)
      # has no minimum in between. Whole numbers again, so no CpG is lost to rounding.
      reads_low <- part_a_w[, first] + range_w$low; reads_high <- part_a_w[, first] + range_w$high
      possible <- which(pmax(abs(part_a_x[, first] + range_x$low), abs(part_a_x[, first] + range_x$high)) >
                          threshold * pmin(reads_low * (slice_total_w - reads_low), reads_high * (slice_total_w - reads_high)))
      if (length(possible) == 0) next
      reads_1 <- part_a_w[possible, first] + part_b_w[possible, , drop = FALSE]
      called <- abs(part_a_x[possible, first] + part_b_x[possible, , drop = FALSE]) >
        threshold * reads_1 * (slice_total_w[possible] - reads_1)                      # CpGs by choices in the second cohort
      storage.mode(called) <- "integer"
      by_pattern <- rowsum(called, slice_pattern[possible])                            # row names: the patterns present
      rows_here <- (first - 1L) * n_second + seq_len(n_second)
      per_labelling[rows_here, ] <- per_labelling[rows_here, ] +
        crossprod(by_pattern, pattern_features[as.integer(rownames(by_pattern)), , drop = FALSE])
      for (scheme in names(in_scheme)) {
        members <- in_scheme[[scheme]][rows_here]
        if (any(members)) {
          times_called[slice[possible], scheme] <- times_called[slice[possible], scheme] +
            rowSums(if (all(members)) called else called[, members, drop = FALSE])
        }
      }
    }
  }
  rm(part_a_w, part_a_x, part_b_w, part_b_x, range_w, range_x, called, reads_1); gc(verbose = FALSE)
  # The real labelling in the enumeration must be the list of the exact rule.
  real_features <- colSums(list_features(eligible[called_real], in_track[called_real, , drop = FALSE]))
  stopifnot(all(per_labelling[real_row, ] == real_features))
  if (mirror_pairs) stopifnot(per_labelling[labelling_row(matrix(exchanged)), "called"] == n_called_exchanged)
  # Cross-product odds ratio of comparison 1 (eligible) for the list of every labelling.
  labelling_odds <- vapply(tracks, function(tr) {
    dml_in <- per_labelling[, paste0("called_in_", tr)]
    dml_out <- per_labelling[, "called"] - dml_in
    other_in <- eligible_in_track[[tr]] - per_labelling[, paste0("called_eligible_in_", tr)]
    other_out <- n_eligible - per_labelling[, "called_eligible"] - other_in
    dml_in * other_out / (dml_out * other_in)
  }, numeric(nrow(per_labelling)))

  # ---- 8b. Comparison 3 for random labellings --------------------------------------------------
  log_msg("[%s] random labellings: comparison 3 under shuffled labels", mode_name)
  # Comparison 3 for the list that one labelling calls. The within-group SD and
  # its deciles are those of the four groups that this labelling defines (label
  # by cohort), as they are for the real labels: a CpG called under shuffled
  # labels differs between the shuffled groups, and that difference must stay
  # out of its SD as the real difference stays out of the SD of a real DML.
  # Called CpGs outside the eligible range have no comparison CpGs and add nothing.
  stopifnot(all(pm_span[is_dml$ploidy] > threshold))
  cand_c <- num_c[candidate, ]; cand_w <- num_w[candidate, ]
  eligible_sum_sq <- rowSums(pm_eligible^2)
  stratum_offset <- (meth_bin[eligible] - 1L) * n_sd_groups
  n_strata_all <- max(meth_bin[eligible]) * n_sd_groups
  track_cell <- 2L - in_track[eligible, ]                    # 1 inside the track, 2 outside
  cand_row_in_eligible <- match(candidate, which(eligible))  # NA for a candidate outside the eligible range
  matched_odds_for <- function(called_candidates, labels) {
    is_called <- logical(n_eligible)
    is_called[cand_row_in_eligible[called_candidates & !is.na(cand_row_in_eligible)]] <- TRUE
    sd_now <- pooled_within_sd(pm_eligible, eligible_sum_sq, group_membership(labels))
    stratum <- stratum_offset + rank_group(sd_now, n_sd_groups)
    vapply(seq_along(tracks), function(j) {
      cell <- (stratum - 1L) * 4L + 2L * (!is_called) + track_cell[, j]      # cells 1 to 4 as in comparison 3
      m <- matrix(tabulate(cell, 4L * n_strata_all), ncol = 4, byrow = TRUE)
      mantel_haenszel(t(m[, 1]), t(m[, 2]), t(m[, 3]), t(m[, 4]))
    }, numeric(1))
  }
  # With the real labels and the list of step 05 this is comparison 3 of the ploidy list.
  if (!isTRUE(all.equal(matched_odds_for(is_dml$ploidy[candidate], real_labels$ploidy),
                        unname(matched_real$ploidy[tracks]), tolerance = 1e-9))) {
    stop("[", mode_name, "] comparison 3 of the ploidy list is not reproduced by the code that serves the shuffled labels",
         call. = FALSE)
  }

  for (scheme in names(label_draws)) {
    # Every labelling of the scheme.
    members <- in_scheme[[scheme]]
    n_members <- as.numeric(sum(members))                       # evaluated labellings; a double, because products below exceed the integer range
    n_labellings <- labellings_per_row * n_members
    every <- function(x) rep(x, each = labellings_per_row)      # one value per labelling, mirror images included
    n_called <- per_labelling[members, "called"]
    enum_note <- paste("every labelling of the scheme (complete enumeration, no random draws);",
                       if (mirror_pairs) "the real labelling and its mirror image are among them" else "the real labelling is among them")
    times_all <- numeric(n_cpg); times_all[candidate] <- times_called[, scheme]
    perm_count <- function(quantity, value) count_row(mode_name, "ploidy", scheme, quantity, value, notes = enum_note)
    counts[[length(counts) + 1]] <- rbind(
      perm_count("permutations", n_labellings),
      perm_count("CpGs called per permutation, median", median(every(n_called))),
      perm_count("CpGs called per permutation, 5th percentile", quantile(every(n_called), 0.05, names = FALSE)),
      perm_count("CpGs called per permutation, 95th percentile", quantile(every(n_called), 0.95, names = FALSE)),
      perm_count("CpGs called per permutation, smallest", min(n_called)),
      perm_count("CpGs called per permutation, largest", max(n_called)),
      perm_count("permutations that call at least as many CpGs as the real labels", labellings_per_row * sum(n_called >= n_called_real)),
      perm_count("share of permutations that call at least as many CpGs as the real labels", mean(n_called >= n_called_real)),
      perm_count("called CpGs outside the eligible range, percent of all called CpGs",
                 100 * (1 - sum(per_labelling[members, "called_eligible"]) / sum(n_called))))
    for (width in bin_widths) {
      bins <- in_bin(width)
      for (b in sort(unique(bins[eligible]))) {
        here_bin <- eligible & bins == b
        counts[[length(counts) + 1]] <- count_row(
          mode_name, "ploidy", scheme, "eligible CpGs called per permutation, percent of the eligible CpGs in the bin",
          100 * sum(times_all[here_bin]) / (sum(here_bin) * n_members), bin_label(b, width), width, enum_note)
      }
    }

    # Random labellings of the scheme: comparison 3, and the lookup in the enumeration.
    draws <- label_draws[[scheme]]
    matched_perm <- matrix(0, n_perm, length(tracks))
    draw_features <- matrix(0, n_perm, ncol(cand_features))
    for (first in seq(1, n_perm, by = perm_chunk)) {
      cols <- first:min(n_perm, first + perm_chunk - 1)
      called <- difference_above_threshold(cand_c, cand_w, draws[, cols, drop = FALSE], threshold) + 0
      for (i in seq_along(cols)) matched_perm[cols[i], ] <- matched_odds_for(called[, i] == 1, draws[, cols[i]])
      draw_features[cols, ] <- crossprod(called, cand_features)
    }
    rm(called)
    if (!all(draw_features == per_labelling[draw_rows[[scheme]], ])) {
      stop("[", mode_name, "] ", scheme, ": the enumeration and the random labellings give different counts", call. = FALSE)
    }
    draw_note <- sprintf(paste("%d random labellings of the scheme, drawn with set.seed(%d); they carry Monte Carlo error (95%% range",
                               "of each percentile from order statistics, standard error of the share), and a difference between the",
                               "two schemes inside those ranges is sampling error"), n_perm, seed_permutation)
    counts[[length(counts) + 1]] <- rbind(
      count_row(mode_name, "ploidy", scheme, "random permutations drawn for comparison 3 (matched odds ratios)", n_perm, notes = draw_note),
      count_row(mode_name, "ploidy", scheme, sprintf("of these, the real labelling%s", if (mirror_pairs) " or its mirror image" else ""),
                sum(draw_rows[[scheme]] == real_row), notes = draw_note))
    for (j in seq_along(tracks)) {
      inside <- eligible & in_track[, j]; outside <- eligible & !in_track[, j]
      odds_j <- labelling_odds[members, j]
      finite <- all(is.finite(odds_j))
      matched_j <- matched_perm[, j]
      matched_share <- mean(matched_j >= matched_real$ploidy[[tracks[j]]])
      mc <- function(q) quantile_mc_range(matched_j, q)
      pct_inside <- 100 * sum(times_all[inside]) / (sum(inside) * n_members)
      pct_outside <- 100 * sum(times_all[outside]) / (sum(outside) * n_members)
      perm_rows[[length(perm_rows) + 1]] <- data.table(
        mode = mode_name, dml_set = "ploidy", permutation_scheme = scheme, track = tracks[j], n_permutations = n_labellings,
        odds_ratio_real = real_odds[[j]], odds_ratio_real_cross_product = labelling_odds[real_row, j],
        perm_median = if (finite) median(every(odds_j)) else NA_real_,
        perm_percentile_2.5 = if (finite) quantile(every(odds_j), 0.025, names = FALSE) else NA_real_,
        perm_percentile_97.5 = if (finite) quantile(every(odds_j), 0.975, names = FALSE) else NA_real_,
        share_perm_at_least_real = if (finite) mean(odds_j >= labelling_odds[real_row, j]) else NA_real_,
        pct_eligible_called_in_track = pct_inside, pct_eligible_called_outside_track = pct_outside,
        ratio_in_to_outside = pct_inside / pct_outside,
        matched_n_permutations = n_perm,
        matched_odds_ratio_real = matched_real$ploidy[[tracks[j]]],
        matched_perm_median = median(matched_j), matched_perm_median_mc_low = mc(0.5)[1], matched_perm_median_mc_high = mc(0.5)[2],
        matched_perm_percentile_2.5 = quantile(matched_j, 0.025, names = FALSE),
        matched_perm_percentile_2.5_mc_low = mc(0.025)[1], matched_perm_percentile_2.5_mc_high = mc(0.025)[2],
        matched_perm_percentile_97.5 = quantile(matched_j, 0.975, names = FALSE),
        matched_perm_percentile_97.5_mc_low = mc(0.975)[1], matched_perm_percentile_97.5_mc_high = mc(0.975)[2],
        matched_share_perm_at_least_real = matched_share,
        matched_share_standard_error = sqrt(matched_share * (1 - matched_share) / n_perm),
        notes = paste(c(sprintf("unmatched columns (cross-product ratios): %s", enum_note),
                        if (!finite) sprintf("%g labellings without a finite odds ratio, no summary", labellings_per_row * sum(!is.finite(odds_j))),
                        sprintf("matched columns: %s", draw_note)), collapse = ". "))
    }
    log_msg("[%s] %s: %g labellings, median %g CpGs called (real labels %d), %g labellings at or above the real number",
            mode_name, scheme, n_labellings, median(every(n_called)), n_called_real, labellings_per_row * sum(n_called >= n_called_real))
  }
  # Heat: no permutation.
  heat_odds <- vapply(tracks, function(tr) fisher_odds(table_counts(is_dml$heat, eligible, in_track[, tr]), interval = FALSE),
                      numeric(1))
  heat_row <- data.table(mode = mode_name, dml_set = "heat", permutation_scheme = "none", track = tracks, n_permutations = 0,
                         odds_ratio_real = heat_odds, matched_n_permutations = 0,
                         matched_odds_ratio_real = unname(matched_real$heat[tracks]), notes = no_heat_permutation)
  perm_rows[[length(perm_rows) + 1]] <- heat_row
  counts[[length(counts) + 1]] <- count_row(mode_name, "heat", "none", "CpGs called per permutation", NA_real_,
                                            notes = no_heat_permutation)

  list(odds = rbindlist(odds), sd = rbindlist(sd_rows), perm = rbindlist(perm_rows, use.names = TRUE, fill = TRUE),
       counts = rbindlist(counts))
}

results <- lapply(config$modes, analyse_mode)
odds_tab  <- rbindlist(lapply(results, `[[`, "odds"))
sd_tab    <- rbindlist(lapply(results, `[[`, "sd"))
perm_tab  <- rbindlist(lapply(results, `[[`, "perm"), use.names = TRUE, fill = TRUE)
count_tab <- rbindlist(lapply(results, `[[`, "counts"))

# ---- Write ---------------------------------------------------------------------------------
round_columns <- function(tab, columns, digits) {
  tab[, (columns) := lapply(.SD, round, digits), .SDcols = columns]
}
round_columns(odds_tab, c("odds_ratio", "ci_low", "ci_high"), 4)
round_columns(sd_tab, c("sd_from", "sd_to", "median_sd_dml", "median_sd_other_eligible", "pct_dml_in_gene_body",
                        "pct_other_eligible_in_gene_body"), 2)
round_columns(perm_tab, setdiff(names(perm_tab), c("mode", "dml_set", "permutation_scheme", "track", "n_permutations",
                                                   "matched_n_permutations", "notes")), 4)
count_tab[value != round(value), value := signif(value, 6)]      # whole numbers (counts) are written in full
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
cat("\nLabel permutation, ploidy, comparison 1 (eligible) over every labelling:\n")
print(perm_tab[dml_set == "ploidy", .(mode, permutation_scheme, track, n_permutations, odds_ratio_real, odds_ratio_real_cross_product,
                                      perm_median, perm_percentile_2.5, perm_percentile_97.5, share_perm_at_least_real,
                                      ratio_in_to_outside)], nrows = 200)
cat("\nLabel permutation, ploidy, comparison 3 over the random labellings (Monte Carlo ranges in the file):\n")
print(perm_tab[dml_set == "ploidy", .(mode, permutation_scheme, track, matched_n_permutations, matched_odds_ratio_real,
                                      matched_perm_median, matched_perm_percentile_2.5, matched_perm_percentile_97.5,
                                      matched_share_perm_at_least_real, matched_share_standard_error)], nrows = 200)
cat("\nMedian within-group SD (percentage points):\n"); print(sd_tab[sd_group == "all"])
cat("\nDifference rule, q-value rule, and CpGs called per permutation:\n")
print(count_tab[mean_methylation_bin == "all", .(mode, dml_set, permutation_scheme, quantity = substr(quantity, 1, 80),
                                                 value = sprintf("%.6g", value))], nrows = 300)
log_msg("feature_odds_matched finished in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))

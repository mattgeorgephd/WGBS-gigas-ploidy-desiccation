# 05_differential_methylation.R
#
# Purpose   Differentially methylated loci (DML) for ploidy and for heat, with
#           methylKit, on the 22 libraries, for the analysis of the run
#           (config$analysis; "legacy" is the paper's, run_config.R):
#             1. import each library's coverage file; keep CpGs with at least
#                10 reads (5 with option cov5) (and, with option corrected, with fewer
#                reads than the library's 99.9th percentile of coverage, taken
#                over its CpGs with at least config$import_min_coverage reads)
#             2. keep the CpGs present in all 22 libraries and normalize coverage
#                between libraries (methylKit's median scaling)
#             3. per CpG, logistic regression of methylation on ploidy with heat
#                as covariate, and on heat with ploidy as covariate; SLIM q-values
#                (option mnchisq, not adopted: with methylKit's overdispersion
#                correction, overdispersion = "MN", test = "Chisq")
#             4. DML = q < the cutoff of the analysis (config$analysis$dml_qvalue:
#                config$dml_qvalue, 0.01, in the paper; 0.05 with option q05) and
#                |difference| > the threshold of the analysis
#                (config$analysis$dml_difference: config$dml_difference, 20
#                points, in the paper; 50 with option d50)
#             5. loci in both DML lists ("shared"; the legacy "interaction-DML")
#             With option ctexcl (not adopted) the CpGs of the C-to-T variant list
#             are removed after the test from the tested CpGs and the DML lists,
#             as in Venkataraman et al. 2022: the fit, its p-values and q-values
#             are those of all CpGs (the fit of the analysis without ctexcl is
#             reused), and only rows are dropped.
#             With option d50 (not adopted) the DML lists are taken at a
#             difference of more than 50 points, as methylKit getMethylDiff with
#             difference = 50 takes them (Venkataraman et al. 2022). The threshold
#             is applied after the test, so the fit of the analysis without d50
#             is reused (its fit id); nothing is fitted for it when that fit is in
#             the cache. The counts at 5 to 50 points of
#             dml_counts_by_threshold_<id>.csv are those of that analysis.
#             With option q05 (not adopted) the DML lists are taken at a q-value
#             below 0.05, as getMethylDiff with qvalue = 0.05 takes them. The
#             q-values are those of the fit (SLIM over all tested CpGs), and the
#             cutoff is applied after the test, so the fit of the analysis without
#             q05 is reused as with d50. The counts of
#             dml_counts_by_threshold_<id>.csv are taken at q < 0.05.
#             With option ctrl (not adopted; issue I-07) the ploidy contrast is
#             tested on the 12 control libraries alone (diploid against triploid
#             controls, no covariate), on the CpGs and with the normalisation of
#             the analysis (all 22 libraries); its meth.diff is triploid minus
#             diploid among the controls and its q-values are SLIM over all
#             tested CpGs. The heat contrast is that of the fit of the analysis
#             without ctrl, reused. The shared set is the overlap of this ploidy
#             list with the heat list.
#             With option inter (not adopted; issue I-02) a third contrast is
#             tested, the ploidy by heat interaction (indicator "triploid and
#             heat-treated", covariates ploidy and heat), with its own DML list
#             by the rule of the analysis; its difference is (triploid heat -
#             triploid control) - (diploid heat - diploid control) from pooled
#             percentages. The ploidy and heat contrasts, and the shared set, are
#             those of the analysis without inter (its fit is reused). Both
#             options use the definitions of analyses A4_ploidy_in_control and
#             A2_interaction of 06_diagnostics/dml_sensitivity.R (section
#             "Options ctrl and inter" below).
#
# Run       Rscript 01_pipeline/05_differential_methylation.R     (after step 01)
#           WGBS_CORES sets the number of forked workers of the model fits
#           (default 1; leave it at 1 on Windows). Time and memory are in the
#           Notes.
#           The paper's analysis is fitted in one piece, as it always was (and so
#           are ctexcl, d50 and q05, which reuse the paper's fit). Every other optional
#           analysis is fitted by sequence: the import and each part of
#           the genome (a linkage group, or a group of small sequences) are saved
#           in the cache as they are done, so a run that is stopped loses only
#           the part it was fitting; the same command run again goes on from the
#           saved parts. Results are identical to a fit in one piece (section
#           "Fit by sequence" below; tests/check_fit_by_sequence.R).
#           An analysis with option ctrl or inter takes the fit of the analysis
#           without the token as that analysis takes it (fitting it first when it
#           is not in the cache) and fits only its new contrast, by sequence,
#           resumable in the same way (tests/check_contrast_variants.R).
#
# Inputs    03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files listed in 03_input_files/gannet_manifest.csv (cache)
# Outputs   04_output/05_differential_methylation/tables/ (the paper run; an optional
#           analysis writes to 04_output/not_adopted/<id>/05_differential_methylation/ with
#           <id> in place of <mode>, and rebuilds not_adopted/comparisons/.../dml_counts_by_mode.csv)
#             library_filtering_<mode>.csv       per-library CpG counts, cutoffs, scaling
#             united_cpg_<mode>.csv.gz           normalized counts at every analysed CpG, all libraries
#                                                (median-scaled and rounded as by methylKit; not raw
#                                                read counts)
#             diffmeth_<contrast>_<mode>.csv.gz  p, q, difference at every analysed CpG
#             dml_<contrast>_<mode>.csv          DML at the thresholds of the analysis (q, difference)
#             dml_shared_<mode>.csv              loci in both DML lists, with both differences
#             united_at_dml_<id>.csv             (optional analysis) the rows of united_cpg_<id>.csv.gz at the
#                                                DML of either list, with united_row; tracked, and read by
#                                                step 13, where the per-CpG tables are not tracked
#             dml_counts_by_threshold_<mode>.csv DML counts at difference thresholds of 5 to 50 points
#             dml_counts_by_mode.csv             headline counts of the analysis (column mode: its id); with
#                                                option inter four columns more, at the end: interaction_dml,
#                                                interaction_hyper, interaction_hypo, interaction_in_shared
#                                                (interaction DML among the shared loci)
#             diffmeth_interaction_<id>.csv.gz   (inter) p, q and the interaction difference (column meth.diff;
#                                                not methylKit's) at every analysed CpG, with the pooled
#                                                percentages of the four groups (pct_diploid_control, ...)
#             dml_interaction_<id>.csv           (inter) the interaction DML, by the rule of the analysis
#             ct_exclusion_summary_<id>.csv      (ctexcl) per contrast: CpGs and DML before and after the
#                                                removal, removed, and kept CpGs and DML that the variant
#                                                screen could not assess
#             ct_excluded_cpgs_<id>.csv.gz       (ctexcl) the tested CpGs removed, with p, q and difference
#                                                per contrast and whether each was a DML
#           Of an optional analysis other than 'corrected', united_cpg_<id>.csv.gz and
#           diffmeth_<contrast>_<id>.csv.gz are not tracked by git (.gitignore); rerun this
#           step for the analysis to make them (half a minute with the fit in the cache).
#           With option ctexcl the step stops when the variant list is not the one
#           step 04 used (md5 in the stamp of the analysis folder).
#           04_output/05_differential_methylation/bed/dml_<set>_<mode>.bed   (inter: also dml_interaction_<id>.bed, 4th column
#             the interaction difference)
#           04_output/05_differential_methylation/cache/methylkit_<fit id>.rds   (not tracked; the fit
#             id is the analysis id without ctexcl, d50 and q05, "legacy" for the paper; for a fit id
#             with ctrl or inter the file holds the new contrast only, p-values and methylKit's
#             differences, and names the fit it was made on, methylkit_<base fit id>.rds)
#           04_output/05_differential_methylation/cache/parts_<fit id>/   (not tracked; the saved
#             import and parts of a fit by sequence while it is under way; deleted when the fit is
#             written)
#           04_output/05_differential_methylation/session_info.txt
#
# Manuscript items: Results 3.2 (CpG and DML counts); data for Figures 2, 3, 4.
# Notes     - Steps 1 and 2 follow methylKit's rules (methRead, filterByCoverage,
#             normalizeCoverage, unite) but are applied one library at a time
#             with data.table, so the import needs about 4 GB instead of 9 GB; see
#             02_R_functions/methylkit_helpers.R. The tests are methylKit's own.
#           - Differences are group 1 minus group 0: triploid minus diploid for
#             ploidy, heat minus control for heat. "hyper" means higher in
#             triploid, or higher in heat. A locus in both lists has two
#             differences; where it is counted as hyper or hypo, the sign is that
#             of its heat difference alone: n_hyper, n_hypo and pct_hyper of the
#             set "shared" in dml_counts_by_threshold_<mode>.csv, and the fourth
#             column of bed/dml_shared_<mode>.bed. Its class (both hyper, both
#             hypo, opposite) is in dml_shared_<mode>.csv and dml_counts_by_mode.csv.
#           - The heat contrast compares two cohorts that were sequenced
#             separately (issue I-01). The shared set is an overlap of two
#             lists, not a test of interaction (issue I-02).
#           - The model has no overdispersion term. Its DML lists reproduce the
#             draft; they are not established differences (issue I-06;
#             06_diagnostics/dml_sensitivity.R).
#           - In united_cpg_<mode>.csv.gz the three columns of a library are
#             rounded separately after scaling, as methylKit does, so coverage
#             can differ by 1 from numCs + numTs (19,811 cells of the legacy
#             table). methylKit's tests and the diagnostics use numCs + numTs.
#             The column scaling_factor of library_filtering_<mode>.csv is
#             written with 15 significant digits. To reproduce the normalized
#             counts from raw counts, compute the factor as max(median_coverage)
#             / median_coverage: the written 1.21428571428571 (17/14) scales a
#             coverage of 133 reads to 161 where the factor itself gives 162.
#           - Thresholds are applied as methylKit applies them (getMethylDiff),
#             in both modes: abs(meth.diff) > threshold, on the difference as
#             methylKit computed it, a double that is not rounded. The DML lists
#             and the counts by threshold are taken this way, and the counts of
#             the draft depend on it: the legacy script reports 5,910 heat loci
#             at 15 points, which holds only because two differences of exactly
#             15 points are held as -15.000000000000002. The tables are written
#             with 15 significant digits, so a count retaken from a written
#             table differs at such ties. Three loci of the heat contrast are
#             concerned, the same in both modes:
#               NC_047559.1:17861107   heat 121/352, control 158/320
#                                      written -15, held -15.000000000000002
#               NC_047562.1:49679999   heat 152/280, control 194/280
#                                      written -15, held -15.000000000000002
#               NW_022994958.1:20934   heat 128/256, control 182/260
#                                      written -20, held -19.999999999999996
#             The first two are counted here at 15 points and are not counted
#             from diffmeth_heat_<mode>.csv.gz (5,908 instead of 5,910 in legacy
#             mode, 5,669 instead of 5,671 in corrected mode). The third is not
#             a DML at 20 points here and would be one from the written table
#             under "at least 20" (2,352 instead of 2,351; 2,274 instead of
#             2,273). With ">" at the configured threshold of 20 points the
#             written tables give exactly the DML lists of this step.
#           - One fit (the paper's, or one optional analysis) took 22 minutes at
#             WGBS_CORES=2 on 2026-11-20 (two contrasts of about 10 minutes
#             each), with a peak of about 9.9 GB summed over the R process and
#             its two workers (shared pages counted in each); an estimated 11
#             minutes at WGBS_CORES=4. The import needs about 4 GB. Each library's kept
#             CpGs are dropped as soon as its object is built, and the heat
#             object is built after the ploidy test; neither changes a number.
#           - Fit by sequence, measured on 2026-11-20 (4 cores, 15 GB, nothing
#             else running; memory is the proportional set size summed over the
#             R process and its workers, sampled every second). At 10x on 4
#             cores: 11 parts per contrast, 15 minutes for the paper's analysis
#             and 16 for mnchisq, the import 2 minutes. At 5x (2,297,236 CpGs; 13
#             parts per contrast: the ten linkage groups, 96,732 to 307,802 CpGs,
#             and three groups of other sequences): the import 3.1 to 3.4
#             minutes with a peak of 8.3 GB, in a forked process that gives its
#             memory back when it ends; the ploidy test of NC_047559.1 (270,178
#             CpGs) 11.1 minutes on 1 core (peak 5.0 GB) and 5.7 minutes on 2
#             cores (peak 6.6 GB). The whole fit at 5x is thus about 3.2 hours on
#             1 core and 1.6 hours on 2 cores, plus a few minutes to assemble
#             and write it.
#           - A fit is kept in the cache folder and reused only while what it
#             was made from is unchanged: the fit id, the coverage settings and
#             the model (overdispersion and test),
#             the sample sheet, the coverage files, the methylKit version and
#             the fitting code (fit_fingerprint() below). Otherwise, and when a
#             cache file cannot be read, the analysis is fitted again; nothing has
#             to be deleted by hand. The thresholds of step 4 are applied after
#             the fit, so changing them does not refit (options d50 and q05 reuse
#             the fit of their analysis without them). With the fit in the
#             cache the step takes about half a minute. Not part of the
#             fingerprint: the version of R and of every package other than
#             methylKit (stats::glm.fit runs inside the tests, data.table reads
#             the files). After an upgrade of R or of those packages, delete
#             the cache folder to fit again.
#           - A fit by sequence (every fit id other than "legacy") has the same
#             fingerprint plus its method and the code of the fit by sequence,
#             so a change of that code refits the optional analyses and not the
#             paper's. Its saved parts carry the first 12 characters of
#             fingerprint_key() of that fingerprint (not of the locale) in their
#             names and are reused only under it; saved files of another key are
#             moved to parts_<fit id>/stale/ (only the latest such set is kept).
#             A fit by sequence also uses a cached fit made in one piece by
#             fit_mode() from the same inputs and code (the two are identical).
#             To resume a stopped fit, run the same command again; to start it
#             afresh, delete cache/parts_<fit id>/. One run at a time per fit id:
#             parts_<fit id>/lock names the process that fits it; a second run
#             stops while that process lives, and takes the lock over after it
#             ended. A process forked by a run that was killed (the import, or a
#             worker of a test) can outlive it and keep its memory: before a
#             rerun, list such processes with pgrep -af 05_differential and stop
#             them.
#           - Options ctrl and inter (section "Options ctrl and inter" below):
#             the fit of the analysis without the token is taken by cached_fit()
#             under its own, unchanged fingerprint; only the new contrast is
#             fitted (about half the time of a fit: an estimated 12 minutes at
#             10x and 45 at 5x on WGBS_CORES=2), and its fingerprint holds the
#             key of that fit's fingerprint, so a base fit made again refits it.
#             None of the code of the paper's fit or of the fit by sequence was
#             changed for them, so every cached fit keeps its fingerprint
#             (tests/check_settings.R compares that code with commit 76a900b).
#           - Known limits, not fixed: when a fork fails (the import or the
#             memory watch), there is no fallback to running in the process; the
#             step stops. On Windows the import runs in the R process, which then
#             keeps its memory during the tests. Memory is logged on Linux only
#             (it is read from /proc); on macOS and Windows the log says "not
#             measured". The memory watch, a forked R process, adds about 0.1 GB.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))

step <- "05_differential_methylation"
t0 <- start_step(step, "differential methylation for ploidy and heat")
lib <- read_libraries()

# ---- Fit one mode ------------------------------------------------------------------
# Returns the filtering record of every library, the united table and the two
# tests. What a fit depends on is used inside this function or listed in
# fit_fingerprint() below; cached_fit() decides whether this function runs.
fit_mode <- function(a) {
  mode <- a$fit_id
  cfg_a <- modifyList(config, list(min_coverage = a$min_coverage))   # 10, or 5 with option cov5
  treat_ploidy <- as.integer(lib$ploidy == "triploid")     # 0 diploid, 1 triploid
  treat_heat   <- as.integer(lib$treatment == "heat")      # 0 control, 1 heat

  # Pass over the libraries: filter each one, keep its surviving CpGs as
  # (position key, counts), and narrow the set of CpGs present in all libraries.
  # The md5 of each file is checked as it is read, so that the checksums in the
  # fingerprint of the fit are those of the files the fit was made from.
  kept <- vector("list", nrow(lib)); records <- vector("list", nrow(lib))
  chr_levels <- character(0); common <- NULL
  for (i in seq_len(nrow(lib))) {
    r <- read_filtered_counts(fetch_gannet_file(lib$cov_file[i], config, verify_md5 = TRUE), a$filter, cfg_a, chr_levels)
    chr_levels <- r$chr_levels
    common <- if (is.null(common)) r$counts$pos_key else common[common %in% r$counts$pos_key]
    kept[[i]] <- r$counts
    records[[i]] <- cbind(library_id = lib$library_id[i], r$record)
    log_msg("[%s] %-10s kept %d CpGs (median coverage %g); in all libraries so far: %d", mode,
            lib$library_id[i], r$record$n_cpg_kept, r$record$median_coverage, length(common))
  }
  filtering <- rbindlist(records)
  # methylKit's median normalization: each library's counts are multiplied by
  # (largest median coverage among libraries) / (its own median coverage). The
  # medians are those of each library's full filtered set, as in methylKit.
  filtering[, scaling_factor := max(median_coverage) / median_coverage]

  pos <- key_to_position(common, chr_levels)
  objects <- vector("list", nrow(lib))
  for (i in seq_len(nrow(lib))) {                 # each library's kept CpGs are dropped as soon as its object is built
    k <- kept[[i]][match(common, pos_key)]
    objects[[i]] <- make_methylraw(pos$chr, pos$start, scale_counts_like_methylkit(k$coverage, k$nC, k$nT, filtering$scaling_factor[i]),
                                   as.character(lib$legacy_sample_index[i]))
    kept[i] <- list(NULL)
  }
  rm(kept, records, k); gc(verbose = FALSE)
  raw_list <- new("methylRawList", objects, treatment = treat_ploidy)
  meth_ploidy <- unite(raw_list, destrand = FALSE)        # same CpGs in every object; sets the row order
  rm(objects, raw_list); gc(verbose = FALSE)
  log_msg("[%s] %d CpGs with at least %dx in all %d libraries", mode, nrow(meth_ploidy),
          cfg_a$min_coverage, nrow(lib))

  # The paper's model is the first call of each pair, with exactly its arguments
  # of the draft; option mnchisq (not adopted) adds the overdispersion correction
  # with the chi-squared test, the only test offered with it.
  log_msg("[%s] fitting ploidy (covariate: heat) on %d core(s)%s", mode, config$cores,
          if (a$overdispersion == "MN") "; overdispersion MN, chi-squared test" else "")
  dm_ploidy <- if (a$overdispersion == "none") {
    calculateDiffMeth(meth_ploidy, covariates = data.frame(heat = lib$treatment),
                      mc.cores = config$cores)
  } else {
    calculateDiffMeth(meth_ploidy, covariates = data.frame(heat = lib$treatment),
                      overdispersion = "MN", test = "Chisq", mc.cores = config$cores)
  }
  meth_heat <- reorganize(meth_ploidy, sample.ids = meth_ploidy@sample.ids, treatment = treat_heat)
  log_msg("[%s] fitting heat (covariate: ploidy)", mode)
  dm_heat <- if (a$overdispersion == "none") {
    calculateDiffMeth(meth_heat, covariates = data.frame(ploidy = lib$ploidy),
                      mc.cores = config$cores)
  } else {
    calculateDiffMeth(meth_heat, covariates = data.frame(ploidy = lib$ploidy),
                      overdispersion = "MN", test = "Chisq", mc.cores = config$cores)
  }
  list(mode = mode, filtering = filtering, united = meth_ploidy,
       diffmeth = list(ploidy = dm_ploidy, heat = dm_heat))
}

# ---- Fit by sequence (every fit id other than "legacy") ---------------------------
# The paper's fit is made in one piece by fit_mode() above, which is left exactly
# as it was; ctexcl reuses it. Every other optional analysis is fitted one part of the genome at a time, and
# each part is saved in the cache as soon as it is done, so that a run stopped
# partway loses only the part it was fitting (the fit at 5x needs more than the 2
# hours of a tracked task on one core). The results are those of a fit in one
# piece, to the last bit:
#   - methylKit's test is made CpG by CpG. For more than two libraries
#     calculateDiffMeth() calls logReg() on each row of the united table: the
#     logistic regression of that CpG's 22 libraries, and with overdispersion =
#     "MN" the overdispersion factor of that CpG, from its own residuals. Nothing
#     is shared between CpGs (only overdispersion = "shrinkMN", not used here,
#     pools them). A p-value that comes out NaN is set to 1 inside logReg(), CpG
#     by CpG. A part is the same call of calculateDiffMeth() on a subset of the
#     rows, so every CpG gets the same p-value and difference.
#   - What methylKit computes over all CpGs is computed over all CpGs. The import,
#     the coverage filters, the median scaling of each library (normalizeCoverage:
#     the medians of each library's filtered CpGs) and the united table are made
#     once, genome-wide, before the table is split (import_united(), a copy of the
#     first half of fit_mode()). The SLIM q-values are computed once, after all
#     parts are in, from the p-values of all CpGs in the order of the united
#     table, with methylKit's own function and arguments (p.adjusted(method =
#     "SLIM"), as calculateDiffMeth() calls it); the q-values that
#     calculateDiffMeth() returns for a part are not used.
#   - The parts are assembled into the methylDiff objects that calculateDiffMeth()
#     builds, row names and slots included. tests/check_fit_by_sequence.R checks
#     with identical() that the result equals the fit in one piece (the paper's
#     analysis and mnchisq at 10x).
# Parts: each of the ten linkage groups (config$linkage_groups) is a part; the
# other sequences (the mitochondrial genome and the unplaced scaffolds) are taken
# in the order of the united table and grouped into parts of at most 50,000 CpGs
# (a sequence is never split). Parts are fitted in the order of their first row
# in the united table, all parts of the ploidy test, then all parts of the heat
# test. Saved in <cache folder of the step>/parts_<fit id>/ (git-ignored), under
# names that begin with the first 12 characters of fingerprint_key() of the fit's
# fingerprint (fit_fingerprint() below):
#   <fp>_united.rds                filtering record and united table (the import)
#   <fp>_<contrast>_<part>.rds     p-values and differences of one part
#   lock                           the process that is fitting (see the Notes)
#   stale/                         the latest set of saved files of another key
# Each file is written under a temporary name (<file>.<pid>.part) and then
# renamed. A rerun reuses every file of the same key that can be read (a broken
# import file is imported again) and fits the rest. When all parts are in,
# cached_fit() writes the assembled fit to methylkit_<fit id>.rds and deletes the
# folder.

# The filtering record and the united table, as fit_mode() makes them. This is a
# copy of the first half of fit_mode(), kept separate so that the paper's
# fit_mode() and its cached fit stay as they were; tests/check_fit_by_sequence.R
# checks that the two give identical tables.
import_united <- function(a) {
  mode <- a$fit_id
  cfg_a <- modifyList(config, list(min_coverage = a$min_coverage))   # 10, or 5 with option cov5
  treat_ploidy <- as.integer(lib$ploidy == "triploid")     # 0 diploid, 1 triploid
  kept <- vector("list", nrow(lib)); records <- vector("list", nrow(lib))
  chr_levels <- character(0); common <- NULL
  for (i in seq_len(nrow(lib))) {
    r <- read_filtered_counts(fetch_gannet_file(lib$cov_file[i], config, verify_md5 = TRUE), a$filter, cfg_a, chr_levels)
    chr_levels <- r$chr_levels
    common <- if (is.null(common)) r$counts$pos_key else common[common %in% r$counts$pos_key]
    kept[[i]] <- r$counts
    records[[i]] <- cbind(library_id = lib$library_id[i], r$record)
    log_msg("[%s] %-10s kept %d CpGs (median coverage %g); in all libraries so far: %d", mode,
            lib$library_id[i], r$record$n_cpg_kept, r$record$median_coverage, length(common))
  }
  filtering <- rbindlist(records)
  filtering[, scaling_factor := max(median_coverage) / median_coverage]
  pos <- key_to_position(common, chr_levels)
  objects <- vector("list", nrow(lib))
  for (i in seq_len(nrow(lib))) {
    k <- kept[[i]][match(common, pos_key)]
    objects[[i]] <- make_methylraw(pos$chr, pos$start, scale_counts_like_methylkit(k$coverage, k$nC, k$nT, filtering$scaling_factor[i]),
                                   as.character(lib$legacy_sample_index[i]))
    kept[i] <- list(NULL)
  }
  rm(kept, records, k); gc(verbose = FALSE)
  raw_list <- new("methylRawList", objects, treatment = treat_ploidy)
  meth_ploidy <- unite(raw_list, destrand = FALSE)
  rm(objects, raw_list); gc(verbose = FALSE)
  log_msg("[%s] %d CpGs with at least %dx in all %d libraries", mode, nrow(meth_ploidy),
          cfg_a$min_coverage, nrow(lib))
  list(filtering = filtering, united = meth_ploidy)
}

# The parts of the united table: a list of row numbers per part, named by part,
# in the order of each part's first row.
sequence_parts <- function(chr, max_other = 50000L) {
  seqs <- unique(chr)
  n <- tabulate(match(chr, seqs), nbins = length(seqs))
  part <- seqs
  own <- seqs %in% config$linkage_groups
  group <- integer(length(seqs)); g <- 0L; size <- Inf
  for (i in which(!own)) {
    if (size + n[i] > max_other) { g <- g + 1L; size <- 0 }
    group[i] <- g; size <- size + n[i]
  }
  for (j in seq_len(g)) {
    s <- seqs[!own & group == j]
    part[!own & group == j] <- if (length(s) == 1) sprintf("other%02d_%s", j, s)
                               else sprintf("other%02d_%s_to_%s", j, s[1], s[length(s)])
  }
  part <- gsub("[^A-Za-z0-9._-]", "_", part)
  row_part <- part[match(chr, seqs)]
  rows <- split(seq_along(chr), factor(row_part, levels = unique(row_part)))
  stopifnot(sum(lengths(rows)) == length(chr), !anyDuplicated(names(rows)))
  rows
}

# The test of one contrast on the rows of one part: the calls of fit_mode(), on a
# methylBase of those rows. Returns the p-values and differences.
test_part <- function(obj, contrast, a) {
  covariates <- if (contrast == "ploidy") data.frame(heat = lib$treatment) else data.frame(ploidy = lib$ploidy)
  dm <- if (a$overdispersion == "none") {
    calculateDiffMeth(obj, covariates = covariates, mc.cores = config$cores)
  } else {
    calculateDiffMeth(obj, covariates = covariates, overdispersion = "MN", test = "Chisq", mc.cores = config$cores)
  }
  d <- getData(dm)
  list(pvalue = d$pvalue, meth.diff = d$meth.diff)
}

# The methylDiff object that calculateDiffMeth() returns for `obj` (the whole
# united table, with the treatment of the contrast), from the p-values and
# differences of all its rows: the same data frame (row names, columns and their
# names), the q-values of methylKit's p.adjusted() with its default method SLIM
# over all p-values, and the same slots.
diffmeth_from_parts <- function(obj, pvalue, meth.diff) {
  stopifnot(length(pvalue) == nrow(obj), length(meth.diff) == nrow(obj), is.double(pvalue), !anyNA(pvalue))
  p_adjusted <- utils::getFromNamespace("p.adjusted", "methylKit")
  subst <- S3Part(obj, strictS3 = TRUE)
  x <- data.frame(subst[, 1:4], pvalue, p_adjusted(pvalue, method = "SLIM"), meth.diff = meth.diff,
                  stringsAsFactors = FALSE)
  colnames(x)[5:7] <- c("pvalue", "qvalue", "meth.diff")
  new("methylDiff", x, sample.ids = obj@sample.ids, assembly = obj@assembly, context = obj@context,
      destranded = obj@destranded, treatment = obj@treatment, resolution = obj@resolution)
}

# saveRDS() under a temporary name that carries the process id (<file>.<pid>.part),
# then rename: a file with the final name is complete, and two processes never
# write the same temporary file.
save_rds_atomic <- function(x, file) {
  tmp <- sprintf("%s.%d.part", file, Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, file)) stop("Could not move ", tmp, " to ", file, ".", call. = FALSE)
  invisible(file)
}

parts_dir <- function(a) file.path(step_cache_dir(step), paste0("parts_", a$fit_id))

# The key of the saved parts: md5 of the fingerprint serialized in format 2
# without its 14-byte header (the header names the version of R that wrote it).
# Format 3, which digest_md5() below uses for the stamp, records the native
# encoding of the session, so its md5 changes with the locale (LC_ALL=C and
# C.UTF-8 differ); this key does not. The fingerprint holds ASCII text only.
fingerprint_key <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  writeBin(serialize(x, NULL, version = 2)[-(1:14)], f)
  unname(tools::md5sum(f))
}

# ---- The lock of a fit by sequence ------------------------------------------------
# parts_<fit id>/lock holds one line: the process id, its start time in clock
# ticks since boot (field 22 of /proc/<pid>/stat; "-" where there is no /proc),
# and the time the lock was taken (UTC). A run that finds the lock of a live
# process stops; the lock of a process that has ended is taken over. A process
# id that was given to a new process is told apart by its start time (Linux).
process_start <- function(pid) {
  s <- tryCatch(suppressWarnings(readLines(sprintf("/proc/%d/stat", as.integer(pid)), warn = FALSE)), error = function(e) character(0))
  if (length(s) == 0) return(NA_character_)
  strsplit(sub("^.*\\) ", "", s[1]), " ", fixed = TRUE)[[1]][20]   # field 22 (fields from the state on, which is field 3)
}
read_lock <- function(lock) {
  l <- tryCatch(suppressWarnings(strsplit(readLines(lock, n = 1, warn = FALSE), " ", fixed = TRUE)[[1]]), error = function(e) character(0))
  list(pid = suppressWarnings(as.integer(l[1])), start = if (length(l) >= 2) l[2] else NA_character_,
       taken = if (length(l) >= 3) paste(l[-(1:2)], collapse = " ") else "unknown")
}
lock_is_live <- function(lock) {
  if (!file.exists(lock)) return(FALSE)
  l <- read_lock(lock)
  if (is.na(l$pid) || l$pid == Sys.getpid()) return(FALSE)
  if (file.exists("/proc/self/stat")) return(identical(process_start(l$pid), l$start))
  isTRUE(tools::pskill(l$pid, 0L))                               # signal 0: does the process exist
}
take_parts_lock <- function(a) {
  dir <- parts_dir(a)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(dir, "lock")
  if (lock_is_live(lock)) {
    l <- read_lock(lock)
    stop("Another run of step 05 is fitting the analysis '", a$fit_id, "' (process ", l$pid, ", lock taken ", l$taken,
         ", in ", lock, "). Let it finish, or stop it and its forked workers (list them with: pgrep -af 05_differential) ",
         "and run the step again; the saved parts are kept.", call. = FALSE)
  }
  if (file.exists(lock)) {
    l <- read_lock(lock)
    log_msg("[%s] taking over the lock of process %s (lock taken %s), which is no longer running", a$fit_id, l$pid, l$taken)
  }
  start <- process_start(Sys.getpid())
  tmp <- sprintf("%s.%d.part", lock, Sys.getpid())
  writeLines(paste(Sys.getpid(), if (is.na(start)) "-" else start, format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")), tmp)
  if (!file.rename(tmp, lock)) stop("Could not write the lock ", lock, ".", call. = FALSE)
  if (!identical(read_lock(lock)$pid, Sys.getpid())) {
    stop("Another run of step 05 took the lock ", lock, " at the same moment; run the step again later.", call. = FALSE)
  }
  lock
}
drop_parts_lock <- function(lock) {
  if (file.exists(lock) && identical(read_lock(lock)$pid, Sys.getpid())) unlink(lock)
  invisible(NULL)
}

# What a run leaves behind when it is stopped, removed once a valid fit of the
# analysis is in the cache: the folder parts_<fit id>/ (unless a live run holds
# its lock) and temporary files <cache file>.<pid>.part of processes that have ended.
remove_leftovers <- function(a, cache_file) {
  dir <- parts_dir(a)
  if (dir.exists(dir)) {
    if (lock_is_live(file.path(dir, "lock"))) {
      log_msg("[%s] %s is locked by a running process and is left as it is", a$fit_id, basename(dir))
    } else {
      unlink(dir, recursive = TRUE)
      log_msg("[%s] the fit is in the cache; leftover %s/ removed", a$fit_id, basename(dir))
    }
  }
  tmp <- list.files(dirname(cache_file), full.names = TRUE)
  tmp <- tmp[startsWith(basename(tmp), paste0(basename(cache_file), ".")) & endsWith(tmp, ".part")]
  pid <- suppressWarnings(as.integer(sub("^.*\\.([0-9]+)\\.part$", "\\1", basename(tmp))))
  dead <- is.na(pid) | !vapply(pid, function(p) !is.na(p) && isTRUE(tools::pskill(p, 0L)), logical(1))
  if (any(dead)) {
    unlink(tmp[dead])
    log_msg("[%s] %d leftover temporary file(s) of %s removed", a$fit_id, sum(dead), basename(cache_file))
  }
  invisible(NULL)
}

# The fit of an optional analysis, part by part, reusing the parts saved by an
# earlier run with the same fingerprint. Returns what fit_mode() returns.
fit_by_sequence <- function(a, fingerprint) {
  mode <- a$fit_id
  fp <- fingerprint_key(fingerprint)
  prefix <- paste0(substr(fp, 1, 12), "_")
  dir <- parts_dir(a)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # The caller holds the lock, so a temporary file here is left by a run that
  # was stopped. Saved files of another fingerprint are moved to stale/, which
  # keeps only the latest such set (they can be inspected, or moved back by hand).
  entries <- setdiff(list.files(dir), c("lock", "stale"))
  leftover <- entries[endsWith(entries, ".part")]
  if (length(leftover) > 0) unlink(file.path(dir, leftover))
  stale <- setdiff(entries, leftover)
  stale <- stale[!startsWith(stale, prefix)]
  if (length(stale) > 0) {
    unlink(file.path(dir, "stale"), recursive = TRUE)
    dir.create(file.path(dir, "stale"))
    file.rename(file.path(dir, stale), file.path(dir, "stale", stale))
    log_msg("[%s] %d saved file(s) of another fingerprint moved to %s/stale/ (an older stale set there deleted)",
            mode, length(stale), basename(dir))
  }
  part_file <- function(name) file.path(dir, paste0(prefix, name, ".rds"))
  saved <- function(file, ok) {
    if (!file.exists(file)) return(NULL)
    x <- tryCatch(suppressWarnings(readRDS(file)), error = function(e) NULL)
    if (is.list(x) && identical(x$fingerprint_md5, fp) && isTRUE(ok(x))) x else NULL
  }
  gb <- function(x) if (is.na(x)) "not measured" else sprintf("%.2f GB", x)

  # The import and the united table, genome-wide.
  file <- part_file("united")
  imp <- saved(file, function(x) is(x$united, "methylBase") && is.data.frame(x$filtering))
  if (is.null(imp)) {
    # Where R can fork, the import runs in a forked process that saves it and
    # exits, and this process reads the saved file: the memory of the import
    # (about 8 GB at 5x) is then given back before the tests start, instead of
    # staying in this process, where every worker of the tests would inherit it.
    t1 <- Sys.time(); watch <- memory_watch_start()
    make_import <- function() { save_rds_atomic(c(import_united(a), list(fingerprint_md5 = fp)), file); TRUE }
    if (.Platform$OS.type == "unix") {
      done <- parallel::mccollect(parallel::mcparallel(make_import()))[[1]]
      if (!isTRUE(done)) stop("The import of the analysis '", mode, "' failed: ",
                              if (inherits(done, "try-error")) conditionMessage(attr(done, "condition")) else "no result (out of memory?)",
                              call. = FALSE)
    } else {
      make_import()
    }
    peak <- memory_watch_stop(watch)
    imp <- tryCatch(readRDS(file), error = function(e) NULL)
    if (!(is.list(imp) && identical(imp$fingerprint_md5, fp) && is(imp$united, "methylBase"))) {
      unlink(file)
      stop("The import just saved as ", file, " cannot be read back; it was deleted, and the next run imports again.", call. = FALSE)
    }
    log_msg("[%s] import saved as %s/%s: %.1f min, peak memory %s", mode, basename(dir), basename(file),
            as.numeric(difftime(Sys.time(), t1, units = "mins")), gb(peak))
  } else {
    log_msg("[%s] reusing the saved import %s/%s: %d CpGs", mode, basename(dir), basename(file), nrow(imp$united))
  }
  united <- imp$united
  rows <- sequence_parts(united$chr)
  log_msg("[%s] fit by sequence: %d parts per contrast (%s), %d CpGs, %d core(s)%s", mode, length(rows),
          paste(sprintf("%s %d", names(rows), lengths(rows)), collapse = ", "), nrow(united), config$cores,
          if (a$overdispersion == "MN") "; overdispersion MN, chi-squared test" else "")

  treat_heat <- as.integer(lib$treatment == "heat")        # 0 control, 1 heat
  results <- list()
  for (contrast in c("ploidy", "heat")) {
    pvalue <- rep(NA_real_, nrow(united)); meth.diff <- pvalue
    for (k in seq_along(rows)) {
      name <- names(rows)[k]; r <- rows[[k]]
      file <- part_file(paste(contrast, name, sep = "_"))
      x <- saved(file, function(x) identical(x$contrast, contrast) && identical(x$part, name) && identical(x$rows, r) &&
                   is.double(x$pvalue) && length(x$pvalue) == length(r) && length(x$meth.diff) == length(r))
      if (is.null(x)) {
        t1 <- Sys.time(); watch <- memory_watch_start()
        obj <- methylKit::select(united, r)
        if (contrast == "heat") obj <- methylKit::reorganize(obj, sample.ids = obj@sample.ids, treatment = treat_heat)
        x <- suppressMessages(test_part(obj, contrast, a))
        rm(obj)
        peak <- memory_watch_stop(watch)
        # A worker that failed (killed for lack of memory, for example) leaves no p-value; such a part is not saved.
        if (!is.double(x$pvalue) || length(x$pvalue) != length(r) || anyNA(x$pvalue) || length(x$meth.diff) != length(r)) {
          stop("The ", contrast, " test of part ", name, " did not return a p-value for every CpG (a worker may have failed); ",
               "run the step again (the saved parts are kept), with fewer cores if memory is short.", call. = FALSE)
        }
        minutes <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
        x <- c(x, list(fingerprint_md5 = fp, contrast = contrast, part = name, rows = r,
                       minutes = minutes, peak_gb = peak, cores = config$cores))
        save_rds_atomic(x, file)
        log_msg("[%s] %s part %d of %d, %s: %d CpGs fitted in %.1f min on %d core(s), peak memory %s; saved",
                mode, contrast, k, length(rows), name, length(r), minutes, config$cores, gb(peak))
      } else {
        log_msg("[%s] %s part %d of %d, %s: %d CpGs, saved part reused (fitted in %.1f min on %d core(s), peak memory %s)",
                mode, contrast, k, length(rows), name, length(r), x$minutes, x$cores, gb(x$peak_gb))
      }
      pvalue[r] <- x$pvalue; meth.diff[r] <- x$meth.diff
    }
    results[[contrast]] <- list(pvalue = pvalue, meth.diff = meth.diff)
  }
  rm(x); invisible(gc(verbose = FALSE))

  # q-values over all CpGs, once, and the methylDiff objects of a fit in one piece.
  log_msg("[%s] all parts in; SLIM q-values over all %d CpGs", mode, nrow(united))
  dm_ploidy <- diffmeth_from_parts(united, results$ploidy$pvalue, results$ploidy$meth.diff)
  meth_heat <- reorganize(united, sample.ids = united@sample.ids, treatment = treat_heat)
  dm_heat <- diffmeth_from_parts(meth_heat, results$heat$pvalue, results$heat$meth.diff)
  rm(meth_heat)
  list(mode = mode, filtering = imp$filtering, united = united,
       diffmeth = list(ploidy = dm_ploidy, heat = dm_heat))
}

# md5 of a fingerprint, for the stamp of the analysis folder (a record only; it
# depends on the locale, see fingerprint_key()).
digest_md5 <- function(x) { f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, version = 3, compress = FALSE); unname(tools::md5sum(f)) }

# ---- What a fit depends on -----------------------------------------------------------
# A cached fit is reused only while all of this is unchanged: the mode and the
# coverage settings, the columns of the sample sheet that the fit uses, the
# coverage files (size and md5 as the manifest lists them; fit_mode() checks
# the files against it), the methylKit version, and the code that builds the
# table of counts and runs the tests. The code is compared as R parses it, so a
# change of comments or layout does not count and any other change does, a
# reworded log message included. The version of R and of packages other than
# methylKit is not compared (see Notes). The thresholds that define a DML are applied
# after the fit (write_mode) and are not part of it; the number of cores is not
# either, because the results do not depend on it.
# A fit by sequence (by_sequence = TRUE, every fit id other than "legacy") adds two parts:
# the method and the code of the fit by sequence, so that a change of that code
# refits the optional analyses and leaves the paper's fingerprint, and its cached
# fit, as they were. The size of the groups of small sequences is in it too (the
# default max_other of sequence_parts()), so changing it refits the optional
# analyses although a fit does not depend on how its CpGs are grouped. Left out
# because they touch no result: the memory watch, the lock and the writing of files.
fit_fingerprint <- function(a, by_sequence = FALSE) {
  cfg_a <- modifyList(config, list(min_coverage = a$min_coverage))
  manifest <- read_gannet_manifest()
  files <- manifest[match(lib$cov_file, manifest$file), ]
  # Every helper that fit_mode() calls, directly or through another helper, must be listed here: an edit of
  # a helper that is not listed does not count as a change of the fitting code.
  fitting_code <- list(fit_mode = fit_mode, read_filtered_counts = read_filtered_counts, read_cov = read_cov,
                       first_line_lost = first_line_lost, upper_coverage_cutoff = upper_coverage_cutoff,
                       cpg_key = cpg_key, key_to_position = key_to_position,
                       scale_counts_like_methylkit = scale_counts_like_methylkit, make_methylraw = make_methylraw)
  fp <- list(mode = a$fit_id,
             coverage_settings = unlist(cfg_a[c("min_coverage", "import_min_coverage", "upper_coverage_percentile")]),
             sample_sheet = paste(lib$library_id, lib$legacy_sample_index, lib$ploidy, lib$treatment, lib$cov_file),
             coverage_files = paste(files$file, files$bytes, files$md5),
             methylKit_version = as.character(utils::packageVersion("methylKit")),
             fitting_code = lapply(fitting_code, function(f) deparse(utils::removeSource(f))),
             model = c(overdispersion = a$overdispersion, test = a$test))
  if (by_sequence) {
    # Every function of the fit by sequence that can change a result (its helpers of the import are listed above).
    by_sequence_code <- list(fit_by_sequence = fit_by_sequence, import_united = import_united, sequence_parts = sequence_parts,
                             test_part = test_part, diffmeth_from_parts = diffmeth_from_parts)
    fp$fit_method <- "by sequence: parts of the united table tested separately, q-values over all CpGs"
    fp$by_sequence_code <- lapply(by_sequence_code, function(f) deparse(utils::removeSource(f)))
  }
  fp
}

# ---- The fit of one mode, from the cache where that is safe ------------------------
# A cache file that cannot be read (a run was stopped while it was written, or
# the file is not a fit of this step) is fitted again, as is a fit made from
# anything other than what fit_fingerprint() returns now. A new fit is written
# under a temporary name and renamed, so a file with the final name is complete.
# The fit id "legacy" (the paper's analysis, and ctexcl, d50 and q05, which reuse it) is
# fitted in one piece (fit_mode()), so methylkit_legacy.rds is always the
# paper's fit; every other fit id by sequence (fit_by_sequence(), under the lock
# of parts_<fit id>/; its saved parts are deleted once the assembled fit is
# written). A fit by sequence also takes a cached fit made in one piece from the
# same inputs and code (the fingerprint of fit_mode()), since the two give
# identical results; the fits in one piece made before the fit by sequence
# existed (mnchisq, corrected, corrected-mnchisq) stay valid that way. When a
# valid fit is found, what a stopped run left (remove_leftovers()) is removed.
# fit_whole and fit_parts are the fitting functions; tests/check_cached_fit.R
# passes stand-ins, the step never does.
cached_fit <- function(a, fit_whole = fit_mode, fit_parts = fit_by_sequence) {
  mode <- a$fit_id
  cache_file <- file.path(step_cache_dir(step), sprintf("methylkit_%s.rds", mode))
  by_sequence <- a$fit_id != "legacy"
  fingerprint <- fit_fingerprint(a, by_sequence)
  if (file.exists(cache_file)) {
    cached <- tryCatch(suppressWarnings(readRDS(cache_file)), error = function(e) NULL)
    readable <- is.list(cached) && all(c("mode", "filtering", "united", "diffmeth") %in% names(cached))
    differs <- function(fp) c(names(fp)[!vapply(names(fp), function(part) identical(cached$fingerprint[[part]], fp[[part]]), logical(1))],
                              setdiff(names(cached$fingerprint), names(fp)))
    changed <- if (readable) differs(fingerprint)
    if (readable && length(changed) == 0) {
      log_msg("[%s] using the cached fit %s: settings, sample sheet, coverage files, methylKit version and fitting code are unchanged",
              mode, basename(cache_file))
      remove_leftovers(a, cache_file)
      return(cached)
    }
    if (readable && by_sequence && length(differs(fit_fingerprint(a, FALSE))) == 0) {
      log_msg("[%s] using the cached fit %s, made in one piece: settings, sample sheet, coverage files, methylKit version and fitting code are unchanged, and a fit by sequence gives identical results",
              mode, basename(cache_file))
      remove_leftovers(a, cache_file)
      return(cached)
    }
    log_msg("[%s] the cached fit %s %s; fitting again", mode, basename(cache_file),
            if (!readable) "cannot be read (cut short, or not a fit of this step)"
            else if (is.null(cached$fingerprint)) "does not say what it was made from (written before this step recorded it)"
            else paste("differs from this run in:", paste(gsub("_", " ", changed), collapse = ", ")))
  }
  if (by_sequence) {
    lock <- take_parts_lock(a)
    on.exit(drop_parts_lock(lock), add = TRUE)
  }
  fit <- c(if (by_sequence) fit_parts(a, fingerprint) else fit_whole(a), list(fingerprint = fingerprint))
  save_rds_atomic(fit, cache_file)
  if (by_sequence) {
    unlink(parts_dir(a), recursive = TRUE)
    log_msg("[%s] fit written to %s; its saved parts deleted", mode, basename(cache_file))
  }
  fit
}

# ---- Options ctrl and inter: a variant of the contrasts -------------------------------
# An analysis with token ctrl or inter (config$analysis$contrast_variant) has the
# fit of the analysis without the token (its base_fit_id: the paper's fit for
# ctrl and inter, the fit of mnchisq for mnchisq-ctrl, ...), taken with
# cached_fit() above exactly as that analysis takes it, and one contrast more:
#   ctrl   ploidy among the 12 control libraries alone: triploid (1) against
#          diploid (0), no covariate (all are controls). It replaces the ploidy
#          contrast of the base fit; the heat contrast is that of the base fit.
#          meth.diff is methylKit's: triploid minus diploid among the controls.
#   inter  the indicator "triploid and heat-treated" tested with ploidy and heat
#          as covariates, all 22 libraries: with both main effects in the null
#          model the full model has one parameter per group, so the test is that
#          of the ploidy by heat interaction (checked in variant_design()). Its
#          meth.diff is NOT methylKit's for this coding (triploid heat minus all
#          other libraries) but (triploid heat - triploid control) - (diploid
#          heat - diploid control) from the pooled, coverage-weighted percentages
#          (methylated over methylated plus unmethylated reads of a group),
#          computed in interaction_groups(); positive ("hyper") means that heat
#          minus control is higher in triploids. It is added as a third contrast,
#          "interaction"; the ploidy and heat contrasts are those of the base fit.
# These are the definitions of analyses A4_ploidy_in_control and A2_interaction
# of 06_diagnostics/dml_sensitivity.R (and, with mnchisq, of their _mn_chisq rows),
# so that their results can be compared. The CpGs are those of the united table
# of the base fit (coverage filter over all 22 libraries) with its normalisation;
# the new contrast is tested on those rows (ctrl: the 12 control columns of
# them), and its SLIM q-values are computed once over all of them.
# The new contrast is fitted by sequence, as fit_by_sequence() fits a contrast
# (sequence_parts(), the same call of calculateDiffMeth() on each part, q-values
# over all CpGs by diffmeth_from_parts()), with the methylKit settings of the
# analysis (overdispersion "MN" and the chi-squared test with mnchisq). Its parts
# are saved in parts_<fit id>/ under the key of its fingerprint, and the result
# (p-values and methylKit's differences of the new contrast only, not a copy of
# the base fit) in methylkit_<fit id>.rds. Its fingerprint (variant_fingerprint())
# holds the key of the fingerprint of the base fit as cached_fit() returned it,
# so a base fit made again (other inputs or fitting code) refits the new contrast
# too; it also holds the library sheet, the design, the model, the methylKit
# version and the code below. Nothing of the paper's fit or of any other fit,
# and none of their fingerprints, changes.

# The libraries, tested indicator and covariates of the new contrast.
variant_design <- function(a) {
  if (a$contrast_variant == "ctrl") {
    keep <- lib$treatment == "control"
    d <- list(contrast = "ploidy", keep = keep, treatment = as.integer(lib$ploidy[keep] == "triploid"), covariates = NULL,
              described = "ploidy among the control libraries: triploid (1) against diploid (0), no covariate")
  } else if (a$contrast_variant == "inter") {
    keep <- rep(TRUE, nrow(lib))
    d <- list(contrast = "interaction", keep = keep, treatment = as.integer(lib$ploidy == "triploid" & lib$treatment == "heat"),
              covariates = data.frame(ploidy = lib$ploidy, heat = lib$treatment),
              described = "triploid and heat-treated (1) against the other libraries (0), covariates ploidy and heat")
  } else {
    stop("variant_design(): the analysis '", a$id, "' has no variant of the contrasts.", call. = FALSE)
  }
  if (sum(d$treatment == 1) < 2 || sum(d$treatment == 0) < 2) stop("Fewer than two libraries in a group of ", d$described, ".", call. = FALSE)
  frame <- data.frame(tested = d$treatment)
  if (!is.null(d$covariates)) frame <- cbind(frame, d$covariates)
  design <- stats::model.matrix(~ ., frame)
  if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) stop("The design of ", d$described, " is not of full rank.", call. = FALSE)
  if (d$contrast == "interaction") {
    # one parameter per group: the coded model spans the model ploidy * heat, so its test is the interaction test
    factorial <- stats::model.matrix(~ ploidy * heat, d$covariates)
    if (ncol(design) != ncol(factorial) || qr(cbind(design, factorial))$rank != ncol(factorial)) {
      stop("The coding of the interaction does not span the model ploidy * heat.", call. = FALSE)
    }
  }
  d
}
# The united table (or a part of it) with the libraries and the tested indicator of the new contrast.
variant_object <- function(obj, d) {
  methylKit::reorganize(obj, sample.ids = obj@sample.ids[d$keep], treatment = d$treatment)
}
# The test of the new contrast on the rows of one part, as test_part() tests the paper's contrasts.
variant_test_part <- function(obj, d, a) {
  dm <- if (is.null(d$covariates)) {
    if (a$overdispersion == "none") calculateDiffMeth(obj, mc.cores = config$cores)
    else calculateDiffMeth(obj, overdispersion = "MN", test = "Chisq", mc.cores = config$cores)
  } else {
    if (a$overdispersion == "none") calculateDiffMeth(obj, covariates = d$covariates, mc.cores = config$cores)
    else calculateDiffMeth(obj, covariates = d$covariates, overdispersion = "MN", test = "Chisq", mc.cores = config$cores)
  }
  x <- getData(dm)
  list(pvalue = x$pvalue, meth.diff = x$meth.diff)
}
# The pooled, coverage-weighted percent methylation of the four groups at every
# CpG of a united table, and the interaction difference (as interaction_contrast()
# of 06_diagnostics/dml_sensitivity.R computes it).
interaction_groups <- function(united) {
  x <- as.data.table(getData(united))
  pct <- function(ploidy, treatment) {
    j <- which(lib$ploidy == ploidy & lib$treatment == treatment)
    m <- rowSums(as.matrix(x[, paste0("numCs", j), with = FALSE]))
    100 * m / (m + rowSums(as.matrix(x[, paste0("numTs", j), with = FALSE])))
  }
  g <- data.table(pct_diploid_control = pct("diploid", "control"), pct_diploid_heat = pct("diploid", "heat"),
                  pct_triploid_control = pct("triploid", "control"), pct_triploid_heat = pct("triploid", "heat"))
  g[, interaction_diff := (pct_triploid_heat - pct_triploid_control) - (pct_diploid_heat - pct_diploid_control)]
  if (anyNA(g$interaction_diff)) stop("A group without reads at a tested CpG: the interaction difference cannot be computed.", call. = FALSE)
  g
}

# The new contrast, part by part, on the united table of the base fit, reusing the
# parts saved by an earlier run with the same fingerprint. Returns its p-values and
# methylKit's differences for every CpG, in the order of the united table.
fit_variant_by_sequence <- function(a, base, fingerprint) {
  mode <- a$fit_id
  fp <- fingerprint_key(fingerprint)
  prefix <- paste0(substr(fp, 1, 12), "_")
  dir <- parts_dir(a)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  entries <- setdiff(list.files(dir), c("lock", "stale"))
  leftover <- entries[endsWith(entries, ".part")]
  if (length(leftover) > 0) unlink(file.path(dir, leftover))
  stale <- setdiff(entries, leftover)
  stale <- stale[!startsWith(stale, prefix)]
  if (length(stale) > 0) {
    unlink(file.path(dir, "stale"), recursive = TRUE)
    dir.create(file.path(dir, "stale"))
    file.rename(file.path(dir, stale), file.path(dir, "stale", stale))
    log_msg("[%s] %d saved file(s) of another fingerprint moved to %s/stale/ (an older stale set there deleted)",
            mode, length(stale), basename(dir))
  }
  gb <- function(x) if (is.na(x)) "not measured" else sprintf("%.2f GB", x)
  d <- variant_design(a)
  united <- base$united
  rows <- sequence_parts(united$chr)
  log_msg("[%s] %s: fit by sequence on the united table of the fit %s, %d parts (%s), %d CpGs, %d libraries, %d core(s)%s",
          mode, d$described, a$base_fit_id, length(rows), paste(sprintf("%s %d", names(rows), lengths(rows)), collapse = ", "),
          nrow(united), sum(d$keep), config$cores, if (a$overdispersion == "MN") "; overdispersion MN, chi-squared test" else "")
  pvalue <- rep(NA_real_, nrow(united)); meth.diff <- pvalue
  for (k in seq_along(rows)) {
    name <- names(rows)[k]; r <- rows[[k]]
    file <- file.path(dir, paste0(prefix, d$contrast, "_", name, ".rds"))
    x <- if (file.exists(file)) tryCatch(suppressWarnings(readRDS(file)), error = function(e) NULL)
    if (!(is.list(x) && identical(x$fingerprint_md5, fp) && identical(x$contrast, d$contrast) && identical(x$part, name) &&
          identical(x$rows, r) && is.double(x$pvalue) && length(x$pvalue) == length(r) && length(x$meth.diff) == length(r))) {
      t1 <- Sys.time(); watch <- memory_watch_start()
      obj <- variant_object(methylKit::select(united, r), d)
      x <- suppressMessages(variant_test_part(obj, d, a))
      rm(obj)
      peak <- memory_watch_stop(watch)
      if (!is.double(x$pvalue) || length(x$pvalue) != length(r) || anyNA(x$pvalue) || length(x$meth.diff) != length(r)) {
        stop("The ", d$contrast, " test of part ", name, " did not return a p-value for every CpG (a worker may have failed); ",
             "run the step again (the saved parts are kept), with fewer cores if memory is short.", call. = FALSE)
      }
      minutes <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
      x <- c(x, list(fingerprint_md5 = fp, contrast = d$contrast, part = name, rows = r, minutes = minutes, peak_gb = peak, cores = config$cores))
      save_rds_atomic(x, file)
      log_msg("[%s] %s part %d of %d, %s: %d CpGs fitted in %.1f min on %d core(s), peak memory %s; saved",
              mode, d$contrast, k, length(rows), name, length(r), minutes, config$cores, gb(peak))
    } else {
      log_msg("[%s] %s part %d of %d, %s: %d CpGs, saved part reused (fitted in %.1f min on %d core(s), peak memory %s)",
              mode, d$contrast, k, length(rows), name, length(r), x$minutes, x$cores, gb(x$peak_gb))
    }
    pvalue[r] <- x$pvalue; meth.diff[r] <- x$meth.diff
  }
  list(mode = mode, contrast = d$contrast, base_fit_id = a$base_fit_id, pvalue = pvalue, meth.diff = meth.diff)
}

# What the new contrast depends on (see the section header): a change of any of it refits the new contrast only.
variant_fingerprint <- function(a, base_fingerprint) {
  d <- variant_design(a)
  variant_code <- list(variant_design = variant_design, variant_object = variant_object, variant_test_part = variant_test_part,
                       fit_variant_by_sequence = fit_variant_by_sequence, sequence_parts = sequence_parts,
                       diffmeth_from_parts = diffmeth_from_parts)
  list(mode = a$fit_id, contrast_variant = a$contrast_variant, base_fit_id = a$base_fit_id,
       base_fingerprint_key = fingerprint_key(base_fingerprint),
       sample_sheet = paste(lib$library_id, lib$legacy_sample_index, lib$ploidy, lib$treatment),
       design = c(contrast = d$contrast, libraries = paste(lib$library_id[d$keep], collapse = ","),
                  tested = paste(d$treatment, collapse = ""),
                  covariates = if (is.null(d$covariates)) "none" else paste(names(d$covariates), collapse = " + ")),
       model = c(overdispersion = a$overdispersion, test = a$test),
       methylKit_version = as.character(utils::packageVersion("methylKit")),
       fit_method = "by sequence: parts of the united table of the base fit tested separately, q-values over all CpGs",
       variant_code = lapply(variant_code, function(f) deparse(utils::removeSource(f))))
}

# The fit of an analysis with ctrl or inter: the base fit from cached_fit(), and
# the new contrast from methylkit_<fit id>.rds where its fingerprint is that of
# this run, else fitted by sequence (under the lock of parts_<fit id>/). Returns
# what cached_fit() returns, with the fingerprint of the new contrast, and for
# inter the four group percentages (interaction_groups). fit_base and fit_parts
# are for tests (tests/check_contrast_variants.R); the step never passes them.
cached_variant_fit <- function(a, fit_base = cached_fit, fit_parts = fit_variant_by_sequence) {
  mode <- a$fit_id
  base <- fit_base(analysis_spec(a$base_fit_id))
  fingerprint <- variant_fingerprint(a, base$fingerprint)
  cache_file <- file.path(step_cache_dir(step), sprintf("methylkit_%s.rds", mode))
  x <- NULL
  if (file.exists(cache_file)) {
    cached <- tryCatch(suppressWarnings(readRDS(cache_file)), error = function(e) NULL)
    readable <- is.list(cached) && all(c("mode", "contrast", "pvalue", "meth.diff", "fingerprint") %in% names(cached)) &&
      length(cached$pvalue) == nrow(base$united)
    changed <- if (readable) c(names(fingerprint)[!vapply(names(fingerprint), function(p) identical(cached$fingerprint[[p]], fingerprint[[p]]), logical(1))],
                               setdiff(names(cached$fingerprint), names(fingerprint)))
    if (readable && length(changed) == 0) {
      log_msg("[%s] using the cached %s contrast %s (on the fit %s): base fit, design, model, methylKit version and code are unchanged",
              mode, cached$contrast, basename(cache_file), a$base_fit_id)
      remove_leftovers(a, cache_file)
      x <- cached
    } else {
      log_msg("[%s] the cached contrast %s %s; fitting again", mode, basename(cache_file),
              if (!readable) "cannot be read (cut short, or not of this analysis)"
              else paste("differs from this run in:", paste(gsub("_", " ", changed), collapse = ", ")))
    }
  }
  if (is.null(x)) {
    lock <- take_parts_lock(a)
    on.exit(drop_parts_lock(lock), add = TRUE)
    x <- c(fit_parts(a, base, fingerprint), list(fingerprint = fingerprint))
    save_rds_atomic(x, cache_file)
    unlink(parts_dir(a), recursive = TRUE)
    log_msg("[%s] %s contrast written to %s; its saved parts deleted", mode, x$contrast, basename(cache_file))
  }
  d <- variant_design(a)
  obj <- variant_object(base$united, d)
  diffmeth <- base$diffmeth
  groups <- NULL
  if (d$contrast == "ploidy") {
    diffmeth$ploidy <- diffmeth_from_parts(obj, x$pvalue, x$meth.diff)
  } else {
    groups <- interaction_groups(base$united)
    diffmeth$interaction <- diffmeth_from_parts(obj, x$pvalue, groups$interaction_diff)
  }
  rm(obj)
  list(mode = mode, filtering = base$filtering, united = base$united, diffmeth = diffmeth,
       interaction_groups = groups, base_fit_id = a$base_fit_id, fingerprint = fingerprint)
}

# ---- Write the outputs of one mode -----------------------------------------------
write_mode <- function(fit, a) {
  mode <- a$id                                  # names and the column mode carry the analysis id
  min_diff <- a$dml_difference                  # 20 points (config$dml_difference), or 50 with option d50
  max_q <- a$dml_qvalue                         # 0.01 (config$dml_qvalue), or 0.05 with option q05
  tab <- function(name) analysis_path(a, step, "tables", name)

  write_table(fit$filtering, tab(sprintf("library_filtering_%s.csv", mode)))

  united <- as.data.table(getData(fit$united))
  tests <- lapply(fit$diffmeth, function(x) as.data.table(getData(x)))
  # Option inter: the third contrast, "interaction", whose meth.diff is the
  # interaction difference (cached_variant_fit()); its four group percentages are
  # written next to it.
  groups <- fit$interaction_groups
  for (contrast in names(tests)) {
    stopifnot(identical(tests[[contrast]]$chr, united$chr), identical(tests[[contrast]]$start, united$start))
    # the DML of the fit on all CpGs, as methylKit's getMethylDiff gives them
    stopifnot(tests[[contrast]][qvalue < max_q & abs(meth.diff) > min_diff, .N] ==
                nrow(dml_table(fit$diffmeth[[contrast]], min_diff, max_q)))
  }
  # Option ctexcl: the listed CpGs leave the tested CpGs after the test; p-values,
  # q-values and differences of the other CpGs are those of the fit on all CpGs.
  if (a$ct_excluded) {
    variants <- ct_variant_table()
    listed <- !is.na(variants[united[, .(chr, start)], on = c("chr", "start"), which = TRUE, mult = "first"])
    assess <- ct_variant_assessability()
    assessable <- assess[united[, .(chr, start)], on = c("chr", "start"), assessable]
    assessable[is.na(assessable)] <- FALSE
    removed <- united[listed, .(chr, start, end)]
    summary <- list()
    for (contrast in names(tests)) {
      x <- tests[[contrast]]
      is_dml <- x$qvalue < max_q & abs(x$meth.diff) > min_diff
      removed[, (paste0(c("pvalue_", "qvalue_", "meth.diff_", "dml_"), contrast)) :=
                list(x$pvalue[listed], x$qvalue[listed], x$meth.diff[listed], is_dml[listed])]
      summary[[contrast]] <- data.table(mode = mode, contrast = contrast, cpg_tested_before = nrow(x), cpg_removed = sum(listed),
                                        cpg_tested_after = sum(!listed), dml_before = sum(is_dml), dml_removed = sum(is_dml & listed),
                                        dml_after = sum(is_dml & !listed),
                                        cpg_kept_not_assessable = sum(!listed & !assessable),
                                        dml_kept_not_assessable = sum(is_dml & !listed & !assessable),
                                        min_crick_reads_assessable = config$ct_variants$min_crick_assessable)
    }
    write_table(rbindlist(summary), tab(sprintf("ct_exclusion_summary_%s.csv", mode)))
    write_table(removed, tab(sprintf("ct_excluded_cpgs_%s.csv.gz", mode)))
    log_msg("[%s] C-to-T exclusion: %d of %d tested CpGs removed; ploidy DML %d -> %d, heat DML %d -> %d; kept CpGs not assessable: %d",
            mode, sum(listed), nrow(united), summary$ploidy$dml_before, summary$ploidy$dml_after,
            summary$heat$dml_before, summary$heat$dml_after, summary$ploidy$cpg_kept_not_assessable)
    united <- united[!listed]
    tests <- lapply(tests, function(x) x[!listed])
    if (!is.null(groups)) groups <- groups[!listed]
  }
  write_table(united, tab(sprintf("united_cpg_%s.csv.gz", mode)))

  dml <- list()
  for (contrast in names(tests)) {
    all_cpg <- copy(tests[[contrast]])
    all_cpg[, united_row := .I]                 # row in united_cpg_<mode>; the legacy row name (renumbered after an exclusion)
    if (contrast == "interaction") {
      stopifnot(!is.null(groups), nrow(groups) == nrow(all_cpg), identical(groups$interaction_diff, all_cpg$meth.diff))
      all_cpg <- cbind(all_cpg, groups[, .(pct_diploid_control, pct_diploid_heat, pct_triploid_control, pct_triploid_heat)])
    }
    write_table(all_cpg, tab(sprintf("diffmeth_%s_%s.csv.gz", contrast, mode)))
    d <- all_cpg[qvalue < max_q & abs(meth.diff) > min_diff]
    stopifnot(all(d$united_row <= nrow(united)))
    d[, direction := fifelse(meth.diff > 0, "hyper", "hypo")]
    write_table(d, tab(sprintf("dml_%s_%s.csv", contrast, mode)))
    dml[[contrast]] <- d
  }

  # Loci in both lists. This is the legacy "interaction-DML" set; it is an
  # overlap of two main-effect lists, not a tested interaction (issue I-02).
  shared <- merge(dml$ploidy, dml$heat, by = c("chr", "start", "end", "strand", "united_row"),
                  suffixes = c("_ploidy", "_heat"))
  shared[, class := fifelse(meth.diff_ploidy > 0 & meth.diff_heat > 0, "both_hyper",
                    fifelse(meth.diff_ploidy < 0 & meth.diff_heat < 0, "both_hypo", "opposite"))]
  setorder(shared, united_row)
  write_table(shared, tab(sprintf("dml_shared_%s.csv", mode)))

  # Optional analysis: the rows of the united table at the DML of either list
  # (the loci in both lists are among them), with their united_row. The per-CpG
  # tables of an optional analysis other than 'corrected' are not tracked by
  # git; this small table is, and step 13 reads it.
  if (!a$paper) {
    at_dml <- sort(unique(c(dml$ploidy$united_row, dml$heat$united_row)))
    write_table(cbind(united_row = at_dml, united[at_dml]), tab(sprintf("united_at_dml_%s.csv", mode)))
  }

  # BED files (0-based start, as in the coverage files): chr, start, end, difference.
  bed_dir <- analysis_dir(a, step, "bed")
  write_bed <- function(d, diff_col, name) {
    fwrite(d[, .(chr, start, end, get(diff_col))], file.path(bed_dir, sprintf("dml_%s_%s.bed", name, mode)),
           sep = "\t", col.names = FALSE, scipen = 50, eol = "\n")   # fwrite would write CRLF on Windows
  }
  write_bed(dml$ploidy, "meth.diff", "ploidy"); write_bed(dml$heat, "meth.diff", "heat")
  write_bed(shared, "meth.diff_heat", "shared")
  if (!is.null(dml$interaction)) write_bed(dml$interaction, "meth.diff", "interaction")

  # DML counts across difference thresholds (Figure 2A, 2B), at the q-value cutoff of the analysis.
  sweep <- rbindlist(lapply(config$dml_difference_sweep, function(th) {
    p <- tests$ploidy[qvalue < max_q & abs(meth.diff) > th]
    h <- tests$heat[qvalue < max_q & abs(meth.diff) > th]
    s <- merge(p, h, by = c("chr", "start"), suffixes = c("_ploidy", "_heat"))
    data.table(mode = mode, difference = th, set = c("ploidy", "heat", "shared"),
               n = c(nrow(p), nrow(h), nrow(s)),
               n_hyper = c(sum(p$meth.diff > 0), sum(h$meth.diff > 0), sum(s$meth.diff_heat > 0)),
               n_hypo  = c(sum(p$meth.diff < 0), sum(h$meth.diff < 0), sum(s$meth.diff_heat < 0)))
  }))
  sweep[, pct_hyper := fifelse(n > 0, 100 * n_hyper / n, NA_real_)]
  write_table(sweep, tab(sprintf("dml_counts_by_threshold_%s.csv", mode)))

  headline <- data.table(mode = mode, cpg_analysed = nrow(united),
             ploidy_dml = nrow(dml$ploidy), ploidy_hyper = sum(dml$ploidy$direction == "hyper"),
             ploidy_hypo = sum(dml$ploidy$direction == "hypo"),
             heat_dml = nrow(dml$heat), heat_hyper = sum(dml$heat$direction == "hyper"),
             heat_hypo = sum(dml$heat$direction == "hypo"),
             shared = nrow(shared), shared_both_hyper = sum(shared$class == "both_hyper"),
             shared_both_hypo = sum(shared$class == "both_hypo"), shared_opposite = sum(shared$class == "opposite"))
  # Option inter only: four columns more, after the others (rebuild_comparison()
  # leaves them empty for the analyses without inter); the table of every other
  # analysis, the paper's included, keeps its columns.
  if (!is.null(dml$interaction)) {
    i <- dml$interaction
    in_shared <- sum(i$united_row %in% shared$united_row)
    headline[, `:=`(interaction_dml = nrow(i), interaction_hyper = sum(i$direction == "hyper"),
                    interaction_hypo = sum(i$direction == "hypo"), interaction_in_shared = in_shared)]
  }
  headline
}

# ---- Run ----------------------------------------------------------------------------
an <- config$analysis
if (an$paper && config$run_kind == "optional") {
  log_msg("the analysis of this optional run is the paper's: its DML are the paper's (04_output/05_differential_methylation/); nothing to do")
  quit(save = "no", status = 0)
}
check_stamp(an, needs = "04")
# Option ctexcl: the variant list must be the one step 04 used (its md5 is in the stamp).
if (an$ct_excluded) {
  if (!ct_variant_list_present()) invisible(ct_variant_table())          # stops with the command that makes the list
  stamped <- read_stamp(an$id)$ct_variant_list_md5
  if (!identical(stamped, ct_variant_list_md5())) {
    stop("The C-to-T variant list (04_output/not_adopted/ct_variants/tables/ct_variant_cpgs.csv.gz, md5 ", ct_variant_list_md5(),
         ") is not the list that step 04 used for the analysis '", an$id, "' (md5 in its stamp: '", stamped, "'). Run step 04 for it again.",
         call. = FALSE)
  }
}
fit <- if (an$contrast_variant == "none") cached_fit(an) else cached_variant_fit(an)     # options ctrl and inter: the base fit and the new contrast
headline <- write_mode(fit, an)
fingerprint_md5 <- digest_md5(fit$fingerprint)
rm(fit); invisible(gc(verbose = FALSE))
write_table(headline, analysis_path(an, step, "tables", "dml_counts_by_mode.csv"))
print(t(headline))
if (config$run_kind == "optional") {
  write_stamp(an, step, fit_fingerprint_md5 = fingerprint_md5)
  rebuild_comparison(step, "dml_counts_by_mode.csv")
}

finish_step(step, t0)

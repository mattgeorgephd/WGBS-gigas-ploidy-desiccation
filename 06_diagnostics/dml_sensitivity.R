# dml_sensitivity.R
#
# Purpose   Sensitivity analyses of the differential methylation calls of step 05.
#           The audit (issues.md) raises several defensible alternatives to the
#           analysis that the draft ran. This script repeats the per-CpG tests
#           under each alternative and reports how the lists of differentially
#           methylated loci (DML) change, so that the authors can decide what the
#           paper should report. It is a diagnostics script: nothing in
#           01_pipeline/ reads its outputs.
#             A0  primary: exactly step 05. Ploidy with heat as covariate, heat
#                 with ploidy as covariate, methylKit defaults. This is the
#                 equivalence check; it must reproduce step 05.
#             A1  overdispersion (issue I-06): the same two contrasts with the
#                 variation among animals estimated ("MN"), once with the F test
#                 and once with the chi-squared test.
#             A2  a tested ploidy by heat interaction (issue I-02), without and
#                 with "MN" and the F test, cross-tabulated with the loci that
#                 are in both primary lists (the legacy "interaction-DML").
#                 Heat is confounded with cohort (issue I-01), so this is a
#                 ploidy by cohort interaction, and must be described as one.
#             A3  both contrasts without the heat libraries that also had the
#                 acute heat shock (issue I-04).
#             A4  ploidy within one cohort (issues I-07 and I-01): diploid
#                 against triploid among the controls alone, and among the
#                 heat-treated libraries alone. The control-only contrast is what
#                 the draft's phrase "under control conditions" describes. Both
#                 are free of the cohort confound.
#             A5  heat within one ploidy: among diploids alone, and among
#                 triploids alone. Still a contrast between cohorts (issue I-01).
#           Each analysis is one row of the registry below and one call of
#           methylKit::calculateDiffMeth. Every analysis runs in mode "legacy";
#           A0, A1 with the F test, and A2 also run in mode "corrected".
#           The size of a list does not show by itself whether groups differ.
#           The script therefore also writes what a list size must be read
#           with: the CpGs above the difference cutoff whatever their q-value,
#           the share of small p-values, the test statistic that a call needs,
#           how many CpGs pass the difference cutoff when the ploidy labels are
#           permuted, and the dispersion of the read counts next to the
#           dispersion that binomial reads would give.
#
# Run       Rscript 06_diagnostics/dml_sensitivity.R
#             Full run, after step 05. One model fit per analysis and mode on
#             every tested CpG: 22 fits. With WGBS_CORES=3 on a 4-core machine
#             that other jobs shared, a fit took 5 to 13 minutes (8 on average),
#             which makes about 3 hours for the 22 fits; expect about 2 hours
#             in all on 4 free cores. Fits are cached, so an interrupted run
#             resumes where it stopped; a cache file that cannot be read is
#             fitted again. A run that finds every fit in the cache takes about
#             10 minutes, most of it for the permutation reference.
#           WGBS_SENSITIVITY_SUBSET=20000 Rscript 06_diagnostics/dml_sensitivity.R
#             Test mode: the same code on a seeded random sample of that many
#             CpGs of each united table, sorted by row. 10 to 15 minutes on 2
#             cores. Its outputs test the code and are not results; they go to
#             a separate folder (below). The variable must hold a whole number:
#             a value that is empty or not a number stops the run.
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           04_output/05_differential_methylation/tables/
#             united_cpg_<mode>.csv.gz              tested CpGs with normalized counts
#             diffmeth_{ploidy,heat}_<mode>.csv.gz  step 05 results (equivalence check)
#             dml_{ploidy,heat,shared}_<mode>.csv   the primary lists and the loci in both
#             library_filtering_<mode>.csv          scaling factor of each library (dispersion reference)
#
# Writes    Everything goes to one output folder, <out>:
#             full run   06_diagnostics/dml_sensitivity/
#             test mode  06_diagnostics/dml_sensitivity/test_subset/   (git-ignored)
#           <out>/analyses.csv                      the registry, with group sizes
#           <out>/equivalence_with_step05.csv       A0 against step 05
#           <out>/dml_counts.csv                    per analysis and mode: DML, overlap with the primary list,
#                                                   and what the size of the list must be read with
#           <out>/dml_<analysis id>_<mode>.csv      the CpGs that pass, one file per analysis and mode
#           <out>/interaction_vs_shared_<mode>.csv  tested interaction against the loci in both primary lists
#           <out>/permutation_reference.csv         CpGs above the difference cutoff under permuted ploidy labels
#           <out>/dispersion_summary.csv            dispersion of the counts at a sample of eligible CpGs,
#                                                   observed and for simulated binomial reads
#           <out>/figure_dml_sensitivity.{png,pdf}  DML per analysis, with the part that is in the primary list
#           <out>/session_info.txt
#           <out>/cache/fit_<analysis id>_<mode>.rds   model fits (git-ignored)
#           Every CSV has a column run_scope that says whether it comes from a
#           full run or a test run. Test mode also writes
#           <out>/TEST_RUN_NOT_RESULTS.txt.
#
# Notes     - The tests are methylKit's (Akalin et al. 2012, Genome Biology
#             13:R87, doi:10.1186/gb-2012-13-10-r87): per CpG, a logistic
#             regression of methylated and unmethylated read counts on the tested
#             indicator and the covariates, compared by deviance with the model
#             without the indicator. That paper does not mention overdispersion;
#             it takes the read coverage summed over a group as the sample size.
#           - What the options do, read from the code of methylKit 1.28.0
#             (function logReg) and confirmed against stats::glm:
#               overdispersion "none"  the deviance difference is referred to a
#                 chi-squared distribution with 1 degree of freedom. The
#                 argument `test` has no effect, so A0 with test "Chisq" is the
#                 default call of step 05.
#               overdispersion "MN"    a dispersion is estimated per CpG as the
#                 Pearson chi-squared statistic of the full model over its
#                 residual degrees of freedom (libraries minus parameters; the
#                 estimator that the methylKit documentation attributes to
#                 McCullagh and Nelder 1989, Generalized Linear Models, 2nd ed.,
#                 Chapman and Hall). Estimates below 1 are set to 1, so the
#                 correction never makes a test less conservative. The deviance
#                 difference is divided by the dispersion and referred to
#                 chi-squared with 1 degree of freedom (test "Chisq"), or to F
#                 with 1 and the residual degrees of freedom (test "F"). The F
#                 distribution is used only where the dispersion is above 1.
#             Differences in percent methylation do not depend on these options.
#           - Which of the two tests. In the benchmark of Wreczycka et al. 2017
#             (Journal of Biotechnology 261:105-115,
#             doi:10.1016/j.jbiotec.2017.08.007; read in the preprint version of
#             August 2017, doi:10.1101/109512), the default test without the
#             correction had the lowest specificity and the best overall
#             accuracy, the F test with the correction gave results similar to
#             DSS, one of the stringent methods, and the chi-squared test with
#             the correction is named as a good compromise between stringent
#             and relaxed methods. Their simulation had three samples per group.
#             The help page of stats::anova.glm (R 4.3.3) states that the F test
#             is the most appropriate one when the dispersion is estimated by
#             moments. Both tests are run here. Under relabelled ploidy labels
#             the chi-squared form does not hold its level on these libraries
#             (issues.md, I-06, item 4), and neither form is calibrated where
#             the dispersion estimate is not above 1 (item 7).
#           - q-values are SLIM q-values (Wang, Tuominen and Tsai 2011,
#             Bioinformatics 27:225-231, doi:10.1093/bioinformatics/btq650),
#             computed within each analysis. SLIM estimates the share of true
#             null hypotheses from all p-values of an analysis, so a q-value
#             depends on the set of CpGs tested. Every analysis therefore uses
#             the CpG set of the united table (at least config$min_coverage
#             reads in all libraries), also when it uses fewer libraries. In
#             methylKit 1.28.0 (function QValuesfun) a SLIM q-value is the
#             Benjamini-Hochberg adjusted p-value multiplied by that share.
#           - A DML is a CpG with q < config$dml_qvalue and an absolute effect
#             size above config$dml_difference percentage points. The effect
#             size is the difference in pooled, coverage-weighted percent
#             methylation between the two groups (methylKit's meth.diff):
#             triploid minus diploid, or heat minus control, among the libraries
#             of the analysis. "hyper" is a positive effect size.
#           - A2 tests the indicator "triploid and heat-treated" with ploidy and
#             heat as covariates. With both main effects in the null model, the
#             test of that indicator is the test of the ploidy by heat
#             interaction: the full model has one parameter per group. The
#             script confirms this against a model written as ploidy * heat.
#             methylKit's meth.diff for this coding is triploid heat minus all
#             other libraries, which is not the interaction. The effect size of
#             A2 is computed here as (triploid heat minus triploid control)
#             minus (diploid heat minus diploid control), from the same pooled
#             percentages; a positive value ("hyper" in the list) means that
#             heat minus control is higher in triploids than in diploids. The
#             test is on the log-odds scale and the effect size on the
#             percentage scale, as for the main effects.
#           - Reading the size of a list. A long list is not by itself evidence
#             of differences between the groups, and a short one is not
#             evidence of their absence. dml_counts.csv holds four things to
#             read a list size with.
#             (1) n_effect_above_cutoff: the CpGs above the difference cutoff
#             whatever their q-value. Where it equals n_dml, the q-value
#             criterion removes nothing and the list is the set of CpGs above
#             the difference cutoff. The difference between two small groups
#             is noisier than the difference between two large ones, so that
#             set is larger by chance in the analyses on fewer libraries (A3,
#             A4, A5) and for the interaction contrast of A2, whose four groups
#             are each about half the size of a main-effect group. The
#             permutation reference (next note) shows how large it is by chance.
#             (2) The share of the eligible CpGs (dml_eligible) with a p-value
#             below each value of p_levels. Where the groups do not differ, a
#             test that holds its level gives about that share. A larger share
#             means differences, or a test that does not hold its level; the
#             test without the correction does not hold it where the dispersion
#             is above 1 (dispersion_summary.csv).
#             (3) statistic_needed_for_a_single_call: the test statistic at the
#             p-value config$dml_qvalue / (number of CpGs tested). A CpG with a
#             smaller p-value has q below the cutoff whatever the other CpGs do
#             (Benjamini and Hochberg 1995, Journal of the Royal Statistical
#             Society B 57:289-300; the share of true nulls that SLIM multiplies
#             by is at most 1). Where k CpGs are called together, a p-value k
#             times as large is enough. For the F test this statistic rises
#             steeply as the residual degrees of freedom fall, so a single CpG
#             can hardly be called from a few animals, and a short corrected
#             list can mean little power and not the absence of differences.
#             (4) median_effect_size_eligible. A value away from 0 means that
#             one group is shifted as a whole at CpGs with intermediate
#             methylation. Many CpGs then pass the difference cutoff in one
#             direction (n_hyper against n_hypo) without standing out from the
#             other CpGs.
#             Methods made for few animals share information between CpGs when
#             they estimate the variation among animals, or test regions: DSS
#             (Feng, Conneely and Wu 2014, Nucleic Acids Research 42:e69,
#             doi:10.1093/nar/gku154; general designs: Park and Wu 2016,
#             Bioinformatics 32:1446-1453, doi:10.1093/bioinformatics/btw026)
#             and dmrseq (Korthauer et al. 2019, Biostatistics 20:367-383,
#             doi:10.1093/biostatistics/kxy007). They are not run here.
#           - Permutation reference (permutation_reference.csv). For each
#             contrast and set of libraries, the label of the contrast is
#             permuted among the libraries of each cohort, group sizes kept;
#             the effect size is recomputed from the pooled counts, and the
#             CpGs above the difference cutoff are counted. No model is fitted.
#             Ploidy is the label permuted for the ploidy contrasts and for the
#             interaction. All labellings are used where at most n_permutations
#             exist, otherwise the observed one and a seeded random sample. The
#             share of the labellings that give at least as many CpGs as the
#             observed labels, the observed labelling counted among them, is a
#             permutation p-value (Phipson and Smyth 2010, Statistical
#             Applications in Genetics and Molecular Biology 9:Article 39,
#             doi:10.2202/1544-6115.1585). Its null hypothesis is that the
#             ploidy labels are exchangeable among the libraries of a cohort;
#             its test statistic is the number of CpGs above the cutoff. A
#             large share means that the list of an analysis whose q-value
#             criterion removes nothing is no longer than random labels give.
#             A small share means that the diploid and the triploid libraries
#             differ (for the interaction: that their difference is not the
#             same in the two cohorts). It does not say which CpGs differ; a
#             shift of one whole group gives a small share as well (see (4)),
#             and so does anything else that separates the diploid from the
#             triploid libraries of a cohort, such as a preparation batch. The
#             heat label cannot be permuted, because every heat-treated library
#             is of one cohort and every control of the other (issue I-01): the
#             heat contrasts have no such reference.
#           - Dispersion (dispersion_summary.csv): the Pearson chi-squared
#             statistic of the model ploidy + heat over its residual degrees of
#             freedom, at a seeded sample of eligible CpGs. The binomial model
#             of the primary analysis assumes a dispersion of 1 for reads as
#             sequenced. The counts that are tested are the normalized counts
#             of step 05: the counts of a library are multiplied by its scaling
#             factor (1 or more) and rounded, which multiplies their variance,
#             relative to binomial counts of the same coverage, by about that
#             factor. The reference for the observed row is therefore the first
#             simulated row: binomial reads drawn at the same CpGs from the
#             coverage of each library before scaling and the proportion fitted
#             by the model, then scaled and rounded as in step 05. The second
#             simulated row leaves the reads as drawn; it shows that the
#             estimator gives about 1 for binomial counts.
#           - The primary lists and the loci in both are read from step 05. A0
#             is fitted first. In a full run it must reproduce the p-values,
#             q-values and lists of step 05 exactly, and the script stops before
#             any other analysis is fitted if it does not. The message names
#             what differs, and equivalence_with_step05.csv holds the
#             comparison in that case too. In test mode the p-values and
#             differences of A0 must equal those of step 05, but q-values are
#             computed on the sample, so the A0 lists can differ slightly from
#             the primary lists.
#           - glm.fit and methylKit raise warnings at some CpGs: a fitted
#             probability numerically 0 or 1, no convergence, a full model that
#             fits no better than the reduced one, a p-value that cannot be
#             computed and is set to 1. These warnings are counted by message
#             and logged, not printed one by one. A CpG at which a group has no
#             methylated read does not always raise one, so the counts are not
#             counts of such CpGs. With more than one core the warnings are
#             raised in the worker processes and do not reach the log. Any
#             other warning stops the run.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit", "ggplot2", "patchwork"))
if (!requireNamespace("R.utils", quietly = TRUE)) stop("Missing R package: R.utils (needed to read .csv.gz files).", call. = FALSE)

t0 <- start_step("dml_sensitivity", "sensitivity analyses of the differential methylation calls")

# ---- Settings -------------------------------------------------------------------------
seed              <- 20261119     # test subset, dispersion sample and simulation, permuted labels, interaction coding check
n_dispersion      <- 5000         # eligible CpGs in the dispersion sample
dispersion_marks  <- c(1, 1.5, 3) # the share of CpGs with a dispersion above each value is reported
p_levels          <- c(0.05, 0.01, 0.001)  # the share of eligible CpGs with a p-value below each level is reported
n_permutations    <- 1000         # labellings in the permutation reference (all of them where fewer exist)
permutation_cells <- 4e6          # CpGs times labellings held at once in the permutation reference (bounds its memory)
n_coding_check    <- 200          # CpGs on which the interaction coding is compared with ploidy * heat
equal_tolerance   <- 1e-6         # largest difference accepted where two computations must agree
                                  # (relative for p- and q-values, percentage points for effect sizes)
min_subset        <- 1000         # smallest test subset accepted (SLIM needs a distribution of p-values)

# ---- Full run or test mode; the one output folder ----------------------------------------
# Every file of this script is written through out_file(), so a test run cannot
# write into the folder of the full run. The variable selects test mode as soon
# as it is set: a value that is empty or not a number stops the run, so that a
# mistyped test never starts the full analysis.
subset_request <- Sys.getenv("WGBS_SENSITIVITY_SUBSET", unset = NA_character_)
test_mode <- !is.na(subset_request)
subset_n  <- NA_integer_
if (test_mode) {
  if (!grepl("^[0-9]+$", subset_request)) {
    stop("WGBS_SENSITIVITY_SUBSET must be a whole number of CpGs, not '", subset_request,
         "'. Unset the variable for the full run.", call. = FALSE)
  }
  subset_n <- suppressWarnings(as.integer(subset_request))
  if (is.na(subset_n) || subset_n < min_subset) {
    stop("WGBS_SENSITIVITY_SUBSET must be between ", min_subset, " and the number of tested CpGs.", call. = FALSE)
  }
}
out_dir <- here("06_diagnostics", "dml_sensitivity")
if (test_mode) out_dir <- file.path(out_dir, "test_subset")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
out_file  <- function(name) file.path(out_dir, name)
scope_text <- if (test_mode) sprintf("TEST RUN: random subset of %d CpGs per mode; not a result", subset_n) else "all tested CpGs"
log_msg("%s; outputs in %s", scope_text, out_dir)
if (test_mode) {
  writeLines(c(sprintf("The files in this folder come from a test run of 06_diagnostics/dml_sensitivity.R on a random sample of %d CpGs per mode (seed %d).",
                       subset_n, seed),
               "They test the code. They are not results. The results of a full run are one folder up."),
             out_file("TEST_RUN_NOT_RESULTS.txt"))
}
step05_file <- function(name) {
  f <- output_path("05_differential_methylation", "tables", name)
  if (!file.exists(f)) stop("Missing ", f, ". Run 01_pipeline/05_differential_methylation.R first.", call. = FALSE)
  f
}
# The sample of a test run must fit into the united table of every mode. This is
# checked here, so that a sample that is too large stops the run before any fit.
if (test_mode) {
  united_rows <- vapply(config$modes, function(m) nrow(fread(step05_file(sprintf("united_cpg_%s.csv.gz", m)), select = 1L)),
                        integer(1))
  if (any(united_rows < subset_n)) {
    stop("WGBS_SENSITIVITY_SUBSET (", subset_n, ") is larger than the number of tested CpGs: ",
         paste(sprintf("%d in mode %s", united_rows, config$modes), collapse = ", "), ".", call. = FALSE)
  }
}

# ---- Libraries, sets of libraries, contrasts ------------------------------------------------
lib <- read_libraries()
if (!all(lib$acute_heat_shock %in% c("yes", "no"))) {
  stop("Column acute_heat_shock of wgbs_libraries.csv must be 'yes' or 'no' in every row.", call. = FALSE)
}
if (is.null(lib$cohort) || anyNA(lib$cohort) || !all(nzchar(lib$cohort))) {
  stop("Column cohort of wgbs_libraries.csv must name the experiment of every library.", call. = FALSE)
}

# Sets of libraries, defined from the library sheet and never from positions.
library_sets <- list(
  all                      = list(label = "all libraries",            keep = rep(TRUE, nrow(lib))),
  without_acute_heat_shock = list(label = "without acute heat shock", keep = lib$acute_heat_shock == "no"),
  control_only             = list(label = "controls only",            keep = lib$treatment == "control"),
  heat_only                = list(label = "heat-treated only",        keep = lib$treatment == "heat"),
  diploid_only             = list(label = "diploids only",            keep = lib$ploidy == "diploid"),
  triploid_only            = list(label = "triploids only",           keep = lib$ploidy == "triploid"))

# The tested indicator of each contrast (1 against 0), what its effect size is,
# and the step 05 list that its DML are compared with. `relabelled` is the label
# that the permutation reference permutes among the libraries of each cohort:
# ploidy for the interaction, because the heat label cannot be permuted.
contrast_defs <- list(
  ploidy = list(indicator = lib$ploidy == "triploid", one = "triploid", zero = "diploid", primary = "ploidy",
                effect = "triploid minus diploid, pooled percent methylation",
                primary_label = "primary ploidy DML (step 05)",
                relabelled = lib$ploidy == "triploid", relabelled_name = "ploidy"),
  heat = list(indicator = lib$treatment == "heat", one = "heat", zero = "control", primary = "heat",
              effect = "heat minus control, pooled percent methylation",
              primary_label = "primary heat DML (step 05)",
              relabelled = lib$treatment == "heat", relabelled_name = "heat"),
  interaction = list(indicator = lib$ploidy == "triploid" & lib$treatment == "heat", one = "triploid heat",
                     zero = "all other libraries", primary = "shared",
                     effect = "(triploid heat minus triploid control) minus (diploid heat minus diploid control), pooled percent methylation",
                     primary_label = "loci in both primary lists (step 05; the legacy interaction-DML)",
                     relabelled = lib$ploidy == "triploid", relabelled_name = "ploidy"))

# Covariates, as step 05 passes them to methylKit.
covariate_columns <- list(ploidy = lib$ploidy, heat = lib$treatment)
covariate_frame <- function(covariate_spec, keep) {
  if (covariate_spec == "none") return(NULL)
  terms <- strsplit(covariate_spec, " + ", fixed = TRUE)[[1]]
  if (!all(terms %in% names(covariate_columns)) || anyDuplicated(terms)) {
    stop("Unknown covariates '", covariate_spec, "'. Use 'none', 'ploidy', 'heat', or 'ploidy + heat'.", call. = FALSE)
  }
  as.data.frame(lapply(covariate_columns[terms], function(v) v[keep]))
}

# ---- Registry of analyses --------------------------------------------------------------------
# One entry per analysis. To add an analysis, add an entry; its id becomes part
# of the output file names.
#   contrast        "ploidy", "heat", or "interaction" (a name of contrast_defs)
#   libraries       a name of library_sets
#   covariates      "none", "ploidy", "heat", or "ploidy + heat"
#   overdispersion  "none" or "MN";  test  "Chisq" or "F" (F needs "MN")
#   modes           "legacy" or both_modes
both_modes <- "legacy; corrected"
analysis <- function(analysis_id, question, contrast, libraries, covariates, overdispersion, test, modes = "legacy") {
  data.table(analysis_id, question, contrast, libraries, covariates, overdispersion, test, modes)
}
registry <- rbindlist(list(
  analysis("A0_ploidy", "Reference: the ploidy DML of step 05, as in the draft.",
           "ploidy", "all", "heat", "none", "Chisq", both_modes),
  analysis("A0_heat", "Reference: the heat DML of step 05, as in the draft.",
           "heat", "all", "ploidy", "none", "Chisq", both_modes),
  analysis("A1_ploidy_mn_f", "Which ploidy DML remain when the variation among animals is estimated and the F test is used? Issue I-06.",
           "ploidy", "all", "heat", "MN", "F", both_modes),
  analysis("A1_heat_mn_f", "Which heat DML remain when the variation among animals is estimated and the F test is used? Issue I-06.",
           "heat", "all", "ploidy", "MN", "F", both_modes),
  analysis("A1_ploidy_mn_chisq", "As A1_ploidy_mn_f with the chi-squared test, the less stringent form of the correction. Issue I-06.",
           "ploidy", "all", "heat", "MN", "Chisq"),
  analysis("A1_heat_mn_chisq", "As A1_heat_mn_f with the chi-squared test, the less stringent form of the correction. Issue I-06.",
           "heat", "all", "ploidy", "MN", "Chisq"),
  analysis("A2_interaction", "Where does the difference between the heat cohort and the control cohort depend on ploidy? The tested interaction. Issues I-02 and I-01.",
           "interaction", "all", "ploidy + heat", "none", "Chisq", both_modes),
  analysis("A2_interaction_mn_f", "The same interaction test with the variation among animals estimated and the F test. Issues I-02 and I-06.",
           "interaction", "all", "ploidy + heat", "MN", "F", both_modes),
  analysis("A3_ploidy_no_acute", "Do the ploidy DML depend on the heat libraries that also had the acute heat shock? Issue I-04.",
           "ploidy", "without_acute_heat_shock", "heat", "none", "Chisq"),
  analysis("A3_heat_no_acute", "Do the heat DML depend on the heat libraries that also had the acute heat shock? Issue I-04.",
           "heat", "without_acute_heat_shock", "ploidy", "none", "Chisq"),
  analysis("A4_ploidy_in_control", "Which CpGs differ between diploids and triploids among the controls alone ('under control conditions')? Issues I-07 and I-01.",
           "ploidy", "control_only", "none", "none", "Chisq"),
  analysis("A4_ploidy_in_control_mn_f", "As A4_ploidy_in_control with the variation among animals estimated and the F test.",
           "ploidy", "control_only", "none", "MN", "F"),
  analysis("A4_ploidy_in_heat", "Which CpGs differ between diploids and triploids among the heat-treated animals alone? Issues I-07 and I-01.",
           "ploidy", "heat_only", "none", "none", "Chisq"),
  analysis("A4_ploidy_in_heat_mn_f", "As A4_ploidy_in_heat with the variation among animals estimated and the F test.",
           "ploidy", "heat_only", "none", "MN", "F"),
  analysis("A5_heat_in_diploid", "Which CpGs differ between the heat cohort and the control cohort among diploids alone? Still a cohort contrast. Issue I-01.",
           "heat", "diploid_only", "none", "none", "Chisq"),
  analysis("A5_heat_in_triploid", "Which CpGs differ between the heat cohort and the control cohort among triploids alone? Still a cohort contrast. Issue I-01.",
           "heat", "triploid_only", "none", "none", "Chisq")))

# The two rows that must be step 05 itself. The equivalence check uses them.
primary_ids <- c(ploidy = "A0_ploidy", heat = "A0_heat")

# Check every row and add what follows from the library sheet: group sizes,
# residual degrees of freedom, the reference distribution of the test statistic,
# and a short label for the figure.
describe_analysis <- function(spec) {
  id <- spec$analysis_id
  bad <- function(...) stop("Registry row ", id, ": ", ..., call. = FALSE)
  if (!grepl("^A[0-9]+_[a-z0-9_]+$", id)) bad("the id must look like 'A7_short_name'.")
  if (!spec$contrast %in% names(contrast_defs)) bad("unknown contrast '", spec$contrast, "'.")
  if (!spec$libraries %in% names(library_sets)) bad("unknown set of libraries '", spec$libraries, "'.")
  if (!spec$overdispersion %in% c("none", "MN") || !spec$test %in% c("Chisq", "F")) bad("unknown overdispersion or test.")
  if (spec$overdispersion == "none" && spec$test != "Chisq") {
    bad("without an overdispersion correction methylKit always applies the chi-squared test; write test 'Chisq'.")
  }
  spec_modes <- strsplit(spec$modes, "; ", fixed = TRUE)[[1]]
  if (!all(spec_modes %in% config$modes)) bad("mode '", paste(setdiff(spec_modes, config$modes), collapse = "', '"), "' is not in config$modes.")

  keep <- library_sets[[spec$libraries]]$keep
  def  <- contrast_defs[[spec$contrast]]
  tested <- as.integer(def$indicator[keep])
  if (sum(tested == 1) < 2 || sum(tested == 0) < 2) bad("fewer than two libraries in a group of the tested indicator.")
  frame <- data.frame(tested = tested)
  covariates <- covariate_frame(spec$covariates, keep)
  if (!is.null(covariates)) frame <- cbind(frame, covariates)
  design <- stats::model.matrix(~ ., frame)
  if (qr(design)$rank < ncol(design)) bad("the design is not of full rank (a covariate is constant or aliased among these libraries).")
  if (nrow(design) - ncol(design) < 1) bad("no residual degrees of freedom.")
  cells <- table(droplevels(lib$group[keep]))
  # The indicator "triploid and heat" tests the interaction only when both main
  # effects are covariates, so that the full model has one parameter per group.
  if (spec$contrast == "interaction" && (spec$covariates != "ploidy + heat" || ncol(design) != length(cells))) {
    bad("the interaction needs covariates 'ploidy + heat' and libraries from all four groups.")
  }
  family <- sub("_.*$", "", id)
  residual_df <- nrow(design) - ncol(design)
  data.table(family = family,
             label = paste0(family, ": ", library_sets[[spec$libraries]]$label,
                            if (spec$overdispersion == "MN") sprintf("; MN, %s test", if (spec$test == "F") "F" else "chi-squared") else ""),
             tested_indicator = sprintf("%s (1) against %s (0)", def$one, def$zero),
             n_libraries = sum(keep), n_indicator_1 = sum(tested == 1), n_indicator_0 = sum(tested == 0),
             group_sizes = paste(sprintf("%s %d", names(cells), as.integer(cells)), collapse = "; "),
             residual_df = residual_df,
             reference_distribution = if (spec$overdispersion == "MN" && spec$test == "F") {
               sprintf("F with 1 and %d degrees of freedom (chi-squared with 1 where the dispersion estimate is not above 1)", residual_df)
             } else "chi-squared with 1 degree of freedom",
             effect_size_definition = def$effect, primary_list = def$primary_label)
}
if (anyDuplicated(registry$analysis_id)) stop("Duplicated analysis id in the registry.", call. = FALSE)
registry <- cbind(registry, rbindlist(lapply(seq_len(nrow(registry)), function(k) describe_analysis(registry[k]))))
if (anyDuplicated(registry[, .(contrast, label)])) stop("Two analyses of one contrast have the same figure label.", call. = FALSE)
for (contrast_i in names(primary_ids)) {
  other <- setdiff(names(primary_ids), contrast_i)
  ok <- registry[analysis_id == primary_ids[[contrast_i]],
                 .N == 1 && contrast == contrast_i && libraries == "all" && covariates == other && overdispersion == "none" && modes == both_modes]
  if (!isTRUE(ok)) stop("Registry row ", primary_ids[[contrast_i]], " must stay the analysis of step 05.", call. = FALSE)
}
write_table(cbind(registry[, .(analysis_id, question, contrast, tested_indicator, libraries, n_libraries, n_indicator_1, n_indicator_0,
                               group_sizes, covariates, overdispersion, test, residual_df, reference_distribution,
                               effect_size_definition, primary_list, modes, figure_label = label)],
                  q_cutoff = config$dml_qvalue, difference_cutoff = config$dml_difference, run_scope = scope_text),
            out_file("analyses.csv"))
log_msg("%d analyses in the registry; DML: q < %g and absolute effect size above %g points", nrow(registry),
        config$dml_qvalue, config$dml_difference)

# ---- Helpers ---------------------------------------------------------------------------------
# Evaluate an expression and count the warnings it raises by message. The per-CpG
# warnings of methylKit and glm.fit are expected and are returned as counts; any
# other warning is an error.
expected_warnings <- c("fitted probabilities numerically 0 or 1 occurred", "algorithm did not converge",
                       "Full model does not deviate from Reduced model", "Replacing NaN with 1")
count_warnings <- function(expr) {
  counts <- new.env()
  value <- withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (!any(vapply(expected_warnings, grepl, logical(1), x = msg, fixed = TRUE))) stop("Unexpected warning: ", msg, call. = FALSE)
    counts[[msg]] <- (if (is.null(counts[[msg]])) 0L else counts[[msg]]) + 1L
    invokeRestart("muffleWarning")
  })
  tally <- unlist(as.list(counts))
  list(value = value, warnings = if (length(tally) > 0) tally[order(names(tally))] else integer(0))
}
describe_warnings <- function(tally, forked = FALSE) {
  if (length(tally) > 0) return(paste(sprintf("%d x '%s'", tally, names(tally)), collapse = "; "))
  if (forked) "not collected with more than one core" else "none"
}
largest_relative_difference <- function(x, y) max(abs(x - y) / pmax(abs(y), .Machine$double.xmin))

# A methylBase object for all libraries, built from the united table the way
# methylKit::unite() builds one. The library order of the table is the row order
# of read_libraries().
united_columns <- c("chr", "start", "end", "strand",
                    as.vector(rbind(paste0("coverage", seq_len(nrow(lib))), paste0("numCs", seq_len(nrow(lib))),
                                    paste0("numTs", seq_len(nrow(lib))))))
make_methylbase <- function(united) {
  first <- seq(5, by = 3, length.out = nrow(lib))
  methods::new("methylBase", as.data.frame(united), sample.ids = lib$library_id, assembly = "Roslin_v1", context = "CpG",
               treatment = as.integer(contrast_defs$ploidy$indicator), coverage.index = first, numCs.index = first + 1,
               numTs.index = first + 2, destranded = FALSE, resolution = "base")
}

# Pooled, coverage-weighted percent methylation of a group of libraries: all
# methylated reads of the group over all its reads (methylKit's effect "wmean").
pooled_percent <- function(counts, columns) {
  if (length(columns) < 1) stop("A group without libraries.", call. = FALSE)
  methylated <- rowSums(counts$numCs[, columns, drop = FALSE])
  100 * methylated / (methylated + rowSums(counts$numTs[, columns, drop = FALSE]))
}
indicator_difference <- function(counts, spec) {
  keep <- library_sets[[spec$libraries]]$keep
  indicator <- contrast_defs[[spec$contrast]]$indicator
  pooled_percent(counts, which(keep & indicator)) - pooled_percent(counts, which(keep & !indicator))
}
# The four group percentages and the interaction contrast.
interaction_contrast <- function(counts, spec) {
  keep <- library_sets[[spec$libraries]]$keep
  pct <- function(ploidy_level, treatment_level) {
    pooled_percent(counts, which(keep & lib$ploidy == ploidy_level & lib$treatment == treatment_level))
  }
  groups <- data.table(pct_diploid_control = pct("diploid", "control"), pct_diploid_heat = pct("diploid", "heat"),
                       pct_triploid_control = pct("triploid", "control"), pct_triploid_heat = pct("triploid", "heat"))
  groups[, effect_size := (pct_triploid_heat - pct_triploid_control) - (pct_diploid_heat - pct_diploid_control)]
  groups
}

# One model fit: pvalue, qvalue and methylKit's meth.diff for every CpG of the
# united table, from the cache when the cached fit was made from the same input
# with the same settings. A cache file that cannot be read (a run was stopped
# while it was written) is fitted again; a fit is written under a temporary name
# and renamed, so that a file with the final name is complete.
fit_analysis <- function(spec, mode_i, base_object, united, input_id) {
  keep <- library_sets[[spec$libraries]]$keep
  tested <- as.integer(contrast_defs[[spec$contrast]]$indicator[keep])
  covariates <- covariate_frame(spec$covariates, keep)
  signature <- c(analysis_id = spec$analysis_id, mode = mode_i, input_id, libraries = paste(lib$library_id[keep], collapse = ","),
                 tested = paste(tested, collapse = ""), covariates = spec$covariates,
                 covariate_values = paste(unlist(lapply(covariates, as.character)), collapse = ","),
                 overdispersion = spec$overdispersion, test = spec$test)
  cache_file <- file.path(cache_dir, sprintf("fit_%s_%s.rds", spec$analysis_id, mode_i))
  if (file.exists(cache_file)) {
    cached <- tryCatch(suppressWarnings(readRDS(cache_file)), error = function(e) NULL)
    readable <- is.list(cached) && is.data.frame(cached$fit) && is.numeric(cached$minutes) &&
      all(c("pvalue", "qvalue", "meth.diff") %in% names(cached$fit))
    if (readable && identical(cached$signature, signature) && nrow(cached$fit) == nrow(united)) {
      log_msg("[%s] %-26s using the cached fit %s (%.1f min when fitted)", mode_i, spec$analysis_id, basename(cache_file), cached$minutes)
      return(cached$fit)
    }
    log_msg("[%s] %-26s the cached fit %s %s; fitting again", mode_i, spec$analysis_id, basename(cache_file),
            if (readable) "was made from other input or settings" else "cannot be read (cut short or not a fit of this script)")
  }
  started <- Sys.time()
  object <- methylKit::reorganize(base_object, sample.ids = lib$library_id[keep], treatment = tested)
  caught <- count_warnings(suppressMessages(
    methylKit::calculateDiffMeth(object, covariates = covariates, overdispersion = spec$overdispersion, test = spec$test,
                                 mc.cores = config$cores)))
  fit <- as.data.table(methylKit::getData(caught$value))
  if (nrow(fit) != nrow(united) || !identical(fit$chr, united$chr) || !identical(fit$start, united$start)) {
    stop(spec$analysis_id, ": methylKit did not return one row per CpG in the order of the united table.", call. = FALSE)
  }
  fit <- fit[, .(pvalue, qvalue, meth.diff)]
  if (anyNA(fit) || any(fit$pvalue < 0 | fit$pvalue > 1) || any(fit$qvalue < 0)) {
    stop(spec$analysis_id, ": missing or impossible p-values, q-values or differences.", call. = FALSE)
  }
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  part_file <- paste0(cache_file, ".part")
  saveRDS(list(signature = signature, fit = fit, minutes = minutes, warnings = caught$warnings), part_file)
  if (!file.rename(part_file, cache_file)) stop("Could not move ", part_file, " to ", cache_file, ".", call. = FALSE)
  log_msg("[%s] %-26s fitted %d CpGs on %d libraries in %.1f min (%d core(s)); warnings: %s", mode_i, spec$analysis_id,
          nrow(fit), sum(keep), minutes, config$cores, describe_warnings(caught$warnings, forked = config$cores > 1))
  fit
}

# The interaction coding against a model written with the interaction term:
# likelihood ratio test of ploidy:heat, binomial, at a sample of eligible CpGs.
check_interaction_coding <- function(fit, counts, eligible, spec) {
  keep <- library_sets[[spec$libraries]]$keep
  set.seed(seed)
  picked <- if (length(eligible) > n_coding_check) sort(sample(eligible, n_coding_check)) else eligible
  if (length(picked) == 0) stop("No eligible CpG to check the interaction coding of ", spec$analysis_id, " on.", call. = FALSE)
  frame <- data.frame(ploidy = lib$ploidy[keep], heat = lib$treatment[keep])
  additive  <- stats::model.matrix(~ ploidy + heat, frame)
  factorial <- stats::model.matrix(~ ploidy * heat, frame)
  p_factorial <- count_warnings(vapply(picked, function(i) {
    y <- counts$numCs[i, keep]; n <- y + counts$numTs[i, keep]
    d0 <- stats::glm.fit(additive, y / n, weights = n, family = stats::binomial())$deviance
    d1 <- stats::glm.fit(factorial, y / n, weights = n, family = stats::binomial())$deviance
    stats::pchisq(d0 - d1, 1, lower.tail = FALSE)
  }, numeric(1)))$value
  worst <- largest_relative_difference(fit$pvalue[picked], p_factorial)
  if (!is.finite(worst) || worst > equal_tolerance) {
    stop(spec$analysis_id, ": the coded test is not the test of ploidy:heat (largest relative difference in p ", signif(worst, 3), ").",
         call. = FALSE)
  }
  log_msg("  interaction coding of %s agrees with ploidy * heat at %d CpGs (largest relative difference in p %.1e)",
          spec$analysis_id, length(picked), worst)
}

# ---- Permutation reference for the effect-size criterion -----------------------------------------
# How many CpGs pass the difference cutoff when the labels carry no information?
# The label of a contrast (contrast_defs: `relabelled`) is permuted among the
# libraries of each cohort, the effect size is recomputed from the pooled counts,
# and the CpGs above config$dml_difference are counted. No model is fitted.

# Labellings of a 0/1 label that keep the number of ones within each stratum:
# all of them where at most `n_max` exist, otherwise the observed labelling and
# n_max - 1 random ones. One column per labelling; the observed one is the first.
labellings <- function(observed, strata, n_max) {
  blocks <- split(seq_along(observed), strata)
  ones <- vapply(blocks, function(b) sum(observed[b]), numeric(1))
  n_possible <- prod(choose(lengths(blocks), ones))
  if (n_possible <= n_max) {
    # Every choice of the positions of the ones within each stratum.
    choices <- Map(function(b, k) utils::combn(length(b), k, FUN = function(i) b[i], simplify = FALSE), blocks, ones)
    grid <- expand.grid(lapply(choices, seq_along))
    every <- vapply(seq_len(nrow(grid)), function(g) {
      z <- integer(length(observed))
      for (s in seq_along(choices)) z[choices[[s]][[grid[g, s]]]] <- 1L
      z
    }, integer(length(observed)))
    is_observed <- colSums(every != observed) == 0
    if (ncol(every) != n_possible || sum(is_observed) != 1) stop("The labellings were not enumerated correctly.", call. = FALSE)
    labels <- cbind(every[, is_observed, drop = FALSE], every[, !is_observed, drop = FALSE])
  } else {
    set.seed(seed)
    labels <- cbind(observed, vapply(seq_len(n_max - 1L), function(draw) {
      z <- observed
      for (b in blocks) z[b] <- observed[b][sample.int(length(b))]
      z
    }, integer(length(observed))))
  }
  list(labels = unname(labels), n_possible = n_possible, all_used = n_possible <= n_max)
}

# Effect size of a contrast for each column of `labels` (0/1 labels of the
# libraries in the columns of `cs` and `ns`, the methylated reads and all reads).
# The arithmetic is that of the observed effect size, operation for operation:
# methylKit's 100 * (m1/n1 - m0/n0) for a difference between two groups, and
# interaction_contrast() for the interaction, in which the labels mark the
# triploids and `is_heat` the libraries of the heat cohort. A difference of
# exactly the cutoff then falls on the same side of it as in the DML rule; in
# floating point, 100 * (0.6 - 0.8) is not 100 * 0.6 - 100 * 0.8.
effect_of_labellings <- function(cs, ns, labels, contrast_i, is_heat) {
  if (contrast_i == "interaction") {
    pct <- function(methylated, reads) 100 * methylated / reads
    c_th <- cs %*% (labels * is_heat);  n_th <- ns %*% (labels * is_heat)       # triploid heat under each labelling
    c_tc <- cs %*% (labels * !is_heat); n_tc <- ns %*% (labels * !is_heat)      # triploid control
    c_heat <- rowSums(cs[, is_heat, drop = FALSE]);     n_heat <- rowSums(ns[, is_heat, drop = FALSE])
    c_control <- rowSums(cs[, !is_heat, drop = FALSE]); n_control <- rowSums(ns[, !is_heat, drop = FALSE])
    (pct(c_th, n_th) - pct(c_tc, n_tc)) - (pct(c_heat - c_th, n_heat - n_th) - pct(c_control - c_tc, n_control - n_tc))
  } else {
    c_one <- cs %*% labels; n_one <- ns %*% labels
    100 * (c_one / n_one - (rowSums(cs) - c_one) / (rowSums(ns) - n_one))
  }
}

# One row of permutation_reference.csv: the analyses `spec_rows` share a contrast
# and a set of libraries, and therefore an effect size and a reference.
# `effect_size` is the effect size that their DML rule uses.
permutation_reference <- function(counts, spec_rows, effect_size, mode_i) {
  spec <- spec_rows[1]
  keep <- library_sets[[spec$libraries]]$keep
  def <- contrast_defs[[spec$contrast]]
  drawn <- labellings(as.integer(def$relabelled[keep]), lib$cohort[keep], n_permutations)
  n_used <- ncol(drawn$labels)
  described <- data.table(mode = mode_i, contrast = spec$contrast, libraries = spec$libraries,
                          analyses = paste(spec_rows$analysis_id, collapse = "; "), effect_size_definition = def$effect,
                          n_cpg_tested = length(effect_size), difference_cutoff = config$dml_difference,
                          n_effect_above_cutoff = sum(abs(effect_size) > config$dml_difference))
  quantile_names <- paste0("n_above_cutoff_", c("min", "q05", "q25", "median", "q75", "q95", "max"))
  if (n_used == 1) {
    # One labelling only: the label is constant within every cohort (heat).
    return(cbind(described, permuted_label = sprintf("none: the %s label does not vary within a cohort", def$relabelled_name),
                 n_labellings_possible = drawn$n_possible, n_labellings_used = NA_integer_, labellings_used = "none",
                 as.data.table(stats::setNames(as.list(rep(NA_real_, length(quantile_names))), quantile_names)),
                 n_labellings_at_least_observed = NA_integer_, share_of_labellings_at_least_observed = NA_real_, run_scope = scope_text))
  }
  cs <- counts$numCs[, keep, drop = FALSE] + 0             # as double, for the matrix products
  ns <- cs + counts$numTs[, keep, drop = FALSE]
  is_heat <- lib$treatment[keep] == "heat"
  # The first labelling is the observed one. It must give the effect size of the
  # DML rule, and its count is the count by that effect size.
  worst <- max(abs(effect_of_labellings(cs, ns, drawn$labels[, 1, drop = FALSE], spec$contrast, is_heat)[, 1] - effect_size))
  if (!is.finite(worst) || worst > equal_tolerance) {
    stop("Permutation reference of ", spec$contrast, ", ", spec$libraries, ": the observed labels do not give the effect size of ",
         spec$analysis_id, " (largest difference ", signif(worst, 3), " points).", call. = FALSE)
  }
  n_above <- integer(n_used)
  n_above[1] <- described$n_effect_above_cutoff
  width <- max(1L, as.integer(permutation_cells %/% nrow(cs)))
  for (from in seq(2L, n_used, by = width)) {
    columns <- from:min(from + width - 1L, n_used)
    effect <- effect_of_labellings(cs, ns, drawn$labels[, columns, drop = FALSE], spec$contrast, is_heat)
    if (anyNA(effect)) stop("Permutation reference of ", spec$contrast, ", ", spec$libraries, ": a group without reads.", call. = FALSE)
    n_above[columns] <- as.integer(colSums(abs(effect) > config$dml_difference))
  }
  quantiles <- stats::quantile(n_above, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), names = FALSE)
  n_at_least <- sum(n_above >= n_above[1])
  log_msg(paste("[%s] permutation reference, %s, %s: %d CpGs above the cutoff with the observed labels;",
                "%s labellings give a median of %g (5th to 95th percentile %g to %g); %d of them give at least %d"),
          mode_i, spec$contrast, library_sets[[spec$libraries]]$label, n_above[1],
          if (drawn$all_used) sprintf("all %d", n_used) else sprintf("%d of %.0f", n_used, drawn$n_possible),
          quantiles[4], quantiles[2], quantiles[6], n_at_least, n_above[1])
  cbind(described, permuted_label = sprintf("%s, among the libraries of each cohort", def$relabelled_name),
        n_labellings_possible = drawn$n_possible, n_labellings_used = n_used,
        labellings_used = if (drawn$all_used) "all" else sprintf("the observed labelling and %d random ones (seed %d)", n_used - 1L, seed),
        as.data.table(stats::setNames(as.list(quantiles), quantile_names)),
        n_labellings_at_least_observed = n_at_least, share_of_labellings_at_least_observed = n_at_least / n_used, run_scope = scope_text)
}

# ---- Dispersion of the counts ---------------------------------------------------------------------
# The Pearson chi-squared statistic of the model ploidy + heat, all libraries,
# over its residual degrees of freedom, for each row of a pair of count matrices.
# This is the moment estimate of a quasi-binomial GLM, and the quantity that "MN"
# estimates in A1 before it sets values below 1 to 1. Returns the estimates, the
# fitted proportions, and the counted warnings.
pearson_dispersion <- function(numCs, numTs, design) {
  caught <- count_warnings(lapply(seq_len(nrow(numCs)), function(i) {
    y <- numCs[i, ]; n <- y + numTs[i, ]
    fit <- stats::glm.fit(design, y / n, weights = n, family = stats::quasibinomial())
    mu <- fit$fitted.values
    list(phi = sum((y - n * mu)^2 / (n * mu * (1 - mu))) / fit$df.residual, mu = mu)
  }))
  list(phi = vapply(caught$value, `[[`, numeric(1), "phi"), mu = do.call(rbind, lapply(caught$value, `[[`, "mu")),
       warnings = caught$warnings)
}

# Scaling factor of each library, as step 05 computes it (the largest median
# coverage among libraries over the library's own), from library_filtering_<mode>.csv.
read_scaling_factors <- function(mode_i) {
  path <- step05_file(sprintf("library_filtering_%s.csv", mode_i))
  filtering <- fread(path)
  if (!all(c("library_id", "median_coverage", "scaling_factor") %in% names(filtering)) || !identical(filtering$library_id, lib$library_id)) {
    stop(basename(path), " does not list the libraries of the library sheet in order, with median_coverage and scaling_factor.", call. = FALSE)
  }
  scaling <- max(filtering$median_coverage) / filtering$median_coverage
  if (anyNA(scaling) || any(scaling < 1) || max(abs(scaling - filtering$scaling_factor)) > equal_tolerance) {
    stop("The scaling factors of ", basename(path), " are not the largest median coverage over each library's median coverage.", call. = FALSE)
  }
  scaling
}

# Rows of dispersion_summary.csv for one mode, at a seeded sample of eligible CpGs:
#   observed    the normalized counts of the united table, which the tests use
#   simulated   binomial reads at the same CpGs: for every library, its coverage
#               before scaling and the proportion fitted by the model. Once
#               passed through the scaling and rounding of step 05, as the
#               observed counts were, and once left as drawn.
# The first simulated row is what the binomial model of the primary analysis
# would give for counts of this kind; it is above 1 because step 05 multiplies
# the counts of a library by its scaling factor.
dispersion_summary <- function(united, counts, eligible, mode_i) {
  set.seed(seed)
  picked <- if (length(eligible) > n_dispersion) sort(sample(eligible, n_dispersion)) else eligible
  if (length(picked) < n_dispersion) log_msg("[%s] dispersion: only %d eligible CpGs, all are used", mode_i, length(picked))
  design <- stats::model.matrix(~ ploidy + treatment, lib)
  one_row <- function(source_described, estimate) {
    phi <- estimate$phi[is.finite(estimate$phi)]
    if (length(phi) == 0) stop("No dispersion could be estimated in mode ", mode_i, " (", source_described, ").", call. = FALSE)
    quartiles <- stats::quantile(phi, c(0.25, 0.5, 0.75), names = FALSE)
    shares <- stats::setNames(lapply(dispersion_marks, function(m) mean(phi > m)),
                              paste0("share_above_", gsub(".", "_", as.character(dispersion_marks), fixed = TRUE)))
    log_msg("[%s] dispersion at %d eligible CpGs, %s: median %.2f (quartiles %.2f to %.2f); warnings: %s", mode_i, length(phi),
            source_described, quartiles[2], quartiles[1], quartiles[3], describe_warnings(estimate$warnings))
    cbind(data.table(mode = mode_i, count_source = source_described, model = "ploidy + heat, all libraries",
                     residual_df = nrow(design) - ncol(design), n_cpg_eligible = length(eligible), n_cpg_sampled = length(picked),
                     n_cpg_estimated = length(phi), dispersion_q25 = quartiles[1], dispersion_median = quartiles[2],
                     dispersion_q75 = quartiles[3]),
          as.data.table(shares), run_scope = scope_text)
  }
  observed <- pearson_dispersion(counts$numCs[picked, , drop = FALSE], counts$numTs[picked, , drop = FALSE], design)

  # Coverage of each library before scaling: the whole number that step 05
  # multiplied by the scaling factor and rounded to give the coverage column.
  scaling <- read_scaling_factors(mode_i)
  coverage <- as.matrix(united[picked, paste0("coverage", seq_len(nrow(lib))), with = FALSE]) + 0
  raw_coverage <- round(sweep(coverage, 2, scaling, "/"))
  if (anyNA(raw_coverage) || any(round(sweep(raw_coverage, 2, scaling, "*")) != coverage)) {
    stop("The coverage columns of united_cpg_", mode_i, ".csv.gz are not whole-number coverages multiplied by the scaling factors of ",
         "library_filtering_", mode_i, ".csv and rounded.", call. = FALSE)
  }
  if (anyNA(observed$mu)) stop("Missing fitted proportions in the dispersion sample of mode ", mode_i, ".", call. = FALSE)
  set.seed(seed)
  raw_cs <- matrix(stats::rbinom(length(raw_coverage), raw_coverage, observed$mu), nrow = nrow(raw_coverage))
  raw_ts <- raw_coverage - raw_cs
  scaled <- lapply(seq_len(nrow(lib)), function(j) scale_counts_like_methylkit(raw_coverage[, j], raw_cs[, j], raw_ts[, j], scaling[j]))
  scaled_cs <- vapply(scaled, `[[`, numeric(nrow(raw_coverage)), "numCs")
  scaled_ts <- vapply(scaled, `[[`, numeric(nrow(raw_coverage)), "numTs")
  rbind(one_row("observed (normalized counts of step 05)", observed),
        one_row("simulated binomial reads, scaled and rounded as in step 05", pearson_dispersion(scaled_cs, scaled_ts, design)),
        one_row("simulated binomial reads, not scaled", pearson_dispersion(raw_cs, raw_ts, design)))
}

# ---- One mode: fits, DML lists, comparisons ---------------------------------------------------
# CSVs written here, into the output folder of this run (out_file):
# equivalence_with_step05.csv (rewritten after the check of each mode),
# dml_<analysis id>_<mode>.csv and interaction_vs_shared_<mode>.csv. Returns the
# rows of the other tables that span both modes.

# Rows of equivalence_with_step05.csv of the modes checked so far in this run.
# The table is written as soon as a mode has been checked, also when the check
# fails, so that the folder never keeps the table of an earlier run next to a run
# that stopped.
equivalence_checked <- new.env()
equivalence_checked$rows <- list()

run_mode <- function(mode_i) {
  specs <- registry[vapply(strsplit(modes, "; ", fixed = TRUE), function(m) mode_i %in% m, logical(1))]
  united_path <- step05_file(sprintf("united_cpg_%s.csv.gz", mode_i))
  united <- fread(united_path)
  if (!identical(names(united), united_columns)) {
    stop(basename(united_path), " does not have the columns chr, start, end, strand and coverage, numCs, numTs for ",
         nrow(lib), " libraries in order.", call. = FALSE)
  }
  cpg_rows <- seq_len(nrow(united))               # row of each CpG in united_cpg_<mode> (united_row in step 05)
  if (test_mode) {
    if (subset_n > nrow(united)) stop("WGBS_SENSITIVITY_SUBSET is larger than the ", nrow(united), " CpGs of mode ", mode_i, ".", call. = FALSE)
    set.seed(seed)
    cpg_rows <- sort(sample.int(nrow(united), subset_n))
    united <- united[cpg_rows]
  }
  log_msg("[%s] %d CpGs, %d analyses", mode_i, nrow(united), nrow(specs))
  input_id <- c(united_md5 = unname(tools::md5sum(united_path)),
                cpg_set = if (test_mode) sprintf("random subset of %d, seed %d", subset_n, seed) else "all",
                n_cpg = as.character(nrow(united)), methylKit = as.character(utils::packageVersion("methylKit")))
  count_columns <- function(prefix) {
    m <- as.matrix(united[, paste0(prefix, seq_len(nrow(lib))), with = FALSE])
    if (anyNA(m) || any(m < 0)) stop("Missing or negative counts in ", basename(united_path), ".", call. = FALSE)
    m
  }
  counts <- list(numCs = count_columns("numCs"), numTs = count_columns("numTs"))
  eligible <- which(dml_eligible(mean_percent_methylation(united), config))

  # The step 05 lists among the CpGs tested here: position in `united` and
  # direction (class, for the loci in both lists).
  read_primary <- function(name, value_column) {
    d <- fread(step05_file(sprintf("dml_%s_%s.csv", name, mode_i)))
    if (!all(c("chr", "start", "united_row", value_column) %in% names(d))) stop("Unexpected columns in dml_", name, "_", mode_i, ".csv.", call. = FALSE)
    if (!test_mode && !all(d$united_row %in% cpg_rows)) stop("dml_", name, "_", mode_i, ".csv holds CpGs that are not in the united table.", call. = FALSE)
    d <- d[united_row %in% cpg_rows]
    d[, position := match(united_row, cpg_rows)]
    if (nrow(d) > 0 && (!identical(d$chr, united$chr[d$position]) || !identical(d$start, united$start[d$position]))) {
      stop("dml_", name, "_", mode_i, ".csv and the united table disagree on positions.", call. = FALSE)
    }
    data.table(position = d$position, value = as.character(d[[value_column]]))
  }
  primary <- list(ploidy = read_primary("ploidy", "direction"), heat = read_primary("heat", "direction"),
                  shared = read_primary("shared", "class"))
  if (!all(primary$shared$value %in% c("both_hyper", "both_hypo", "opposite"))) stop("Unexpected class in dml_shared_", mode_i, ".csv.", call. = FALSE)

  # One analysis: the fit (or its cached copy), the effect size, and the calls.
  base_object <- make_methylbase(united)
  analyse <- function(spec) {
    res <- fit_analysis(spec, mode_i, base_object, united, input_id)
    if (spec$contrast == "interaction") {
      if (spec$overdispersion == "none") check_interaction_coding(res, counts, eligible, spec)
      res <- cbind(res[, .(pvalue, qvalue)], interaction_contrast(counts, spec))
    } else {
      # methylKit's difference is used as it is, so that A0 is step 05 to the
      # last digit; the pooled percentages computed here must agree with it.
      if (max(abs(indicator_difference(counts, spec) - res$meth.diff)) > equal_tolerance) {
        stop(spec$analysis_id, ": the pooled group difference computed here is not methylKit's meth.diff.", call. = FALSE)
      }
      res <- res[, .(pvalue, qvalue, effect_size = meth.diff)]
    }
    if (anyNA(res$effect_size)) stop(spec$analysis_id, ": missing effect sizes.", call. = FALSE)
    res[, significant := qvalue < config$dml_qvalue]
    res[, is_dml := significant & abs(effect_size) > config$dml_difference]
    res
  }

  # A0 first, checked against step 05 before any other analysis is fitted, so
  # that a run that cannot reproduce step 05 stops early.
  # CSV: equivalence_with_step05.csv in the output folder of this run, written
  # before the checks stop the run.
  results <- list()
  is_primary <- specs$analysis_id %in% primary_ids
  for (k in which(is_primary)) results[[specs$analysis_id[k]]] <- analyse(specs[k])
  equivalence <- rbindlist(lapply(names(primary_ids), function(contrast_i) {
    res <- results[[primary_ids[[contrast_i]]]]
    ref <- fread(step05_file(sprintf("diffmeth_%s_%s.csv.gz", contrast_i, mode_i)),
                 select = c("chr", "start", "pvalue", "qvalue", "meth.diff", "united_row"))[cpg_rows]
    if (!identical(ref$united_row, cpg_rows) || !identical(ref$chr, united$chr) || !identical(ref$start, united$start)) {
      stop("diffmeth_", contrast_i, "_", mode_i, ".csv.gz does not list the CpGs of the united table in order.", call. = FALSE)
    }
    a0 <- which(res$is_dml)
    row <- data.table(mode = mode_i, contrast = contrast_i, analysis_id = primary_ids[[contrast_i]], n_cpg_compared = nrow(res),
                      largest_relative_difference_pvalue = largest_relative_difference(res$pvalue, ref$pvalue),
                      largest_absolute_difference_effect = max(abs(res$effect_size - ref$meth.diff)),
                      largest_relative_difference_qvalue = if (test_mode) NA_real_ else largest_relative_difference(res$qvalue, ref$qvalue),
                      n_dml_step05 = nrow(primary[[contrast_i]]), n_dml_a0 = length(a0),
                      n_dml_in_both = sum(a0 %in% primary[[contrast_i]]$position),
                      dml_list_identical = if (test_mode) NA else
                        setequal(a0, primary[[contrast_i]]$position) &&
                        identical(fifelse(res$effect_size[primary[[contrast_i]]$position] > 0, "hyper", "hypo"), primary[[contrast_i]]$value))
    row[, p_values_and_effects_equal := largest_relative_difference_pvalue <= equal_tolerance &
          largest_absolute_difference_effect <= equal_tolerance]
    row[, q_values_equal := if (test_mode) NA else largest_relative_difference_qvalue <= equal_tolerance]
    row
  }))
  # The loci in both A0 lists, and their classes by the rule of step 05.
  a0_ploidy <- results[[primary_ids[["ploidy"]]]]; a0_heat <- results[[primary_ids[["heat"]]]]
  a0_class <- fifelse(a0_ploidy$effect_size > 0 & a0_heat$effect_size > 0, "both_hyper",
                      fifelse(a0_ploidy$effect_size < 0 & a0_heat$effect_size < 0, "both_hypo", "opposite"))
  equivalence[, shared_loci_identical := if (test_mode) NA else
    setequal(which(a0_ploidy$is_dml & a0_heat$is_dml), primary$shared$position) &&
    identical(a0_class[primary$shared$position], primary$shared$value)]
  equivalence[, run_scope := scope_text]
  equivalence_checked$rows[[mode_i]] <- equivalence
  write_table(rbindlist(equivalence_checked$rows), out_file("equivalence_with_step05.csv"))
  log_msg("[%s] A0 against step 05: largest relative difference in p %.1e, largest difference in effect size %.1e points", mode_i,
          max(equivalence$largest_relative_difference_pvalue), max(equivalence$largest_absolute_difference_effect))
  # Every check that failed, in words. In test mode the checks of q-values and
  # lists are not made (their columns are NA).
  differing <- unlist(lapply(seq_len(nrow(equivalence)), function(k) {
    eq <- as.list(equivalence[k])
    c(if (!isTRUE(eq$p_values_and_effects_equal)) {
        sprintf("%s, p-values or differences (largest relative difference in p %.3g, largest difference in effect size %.3g points)",
                eq$contrast, eq$largest_relative_difference_pvalue, eq$largest_absolute_difference_effect)
      },
      if (!test_mode && !isTRUE(eq$q_values_equal)) {
        sprintf("%s, q-values (largest relative difference %.3g)", eq$contrast, eq$largest_relative_difference_qvalue)
      },
      if (!test_mode && !isTRUE(eq$dml_list_identical)) {
        sprintf("%s, DML list (%d in step 05, %d here, %d in both; if the three are equal, a direction differs)",
                eq$contrast, eq$n_dml_step05, eq$n_dml_a0, eq$n_dml_in_both)
      })
  }))
  if (!test_mode && !isTRUE(all(equivalence$shared_loci_identical))) differing <- c(differing, "the loci in both lists, or their classes")
  if (length(differing) > 0) {
    stop("A0 does not reproduce step 05 in mode ", mode_i, ". What differs: ", paste(differing, collapse = "; "), ". ",
         "Values must agree within ", equal_tolerance, "; a DML here has q < ", config$dml_qvalue, " and an absolute difference above ",
         config$dml_difference, " points (run_config.R). If run_config.R, the library sheet or the methylKit version changed after ",
         "step 05 was run, run step 05 again. The comparison is in ", out_file("equivalence_with_step05.csv"), ".", call. = FALSE)
  }
  if (test_mode) {
    log_msg("[%s] test mode: q-values are computed on the sample, so the A0 lists are not expected to equal the primary lists exactly (%s)",
            mode_i, paste(sprintf("%s: %d of %d primary DML recovered, %d called", equivalence$contrast, equivalence$n_dml_in_both,
                                  equivalence$n_dml_step05, equivalence$n_dml_a0), collapse = "; "))
  } else {
    log_msg("[%s] A0 reproduces the q-values, both DML lists and the loci in both lists of step 05", mode_i)
  }

  # The other analyses.
  for (k in which(!is_primary)) results[[specs$analysis_id[k]]] <- analyse(specs[k])
  rm(base_object); gc(verbose = FALSE)

  # DML list of each analysis (CSV: dml_<analysis id>_<mode>.csv in the output
  # folder of this run) and one row of dml_counts.csv: the size of the list, its
  # overlap with the primary list of its contrast, and what is needed to read
  # the size: the CpGs above the difference cutoff whatever their q-value, the
  # median effect size and the share of small p-values among the eligible CpGs,
  # and the p-value below which, and the test statistic above which, a CpG has
  # q below the cutoff whatever the other CpGs do.
  counts_rows <- rbindlist(lapply(seq_len(nrow(specs)), function(k) {
    spec <- specs[k]
    res <- results[[spec$analysis_id]]
    reference <- primary[[contrast_defs[[spec$contrast]]$primary]]
    listed <- which(res$is_dml)
    dml <- cbind(united[listed, .(chr, start)], res[listed, .(pvalue, qvalue, effect_size)])
    dml[, direction := fifelse(effect_size > 0, "hyper", "hypo")]
    group_columns <- grep("^pct_", names(res), value = TRUE)          # interaction analyses only
    if (length(group_columns) > 0) dml <- cbind(dml, res[listed, group_columns, with = FALSE])
    dml[, run_scope := scope_text]
    write_table(dml, out_file(sprintf("dml_%s_%s.csv", spec$analysis_id, mode_i)))
    in_reference <- listed %in% reference$position
    same_direction <- if (spec$contrast == "interaction") NA_integer_ else
      sum(dml$direction[in_reference] == reference$value[match(listed[in_reference], reference$position)])
    p_shares <- stats::setNames(lapply(p_levels, function(level) if (length(eligible) > 0) mean(res$pvalue[eligible] < level) else NA_real_),
                                paste0("share_eligible_p_below_",
                                       gsub(".", "_", vapply(p_levels, format, character(1), scientific = FALSE), fixed = TRUE)))
    p_single <- config$dml_qvalue / nrow(res)
    cbind(data.table(analysis_id = spec$analysis_id, mode = mode_i, contrast = spec$contrast, label = spec$label,
                     n_cpg_tested = nrow(res), n_effect_above_cutoff = sum(abs(res$effect_size) > config$dml_difference),
                     n_significant = sum(res$significant), n_dml = length(listed),
                     n_hyper = sum(dml$direction == "hyper"), n_hypo = sum(dml$direction == "hypo"),
                     primary_list = spec$primary_list, n_primary = nrow(reference), n_in_primary = sum(in_reference),
                     n_in_primary_same_direction = same_direction,
                     share_of_primary_recovered = if (nrow(reference) > 0) sum(in_reference) / nrow(reference) else NA_real_,
                     share_of_list_in_primary = if (length(listed) > 0) sum(in_reference) / length(listed) else NA_real_,
                     n_cpg_eligible = length(eligible),
                     median_effect_size_eligible = if (length(eligible) > 0) stats::median(res$effect_size[eligible]) else NA_real_),
          as.data.table(p_shares),
          data.table(p_needed_for_a_single_call = p_single,
                     statistic_needed_for_a_single_call = if (spec$overdispersion == "MN" && spec$test == "F") {
                       stats::qf(p_single, 1, spec$residual_df, lower.tail = FALSE)
                     } else stats::qchisq(p_single, 1, lower.tail = FALSE),
                     run_scope = scope_text))
  }))

  # Permutation reference of each effect size: one row per contrast and set of
  # libraries (CSV: permutation_reference.csv, written with the tables that span
  # the modes). Analyses that differ only in the test share their effect size.
  effect_sets <- unique(specs[, .(contrast, libraries)])
  permutation_rows <- rbindlist(lapply(seq_len(nrow(effect_sets)), function(k) {
    same <- specs[specs$contrast == effect_sets$contrast[k] & specs$libraries == effect_sets$libraries[k]]
    effect_size <- results[[same$analysis_id[1]]]$effect_size
    for (id in same$analysis_id[-1]) {
      if (!identical(results[[id]]$effect_size, effect_size)) stop(id, " and ", same$analysis_id[1], " do not share their effect size.", call. = FALSE)
    }
    permutation_reference(counts, same, effect_size, mode_i)
  }))

  # Tested interaction against the loci in both primary lists.
  # CSV: interaction_vs_shared_<mode>.csv in the output folder of this run.
  shared_class <- rep(NA_character_, nrow(united))
  shared_class[primary$shared$position] <- primary$shared$value
  interaction_ids <- specs[contrast == "interaction", analysis_id]
  if (length(interaction_ids) > 0) {
    cross <- rbindlist(lapply(interaction_ids, function(id) {
      res <- results[[id]]
      one_group <- function(cpg_group, members) {
        data.table(analysis_id = id, mode = mode_i, cpg_group = cpg_group, n_cpg = sum(members),
                   n_interaction_significant = sum(res$significant[members]), n_interaction_dml = sum(res$is_dml[members]))
      }
      rbind(one_group("in both primary lists, both_hyper", !is.na(shared_class) & shared_class == "both_hyper"),
            one_group("in both primary lists, both_hypo", !is.na(shared_class) & shared_class == "both_hypo"),
            one_group("in both primary lists, opposite", !is.na(shared_class) & shared_class == "opposite"),
            one_group("in both primary lists, all classes", !is.na(shared_class)),
            one_group("not in both primary lists", is.na(shared_class)),
            one_group("all tested CpGs", rep(TRUE, nrow(res))))
    }))
    cross[, share_interaction_significant := fifelse(n_cpg > 0, n_interaction_significant / n_cpg, NA_real_)]
    cross[, share_interaction_dml := fifelse(n_cpg > 0, n_interaction_dml / n_cpg, NA_real_)]
    cross[, `:=`(q_cutoff = config$dml_qvalue, difference_cutoff = config$dml_difference, run_scope = scope_text)]
    write_table(cross, out_file(sprintf("interaction_vs_shared_%s.csv", mode_i)))
  }

  list(counts = counts_rows, permutation = permutation_rows, dispersion = dispersion_summary(united, counts, eligible, mode_i))
}

# ---- Run every mode; tables that span the modes -------------------------------------------------
# CSVs written here, into the output folder of this run (out_file):
# dml_counts.csv, permutation_reference.csv, dispersion_summary.csv.
by_mode <- lapply(stats::setNames(config$modes, config$modes), run_mode)
dml_counts  <- rbindlist(lapply(by_mode, `[[`, "counts"))
permutation <- rbindlist(lapply(by_mode, `[[`, "permutation"))
dispersion  <- rbindlist(lapply(by_mode, `[[`, "dispersion"))
write_table(dml_counts[, !"label"], out_file("dml_counts.csv"))
write_table(permutation, out_file("permutation_reference.csv"))
write_table(dispersion, out_file("dispersion_summary.csv"))
print(dml_counts[, .(analysis_id, mode, n_cpg_tested, n_effect_above_cutoff, n_significant, n_dml, n_hyper, n_hypo, n_primary, n_in_primary)])

# Files of an earlier run that this run did not write (an analysis that was
# renamed or removed) would be mistaken for current results.
written <- c(sprintf("dml_%s_%s.csv", dml_counts$analysis_id, dml_counts$mode),
             sprintf("interaction_vs_shared_%s.csv", config$modes))
stale <- setdiff(list.files(out_dir, pattern = "^(dml_A|interaction_vs_shared_).*\\.csv$"), written)
if (length(stale) > 0) log_msg("NOTE: not written by this run, delete if no longer wanted: %s", paste(stale, collapse = ", "))

# ---- Figure: DML per analysis, with the part that is in the primary list ---------------------
# One panel per contrast, one column per mode. The colour of a contrast is the
# colour of its DML set in the pipeline figures (config$colours$dml_set).
panel_titles <- c(ploidy = "Ploidy: triploid against diploid",
                  heat = "Heat: heat against control (a contrast between cohorts, issue I-01)",
                  interaction = "Ploidy by heat interaction, tested (a ploidy by cohort interaction, issue I-01)")
reference_phrases <- c(ploidy = "the primary ploidy DML", heat = "the primary heat DML", interaction = "the loci in both primary lists")
thousands <- function(v) format(v, big.mark = ",", trim = TRUE, scientific = FALSE)
plot_contrast <- function(contrast_i, x_title) {
  d <- dml_counts[contrast == contrast_i]
  d[, label := factor(label, levels = rev(unique(label)))]
  d[, mode := factor(mode, levels = config$modes)]
  parts <- rbind(d[, .(label, mode, part = "in_primary", n = n_in_primary)],
                 d[, .(label, mode, part = "not_in_primary", n = n_dml - n_in_primary)])
  parts[, part := factor(part, levels = c("in_primary", "not_in_primary"))]
  ggplot(parts, aes(n, label)) +
    geom_vline(data = unique(d[, .(mode, n_primary)]), aes(xintercept = n_primary), colour = "grey65", linewidth = 0.35) +
    geom_col(aes(fill = part), width = 0.62, colour = "white", linewidth = 0.5, position = position_stack(reverse = TRUE)) +
    geom_label(data = d, aes(n_dml, label, label = thousands(n_dml)), hjust = 0, size = 3.1, colour = "black", fill = "white",
               label.size = 0, label.r = grid::unit(0, "pt"), label.padding = grid::unit(2.5, "pt")) +
    facet_grid(. ~ mode, labeller = labeller(mode = function(v) paste("mode:", v))) +
    scale_fill_manual(values = c(in_primary = config$colours$dml_set[[contrast_defs[[contrast_i]]$primary]], not_in_primary = "grey55"),
                      labels = c(in_primary = paste("also among", reference_phrases[[contrast_i]]), not_in_primary = "not among them"),
                      name = NULL, drop = FALSE) +
    scale_x_continuous(labels = thousands, expand = expansion(mult = c(0, 0.14))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = x_title, y = NULL, title = panel_titles[[contrast_i]]) +
    theme_wgbs(11) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), legend.key.size = grid::unit(9, "pt"), panel.spacing.x = grid::unit(18, "pt"),
          plot.title = element_text(size = 11, face = "bold"), plot.title.position = "plot")
}
shown <- intersect(names(panel_titles), unique(dml_counts$contrast))
if (!setequal(shown, unique(dml_counts$contrast))) stop("A contrast of the registry has no panel title.", call. = FALSE)
x_titles <- c(rep("", length(shown) - 1),
              sprintf("DML (q < %g, absolute effect size above %g percentage points)", config$dml_qvalue, config$dml_difference))
panels <- Map(plot_contrast, shown, x_titles)
bars <- vapply(shown, function(contrast_i) dml_counts[contrast == contrast_i, uniqueN(label)], integer(1))
figure <- wrap_plots(panels, ncol = 1, heights = bars + 1) +
  plot_annotation(
    title = "DML under each alternative analysis",
    subtitle = if (test_mode) scope_text else NULL,
    caption = paste0("Vertical line: size of the step 05 list that the bars are compared with. ",
                     "A row without a number was not run in that mode.\n",
                     "The size of a list does not show by itself whether groups differ: ",
                     "read it with dml_counts.csv and permutation_reference.csv."),
    theme = theme(plot.title = element_text(face = "bold", size = 13), plot.caption = element_text(hjust = 0, size = 9),
                  plot.subtitle = element_text(face = "bold", size = 11)))
save_figure(figure, out_file("figure_dml_sensitivity"), width = 11, height = 2.8 + 0.36 * sum(bars + 1))

# ---- Session information ---------------------------------------------------------------------
writeLines(capture.output(print(utils::sessionInfo())), out_file("session_info.txt"))
log_msg("dml_sensitivity finished in %.1f min; outputs in %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir)

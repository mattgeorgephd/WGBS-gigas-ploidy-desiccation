# run_config.R
#
# Single control surface for the pipeline. Every script in 01_pipeline/ sources
# this file first. Change analysis settings here, never inside a step.
#
# The settings of the paper are the defaults. Optional analyses, which the paper
# did NOT adopt, are off by default and are set in the list `optional` below or
# with the switches of run_pipeline.R; their outputs go to folders named
# not_adopted (04_output/not_adopted/README.md). The issue numbers refer to
# issues.md.
#
# The file has three parts: the settings, the checks of the settings, and the
# number formatting of the R session. The second part checks the settings of the
# analysis: the optional analyses and the environment variables that choose them,
# the frozen settings of the paper, the cores, the thresholds and other numbers
# of steps 03 to 10, and the measure, threshold, counts, tolerances and seeds of
# steps 11 and 12. Such a setting that cannot work stops every step there,
# before the step reads or writes anything. Not checked: the cache folder, the
# facts about the draft and the legacy files (figure numbers, paths of the
# REVIGO exports, the qPCR design and gene names other than the reference gene),
# the sizes, colours and label placement of the figures of step 11, the plot
# colours, the parameters of the C-to-T variant rule (checked by the variant
# step), and an upper bound on WGBS_CORES.

config <- list(

  # ---- Inputs on gannet ------------------------------------------------------
  # Files listed in 03_input_files/gannet_manifest.csv are downloaded here.
  # Default: 03_input_files/gannet_cache/ (git-ignored). Set the environment
  # variable WGBS_CACHE_DIR to use another location, for example a larger drive.
  cache_dir          = Sys.getenv("WGBS_CACHE_DIR", unset = ""),
  download_timeout_s = 7200,

  # ---- Optional analyses: NOT ADOPTED, all off by default ----------------------
  # With every setting below at its default the run produces the results of the
  # paper and nothing else. Settings (a) to (f) together define ONE optional
  # analysis per run, named by its id (04_output/not_adopted/<id>/; the tokens
  # are joined by "-" in the order corrected, cov5, mnchisq, ctexcl, d50, q05);
  # (g) and (h) add steps. The switches of run_pipeline.R turn the same settings on
  # for one run. A step started alone reads WGBS_ANALYSIS, WGBS_GENE_LEVEL_GO and
  # WGBS_DIAGNOSTICS, which run_pipeline.R sets for its steps.
  optional = list(
    coverage_filter_of_methods = FALSE,         # (a) token "corrected": CpGs at or above the 99.9th percentile of
                                                #     coverage in a library removed; every line of each file read (I-05, I-38)
    coverage_standard          = 10L,           # (b) reads per CpG in each of the 22 libraries: 10 (the paper) or
                                                #     5 (token "cov5"; Venkataraman et al. 2022)
    overdispersion_chisq       = FALSE,         # (c) token "mnchisq": methylKit overdispersion = "MN", test = "Chisq" (I-06);
                                                #     no other test is offered
    exclude_ct_variants        = FALSE,         # (d) token "ctexcl": CpGs at C-to-T variants removed from the CpG sets of
                                                #     step 04 and, after the test, from the tested CpGs and DML (I-08);
                                                #     needs the list of 01_pipeline/optional/ct_variant_loci.R
    dml_difference_standard    = 20L,           # (e) absolute difference a DML must exceed, percentage points: 20 (the
                                                #     paper) or 50 (token "d50"; Venkataraman et al. 2022, getMethylDiff
                                                #     difference = 50). Applied after the test: the fit of the analysis
                                                #     without d50 is reused (author's decision of 2026-11-21)
    dml_qvalue_standard        = 0.01,          # (f) SLIM q-value a DML must be below: 0.01 (the paper) or 0.05
                                                #     (token "q05"; getMethylDiff qvalue = 0.05). Applied after the
                                                #     test: the fit of the analysis without q05 is reused (author's
                                                #     request of 2026-11-22)
    gene_level_go              = FALSE,         # (g) step 10 and the gene-level part of step 11 (I-41)
    diagnostics                = character(0)   # (h) "all", or names of scripts of 06_diagnostics/ without .R
  ),

  # ---- C-to-T variant rule (option (d)) ------------------------------------------
  # The rule of BS-Snper as Venkataraman et al. 2022 applied it, computed in R by
  # 01_pipeline/optional/ct_variant_loci.R (README there). Fixed values, not
  # toggles (decision 10.8 of the settings specification): each other value
  # would need its own analysis token.
  ct_variants = list(
    # BS-Snper rule, 2022 values (Venkataraman et al. 2022: defaults, mincover 5)
    minhetfreq = 0.1, minhomfreq = 0.85, minquali = 15L, mincover = 5L, minread2 = 2L,
    errorate = 0.02, mapvalue = 20L,
    # how the rule is applied, as in 2022 (fixed, not toggles)
    cigar = "full",            # author's decision of 2026-11-20 (the one difference from 2022)
    site = "C",                # the C of the CpG on the plus strand only
    filter = "any",            # PASS and Low records both list a CpG
    union = "pooled_all_and_each_library",
    # assessability: a CT call needs at least 3 Crick T reads (BS-Snper post-rule)
    min_crick_assessable = 3L,
    # computing
    window_bp = 500000L, max_reads_per_batch = 120000L,
    download_streams = 4L, download_range_mb = 256, download_retries = 10L,
    download_range_timeout_s = 1800L, user_agent = "wgbs-pipeline/1.0 (R; curl)"
  ),

  # ---- Step 03: global methylation ------------------------------------------
  elisa_input_ng          = 25,   # DNA mass per well assumed in the kit formula (issue I-12)
  elisa_curve_max_percent = 1,    # standards above this are excluded from the fit as saturating
  # Scale for the models reported as primary. "arcsine" is the manuscript's
  # choice; results on all three scales are always written (issue I-10).
  global_methylation_scale = "arcsine",   # "raw", "arcsine", or "log"

  # ---- Steps 04 and 05: coverage and differential methylation ---------------
  # The five settings marked FROZEN are those of the paper and cannot be changed
  # here (the checks below stop): the values of the manuscript map read them. An
  # optional analysis is chosen with `optional` above or a switch of run_pipeline.R.
  coverage_levels           = c(1, 3, 5, 10, 20),  # thresholds for the CpG count table; 1, 5 and 10 must be among
                                                   # them (step 13 reads the counts at 1 and at the standard)
  min_coverage              = 10,    # FROZEN. Reads required per library for differential methylation in the
                                     # paper; option (b) sets 5 for an optional analysis
  import_min_coverage       = 2,     # FROZEN. Rows with fewer reads are dropped on import in step 05
                                     # (legacy methRead(mincov = 2); sets the base of the percentile)
  upper_coverage_percentile = 99.9,  # FROZEN. Used by option (a) only
  dml_qvalue                = 0.01,  # FROZEN. SLIM-adjusted q-value cutoff (a DML has q below it), in the paper;
                                     # option (f) sets 0.05 for an optional analysis (config$analysis$dml_qvalue)
  dml_difference            = 20,    # FROZEN. A DML needs an absolute methylation difference of more than this, percentage points
                                     # (strict >, as in methylKit getMethylDiff), in the paper; option (e) sets 50 for an
                                     # optional analysis (config$analysis$dml_difference)
  dml_difference_sweep      = seq(5, 50, by = 5),  # thresholds for the sensitivity curve (Figure 2A, 2B)
  # Cores for the per-CpG model fits. Results do not depend on it. The default is
  # 1; set the environment variable WGBS_CORES to a whole number of 1 or more to
  # use more. The workers are forked, which Windows cannot do: there the fits run
  # on 1 core whatever the variable says. The text of the variable is checked and
  # turned into a number in the checks below.
  cores                     = Sys.getenv("WGBS_CORES", unset = "1"),

  # ---- Steps 09 to 11: GO enrichment ------------------------------------------
  go_pvalue      = 0.01,   # Fisher p-value cutoff named in the Methods
  go_legacy_topn = 10,     # the legacy "enriched" set was the top 10 terms per ontology (issue I-21)
  # Ontology release of the legacy enrichment (step 09). The nine archived tables
  # behind the draft are reproduced exactly with this release and with none of the
  # eight other releases tested, 3.7.0 to 3.18.0 (issue I-37). It is installed
  # into a side library in the download cache; other steps use the GO.db of the
  # R installation.
  legacy_go_db = list(
    version = "3.13.0",
    url     = "https://bioconductor.org/packages/3.13/data/annotation/src/contrib/GO.db_3.13.0.tar.gz",
    md5     = "796f6a26af270310d8a1f7c9c3913330"
  ),
  # Gene-level enrichment (step 10)
  go_min_genes = 5,        # a term is tested when at least this many universe genes carry it
  go_fdr       = 0.05,     # Benjamini-Hochberg threshold for calling a term enriched

  # ---- Step 11: semantic similarity of GO terms ---------------------------------
  # rrvgo follows the approach of the REVIGO web server in R; it is not the same
  # algorithm. The organism database supplies the information content of each term.
  semantic_similarity = list(
    measure     = "Rel",           # the SimRel measure named in the Methods
    threshold   = 0.7,             # height at which rrvgo cuts the complete-linkage tree of 1 - similarity:
                                   # terms of one cluster are at least 1 - threshold = 0.3 similar. This is
                                   # rrvgo's default. It is NOT REVIGO's "allowed similarity": REVIGO's 0.7
                                   # corresponds to 0.3 here (issues.md, I-22).
    organism_db = "org.Hs.eg.db",
    n_labels    = 15,              # cluster representatives labelled per plot, chosen by score; at most the
                                   # number of colours in cluster_palette below

    # Facts about the draft and about the legacy files.
    figure_numbers  = c(ploidy = 6L, heat = 7L, shared = 8L),   # figure of the draft for each DML set (manuscript)
    figure_ontology = "BP",                                     # the three figures show biological process only (their captions)
    value_tolerance = 0.025,                                    # a REVIGO value "equals" log10 p within this (p rounded to two digits)
    # Plot scripts exported by REVIGO, with their data embedded, and the Fisher tables of another analysis (%s is the
    # ontology). Paths relative to 03_input_files/; the files are copies of files in 08_archive/, and the README of
    # 03_input_files/annotation/ says where each came from. The names of the exports go into the names of the redrawn figures.
    revigo_exports  = c(code_revigo      = "annotation/legacy/revigo_exports/code_Revigo.R",
                        code_old_revigo  = "annotation/legacy/revigo_exports/code_OLD_Revigo.R",
                        revigo_ploidy_bp = "annotation/legacy/revigo_exports/Revigo-PLOIDY-BP.R"),
    other_fisher    = "annotation/legacy/fisher_other_analysis/ploidy-%s-FisherTestResults.csv",

    # Clusters and tables.
    tie_tolerance    = 1e-9,   # uniqueness values closer than this count as equal when a cluster is named
    order_check_n    = 50L,    # random orders of the terms in which each list is clustered again
    order_check_seed = 11L,    # seed of those orders
    gzip_above_mb    = 5,      # a similarity matrix larger than this (MB) is written as .csv.gz

    # Figures. Sizes in inches unless stated.
    cluster_palette  = c("#2364a0", "#69a500", "#aa2300", "#a082f0", "#0000fa", "#d7737d", "#009bb4", "#87196e", "#006437", "#5f19b9",
                         "#9b6e00", "#d20073", "#4664ff", "#ff46c3", "#af37cd"),   # one colour per labelled cluster, by cluster number
                                                                                   # (see the Notes of step 11)
    grey_other       = "grey72",                # every cluster that is not labelled
    revigo_colours   = c(low = "#b30000", high = "#f7d038"),   # red for low values, yellow for high, as in the archived REVIGO screenshots
    tint_share       = 0.16,                    # share of the cluster colour in the fill of its label; the rest is white
    max_circle       = 9,                       # ggplot2 point size of the circle of the largest cluster of a figure
    circle_stroke    = 0.6,                     # ggplot2 stroke of the circles of the labelled clusters
    panel_margin     = 0.14,                    # empty margin around the terms, as a share of the larger axis range
    panel_in         = 6.5,                     # side of a panel as saved (measured: 6.7 to 6.8); sizes the label boxes when labels are placed
    figure_size      = c(15, 9.6),              # width and height of Figures 6 to 8 (two panels)
    gene_figure_size = c(7.6, 9.8),             # width and height of a gene-level figure (one panel)
    export_size      = c(10.5, 8),              # width and height of a redrawn REVIGO export
    export_panel_in  = c(8, 6.6),               # width and height of its panel as saved (measured: 8.3 to 8.8, and 6.9)
    export_margin    = c(0.08, 0.06),           # empty margin around its bubbles, as a share of the range of each axis
    export_bubble    = c(1.2, 9),               # ggplot2 point size of its smallest and largest bubble
    export_stroke    = 0.4,                     # ggplot2 stroke of its bubbles
    label_seed       = 11L,                     # seed of ggrepel
    label_mm         = c(cluster = 2.5, export = 2.4),   # font size of the labels, mm
    label_width      = 28L,                     # characters per line in a label
    subtitle_width   = 90L,                     # characters per line in a panel subtitle
    caption_per_in   = 16,                      # characters of caption per inch of figure width

    # Label placement by rule (place_labels() of step 11). Lengths and costs in inches.
    label_angles     = seq(0, 350, by = 10) * pi / 180,   # directions tried around a circle
    label_rings      = c(0.08, 0.15, 0.22, 0.3, 0.4, 0.5, 0.62, 0.75, 0.9, 1.05, 1.2, 1.4, 1.6, 1.8, 2, 2.2),   # distances tried between the
                                                # circle and the label box
    label_gap        = 0.05,                    # clearance kept around circles and label boxes
    crowd_radius     = 1,                       # labelled circles within this distance make a place crowded (one of the orders of placement)
    charge_circle    = 0.33,                    # cost of an unlabelled circle under a label box
    charge_leader    = 0.5,                     # cost of a labelled circle or a label under a leader line, and of two leader lines that cross
    charge_turn      = 0.06,                    # cost of pointing towards the middle of the labelled circles and not away from it
    leader_steps     = seq(0.05, 0.95, by = 0.05),   # places along a leader line that are tested
    label_padding    = 0.15,                    # padding inside a label box, in lines of 14.4 points (ggrepel label.padding)
    label_box_em     = c(char = 0.57, first_line = 0.75, next_line = 1.08)   # box per character, for one line and per further line, in font
                                                # sizes; measured on PNG and PDF: 0.40 to 0.56, 0.71 to 0.73, 1.08
  ),

  # ---- Step 12: qPCR ------------------------------------------------------------
  qpcr = list(
    reference_gene   = "Actin",    # the only reference gene that was measured
    max_cycles       = 40,         # a well with no Cq by this cycle is a non-detect
    legacy_nondetect = 45,         # the Cq the legacy scripts substituted for a non-detect
    mito_gene        = "COX1",     # measured without a reference gene; never normalized in step 12
    lod_detected     = 0.95,       # a gene is "quantifiable" when at least this share of its wells is detected
    min_detected     = 0.5,        # a gene is not analysed when a group has a smaller share of detected wells
    failed_fence_iqr = 3,          # reference Cq above Q3 + this many interquartile ranges marks a failed sample
    conf_level       = 0.95,       # confidence level of the intervals
    alpha            = 0.05,       # unadjusted p value used in the verdicts
    fdr              = 0.05,       # Benjamini-Hochberg q value used in the verdicts
    cq_tolerance     = 0.0005,     # two Cq values agree when they differ by no more than half of the third decimal
    same_p_tolerance = 1e-6,       # two p values closer than this, relative to the larger one, are one result
    size_check_runs  = 2000,       # simulated data sets per share of censored animals, check of the censored test
    size_check_seed  = 20261119,   # seed of that simulation
    # Facts about the files, from 03_input_files/qpcr/README.md. Step 12 stops if the data disagree.
    design = list(wells_per_animal = 2, animals_per_group = 8, acute_per_group = 2, nondetect_text = "N/A")
  ),

  # ---- Plot colours ----------------------------------------------------------
  colours = list(
    ploidy    = c(diploid = "royalblue1", triploid = "orangered1"),
    treatment = c(control = "grey35", heat = "firebrick2"),
    dml_set   = c(ploidy = "royalblue2", heat = "orangered2", shared = "seagreen3")
  )
)

# Sequence names of the ten linkage groups of the Roslin v1 assembly, in order.
# Everything else in the coverage files is an unplaced scaffold (NW_...) or the
# mitochondrial genome (NC_001276.1).
config$linkage_groups <- sprintf("NC_0475%02d.1", 59:68)
config$mito_sequence  <- "NC_001276.1"

# ---- Checks of the settings ------------------------------------------------------
# Every step sources this file, so every step passes through this block. A value
# that cannot work stops the step here with a message that names the setting;
# without the block such a value fails late (WGBS_CORES after step 05 has
# imported the 22 libraries) or not at all. The checks do not change a value,
# except that cores becomes a number. They add three derived values:
#   config$analysis  the analysis of the run (analysis_spec() of
#                    02_R_functions/analyses.R): id, filter, min_coverage,
#                    overdispersion, test, ct_excluded, dml_difference,
#                    dml_qvalue, fit_id
#   config$addons    gene_level_go (TRUE/FALSE) and diagnostics (names)
#   config$run_kind  "paper" when the analysis is the paper's and no add-on is on,
#                    else "optional"
config <- local({
  cfg <- config
  bad <- function(...) stop("run_config.R: ", ..., call. = FALSE)
  shown <- function(x) paste(deparse(x), collapse = "")
  # `name` is a setting as written after "config$"; `ok` is what must hold for its value x.
  check <- function(name, must, ok) {
    x <- eval(parse(text = paste0("cfg$", name)))
    if (!isTRUE(ok(x))) bad("config$", name, " must be ", must, "; it is ", shown(x), ".")
  }
  numbers <- function(x) is.numeric(x) && length(x) >= 1 && !anyNA(x) && all(is.finite(x))
  number  <- function(x) numbers(x) && length(x) == 1
  whole   <- function(x) number(x) && x == round(x) && x >= 1
  share   <- function(x) number(x) && x > 0 && x <= 1                 # a cutoff on a p-value or q-value
  text    <- function(x) is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
  seed    <- function(x) number(x) && x == round(x) && abs(x) <= .Machine$integer.max   # what set.seed() takes
  at_least_0 <- function(x) number(x) && x >= 0                       # a tolerance or a size limit

  # Optional analyses (settings (a) to (h)) and the environment variables that
  # choose them for a step started alone. The analysis id is defined by
  # analysis_spec() of 02_R_functions/analyses.R.
  source(here::here("02_R_functions", "analyses.R"), local = TRUE)
  flag <- function(x) is.logical(x) && length(x) == 1 && !is.na(x)
  check("optional$coverage_filter_of_methods", "TRUE or FALSE", flag)
  check("optional$overdispersion_chisq", "TRUE or FALSE", flag)
  check("optional$exclude_ct_variants", "TRUE or FALSE", flag)
  check("optional$gene_level_go", "TRUE or FALSE", flag)
  check("optional$coverage_standard", "10 (the paper) or 5", function(x) number(x) && x %in% c(5, 10))
  check("optional$dml_difference_standard", "20 (the paper) or 50", function(x) number(x) && x %in% c(20, 50))
  check("optional$dml_qvalue_standard", "0.01 (the paper) or 0.05", function(x) number(x) && x %in% c(0.01, 0.05))
  diag_table <- here::here("06_diagnostics", "diagnostics.csv")
  diag_names <- if (file.exists(diag_table)) utils::read.csv(diag_table, stringsAsFactors = FALSE)$script else character(0)
  diag_ok <- function(x) is.character(x) && !anyNA(x) && (length(x) == 0 || identical(x, "all") || (all(x %in% diag_names) && !anyDuplicated(x)))
  check("optional$diagnostics", paste0('character(0), "all", or names of scripts listed in 06_diagnostics/diagnostics.csv (',
                                       paste(diag_names, collapse = ", "), ")"), diag_ok)
  o <- cfg$optional
  id <- analysis_id(c(if (o$coverage_filter_of_methods) "corrected", if (o$coverage_standard == 5) "cov5",
                      if (o$overdispersion_chisq) "mnchisq", if (o$exclude_ct_variants) "ctexcl",
                      if (o$dml_difference_standard == 50) "d50", if (o$dml_qvalue_standard == 0.05) "q05"))
  gene_level <- o$gene_level_go
  diagnostics <- o$diagnostics
  # An empty variable is an error, not "not set": it is what a job script leaves
  # behind when the value it was meant to pass is missing.
  env <- function(name) { v <- Sys.getenv(name, unset = NA); if (!is.na(v) && !nzchar(trimws(v))) bad("the environment variable ", name, " is set but empty; unset it or give a value."); v }
  v <- env("WGBS_ANALYSIS")
  if (!is.na(v)) {
    id <- tryCatch(analysis_spec(trimws(v))$id, error = function(e) bad("the environment variable WGBS_ANALYSIS: ", conditionMessage(e)))
  }
  v <- env("WGBS_GENE_LEVEL_GO")
  if (!is.na(v)) {
    if (!trimws(v) %in% c("TRUE", "FALSE")) bad("the environment variable WGBS_GENE_LEVEL_GO must be TRUE or FALSE; it is '", v, "'.")
    gene_level <- trimws(v) == "TRUE"
  }
  v <- env("WGBS_DIAGNOSTICS")
  if (!is.na(v)) {
    d <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
    diagnostics <- if (identical(d, "none")) character(0) else d
    if (!diag_ok(diagnostics)) bad("the environment variable WGBS_DIAGNOSTICS must be none, all, or names of scripts listed in ",
                                   "06_diagnostics/diagnostics.csv joined by commas; it is '", v, "'.")
  }
  cfg$analysis <- analysis_spec(id)
  cfg$addons <- list(gene_level_go = gene_level, diagnostics = diagnostics)
  cfg$run_kind <- if (cfg$analysis$paper && !gene_level && length(diagnostics) == 0) "paper" else "optional"

  # The settings of the paper are frozen (see the settings above).
  frozen <- c(min_coverage = 10, import_min_coverage = 2, upper_coverage_percentile = 99.9, dml_qvalue = 0.01, dml_difference = 20)
  for (nm in names(frozen)) {
    if (!isTRUE(all.equal(cfg[[nm]], frozen[[nm]])) || length(cfg[[nm]]) != 1) {
      bad("config$", nm, " is a setting of the paper and must stay ", frozen[[nm]], " (it is ", shown(cfg[[nm]]),
          "). An optional analysis is chosen with config$optional or a switch of run_pipeline.R.")
    }
  }

  # Cores: the text of WGBS_CORES, or "1" when the variable is not set. An empty
  # variable is an error, not the default: it is what a job script leaves behind
  # when the value it was meant to pass is missing.
  cores_as_set <- paste(cfg$cores, collapse = " ")     # one text, whatever the setting holds
  cores <- trimws(cores_as_set)
  if (!grepl("^[0-9]+$", cores) || as.numeric(cores) < 1 || as.numeric(cores) > .Machine$integer.max) {
    bad("the environment variable WGBS_CORES must be a whole number of 1 or more, or not be set at all (1 core); it is '",
        cores_as_set, "'.")
  }
  cfg$cores <- as.integer(cores)
  if (.Platform$OS.type == "windows" && cfg$cores > 1L) {
    message("WGBS_CORES=", cfg$cores, " is not used on Windows, where the model fits cannot fork workers; using 1 core.")
    cfg$cores <- 1L
  }

  # Thresholds and the other numbers of the analysis.
  check("download_timeout_s", "a number of seconds above 0", function(x) number(x) && x > 0)
  check("elisa_input_ng", "a number above 0", function(x) number(x) && x > 0)
  check("elisa_curve_max_percent", "a number above 0", function(x) number(x) && x > 0)
  check("global_methylation_scale", '"raw", "arcsine" or "log"', function(x) text(x) && x %in% c("raw", "arcsine", "log"))
  check("min_coverage", "a whole number of 1 or more", whole)
  check("import_min_coverage", "a whole number from 1 to config$min_coverage", function(x) whole(x) && x <= cfg$min_coverage)
  check("coverage_levels", "whole numbers of 1 or more in rising order, with 1, 5 and 10 among them",
        function(x) numbers(x) && all(x == round(x) & x >= 1) && all(diff(x) > 0) && all(c(1, 5, 10) %in% x))
  # A percentile typed as a share (0.999 for 99.9) would remove nearly every CpG of option (a) without an error.
  check("upper_coverage_percentile", "a percentile above 50 and at most 100 (99.9 for the 99.9th percentile, not 0.999)",
        function(x) number(x) && x > 50 && x <= 100)
  check("dml_qvalue", "a number above 0 and at most 1", share)
  # A difference typed as a share (0.2 for 20 points) would give DML lists many times longer without an error.
  points <- function(x) x == 0 | (x >= 1 & x < 100)
  check("dml_difference", "0 or a number of percentage points from 1 to below 100 (20 for 20 points, not 0.2)",
        function(x) number(x) && points(x))
  check("dml_difference_sweep", "numbers of percentage points, each 0 or from 1 to below 100, in rising order",
        function(x) numbers(x) && all(points(x)) && all(diff(x) > 0))
  check("go_pvalue", "a number above 0 and at most 1", share)
  check("go_fdr", "a number above 0 and at most 1", share)
  check("go_legacy_topn", "a whole number of 1 or more", whole)
  check("go_min_genes", "a whole number of 1 or more", whole)
  check("legacy_go_db$version", "the version of a GO.db release, for example \"3.13.0\"", function(x) text(x) && grepl("^[0-9]+([.][0-9]+)+$", x))
  check("legacy_go_db$url", paste0("the address of the source package of that release, ending in GO.db_", cfg$legacy_go_db$version, ".tar.gz"),
        function(x) text(x) && grepl("^https?://", x) && endsWith(x, paste0("/GO.db_", cfg$legacy_go_db$version, ".tar.gz")))
  check("legacy_go_db$md5", "the md5 checksum of that package, 32 hexadecimal digits", function(x) text(x) && grepl("^[0-9a-f]{32}$", x))
  check("semantic_similarity$measure", "a measure of GOSemSim that rrvgo offers: \"Rel\", \"Resnik\", \"Lin\", \"Jiang\" or \"Wang\"",
        function(x) text(x) && x %in% c("Rel", "Resnik", "Lin", "Jiang", "Wang"))
  check("semantic_similarity$threshold", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("semantic_similarity$organism_db", "the name of an organism annotation package, for example \"org.Hs.eg.db\"", text)
  check("semantic_similarity$n_labels",
        paste("a whole number from 1 to", length(cfg$semantic_similarity$cluster_palette), "(the number of colours in cluster_palette)"),
        function(x) whole(x) && x <= length(cfg$semantic_similarity$cluster_palette))
  check("semantic_similarity$order_check_n", "a whole number of 1 or more", whole)
  check("semantic_similarity$order_check_seed", "a whole number", seed)
  check("semantic_similarity$label_seed", "a whole number", seed)
  check("semantic_similarity$value_tolerance", "a number of 0 or more", at_least_0)
  check("semantic_similarity$tie_tolerance", "a number of 0 or more", at_least_0)
  check("semantic_similarity$gzip_above_mb", "a number of megabytes, 0 or more", at_least_0)
  check("qpcr$reference_gene", "the name of the reference gene", text)
  check("qpcr$max_cycles", "a number of cycles above 0", function(x) number(x) && x > 0)
  check("qpcr$legacy_nondetect", "a number of cycles above 0", function(x) number(x) && x > 0)
  check("qpcr$lod_detected", "a share of the wells, from 0 to 1", function(x) number(x) && x >= 0 && x <= 1)
  check("qpcr$min_detected", "a share of the wells, from 0 to 1", function(x) number(x) && x >= 0 && x <= 1)
  check("qpcr$conf_level", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("qpcr$alpha", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("qpcr$fdr", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("qpcr$size_check_runs", "a whole number of 1 or more", whole)
  check("qpcr$size_check_seed", "a whole number", seed)
  check("qpcr$failed_fence_iqr", "a number of interquartile ranges, 0 or more", at_least_0)
  check("qpcr$cq_tolerance", "a number of cycles, 0 or more", at_least_0)
  check("qpcr$same_p_tolerance", "a number of 0 or more", at_least_0)
  cfg
})

# ---- Number formatting of the R session ------------------------------------------
# A personal R profile can change how R turns numbers into text, most often with
# options(scipen = 999), and with that the bytes of every table a script writes
# with fwrite() or builds with paste() and format(). Tables are compared byte for
# byte after a rerun, so the options that decide the text of a number or of a
# logical value are set to their defaults here, for the session of every script
# that sources this file. write_table() fixes its own arguments as well
# (02_R_functions/logging.R); the scripts in 06_diagnostics/ call fwrite()
# directly and rely on this line.
options(scipen = 0, digits = 7, OutDec = ".", datatable.logical01 = FALSE)

# run_config.R
#
# Single control surface for the pipeline. Every script in 01_pipeline/ sources
# this file first. Change analysis settings here, never inside a step.
#
# Where a setting has a "legacy" and a "corrected" value, legacy is what was run
# for the manuscript draft of November 2026 and corrected is what the Methods
# describe. The issue numbers refer to issues.md.
#
# The file has three parts: the settings, the checks of the settings, and the
# number formatting of the R session. A setting that cannot work stops every
# step in the second part, before the step reads or writes anything.

config <- list(

  # ---- Inputs on gannet ------------------------------------------------------
  # Files listed in 03_input_files/gannet_manifest.csv are downloaded here.
  # Default: 03_input_files/gannet_cache/ (git-ignored). Set the environment
  # variable WGBS_CACHE_DIR to use another location, for example a larger drive.
  cache_dir          = Sys.getenv("WGBS_CACHE_DIR", unset = ""),
  download_timeout_s = 7200,

  # ---- Analysis modes --------------------------------------------------------
  # One or both of "legacy" and "corrected". Steps 04 to 08 and step 10 write one
  # set of outputs per mode listed here, with the mode in each file name.
  #   legacy     what the legacy scripts did: no upper coverage filter (the call
  #              misnamed the argument, so the filter never ran; issue I-05), and
  #              the first line of every coverage file lost on import (methRead()
  #              was left at header = TRUE; issue I-38)
  #   corrected  CpGs at or above the 99.9th percentile of coverage in a library
  #              are removed from that library, as the Methods describe, and
  #              every line of every coverage file is read
  # Legacy mode is the analysis of the draft. Steps 09, 11 and 13 repeat that
  # analysis or compare the draft with it, and run only when "legacy" is listed:
  # without it each of them says so and ends with status 0, run_pipeline.R leaves
  # out step 10 as well (its retest of the terms of the draft reads the
  # legacy-mode tables), and annotation_summary.csv of step 08 holds no legacy
  # row. Without "corrected", step 13 leaves the corrected-mode columns of its
  # tables empty. Steps 01 to 03 and 12 do not depend on the mode.
  modes = c("legacy", "corrected"),

  # ---- Step 03: global methylation ------------------------------------------
  elisa_input_ng          = 25,   # DNA mass per well assumed in the kit formula (issue I-12)
  elisa_curve_max_percent = 1,    # standards above this are excluded from the fit as saturating
  # Scale for the models reported as primary. "arcsine" is the manuscript's
  # choice; results on all three scales are always written (issue I-10).
  global_methylation_scale = "arcsine",   # "raw", "arcsine", or "log"

  # ---- Steps 04 and 05: coverage and differential methylation ---------------
  coverage_levels           = c(1, 3, 5, 10, 20),  # thresholds for the CpG count table; 1 and min_coverage
                                                   # must be among them (step 13 reads the counts at both)
  min_coverage              = 10,    # reads required per library for differential methylation
  import_min_coverage       = 2,     # rows with fewer reads are dropped on import in step 05
                                     # (legacy methRead(mincov = 2); sets the base of the percentile)
  upper_coverage_percentile = 99.9,  # used in corrected mode only
  dml_qvalue                = 0.01,  # SLIM-adjusted q-value cutoff
  dml_difference            = 20,    # minimum absolute methylation difference, percentage points
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
# imported the 22 libraries) or not at all (an unknown mode writes tables that
# no later step reads). The checks do not change a value, with two exceptions:
# the modes are put in the order legacy, corrected, and cores becomes a number.
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

  # Modes. Tables that span the modes list legacy first, whatever the order typed.
  all_modes <- c("legacy", "corrected")
  check("modes", 'one or both of "legacy" and "corrected", each named once',
        function(x) is.character(x) && length(x) >= 1 && !anyNA(x) && all(x %in% all_modes) && !anyDuplicated(x))
  cfg$modes <- intersect(all_modes, cfg$modes)

  # Cores: the text of WGBS_CORES, or "1" when the variable is not set. An empty
  # variable is an error, not the default: it is what a job script leaves behind
  # when the value it was meant to pass is missing.
  cores <- trimws(cfg$cores)
  if (!grepl("^[0-9]+$", cores) || as.numeric(cores) < 1 || as.numeric(cores) > .Machine$integer.max) {
    bad("the environment variable WGBS_CORES must be a whole number of 1 or more, or not be set at all (1 core); it is '",
        cfg$cores, "'.")
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
  check("coverage_levels", "whole numbers of 1 or more in rising order, with 1 and config$min_coverage among them",
        function(x) numbers(x) && all(x == round(x) & x >= 1) && all(diff(x) > 0) && all(c(1, cfg$min_coverage) %in% x))
  check("upper_coverage_percentile", "a percentile above 0 and at most 100", function(x) number(x) && x > 0 && x <= 100)
  check("dml_qvalue", "a number above 0 and at most 1", share)
  check("dml_difference", "a number of percentage points from 0 to below 100", function(x) number(x) && x >= 0 && x < 100)
  check("dml_difference_sweep", "numbers of percentage points from 0 to below 100, in rising order",
        function(x) numbers(x) && all(x >= 0 & x < 100) && all(diff(x) > 0))
  check("go_pvalue", "a number above 0 and at most 1", share)
  check("go_fdr", "a number above 0 and at most 1", share)
  check("go_legacy_topn", "a whole number of 1 or more", whole)
  check("go_min_genes", "a whole number of 1 or more", whole)
  check("legacy_go_db$version", "the version of a GO.db release, for example \"3.13.0\"", function(x) text(x) && grepl("^[0-9]+([.][0-9]+)+$", x))
  check("legacy_go_db$url", "the address of the source package of that release", function(x) text(x) && grepl("^https?://.+[.]tar[.]gz$", x))
  check("legacy_go_db$md5", "the md5 checksum of that package, 32 hexadecimal digits", function(x) text(x) && grepl("^[0-9a-f]{32}$", x))
  check("semantic_similarity$measure", "a measure of GOSemSim that rrvgo offers: \"Rel\", \"Resnik\", \"Lin\", \"Jiang\" or \"Wang\"",
        function(x) text(x) && x %in% c("Rel", "Resnik", "Lin", "Jiang", "Wang"))
  check("semantic_similarity$threshold", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("semantic_similarity$organism_db", "the name of an organism annotation package, for example \"org.Hs.eg.db\"", text)
  check("semantic_similarity$n_labels",
        paste("a whole number from 1 to", length(cfg$semantic_similarity$cluster_palette), "(the number of colours in cluster_palette)"),
        function(x) whole(x) && x <= length(cfg$semantic_similarity$cluster_palette))
  check("semantic_similarity$order_check_n", "a whole number of 1 or more", whole)
  check("qpcr$reference_gene", "the name of the reference gene", text)
  check("qpcr$max_cycles", "a number of cycles above 0", function(x) number(x) && x > 0)
  check("qpcr$legacy_nondetect", "a number of cycles above 0", function(x) number(x) && x > 0)
  check("qpcr$lod_detected", "a share of the wells, from 0 to 1", function(x) number(x) && x >= 0 && x <= 1)
  check("qpcr$min_detected", "a share of the wells, from 0 to 1", function(x) number(x) && x >= 0 && x <= 1)
  check("qpcr$conf_level", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("qpcr$alpha", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("qpcr$fdr", "a number above 0 and below 1", function(x) number(x) && x > 0 && x < 1)
  check("qpcr$size_check_runs", "a whole number of 1 or more", whole)
  cfg
})

# ---- Number formatting of the R session ------------------------------------------
# A personal R profile can change how R turns numbers into text, most often with
# options(scipen = 999), and with that the bytes of every table a script writes
# with fwrite() or builds with paste() and format(). Tables are compared byte for
# byte after a rerun, so the three options that decide the text of a number are
# set to the defaults of R here, for the session of every script that sources
# this file. write_table() fixes its own arguments as well (02_R_functions/logging.R).
options(scipen = 0, digits = 7, OutDec = ".")

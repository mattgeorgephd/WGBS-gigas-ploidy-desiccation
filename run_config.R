# run_config.R
#
# Single control surface for the pipeline. Every script in 01_pipeline/ sources
# this file first. Change analysis settings here, never inside a step.
#
# Where a setting has a "legacy" and a "corrected" value, legacy is what was run
# for the manuscript draft of November 2026 and corrected is what the Methods
# describe or what the audit recommends. The issue numbers refer to issues.md.

config <- list(

  # ---- Inputs on gannet ------------------------------------------------------
  # Files listed in 03_input_files/gannet_manifest.csv are downloaded here.
  # Default: 03_input_files/gannet_cache/ (git-ignored). Set the environment
  # variable WGBS_CACHE_DIR to use another location, for example a larger drive.
  cache_dir          = Sys.getenv("WGBS_CACHE_DIR", unset = ""),
  download_timeout_s = 7200,

  # ---- Analysis modes --------------------------------------------------------
  # Steps 04 onward write one set of outputs per mode listed here, with the mode
  # in each file name.
  #   legacy     no upper coverage filter (the legacy call misnamed the argument,
  #              so the filter never ran; issue I-05)
  #   corrected  CpGs at or above the 99.9th percentile of coverage in a library
  #              are removed from that library, as the Methods describe
  modes = c("legacy", "corrected"),

  # ---- Step 03: global methylation ------------------------------------------
  elisa_input_ng          = 25,   # DNA mass per well assumed in the kit formula (issue I-12)
  elisa_curve_max_percent = 1,    # standards above this are excluded from the fit as saturating
  # Scale for the models reported as primary. "arcsine" is the manuscript's
  # choice; results on all three scales are always written (issue I-10).
  global_methylation_scale = "arcsine",   # "raw", "arcsine", or "log"

  # ---- Steps 04 and 05: coverage and differential methylation ---------------
  coverage_levels           = c(1, 3, 5, 10, 20),  # thresholds for the CpG count table
  min_coverage              = 10,    # reads required per library for differential methylation
  import_min_coverage       = 2,     # rows with fewer reads are dropped on import in step 05
                                     # (legacy methRead(mincov = 2); sets the base of the percentile)
  upper_coverage_percentile = 99.9,  # used in corrected mode only
  dml_qvalue                = 0.01,  # SLIM-adjusted q-value cutoff
  dml_difference            = 20,    # minimum absolute methylation difference, percentage points
  dml_difference_sweep      = seq(5, 50, by = 5),  # thresholds for the sensitivity curve (Figure 2A, 2B)
  cores                     = 1,     # cores for the per-CpG model fits; results do not depend on it

  # ---- Step 09: GO enrichment ------------------------------------------------
  go_pvalue      = 0.01,   # Fisher p-value cutoff named in the Methods
  go_legacy_topn = 10,     # the legacy "enriched" set was the top 10 terms per ontology (issue I-21)

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

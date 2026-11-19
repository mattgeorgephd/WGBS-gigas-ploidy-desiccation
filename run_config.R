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
  #   legacy     what the legacy scripts did: no upper coverage filter (the call
  #              misnamed the argument, so the filter never ran; issue I-05), and
  #              the first line of every coverage file lost on import (methRead()
  #              was left at header = TRUE; issue I-38)
  #   corrected  CpGs at or above the 99.9th percentile of coverage in a library
  #              are removed from that library, as the Methods describe, and
  #              every line of every coverage file is read
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
  # Cores for the per-CpG model fits (forked workers; leave at 1 on Windows). Results
  # do not depend on it. Set the environment variable WGBS_CORES to override.
  cores                     = max(1L, as.integer(Sys.getenv("WGBS_CORES", unset = "1"))),

  # ---- Steps 09 to 11: GO enrichment ------------------------------------------
  go_pvalue      = 0.01,   # Fisher p-value cutoff named in the Methods
  go_legacy_topn = 10,     # the legacy "enriched" set was the top 10 terms per ontology (issue I-21)
  # Ontology release of the legacy enrichment (step 09). The archived tables are
  # reproduced exactly with this release and with no later one. It is installed
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

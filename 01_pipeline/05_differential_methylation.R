# 05_differential_methylation.R
#
# Purpose   Differentially methylated loci (DML) for ploidy and for heat, with
#           methylKit, on the 22 libraries. For each analysis mode in
#           config$modes (see run_config.R):
#             1. import each library's coverage file; keep CpGs with at least
#                config$min_coverage reads (and, in corrected mode, below the
#                library's 99.9th percentile of coverage)
#             2. keep the CpGs present in all 22 libraries and normalize coverage
#                between libraries (methylKit's median scaling)
#             3. per CpG, logistic regression of methylation on ploidy with heat
#                as covariate, and on heat with ploidy as covariate; SLIM q-values
#             4. DML = q < config$dml_qvalue and |difference| > config$dml_difference
#             5. loci in both DML lists ("shared"; the legacy "interaction-DML")
#
# Inputs    03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files listed in 03_input_files/gannet_manifest.csv (cache)
# Outputs   04_output/05_differential_methylation/tables/
#             library_filtering_<mode>.csv       per-library CpG counts, cutoffs, scaling
#             united_cpg_<mode>.csv.gz           read counts at every analysed CpG, all libraries
#             diffmeth_<contrast>_<mode>.csv.gz  p, q, difference at every analysed CpG
#             dml_<contrast>_<mode>.csv          DML at the configured thresholds
#             dml_shared_<mode>.csv              loci in both DML lists, with both differences
#             dml_counts_by_threshold_<mode>.csv DML counts at 5% to 50% difference
#             dml_counts_by_mode.csv             headline counts, modes side by side
#           04_output/05_differential_methylation/bed/dml_<set>_<mode>.bed
#           04_output/05_differential_methylation/cache/methylkit_<mode>.rds   (not tracked)
#
# Manuscript items: Results 3.2 (CpG and DML counts); data for Figures 2, 3, 4.
# Notes     - Steps 1 and 2 follow methylKit's rules (methRead, filterByCoverage,
#             normalizeCoverage, unite) but are applied one library at a time
#             with data.table, so the step needs about 3 GB instead of 9 GB; see
#             02_R_functions/methylkit_helpers.R. The tests are methylKit's own.
#           - Differences are group 1 minus group 0: triploid minus diploid for
#             ploidy, heat minus control for heat. "hyper" means higher in
#             triploid, or higher in heat.
#           - The heat contrast compares two cohorts that were sequenced
#             separately (issue I-01). The shared set is an overlap of two
#             lists, not a test of interaction (issue I-02).
#           - Fitting takes about 20 minutes per contrast on one core. Results
#             are cached; delete the cache file to refit.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))

step <- "05_differential_methylation"
t0 <- start_step(step, "differential methylation for ploidy and heat")
lib <- read_libraries()
treat_ploidy <- as.integer(lib$ploidy == "triploid")     # 0 diploid, 1 triploid
treat_heat   <- as.integer(lib$treatment == "heat")      # 0 control, 1 heat

# ---- Fit one mode (or load it from the cache) ------------------------------------
fit_mode <- function(mode) {
  cache_file <- file.path(step_cache_dir(step), sprintf("methylkit_%s.rds", mode))
  if (file.exists(cache_file)) {
    log_msg("[%s] loading cached fit %s", mode, basename(cache_file))
    return(readRDS(cache_file))
  }
  # Pass over the libraries: filter each one, keep its surviving CpGs as
  # (position key, counts), and narrow the set of CpGs present in all libraries.
  kept <- vector("list", nrow(lib)); records <- vector("list", nrow(lib))
  chr_levels <- character(0); common <- NULL
  for (i in seq_len(nrow(lib))) {
    r <- read_filtered_counts(fetch_gannet_file(lib$cov_file[i], config), mode, config, chr_levels)
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
  objects <- lapply(seq_len(nrow(lib)), function(i) {
    k <- kept[[i]][match(common, pos_key)]
    make_methylraw(pos$chr, pos$start, scale_counts_like_methylkit(k$coverage, k$nC, k$nT, filtering$scaling_factor[i]),
                   as.character(lib$legacy_sample_index[i]))
  })
  rm(kept, records); gc(verbose = FALSE)
  raw_list <- new("methylRawList", objects, treatment = treat_ploidy)
  meth_ploidy <- unite(raw_list, destrand = FALSE)        # same CpGs in every object; sets the row order
  rm(objects, raw_list); gc(verbose = FALSE)
  meth_heat <- reorganize(meth_ploidy, sample.ids = meth_ploidy@sample.ids, treatment = treat_heat)
  log_msg("[%s] %d CpGs with at least %dx in all %d libraries", mode, nrow(meth_ploidy),
          config$min_coverage, nrow(lib))

  log_msg("[%s] fitting ploidy (covariate: heat) on %d core(s)", mode, config$cores)
  dm_ploidy <- calculateDiffMeth(meth_ploidy, covariates = data.frame(heat = lib$treatment),
                                 mc.cores = config$cores)
  log_msg("[%s] fitting heat (covariate: ploidy)", mode)
  dm_heat <- calculateDiffMeth(meth_heat, covariates = data.frame(ploidy = lib$ploidy),
                               mc.cores = config$cores)
  fit <- list(mode = mode, filtering = filtering, united = meth_ploidy,
              diffmeth = list(ploidy = dm_ploidy, heat = dm_heat))
  saveRDS(fit, cache_file)
  fit
}

# ---- Write the outputs of one mode -----------------------------------------------
write_mode <- function(fit) {
  mode <- fit$mode
  tab <- function(name) output_path(step, "tables", name)

  write_table(fit$filtering, tab(sprintf("library_filtering_%s.csv", mode)))

  united <- as.data.table(getData(fit$united))
  write_table(united, tab(sprintf("united_cpg_%s.csv.gz", mode)))

  dml <- list()
  for (contrast in names(fit$diffmeth)) {
    all_cpg <- as.data.table(getData(fit$diffmeth[[contrast]]))
    stopifnot(identical(all_cpg$chr, united$chr), identical(all_cpg$start, united$start))
    all_cpg[, united_row := .I]                 # row in united_cpg_<mode>; the legacy row name
    write_table(all_cpg, tab(sprintf("diffmeth_%s_%s.csv.gz", contrast, mode)))
    d <- all_cpg[qvalue < config$dml_qvalue & abs(meth.diff) > config$dml_difference]
    stopifnot(nrow(d) == nrow(dml_table(fit$diffmeth[[contrast]], config$dml_difference, config$dml_qvalue)))
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

  # BED files (0-based start, as in the coverage files): chr, start, end, difference.
  bed_dir <- output_dir(step, "bed")
  write_bed <- function(d, diff_col, name) {
    fwrite(d[, .(chr, start, end, get(diff_col))], file.path(bed_dir, sprintf("dml_%s_%s.bed", name, mode)),
           sep = "\t", col.names = FALSE, scipen = 50)
  }
  write_bed(dml$ploidy, "meth.diff", "ploidy"); write_bed(dml$heat, "meth.diff", "heat")
  write_bed(shared, "meth.diff_heat", "shared")

  # DML counts across difference thresholds (Figure 2A, 2B).
  sweep <- rbindlist(lapply(config$dml_difference_sweep, function(th) {
    p <- as.data.table(getData(fit$diffmeth$ploidy))[qvalue < config$dml_qvalue & abs(meth.diff) > th]
    h <- as.data.table(getData(fit$diffmeth$heat))[qvalue < config$dml_qvalue & abs(meth.diff) > th]
    s <- merge(p, h, by = c("chr", "start"), suffixes = c("_ploidy", "_heat"))
    data.table(mode = mode, difference = th, set = c("ploidy", "heat", "shared"),
               n = c(nrow(p), nrow(h), nrow(s)),
               n_hyper = c(sum(p$meth.diff > 0), sum(h$meth.diff > 0), sum(s$meth.diff_heat > 0)),
               n_hypo  = c(sum(p$meth.diff < 0), sum(h$meth.diff < 0), sum(s$meth.diff_heat < 0)))
  }))
  sweep[, pct_hyper := fifelse(n > 0, 100 * n_hyper / n, NA_real_)]
  write_table(sweep, tab(sprintf("dml_counts_by_threshold_%s.csv", mode)))

  data.table(mode = mode, cpg_analysed = nrow(united),
             ploidy_dml = nrow(dml$ploidy), ploidy_hyper = sum(dml$ploidy$direction == "hyper"),
             ploidy_hypo = sum(dml$ploidy$direction == "hypo"),
             heat_dml = nrow(dml$heat), heat_hyper = sum(dml$heat$direction == "hyper"),
             heat_hypo = sum(dml$heat$direction == "hypo"),
             shared = nrow(shared), shared_both_hyper = sum(shared$class == "both_hyper"),
             shared_both_hypo = sum(shared$class == "both_hypo"), shared_opposite = sum(shared$class == "opposite"))
}

# ---- Run ----------------------------------------------------------------------------
headline <- rbindlist(lapply(config$modes, function(mode) {
  fit <- fit_mode(mode)
  out <- write_mode(fit)
  rm(fit); gc(verbose = FALSE)
  out
}))
write_table(headline, output_path(step, "tables", "dml_counts_by_mode.csv"))
print(t(headline))

finish_step(step, t0)

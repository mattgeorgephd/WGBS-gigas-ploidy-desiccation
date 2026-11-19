# 05_differential_methylation.R
#
# Purpose   Differentially methylated loci (DML) for ploidy and for heat, with
#           methylKit, on the 22 libraries. For each analysis mode in
#           config$modes (see run_config.R):
#             1. import each library's coverage file; keep CpGs with at least
#                config$min_coverage reads (and, in corrected mode, with fewer
#                reads than the library's 99.9th percentile of coverage, taken
#                over its CpGs with at least config$import_min_coverage reads)
#             2. keep the CpGs present in all 22 libraries and normalize coverage
#                between libraries (methylKit's median scaling)
#             3. per CpG, logistic regression of methylation on ploidy with heat
#                as covariate, and on heat with ploidy as covariate; SLIM q-values
#             4. DML = q < config$dml_qvalue and |difference| > config$dml_difference
#             5. loci in both DML lists ("shared"; the legacy "interaction-DML")
#
# Run       Rscript 01_pipeline/05_differential_methylation.R     (after step 01)
#           WGBS_CORES sets the number of forked workers of the model fits
#           (default 1; leave it at 1 on Windows). Time and memory are in the
#           Notes.
#
# Inputs    03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files listed in 03_input_files/gannet_manifest.csv (cache)
# Outputs   04_output/05_differential_methylation/tables/
#             library_filtering_<mode>.csv       per-library CpG counts, cutoffs, scaling
#             united_cpg_<mode>.csv.gz           normalized counts at every analysed CpG, all libraries
#                                                (median-scaled and rounded as by methylKit; not raw
#                                                read counts)
#             diffmeth_<contrast>_<mode>.csv.gz  p, q, difference at every analysed CpG
#             dml_<contrast>_<mode>.csv          DML at the configured thresholds
#             dml_shared_<mode>.csv              loci in both DML lists, with both differences
#             dml_counts_by_threshold_<mode>.csv DML counts at difference thresholds of 5 to 50 points
#             dml_counts_by_mode.csv             headline counts, modes side by side
#           04_output/05_differential_methylation/bed/dml_<set>_<mode>.bed
#           04_output/05_differential_methylation/cache/methylkit_<mode>.rds   (not tracked)
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
#           - Fitting both modes takes about 18 minutes at WGBS_CORES=4, 37
#             minutes at WGBS_CORES=2 and about an hour on one core. The import
#             needs about 4 GB. The tests fork their workers: at WGBS_CORES=2 the
#             step needs about 6.5 GB (the R process and its two workers, shared
#             pages counted once), at WGBS_CORES=4 about 10 GB.
#           - A fit is kept in the cache folder and reused only while what it
#             was made from is unchanged: the mode and the coverage settings,
#             the sample sheet, the coverage files, the methylKit version and
#             the fitting code (fit_fingerprint() below). Otherwise, and when a
#             cache file cannot be read, the mode is fitted again; nothing has
#             to be deleted by hand. The thresholds of step 4 are applied after
#             the fit, so changing them does not refit. With both fits in the
#             cache the step takes about half a minute.

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
fit_mode <- function(mode) {
  treat_ploidy <- as.integer(lib$ploidy == "triploid")     # 0 diploid, 1 triploid
  treat_heat   <- as.integer(lib$treatment == "heat")      # 0 control, 1 heat

  # Pass over the libraries: filter each one, keep its surviving CpGs as
  # (position key, counts), and narrow the set of CpGs present in all libraries.
  # The md5 of each file is checked as it is read, so that the checksums in the
  # fingerprint of the fit are those of the files the fit was made from.
  kept <- vector("list", nrow(lib)); records <- vector("list", nrow(lib))
  chr_levels <- character(0); common <- NULL
  for (i in seq_len(nrow(lib))) {
    r <- read_filtered_counts(fetch_gannet_file(lib$cov_file[i], config, verify_md5 = TRUE), mode, config, chr_levels)
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
  list(mode = mode, filtering = filtering, united = meth_ploidy,
       diffmeth = list(ploidy = dm_ploidy, heat = dm_heat))
}

# ---- What a fit depends on -----------------------------------------------------------
# A cached fit is reused only while all of this is unchanged: the mode and the
# coverage settings, the columns of the sample sheet that the fit uses, the
# coverage files (size and md5 as the manifest lists them; fit_mode() checks
# the files against it), the methylKit version, and the code that builds the
# table of counts and runs the tests. The code is compared as R parses it, so a
# change of comments or layout does not count and any other change does, a
# reworded log message included. The thresholds that define a DML are applied
# after the fit (write_mode) and are not part of it; the number of cores is not
# either, because the results do not depend on it.
fit_fingerprint <- function(mode) {
  manifest <- read_gannet_manifest()
  files <- manifest[match(lib$cov_file, manifest$file), ]
  fitting_code <- list(fit_mode = fit_mode, read_filtered_counts = read_filtered_counts, read_cov = read_cov,
                       first_line_lost = first_line_lost, upper_coverage_cutoff = upper_coverage_cutoff,
                       cpg_key = cpg_key, key_to_position = key_to_position,
                       scale_counts_like_methylkit = scale_counts_like_methylkit, make_methylraw = make_methylraw)
  list(mode = mode,
       coverage_settings = unlist(config[c("min_coverage", "import_min_coverage", "upper_coverage_percentile")]),
       sample_sheet = paste(lib$library_id, lib$legacy_sample_index, lib$ploidy, lib$treatment, lib$cov_file),
       coverage_files = paste(files$file, files$bytes, files$md5),
       methylKit_version = as.character(utils::packageVersion("methylKit")),
       fitting_code = lapply(fitting_code, function(f) deparse(utils::removeSource(f))))
}

# ---- The fit of one mode, from the cache where that is safe ------------------------
# A cache file that cannot be read (a run was stopped while it was written, or
# the file is not a fit of this step) is fitted again, as is a fit made from
# anything other than what fit_fingerprint() returns now. A new fit is written
# under a temporary name and renamed, so a file with the final name is complete.
cached_fit <- function(mode) {
  cache_file <- file.path(step_cache_dir(step), sprintf("methylkit_%s.rds", mode))
  fingerprint <- fit_fingerprint(mode)
  if (file.exists(cache_file)) {
    cached <- tryCatch(suppressWarnings(readRDS(cache_file)), error = function(e) NULL)
    readable <- is.list(cached) && all(c("mode", "filtering", "united", "diffmeth") %in% names(cached))
    changed <- if (readable) names(fingerprint)[!vapply(names(fingerprint), function(part) {
      identical(cached$fingerprint[[part]], fingerprint[[part]])
    }, logical(1))]
    if (readable && length(changed) == 0) {
      log_msg("[%s] using the cached fit %s: settings, sample sheet, coverage files, methylKit version and fitting code are unchanged",
              mode, basename(cache_file))
      return(cached)
    }
    log_msg("[%s] the cached fit %s %s; fitting again", mode, basename(cache_file),
            if (!readable) "cannot be read (cut short, or not a fit of this step)"
            else if (is.null(cached$fingerprint)) "does not say what it was made from (written before this step recorded it)"
            else paste("differs from this run in:", paste(gsub("_", " ", changed), collapse = ", ")))
  }
  fit <- c(fit_mode(mode), list(fingerprint = fingerprint))
  part_file <- paste0(cache_file, ".part")
  saveRDS(fit, part_file)
  if (!file.rename(part_file, cache_file)) stop("Could not move ", part_file, " to ", cache_file, ".", call. = FALSE)
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
           sep = "\t", col.names = FALSE, scipen = 50, eol = "\n")   # fwrite would write CRLF on Windows
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
  fit <- cached_fit(mode)
  out <- write_mode(fit)
  rm(fit); gc(verbose = FALSE)
  out
}))
write_table(headline, output_path(step, "tables", "dml_counts_by_mode.csv"))
print(t(headline))

finish_step(step, t0)

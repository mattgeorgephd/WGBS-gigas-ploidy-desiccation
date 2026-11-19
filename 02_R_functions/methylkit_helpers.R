# methylkit_helpers.R
# Building the methylKit input with little memory.
#
# methylKit's own route (methRead, filterByCoverage, normalizeCoverage, unite on
# 22 full libraries) holds every library's full table in memory at once, about
# 9 GB here. The functions below apply the same rules one library at a time and
# keep only the CpGs that survive in all libraries. The arithmetic is methylKit's;
# 06_diagnostics/check_methylkit_equivalence.R confirms the result is identical.

# Read one library and apply the coverage filters.
#   import   rows with fewer than config$import_min_coverage reads are dropped
#            (methRead(mincov = ...)); the upper cutoff is the percentile of
#            coverage among the remaining rows (filterByCoverage(hi.perc = ...)).
#            In mode "legacy" the first line of the file is lost as well, as it
#            was in the legacy import (methRead(header = TRUE); issue I-38). The
#            cutoff is computed from the whole file in both modes, so that the
#            two modes report the same cutoff for a library.
#   lower    coverage >= config$min_coverage            (lo.count)
#   upper    coverage <  cutoff, in mode "corrected" only. The legacy scripts
#            passed the upper bound as `high.perc`, which methylKit ignores (its
#            argument is `hi.perc`), so no upper filter ran (issue I-05).
# Returns the kept rows (pos_key, coverage, nC, nT) and a one-row record.
read_filtered_counts <- function(path, mode, config, chr_levels) {
  if (!mode %in% c("legacy", "corrected")) stop("Unknown mode: ", mode)
  d <- read_cov(path)
  d[, coverage := nC + nT]
  cutoff <- upper_coverage_cutoff(d$coverage, config$import_min_coverage, config$upper_coverage_percentile)
  if (first_line_lost(mode)) d <- d[-1L]
  d <- d[coverage >= config$import_min_coverage]
  n_imported <- nrow(d)
  d <- d[coverage >= config$min_coverage]
  n_lower <- nrow(d)
  if (mode == "corrected") d <- d[coverage < cutoff]
  record <- data.frame(
    n_cpg_imported = n_imported, n_cpg_min_coverage = n_lower,
    first_line_lost = first_line_lost(mode),
    upper_cutoff = cutoff, upper_filter_applied = mode == "corrected",
    n_removed_by_upper = n_lower - nrow(d), n_cpg_kept = nrow(d),
    median_coverage = stats::median(d$coverage), max_coverage_kept = max(d$coverage))
  chr_levels <- union(chr_levels, unique(d$chr))
  list(counts = data.table::data.table(pos_key = cpg_key(d$chr, d$start, chr_levels),
                                       coverage = d$coverage, nC = d$nC, nT = d$nT),
       record = record, chr_levels = chr_levels)
}

# methylKit::normalizeCoverage(method = "median") for one library: counts are
# multiplied by `scaling_factor` (largest median coverage among libraries divided
# by this library's median) and rounded, exactly as methylKit does it.
scale_counts_like_methylkit <- function(coverage, nC, nT, scaling_factor) {
  frac_c <- nC / coverage; frac_t <- nT / coverage
  new_cov <- round(scaling_factor * coverage)
  list(coverage = new_cov, numCs = round(new_cov * frac_c), numTs = round(new_cov * frac_t))
}

# One methylRaw object from position and count vectors.
make_methylraw <- function(chr, start, counts, sample_id) {
  df <- data.frame(chr = chr, start = start, end = start + 2L, strand = "*",
                   coverage = counts$coverage, numCs = counts$numCs, numTs = counts$numTs,
                   stringsAsFactors = FALSE)
  methods::new("methylRaw", df, sample.id = sample_id, assembly = "Roslin_v1",
               context = "CpG", resolution = "base")
}

# Loci passing q-value and difference thresholds, as a data.table.
dml_table <- function(diffmeth, difference, qvalue) {
  data.table::as.data.table(methylKit::getData(
    methylKit::getMethylDiff(diffmeth, difference = difference, qvalue = qvalue)))
}

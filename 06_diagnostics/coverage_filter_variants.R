# coverage_filter_variants.R
#
# Purpose   How far the CpG and DML counts depend on the way the upper coverage
#           filter of the Methods is read (issue I-05). The Methods say that CpGs
#           "in the 99.9th percentile of coverage" were removed. The legacy call
#           misspelled the argument, so no upper filter ran (legacy mode).
#           Corrected mode applies methylKit's filter as the call would have
#           with the argument spelled correctly: one cutoff per library, the
#           percentile taken over the CpGs with at least
#           config$import_min_coverage reads, and a CpG removed from a library
#           at or above its cutoff. The Methods fix none of these choices. This
#           script
#             a. computes the cutoff of each library with the percentile taken
#                over the CpGs with at least 1, 2 and 10 reads, and one cutoff
#                over all libraries pooled;
#             b. applies nine filter rules to the CpGs tested in legacy mode,
#                and counts the CpGs that each rule keeps and the legacy DML
#                (ploidy, heat, loci in both lists) at kept CpGs. The rule of
#                corrected mode must give the CpGs and the DML lists of step 05;
#                the script stops if it does not;
#             c. describes the CpGs that the rule of corrected mode removes: in
#                how many libraries they reach the cutoff, how many lie in the
#                repeat track and on the mitochondrial sequence, and what share
#                lies in runs of removed CpGs;
#             d. gives the evidence on whether the DML counts of b are the
#                lists that step 05 would report under each rule: how many CpGs
#                above the difference threshold fail the q-value rule in the
#                tables of step 05, whether a rule changes the scaling of a
#                library, and whether q-values estimated again on the kept CpGs
#                give the same lists.
#
# Run       Rscript 06_diagnostics/coverage_filter_variants.R      (about 7 minutes; 1.6 GB of memory; one core)
#           Needs the 22 coverage files in the cache (step 01).
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files listed in 03_input_files/gannet_manifest.csv (cache)
#           03_input_files/genome_features/cgigas_uk_roslin_v1_rm.te.bed.gz     the repeat track
#           04_output/05_differential_methylation/tables/, both modes:
#             united_cpg_<mode>.csv.gz, diffmeth_{ploidy,heat}_<mode>.csv.gz,
#             dml_{ploidy,heat,shared}_<mode>.csv, library_filtering_<mode>.csv
# Writes    06_diagnostics/coverage_filter_variants.csv       one row per filter rule
#           06_diagnostics/coverage_filter_cutoffs.csv        one row per library
#           06_diagnostics/coverage_filter_removed_cpgs.csv   one row per set of CpGs
#           06_diagnostics/coverage_filter_qvalue_rule.csv    one row per set of CpGs and contrast
#
# Notes     - What methylKit does (read in the source of methylKit 1.28.0):
#             filterByCoverage(hi.perc = p) computes for each library
#             quantile(coverage, p / 100), R's default quantile (type 7), over
#             the rows that the import kept, and keeps the rows with coverage
#             below that value. Which rows the import kept is set by
#             methRead(mincov = ...). normalizeCoverage() then scales the counts
#             of a library by the largest median coverage among the libraries
#             divided by its own median, and unite() keeps the CpGs present in
#             every library. A CpG is therefore lost when a single library is at
#             or above its cutoff there.
#           - "Base" is the lowest read count of the CpGs over which the
#             percentile is taken. Base 2 is config$import_min_coverage, the
#             mincov of the legacy DML script, and the base of corrected mode.
#             Base 1 is the mincov of the legacy coverage script. Base 10 is the
#             default mincov of methRead.
#           - The cutoffs come from a histogram of coverage over every line of a
#             file, with the arithmetic of stats::quantile(). The script stops
#             unless each is identical to the value that the pipeline's
#             upper_coverage_cutoff() gives for the same file, and unless the
#             cutoffs of base 2 are those in library_filtering_corrected.csv.
#             The pooled cutoff is taken from the summed histograms of all
#             libraries (base 2) and applied to every library. In the cutoff
#             table, n_cpg_base_<b> is the number of CpGs over which the
#             percentile is taken, and n_removed_<...> is the number of CpGs
#             with at least config$min_coverage reads at or above the cutoff:
#             the rows that the filter takes from the library.
#           - Rules. "per library": each library has its own cutoff. "above the
#             cutoff only": a CpG with exactly the cutoff number of reads stays.
#             The last three rules remove a CpG only when at least 2 libraries,
#             more than half of the libraries, or all libraries are at or above
#             their cutoffs. methylKit has no such rule, because its filter acts
#             within one library. Here these three rules drop CpGs from the
#             united table of legacy mode, after the libraries were scaled.
#           - Column cpg_removed counts the CpGs of legacy mode that a rule
#             drops. ploidy_dml, heat_dml and dml_in_both_lists count the DML of
#             legacy mode at kept CpGs. They are the lists that step 05 would
#             give under the rule if two things hold, and the table reports both.
#             (1) libraries_with_changed_median: the number of libraries whose
#             median coverage, over the CpGs that the library keeps under the
#             rule, differs from its median in legacy mode. Where it is 0 the
#             scaling factors, the normalized counts, and with them the p-value
#             and the difference of every kept CpG are those of legacy mode.
#             For corrected mode this is also checked in the tables of step 05:
#             the script stops unless the p-values and the differences of the
#             CpGs of corrected mode are identical in both modes.
#             (2) lists_same_with_reestimated_q: q-values depend on the set of
#             CpGs tested. They are estimated again on the kept CpGs with the
#             function that calculateDiffMeth applies (methylKit's p.adjusted,
#             method "SLIM"), from the p-values in the legacy table. TRUE where
#             the DML lists that follow are the legacy lists at kept CpGs, locus
#             for locus. The script stops unless this estimate reproduces the
#             q-values of step 05 in both modes to a relative difference below
#             q_tolerance.
#           - The q-value table gives, for each set of CpGs and contrast, the
#             CpGs above the difference threshold (n_above_difference), those of
#             them that fail the q-value rule (n_failing_q), and how close the
#             q-values come to the cutoff from either side (largest_q_of_dml,
#             smallest_q_failing; three significant digits).
#           - A run is a set of removed CpGs on one sequence in which neighbours
#             are no more than run_gap_bp apart, the rule by which step 07
#             groups DML into regions. Kept CpGs may lie between them. The share
#             of removed CpGs in runs of run_min_cpgs (10) or more is set
#             against the same share for random sets of as many tested CpGs:
#             n_draws (1,000) sets, each drawn without replacement, after
#             set.seed(seed). Nothing else in the script is random.
#           - First line of a file. Legacy mode lost it (issue I-38). As in the
#             pipeline, the cutoff of a library is computed from every line in
#             both modes. The tested CpGs are looked up in each file as legacy
#             mode read it. The script also recounts the CpGs with at least
#             config$min_coverage reads in all libraries, once with line 1 lost
#             and once with every line read, and stops unless both recounts
#             give exactly the CpGs of united_cpg_legacy.csv.gz.
#           - The repeat track is the transposable element track of step 07
#             (RepeatMasker annotation), with the overlap rule of step 07.
#             Unplaced scaffolds are all sequences other than the ten linkage
#             groups and the mitochondrial sequence (run_config.R).
#           - Raw coverage of the tested CpGs is held in memory only; no table
#             with one row per CpG is written.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges", "methylKit"))
setDTthreads(config$cores)             # results do not depend on the number of threads

# ---- Settings --------------------------------------------------------------------------
percentile   <- config$upper_coverage_percentile
# Lowest read count of the CpGs over which the percentile is taken ("base"):
# every CpG with reads (mincov = 1 of the legacy coverage script), the import
# threshold of the pipeline (mincov = 2 of the legacy DML script; the base of
# corrected mode), and the default mincov of methRead in the installed methylKit.
base_import  <- as.integer(config$import_min_coverage)
bases        <- c(1L, base_import, as.integer(formals(methylKit::methRead)$mincov))
run_gap_bp   <- 1000L                  # removed CpGs no further apart than this form one run (the gap of step 07's DML regions)
run_min_cpgs <- 10L                    # the run statistic counts removed CpGs in runs of at least this many
n_draws      <- 1000L                  # random sets of tested CpGs that give the reference for the run statistic
seed         <- 20261119L              # for those draws; nothing else in this script is random
q_tolerance  <- 1e-9                   # largest relative difference accepted between a q-value estimated here and the one step 05 stored
contrast_names <- c(ploidy = "ploidy", heat = "heat")
stopifnot(!anyDuplicated(bases))

out_file <- function(name) here("06_diagnostics", paste0("coverage_filter_", name, ".csv"))
tab05 <- function(file) output_path("05_differential_methylation", "tables", file)
# Stop with a message unless a statement about the inputs holds.
must <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)

# R's default quantile (type 7) from a histogram: counts[v] is the number of
# CpGs with v reads. The arithmetic is that of stats::quantile(), so the value
# is the one quantile() returns for the CpGs themselves.
hist_quantile <- function(counts, prob) {
  n <- sum(as.numeric(counts))
  must(n >= 1, "A histogram without CpGs has no quantile.")
  index <- 1 + (n - 1) * prob
  lo <- floor(index); hi <- ceiling(index)
  cumulative <- cumsum(as.numeric(counts))
  value_lo <- which(cumulative >= lo)[1L]; value_hi <- which(cumulative >= hi)[1L]
  h <- index - lo
  if (h > 0 && value_hi != value_lo) (1 - h) * value_lo + h * value_hi else as.numeric(value_lo)
}
# Upper cutoff from a histogram: the percentile over the CpGs with at least base_reads reads.
cutoff_of <- function(counts, base_reads) hist_quantile(replace(counts, seq_along(counts) < base_reads, 0L), percentile / 100)

# ---- 1. Tables of step 05 ------------------------------------------------------------------
lib <- read_libraries(); n_lib <- nrow(lib)
tested <- fread(tab05("united_cpg_legacy.csv.gz"), select = c("chr", "start"))     # the CpGs tested in legacy mode
n_tested <- nrow(tested)
seq_levels <- unique(tested$chr)
tested_key <- cpg_key(tested$chr, tested$start, seq_levels)
must(!anyDuplicated(tested_key), "united_cpg_legacy.csv.gz lists a CpG twice.")

# Per mode: the row of each of its CpGs in the legacy table, p-value, q-value and
# difference of every CpG for both contrasts, and the three DML lists as row
# numbers of the legacy table.
read_mode <- function(mode_name) {
  position <- if (mode_name == "legacy") tested else fread(tab05(sprintf("united_cpg_%s.csv.gz", mode_name)), select = c("chr", "start"))
  legacy_row <- match(cpg_key(position$chr, position$start, seq_levels), tested_key)
  must(!anyNA(legacy_row) && !is.unsorted(legacy_row, strictly = TRUE),
       "The CpGs of ", mode_name, " mode are not CpGs of legacy mode in the same order; an upper filter can only remove CpGs.")
  stats <- lapply(contrast_names, function(contrast_name) {
    d <- fread(tab05(sprintf("diffmeth_%s_%s.csv.gz", contrast_name, mode_name)), select = c("chr", "start", "pvalue", "qvalue", "meth.diff"))
    must(identical(d$chr, position$chr) && identical(d$start, position$start) && !anyNA(d$pvalue) && !anyNA(d$qvalue) && !anyNA(d$meth.diff),
         "diffmeth_", contrast_name, "_", mode_name, ".csv.gz does not hold one complete row per CpG of the united table.")
    list(pvalue = d$pvalue, qvalue = d$qvalue, difference = d$meth.diff)
  })
  dml <- lapply(c(contrast_names, shared = "shared"), function(set_name) {
    d <- fread(tab05(sprintf("dml_%s_%s.csv", set_name, mode_name)), select = c("chr", "start"))
    rows <- match(cpg_key(d$chr, d$start, seq_levels), tested_key)
    must(!anyNA(rows) && !anyDuplicated(rows), "dml_", set_name, "_", mode_name, ".csv holds a locus that is not a tested CpG, or a locus twice.")
    sort(rows)
  })
  for (contrast_name in contrast_names) {   # the lists are the two rules of step 05 applied to its table of all CpGs
    s <- stats[[contrast_name]]
    must(identical(dml[[contrast_name]], legacy_row[s$qvalue < config$dml_qvalue & abs(s$difference) > config$dml_difference]),
         "dml_", contrast_name, "_", mode_name, ".csv is not the q-value rule and the difference rule applied to the diffmeth table.")
  }
  must(identical(dml$shared, intersect(dml$ploidy, dml$heat)), "dml_shared_", mode_name, ".csv is not the set of loci in both DML lists.")
  list(legacy_row = legacy_row, stats = stats, dml = dml)
}
step05 <- lapply(c(legacy = "legacy", corrected = "corrected"), read_mode)
record <- lapply(c(legacy = "legacy", corrected = "corrected"), function(mode_name) fread(tab05(sprintf("library_filtering_%s.csv", mode_name))))
must(identical(record$legacy$library_id, lib$library_id) && identical(record$corrected$library_id, lib$library_id),
     "The library_filtering tables of step 05 are not in the order of the sample sheet.")
# At the CpGs of both modes the p-values and the differences must be the same:
# the filter of corrected mode did not change the normalized counts.
for (contrast_name in contrast_names) {
  s_legacy <- step05$legacy$stats[[contrast_name]]; s_corrected <- step05$corrected$stats[[contrast_name]]; rows <- step05$corrected$legacy_row
  must(identical(s_legacy$pvalue[rows], s_corrected$pvalue) && identical(s_legacy$difference[rows], s_corrected$difference),
       "The ", contrast_name, " p-values or differences of a CpG differ between the modes of step 05; the DML of legacy mode at kept CpGs ",
       "are then not the DML of another filter rule.")
}

# ---- 2. One pass over the coverage files -----------------------------------------------------
# Per library: the histogram of coverage over every line, the coverage of the
# tested CpGs, and the CpGs with at least config$min_coverage reads in every
# library read so far (once with every line read, once with line 1 lost).
log_msg("reading %d coverage files", n_lib)
manifest <- read_gannet_manifest()
histograms <- vector("list", n_lib)                 # histograms[[i]][v]: CpGs of library i with v reads
first_line_reads <- integer(n_lib)                  # coverage of the CpG on line 1 of each file
reads_tested <- matrix(NA_integer_, n_tested, n_lib, dimnames = list(NULL, lib$library_id))
in_all_complete <- in_all_legacy <- numeric(0)      # position keys
for (lib_i in seq_len(n_lib)) {
  cov_table <- read_cov(fetch_gannet_file(lib$cov_file[lib_i], config, manifest))
  reads <- cov_table$nC + cov_table$nT
  must(length(reads) > 0L && !anyNA(reads) && min(reads) >= 1L,
       "The coverage file of ", lib$library_id[lib_i], " holds a line without reads or without counts.")
  seq_levels <- union(seq_levels, unique(cov_table$chr))            # indices of earlier names never change
  key <- cpg_key(cov_table$chr, cov_table$start, seq_levels)
  rm(cov_table)
  histograms[[lib_i]] <- tabulate(reads)
  first_line_reads[lib_i] <- reads[1L]
  for (base_reads in bases)
    must(identical(cutoff_of(histograms[[lib_i]], base_reads), as.numeric(upper_coverage_cutoff(reads, base_reads, percentile))),
         "The cutoff from the histogram of ", lib$library_id[lib_i], " (base ", base_reads, ") is not the value of upper_coverage_cutoff().")
  enough <- reads >= config$min_coverage
  in_all_complete <- if (lib_i == 1L) key[enough] else in_all_complete[in_all_complete %in% key[enough]]
  if (first_line_lost("legacy")) enough[1L] <- FALSE                # legacy mode never saw line 1
  in_all_legacy <- if (lib_i == 1L) key[enough] else in_all_legacy[in_all_legacy %in% key[enough]]
  file_row <- match(tested_key, key)
  must(!anyNA(file_row) && all(enough[file_row]),
       "A CpG of united_cpg_legacy.csv.gz has fewer than ", config$min_coverage, " reads in the coverage file of ", lib$library_id[lib_i], ".")
  reads_tested[, lib_i] <- reads[file_row]
  log_msg("%-10s %9d CpGs with reads; in all libraries so far with at least %d reads: %d", lib$library_id[lib_i], length(reads),
          config$min_coverage, length(in_all_legacy))
  rm(reads, key, enough, file_row); invisible(gc(verbose = FALSE))
}
must(setequal(in_all_legacy, tested_key),
     "The CpGs with at least ", config$min_coverage, " reads in all libraries, recounted from the coverage files, are not the CpGs of united_cpg_legacy.csv.gz.")
must(setequal(in_all_complete, in_all_legacy),
     "A CpG on the first line of a coverage file has at least ", config$min_coverage, " reads in all libraries; the rules of this script would then ",
     "start from different CpGs in a complete read.")
n_recounted <- length(in_all_legacy)
rm(in_all_complete, in_all_legacy)

# ---- 3. Cutoffs (a) ----------------------------------------------------------------------------
max_reads <- max(lengths(histograms))
counts <- vapply(histograms, function(h) c(h, integer(max_reads - length(h))), integer(max_reads))   # read count by library
reads_axis <- seq_len(max_reads)
cutoff_library <- vapply(bases, function(base_reads) apply(counts, 2L, cutoff_of, base_reads = base_reads), numeric(n_lib))
dimnames(cutoff_library) <- list(lib$library_id, as.character(bases))
cutoff_pooled <- cutoff_of(rowSums(counts), base_import)
log_msg("pooled cutoff over all libraries (base %d): %g reads", base_import, cutoff_pooled)

# CpGs that a library keeps: at least config$min_coverage reads and fewer reads
# than its upper bound (no more than the bound where above_only). Returns their
# number and their median coverage, the quantity that sets the scaling factor.
library_kept <- function(count_matrix, upper, above_only = FALSE) {
  vapply(seq_len(n_lib), function(lib_i) {
    inside <- reads_axis >= config$min_coverage & (if (above_only) reads_axis <= upper[lib_i] else reads_axis < upper[lib_i])
    kept <- replace(count_matrix[, lib_i], !inside, 0L)
    c(n = sum(kept), median = hist_quantile(kept, 0.5))
  }, numeric(2))
}
no_bound <- rep(Inf, n_lib)
counts_legacy <- counts                                              # the files as legacy mode read them
if (first_line_lost("legacy")) {
  first_cell <- cbind(first_line_reads, seq_len(n_lib))
  counts_legacy[first_cell] <- counts_legacy[first_cell] - 1L
}
kept_legacy    <- library_kept(counts_legacy, no_bound)                                    # legacy mode: lower threshold only, line 1 lost
kept_complete  <- library_kept(counts, no_bound)                                           # lower threshold only, every line read
kept_corrected <- library_kept(counts, cutoff_library[, as.character(base_import)])        # corrected mode
must(all(kept_legacy["n", ] == record$legacy$n_cpg_min_coverage) && all(kept_legacy["median", ] == record$legacy$median_coverage),
     "The histograms do not give the CpG numbers and median coverages of library_filtering_legacy.csv.")
must(isTRUE(all.equal(unname(cutoff_library[, as.character(base_import)]), record$corrected$upper_cutoff, tolerance = 1e-12)) &&
       all(kept_complete["n", ] == record$corrected$n_cpg_min_coverage) && all(kept_corrected["n", ] == record$corrected$n_cpg_kept) &&
       all(kept_complete["n", ] - kept_corrected["n", ] == record$corrected$n_removed_by_upper) &&
       all(kept_corrected["median", ] == record$corrected$median_coverage),
     "The histograms do not give the cutoffs, CpG numbers and median coverages of library_filtering_corrected.csv.")

n_removed_from_library <- function(cutoffs) vapply(seq_len(n_lib), function(lib_i)
  sum(counts[reads_axis >= config$min_coverage & reads_axis >= cutoffs[lib_i], lib_i]), integer(1))
cutoff_table <- data.table(library_id = lib$library_id, cohort = lib$cohort, ploidy = as.character(lib$ploidy), treatment = as.character(lib$treatment))
for (base_reads in bases) set(cutoff_table, j = sprintf("n_cpg_base_%d", base_reads), value = as.integer(colSums(counts[reads_axis >= base_reads, , drop = FALSE])))
for (base_reads in bases) set(cutoff_table, j = sprintf("cutoff_base_%d", base_reads), value = cutoff_library[, as.character(base_reads)])
set(cutoff_table, j = "cutoff_pooled", value = cutoff_pooled)
for (base_reads in bases) set(cutoff_table, j = sprintf("n_removed_base_%d", base_reads), value = n_removed_from_library(cutoff_library[, as.character(base_reads)]))
set(cutoff_table, j = "n_removed_pooled", value = n_removed_from_library(rep(cutoff_pooled, n_lib)))

# ---- 4. Filter rules (b) -----------------------------------------------------------------------
majority <- n_lib %/% 2L + 1L                                        # more than half of the libraries
# One row per rule: where the cutoff comes from, the base of its percentile,
# whether a CpG at the cutoff stays, how many libraries must reach the cutoff
# for a CpG to be removed, and the mode of step 05 that applies the rule.
new_rule <- function(label, cutoff, percentile_base = base_import, above_only = FALSE, min_libraries = 1L, step05_mode = NA_character_)
  data.table(rule = label, cutoff = cutoff, percentile_base = percentile_base, above_only = above_only, min_libraries = min_libraries,
             step05_mode = step05_mode)
rules <- rbind(
  new_rule("no upper filter", "none", percentile_base = NA_integer_, above_only = NA, min_libraries = NA_integer_, step05_mode = "legacy"),
  new_rule(sprintf("per library, base %d", base_import), "per library", step05_mode = "corrected"),
  new_rule(sprintf("per library, base %d, above the cutoff only", base_import), "per library", above_only = TRUE),
  rbindlist(lapply(setdiff(bases, base_import), function(base_reads)
    new_rule(sprintf("per library, base %d", base_reads), "per library", percentile_base = base_reads))),
  new_rule(sprintf("pooled, base %d", base_import), "pooled"),
  new_rule(sprintf("per library, base %d, at least 2 libraries", base_import), "per library", min_libraries = 2L),
  new_rule(sprintf("per library, base %d, more than half of the libraries", base_import), "per library", min_libraries = majority),
  new_rule(sprintf("per library, base %d, all %d libraries", base_import, n_lib), "per library", min_libraries = n_lib))
n_rules <- nrow(rules)
stopifnot(!anyDuplicated(rules$rule))
rule_cutoffs <- function(rule_i) switch(rules$cutoff[rule_i], none = no_bound, pooled = rep(cutoff_pooled, n_lib),
                                        `per library` = unname(cutoff_library[, as.character(rules$percentile_base[rule_i])]))
# Number of libraries in which each tested CpG is at or above the cutoff (above it where above_only).
libraries_at_cutoff <- function(cutoffs, above_only = FALSE) {
  n <- integer(n_tested)
  for (lib_i in seq_len(n_lib)) n <- n + if (above_only) reads_tested[, lib_i] > cutoffs[lib_i] else reads_tested[, lib_i] >= cutoffs[lib_i]
  n
}
keep <- lapply(seq_len(n_rules), function(rule_i) {
  if (rules$cutoff[rule_i] == "none") return(rep(TRUE, n_tested))
  libraries_at_cutoff(rule_cutoffs(rule_i), rules$above_only[rule_i]) < rules$min_libraries[rule_i]
})
# Median coverage of each library over the CpGs it keeps under the rule, every
# line read. Only a rule that removes a CpG as soon as one library is at its
# cutoff filters within a library (methylKit's filterByCoverage); the other
# rules leave the libraries as legacy mode filtered them.
changed_median <- vapply(seq_len(n_rules), function(rule_i) {
  within_library <- rules$cutoff[rule_i] != "none" && rules$min_libraries[rule_i] == 1L
  kept <- if (within_library) library_kept(counts, rule_cutoffs(rule_i), rules$above_only[rule_i]) else kept_complete
  sum(kept["median", ] != kept_legacy["median", ])
}, integer(1))

dml_legacy <- step05$legacy$dml
variants <- data.table(
  rule = rules$rule, cutoff = rules$cutoff, percentile_base_reads = rules$percentile_base,
  removed_when = fifelse(rules$cutoff == "none", NA_character_, fifelse(rules$above_only, "above the cutoff", "at or above the cutoff")),
  in_at_least_n_libraries = rules$min_libraries, step05_mode = rules$step05_mode,
  cpg_kept = vapply(keep, sum, integer(1)), cpg_removed = n_tested - vapply(keep, sum, integer(1)),
  pct_removed = round(100 * (1 - vapply(keep, mean, numeric(1))), 2),
  ploidy_dml = vapply(keep, function(k) sum(k[dml_legacy$ploidy]), integer(1)),
  heat_dml = vapply(keep, function(k) sum(k[dml_legacy$heat]), integer(1)),
  dml_in_both_lists = vapply(keep, function(k) sum(k[dml_legacy$shared]), integer(1)),
  libraries_with_changed_median = changed_median)

# The rule of each mode of step 05 must give that mode's CpGs and DML lists.
for (rule_i in which(!is.na(rules$step05_mode))) {
  mode_name <- rules$step05_mode[rule_i]; kept_rows <- which(keep[[rule_i]])
  must(identical(kept_rows, step05[[mode_name]]$legacy_row),
       "Rule '", rules$rule[rule_i], "' does not give the CpGs of united_cpg_", mode_name, ".csv.gz.")
  for (set_name in names(dml_legacy))
    must(identical(dml_legacy[[set_name]][keep[[rule_i]][dml_legacy[[set_name]]]], step05[[mode_name]]$dml[[set_name]]),
         "The legacy ", set_name, " DML at the CpGs kept by rule '", rules$rule[rule_i], "' are not the list dml_", set_name, "_", mode_name, ".csv.")
}
must(n_recounted == variants$cpg_kept[rules$cutoff == "none"], "The recount from the coverage files and the rule without a filter disagree.")

# ---- 5. The CpGs that corrected mode removes (c) ------------------------------------------------
corrected_rule <- which(rules$step05_mode == "corrected")
n_at_cutoff <- libraries_at_cutoff(rule_cutoffs(corrected_rule))
removed <- n_at_cutoff >= 1L
must(identical(removed, !keep[[corrected_rule]]), "The removed CpGs are not those of the rule of corrected mode.")
n_removed <- sum(removed)
in_repeat <- overlaps_track(cpg_ranges(tested$chr, tested$start), read_feature_track("TE"))
is_mito <- tested$chr == config$mito_sequence
is_unplaced <- !is_mito & !(tested$chr %in% config$linkage_groups)

# Runs: the tested CpGs in genome order; a selected CpG starts a new run when it
# is the first of its sequence or lies more than run_gap_bp after the previous
# selected CpG. Returns the size of the run of each selected CpG.
genome_order <- order(match(tested$chr, unique(tested$chr)), tested$start)
sequence_sorted <- tested$chr[genome_order]; start_sorted <- tested$start[genome_order]
run_size <- function(selected) {                    # selected: logical over the tested CpGs in genome order
  if (!any(selected)) return(integer(0))
  s <- sequence_sorted[selected]; p <- start_sorted[selected]
  run <- cumsum(c(TRUE, s[-1L] != s[-length(s)] | diff(p) > run_gap_bp))
  tabulate(run)[run]
}
in_long_run <- logical(n_tested)
in_long_run[genome_order[removed[genome_order]]] <- run_size(removed[genome_order]) >= run_min_cpgs
set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
random_pct <- replicate(n_draws, {
  selected <- logical(n_tested); selected[sample.int(n_tested, n_removed)] <- TRUE
  100 * mean(run_size(selected) >= run_min_cpgs)
})

# One row per set of CpGs. The run columns are filled for sets of removed CpGs only.
describe_set <- function(label, in_set, within_removed) data.table(
  cpg_set = label, n_cpg = sum(in_set), pct_of_removed = if (within_removed) round(100 * sum(in_set) / n_removed, 2) else NA_real_,
  n_in_repeat_track = sum(in_set & in_repeat), pct_in_repeat_track = round(100 * sum(in_set & in_repeat) / sum(in_set), 2),
  n_mitochondrial = sum(in_set & is_mito), n_unplaced_scaffold = sum(in_set & is_unplaced),
  n_in_long_runs = if (within_removed) sum(in_set & in_long_run) else NA_integer_,
  pct_in_long_runs = if (within_removed) round(100 * sum(in_set & in_long_run) / sum(in_set), 2) else NA_real_)
removed_table <- rbind(
  describe_set("tested in legacy mode", rep(TRUE, n_tested), FALSE),
  describe_set("kept by the rule of corrected mode", !removed, FALSE),
  describe_set("removed by the rule of corrected mode", removed, TRUE),
  describe_set("removed: at or above the cutoff in exactly 1 library", n_at_cutoff == 1L, TRUE),
  describe_set(sprintf("removed: at or above the cutoff in 2 to %d libraries", n_lib - 1L), n_at_cutoff > 1L & n_at_cutoff < n_lib, TRUE),
  describe_set(sprintf("removed: at or above the cutoff in all %d libraries", n_lib), n_at_cutoff == n_lib, TRUE))
# The reference for the run statistic belongs to the row of all removed CpGs.
all_removed <- removed_table$cpg_set == "removed by the rule of corrected mode"
set(removed_table, j = "random_sets_pct_in_long_runs_median", value = fifelse(all_removed, round(median(random_pct), 2), NA_real_))
set(removed_table, j = "random_sets_pct_in_long_runs_max", value = fifelse(all_removed, round(max(random_pct), 2), NA_real_))
set(removed_table, j = "random_sets_drawn", value = fifelse(all_removed, n_draws, NA_integer_))
setnames(removed_table, names(removed_table), sub("long_runs", sprintf("runs_of_%d_or_more", run_min_cpgs), names(removed_table), fixed = TRUE))

# ---- 6. The q-value rule (d) ----------------------------------------------------------------------
# One row per CpG set and contrast: the CpGs above the difference threshold, how
# many of them fail the q-value rule, and whether the resulting list is the
# legacy list at the CpGs of the set.
q_row <- function(cpg_set, q_values, contrast_name, legacy_rows, difference, q) {
  above <- abs(difference) > config$dml_difference
  passes <- q < config$dml_qvalue
  legacy_list <- dml_legacy[[contrast_name]]
  data.table(cpg_set = cpg_set, q_values = q_values, contrast = contrast_name, n_cpg = length(q),
             n_above_difference = sum(above), n_failing_q = sum(above & !passes), n_dml = sum(above & passes),
             largest_q_of_dml = if (any(above & passes)) signif(max(q[above & passes]), 3) else NA_real_,
             smallest_q_failing = if (any(above & !passes)) signif(min(q[above & !passes]), 3) else NA_real_,
             list_is_legacy_dml_at_these_cpgs = identical(legacy_rows[above & passes], legacy_list[legacy_list %in% legacy_rows]))
}
q_stored <- rbindlist(lapply(names(step05), function(mode_name) rbindlist(lapply(contrast_names, function(contrast_name) {
  s <- step05[[mode_name]]$stats[[contrast_name]]
  q_row(paste(mode_name, "mode"), "as stored by step 05", contrast_name, step05[[mode_name]]$legacy_row, s$difference, s$qvalue)
}))))

# q-values of each rule: methylKit's SLIM estimate on the p-values of the kept CpGs.
slim_q <- getFromNamespace("p.adjusted", "methylKit")
log_msg("estimating q-values again on the CpGs kept by each of %d rules, for %d contrasts", n_rules, length(contrast_names))
q_estimated <- rbindlist(lapply(seq_len(n_rules), function(rule_i) rbindlist(lapply(contrast_names, function(contrast_name) {
  s <- step05$legacy$stats[[contrast_name]]; kept_rows <- which(keep[[rule_i]])
  q_here <- slim_q(s$pvalue[kept_rows], method = "SLIM")
  mode_name <- rules$step05_mode[rule_i]
  if (!is.na(mode_name)) {                         # the estimate must reproduce what step 05 stored for this set of CpGs
    q_step05 <- step05[[mode_name]]$stats[[contrast_name]]$qvalue
    must(all(abs(q_here - q_step05) <= q_tolerance * q_step05),
         "The q-values estimated here for rule '", rules$rule[rule_i], "' (", contrast_name, ") are not those of step 05, ", mode_name, " mode.")
  }
  log_msg("  %-55s %-6s done", rules$rule[rule_i], contrast_name)
  q_row(rules$rule[rule_i], "estimated again on the kept CpGs", contrast_name, kept_rows, s$difference[kept_rows], q_here)
}))))
q_table <- rbind(q_stored, q_estimated)
set(variants, j = "lists_same_with_reestimated_q", value = vapply(rules$rule, function(label)
  all(q_estimated$list_is_legacy_dml_at_these_cpgs[q_estimated$cpg_set == label]), logical(1), USE.NAMES = FALSE))

# ---- 7. Write and print -------------------------------------------------------------------------------
write_table(variants, out_file("variants"))
write_table(cutoff_table, out_file("cutoffs"))
write_table(removed_table, out_file("removed_cpgs"))
write_table(q_table, out_file("qvalue_rule"))

options(width = 250)
cat("\nFilter rules: CpGs kept and legacy DML at kept CpGs\n")
print(variants[, .(rule, cpg_kept, cpg_removed, pct_removed, ploidy_dml, heat_dml, dml_in_both_lists, libraries_with_changed_median, lists_same_with_reestimated_q)])
cat("\nCutoffs per library\n")
print(cutoff_table[, c("library_id", grep("^cutoff_|^n_removed_", names(cutoff_table), value = TRUE)), with = FALSE])
cat("\nCpGs removed by the rule of corrected mode\n")
print(removed_table)
cat("\nThe q-value rule\n")
print(q_table)

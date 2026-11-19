# check_methylkit_equivalence.R
#
# Purpose   Step 05 does not build its table of tested CpGs with methylKit's own
#           import. methRead, filterByCoverage, normalizeCoverage and unite hold
#           all 22 libraries in memory at once (about 9 GB), so step 05 applies
#           the same rules one library at a time with the functions of
#           02_R_functions/methylkit_helpers.R and then calls methylKit's test.
#           This script shows that the helper route gives what methylKit's own
#           route gives, in both modes of run_config.R (issues I-05 and I-38).
#             1. Extract. A fixed set of regions, and the first line of every
#                file, are cut from the 22 coverage files into temporary files
#                in the original format.
#             2. Full libraries. For a few libraries the part of both routes
#                that acts on one library (import and filter) is run on the
#                whole coverage file and compared row by row, and with the
#                record that step 05 wrote for the library.
#             3. Both routes on the extract, per mode. methylKit's route:
#                methRead, filterByCoverage, normalizeCoverage, unite, then
#                calculateDiffMeth for ploidy with heat as covariate and for
#                heat with ploidy as covariate. Helper route:
#                read_filtered_counts, scale_counts_like_methylkit,
#                make_methylraw and unite, chained as in step 05, then the same
#                two tests. Compared: rows imported, upper cutoff, rows kept,
#                median coverage and scaling factor of every library; the CpGs
#                and every count of the united table; difference, p-value and
#                q-value of every CpG for both contrasts; the DML lists.
#             4. The legacy call in methylKit's own objects: what methRead loses
#                with its default header = TRUE (I-38); what the misspelled
#                upper filter of the legacy call does (I-05); and what changes
#                when the covariate of the test is coded as the legacy script
#                coded it (a character vector) and not as step 05 codes it (a
#                factor).
#             5. The tables of step 05. Every count that step 05 wrote must be
#                the count of the coverage file scaled by
#                scale_counts_like_methylkit, at every CpG of the table and in
#                every library. Inside the regions of the extract the CpGs of
#                the table must be the CpGs that pass the rule of the mode in
#                every library. This ties the routes compared here to the
#                tables that the pipeline uses.
#             6. Rounding. methylKit::normalizeCoverage against
#                scale_counts_like_methylkit on synthetic counts built to hold
#                many values exactly halfway between two whole numbers.
#             7. The order of operations of step 04, which scales counts for
#                its correlations and PCA as round(new_cov * nc / cv) where
#                methylKit computes round(new_cov * (nc / cv)): how many values
#                of the extract differ, and by how much.
#           One row is written per comparison. If a comparison that must agree
#           does not, the tables are still written and the script then stops
#           with an error that lists the comparisons that failed.
#
# Run       Rscript 06_diagnostics/check_methylkit_equivalence.R      (about 9 minutes; 2 GB of memory; one core)
#           Needs the 22 coverage files in the cache (step 01) and the tables of step 05.
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files listed in 03_input_files/gannet_manifest.csv (cache)
#           04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           04_output/05_differential_methylation/tables/library_filtering_<mode>.csv
# Writes    06_diagnostics/check_methylkit_equivalence.csv             one row per comparison
#           06_diagnostics/check_methylkit_equivalence_libraries.csv   one row per library, mode and scope
#           06_diagnostics/check_methylkit_equivalence_extract.csv     one row per region of the extract
#
# Notes     - What methylKit does (read in the source of methylKit 1.28.0).
#             methRead(pipeline = "bismarkCoverage") reads a file with
#             data.table::fread, by default with header = TRUE, and keeps the
#             rows with at least mincov reads. filterByCoverage(lo.count,
#             hi.perc = p) computes quantile(coverage, p / 100) over the rows
#             that the import kept and keeps the rows with at least lo.count
#             reads and fewer reads than that value. Arguments it does not know
#             (high.perc) are ignored without a message. normalizeCoverage()
#             multiplies the counts of a library by the largest median coverage
#             among the libraries divided by its own median: coverage becomes
#             round(factor * coverage), and each count round(new coverage *
#             (count / coverage)). unite() keeps the CpGs present in every
#             library and sorts them by sequence name and position.
#           - What methylKit does not report. filterByCoverage returns the rows
#             it kept and not its cutoff; normalizeCoverage returns the scaled
#             counts and not its factors. The cutoff and the factor that stand
#             beside those of the helpers (columns
#             upper_cutoff_from_methylkit_rows and
#             scaling_factor_from_methylkit_medians of the library table, and
#             the two comparisons named after them) are the expressions of
#             methylKit's source evaluated in this script: the percentile of
#             the rows that methRead imported, and the largest median among the
#             libraries divided by the library's own median, from the rows that
#             filterByCoverage kept. That methylKit applied these values is
#             shown by other comparisons: the rows it kept, the rows it
#             removed, and its normalized counts.
#           - How a row is judged (column rule). "identical": whole numbers,
#             positions and sets; any difference fails. "tolerance": floating
#             point values (cutoffs, scaling factors, differences, p-values,
#             q-values); the row fails when the largest relative difference,
#             |a - b| / max(|a|, |b|), exceeds float_tolerance (1e-12). Column
#             n_different counts the values that are not bit for bit equal, so
#             a row "within tolerance" is visible as such. "recorded": a
#             difference that is expected; it is put on record and not judged.
#           - The extract. Regions are listed under Settings with the reason
#             for each. Every line of a region is taken from every library, in
#             the order of the file, and the first line of a file stays its
#             first line. The position of the first line of each file is taken
#             from every library, so that a first-line CpG can enter the united
#             table where the data allow it. The percentile cutoffs of the
#             extract are those of the extract files, not of the full files;
#             both routes see the same files. The extract table gives, for each
#             region, what the reason rests on (lines per library, median reads
#             of the region relative to the median of the whole coverage file,
#             lines that the upper filter removes, files that begin there) and
#             the share of the CpGs tested by step 05 that lie in it.
#           - The helper route of part 3 repeats the chain of helper calls of
#             fit_mode() in 01_pipeline/05_differential_methylation.R, because
#             step 05 does not offer that chain as a function. If the chain of
#             step 05 changes, this script must follow. Part 5 fails when a
#             table that step 05 wrote no longer holds the raw counts scaled by
#             the helper, or holds, inside the regions, other CpGs than the
#             rule of its mode keeps.
#           - Both routes run the two tests with the covariates as step 05
#             passes them (factors), so the routes differ only in how the table
#             of counts is built. The legacy script passed character vectors
#             (yes / no, D / T). methylKit turns a factor covariate into its
#             integer codes and a character covariate into an indicator. For
#             two levels that is the same model, fitted by different
#             arithmetic: the methylation difference does not depend on the
#             covariate, and the p-values need not be the same bit for bit.
#             Part 4 therefore fits methylKit's route once more per contrast
#             with the legacy coding and sets it against the helper route with
#             the coding of step 05. The methylation differences and the DML
#             lists must agree; the differences between the p-values and
#             between the q-values are recorded. On the full data the DML lists
#             of legacy mode are compared with the lists of the draft in
#             check_legacy_outputs.R.
#           - glm.fit warns at CpGs where a fitted proportion is numerically 0
#             or 1. Both routes fit the same counts, so the warnings do not
#             bear on the comparison. They are counted by message and written
#             to the log, not printed one by one.
#           - getMethylDiff of methylKit 1.28.0 warns ("no non-missing
#             arguments to max; returning -Inf") when no CpG passes its
#             thresholds, because its row selection takes max() of an empty
#             index. The selection is then empty, as it should be. The warnings
#             of the DML selections are counted by message and logged for each
#             mode, beside the number of selections that returned no CpG.
#           - Full libraries (part 2): the first library of each of the four
#             groups in the order of the sample sheet, every library whose
#             first line has at least config$min_coverage reads (the first-line
#             rule then changes the rows kept), and every library whose upper
#             cutoff in the record of step 05 is not a whole number. The part
#             runs before the routes on the extract because methRead on a full
#             file needs the most memory of all parts; the peak of the run is
#             reached here. methylKit's normalizeCoverage and unite are not run
#             on full libraries: they need all 22 at once. Their arithmetic is
#             checked on the extract (all 22 libraries), on the synthetic
#             counts, and against the tables of step 05 (part 5).
#           - Part 5 does not use methylKit. Counts: at every CpG of
#             united_cpg_<mode>.csv.gz and in every library, the three counts
#             of the table must be the counts of the line of the whole coverage
#             file, scaled by the helper with the factor max(median) / median
#             of the medians that step 05 recorded. This is compared during the
#             pass of part 1, while a whole file is in memory, and once more
#             inside the regions on the lines of the extract files. CpGs:
#             inside the regions, where the extract holds every line of every
#             library, the CpGs of the table must be the CpGs with at least
#             config$min_coverage reads in every library and, in corrected
#             mode, fewer reads than the cutoff that step 05 recorded for the
#             whole file (in legacy mode, not on the first line of a file).
#             The CpGs of the table inside the regions are also set against
#             the united CpGs of the extract. In legacy mode the two sets must
#             be equal. In corrected mode they need not be, because the
#             cutoffs of the extract are not those of the whole files; the
#             numbers are recorded. Which CpGs each mode keeps on the whole
#             files is checked in coverage_filter_variants.R, which derives
#             both sets from the 22 files and stops unless they are the CpGs
#             of step 05.
#           - Synthetic counts (part 6): synthetic_sets (20) sets of four
#             libraries with synthetic_cpgs (30,001) CpGs each, drawn after
#             set.seed(seed). Per library: coverage from a negative binomial
#             distribution, shifted so that its median is a chosen value; the
#             methylated fraction of each CpG from a beta distribution; the
#             methylated count from a binomial distribution. That is three
#             draws per CpG, plus two draws per set for the medians. In each set
#             the medians are chosen so that the scaling factors are 1, 1.5,
#             1.25 and one arbitrary ratio. With the factor 1.5 every odd
#             coverage scales to a value exactly halfway between two whole
#             numbers. The same counts are also scaled in the order of
#             operations of step 04, to show that the test detects a different
#             order.
#           - Part 7 selects the rows of each library by the rule of step 04
#             (its function covered()), at 1 read and at config$min_coverage
#             reads, the two levels at which step 04 computes correlations. At
#             config$min_coverage the selection must give the united CpGs and
#             the scaling factors of the helper route; both are compared and
#             judged like the comparisons of the other parts.
#           - In the library table, lines_read is the number of lines of the
#             file that both routes read: the extract file in scope "extract",
#             the whole coverage file in scope "full library".
#           - Nothing is written with one row per CpG. The extract files are
#             temporary and are deleted at the end.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "methylKit"))
setDTthreads(config$cores)             # results do not depend on the number of threads
t_start <- Sys.time()

# ---- Settings --------------------------------------------------------------------------
float_tolerance <- 1e-12               # largest relative difference accepted between two floating-point values
seed            <- 20261119L           # for the synthetic counts of part 6; nothing else in this script is random
synthetic_sets  <- 20L                 # sets of four synthetic libraries; the libraries of a set are normalized together
synthetic_cpgs  <- 30001L              # CpGs per synthetic library; odd, so that the median coverage is one of the counts
synthetic_top   <- c(15L, 30L, 45L, 60L)   # the largest median coverage of a set is drawn from these values
near_half       <- 1e-9                # a value this close to a half, and not exactly a half, is counted as a near tie

# Regions of the extract. `from` and `to` are 0-based start coordinates as in
# the coverage files; a region holds the CpGs with from <= start < to. Inf means
# the end of the sequence.
locus_start <- 39633323L               # first DML CpG of the high-coverage region of issue I-09 (locus_nc047563_summary.csv)
locus_flank <- 150000L
regions <- data.table(
  region   = c("first 6 Mb of linkage group NC_047559.1", "300 kb around the locus of issue I-09", "mitochondrial sequence",
               "unplaced scaffold NW_022994810.1", "unplaced scaffold NW_022994829.1", "unplaced scaffold NW_022994852.1",
               "unplaced scaffold NW_022994895.1", "unplaced scaffold NW_022994818.1"),
  sequence = c("NC_047559.1", "NC_047563.1", config$mito_sequence,
               "NW_022994810.1", "NW_022994829.1", "NW_022994852.1", "NW_022994895.1", "NW_022994818.1"),
  from     = c(0, locus_start - locus_flank, 0, 0, 0, 0, 0, 0),
  to       = c(6e6, locus_start + locus_flank, Inf, Inf, Inf, Inf, Inf, Inf),
  purpose  = c("ordinary coverage; most CpGs of the extract; several files begin here",
               "a short block with coverage far above the library median in some libraries, inside ordinary sequence",
               "coverage far above the library median in every library",
               "coverage several times the library median along the whole scaffold",
               "a scaffold with many tested CpGs; scaffold names sort after the linkage groups",
               "a scaffold with many tested CpGs; scaffold names sort after the linkage groups",
               "a scaffold on which a file begins",
               "a sequence that some libraries do not have"))
stopifnot(!anyDuplicated(regions$sequence), all(regions$from < regions$to), all(regions$sequence[1:2] %in% config$linkage_groups))

out_file <- function(suffix = "") here("06_diagnostics", paste0("check_methylkit_equivalence", suffix, ".csv"))
tab05 <- function(file) output_path("05_differential_methylation", "tables", file)
# Stop with a message unless a statement about the inputs holds.
must <- function(ok, ...) if (!isTRUE(ok)) stop(..., call. = FALSE)
is_half <- function(x) x %% 1 == 0.5                                   # exactly halfway between two whole numbers
is_near_half <- function(x) abs(x %% 1 - 0.5) < near_half & !is_half(x)
count_text <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# ---- Recording comparisons --------------------------------------------------------------
comparisons <- list()
# Add one row to the result table. `rule`: "identical" (any difference fails),
# "tolerance" (fails above float_tolerance), "recorded" (not judged). A row
# without items is marked "nothing to compare" and does not count as agreement.
record <- function(part, mode_name, quantity, item, n_compared, n_different, largest_absolute = NA_real_,
                   largest_relative = NA_real_, rule = "identical", note = "") {
  must(rule %in% c("identical", "tolerance", "recorded"), "Unknown rule: ", rule)
  must(length(n_different) == 1 && !is.na(n_different), "A comparison without a count of differences: ", quantity)
  result <- if (rule == "recorded") "recorded"
            else if (n_compared == 0) "nothing to compare"
            else if (n_different == 0) "identical"
            else if (rule == "tolerance" && is.finite(largest_relative) && largest_relative <= float_tolerance) "within tolerance"
            else "DIFFERENT"
  comparisons[[length(comparisons) + 1L]] <<- data.table(
    part = part, mode = mode_name, quantity = quantity, item = item,
    n_compared = as.numeric(n_compared), n_different = as.numeric(n_different),
    largest_absolute_difference = as.numeric(largest_absolute), largest_relative_difference = as.numeric(largest_relative),
    rule = rule, result = result, note = note)
  invisible(result)
}
# Two vectors compared item by item: number of items, number that are not equal,
# and for numbers the largest absolute and the largest relative difference.
compare_values <- function(a, b) {
  must(length(a) == length(b) && !anyNA(a) && !anyNA(b), "compare_values needs two complete vectors of one length.")
  out <- list(n = length(a), n_different = sum(a != b), largest_absolute = NA_real_, largest_relative = NA_real_)
  if (is.numeric(a) && is.numeric(b)) {
    gap <- abs(a - b); size <- pmax(abs(a), abs(b))
    out$largest_absolute <- if (length(gap) > 0) max(gap) else 0
    out$largest_relative <- if (length(gap) > 0) max(ifelse(size > 0, gap / size, 0)) else 0
  }
  out
}
record_values <- function(part, mode_name, quantity, item, a, b, rule = "identical", note = "") {
  v <- compare_values(a, b)
  record(part, mode_name, quantity, item, v$n, v$n_different, v$largest_absolute, v$largest_relative, rule, note)
}
# Two sets: size of the union, and the number of members that are in one set only.
record_sets <- function(part, mode_name, quantity, item, a, b, rule = "identical", note = "") {
  record(part, mode_name, quantity, item, length(union(a, b)), length(setdiff(a, b)) + length(setdiff(b, a)), rule = rule, note = note)
}

# ---- Inputs --------------------------------------------------------------------------------
lib <- read_libraries(); n_lib <- nrow(lib)
treat_ploidy <- as.integer(lib$ploidy == "triploid")     # 0 diploid, 1 triploid, as in step 05
treat_heat   <- as.integer(lib$treatment == "heat")      # 0 control, 1 heat
sample_ids   <- as.character(lib$legacy_sample_index)
manifest <- read_gannet_manifest()
cov_file <- vapply(lib$cov_file, fetch_gannet_file, character(1), config = config, manifest = manifest, USE.NAMES = FALSE)

# What step 05 recorded for each library.
records05 <- lapply(setNames(config$modes, config$modes), function(mode_name) {
  recorded <- fread(tab05(sprintf("library_filtering_%s.csv", mode_name)))
  must(identical(recorded$library_id, lib$library_id), "library_filtering_", mode_name, ".csv does not list the libraries of the sample sheet in its order.")
  recorded
})
# The scaling factor as step 05 computes it. The factor written in the table is
# rounded to 15 digits; the medians are written exactly.
factor05 <- lapply(records05, function(recorded) max(recorded$median_coverage) / recorded$median_coverage)

# The tables of tested CpGs that step 05 wrote, for part 5. `keys05` numbers
# their CpGs over the sequences of the tables; `all_keys05` holds every CpG of
# the tables once, and `index05` the place of each CpG of a table in it.
tables05 <- lapply(setNames(config$modes, config$modes), function(mode_name) {
  united05 <- fread(tab05(sprintf("united_cpg_%s.csv.gz", mode_name)))
  must(ncol(united05) == 4 + 3 * n_lib &&
         all(paste0(rep(c("coverage", "numCs", "numTs"), n_lib), rep(seq_len(n_lib), each = 3)) == names(united05)[-(1:4)]),
       "united_cpg_", mode_name, ".csv.gz does not hold coverage, numCs and numTs of the ", n_lib, " libraries.")
  must(!anyDuplicated(united05, by = c("chr", "start")), "united_cpg_", mode_name, ".csv.gz holds a CpG twice.")
  united05
})
levels05 <- unique(unlist(lapply(tables05, function(united05) united05$chr), use.names = FALSE))
keys05 <- lapply(tables05, function(united05) cpg_key(united05$chr, united05$start, levels05))
all_keys05 <- unique(unlist(keys05, use.names = FALSE))
index05 <- lapply(keys05, match, table = all_keys05)
# One library of a table of step 05 against the lines of a coverage file.
# `written`: the table; `n_methylated`, `n_unmethylated`: the two counts on the
# line of each CpG of the table. Returns the number of values compared, the
# number that are not the counts of the line scaled as step 05 scales them, and
# the number that lie exactly halfway between two whole numbers before rounding.
table_against_lines <- function(written, lib_i, n_methylated, n_unmethylated, scaling_factor) {
  reads <- n_methylated + n_unmethylated
  expected <- scale_counts_like_methylkit(reads, n_methylated, n_unmethylated, scaling_factor)
  unrounded <- c(scaling_factor * reads, expected$coverage * (n_methylated / reads), expected$coverage * (n_unmethylated / reads))
  c(values = 3 * length(reads),
    different = sum(vapply(c("coverage", "numCs", "numTs"), function(what) sum(expected[[what]] != written[[paste0(what, lib_i)]]), numeric(1))),
    half = sum(is_half(unrounded)))
}

# All six columns of a coverage file. read_cov() of the pipeline drops the end
# coordinate and the percentage; here they are needed to write lines back
# unchanged, so the percentage is kept as text.
read_whole_lines <- function(path, n_rows = Inf) {
  fread(path, header = FALSE, sep = "\t", nrows = n_rows, showProgress = FALSE,
        colClasses = list(character = c(1, 4), integer = c(2, 3, 5, 6)),
        col.names = c("chr", "start", "end", "percent", "nC", "nT"))
}

# The first line of every file. Its position becomes a region of one CpG, taken
# from every library.
first_line <- rbindlist(lapply(cov_file, read_whole_lines, n_rows = 1L))
first_line[, reads := nC + nT]
must(nrow(first_line) == n_lib && !anyNA(first_line), "A coverage file has no readable first line.")
first_positions <- unique(first_line[, .(sequence = chr, from = as.numeric(start))])
setorder(first_positions, sequence, from)
n_interval_regions <- nrow(regions)
regions <- rbind(regions, first_positions[, .(region = "position of the first line of a file", sequence, from, to = from + 1,
                                              purpose = "the first-line rule of legacy mode")])
extract_levels <- unique(regions$sequence)               # sequence names of the extract, for cpg_key()
# Rows of the extract table: one per region above; the positions of first lines
# that lie outside those regions form one more row.
region_group <- c(seq_len(n_interval_regions), rep(n_interval_regions + 1L, nrow(regions) - n_interval_regions))
n_groups <- n_interval_regions + 1L

part_extract   <- sprintf("extract, all %d libraries", n_lib)
part_coding    <- sprintf("legacy coding of the covariate, extract, all %d libraries", n_lib)
part_tables    <- "tables of step 05"
part_synthetic <- sprintf("synthetic counts (%d libraries of %s CpGs)", 4L * synthetic_sets, count_text(synthetic_cpgs))
part_order     <- function(level) sprintf("order of operations of step 04, CpGs with at least %g read%s in all libraries", level, if (level == 1) "" else "s")
levels_step04  <- unique(c(1, config$min_coverage))      # the two levels at which step 04 computes correlations

# ---- 1. Extract ---------------------------------------------------------------------------
# Index of the first region of `regions` that holds each position; 0 where none does.
region_of <- function(chr, start) {
  found <- integer(length(chr))
  sequence_index <- chmatch(chr, extract_levels, nomatch = 0L)
  for (s in seq_along(extract_levels)) {
    rows <- which(sequence_index == s)
    if (length(rows) == 0L) next
    start_s <- start[rows]; region_s <- integer(length(rows))
    for (r in rev(which(regions$sequence == extract_levels[s]))) region_s[start_s >= regions$from[r] & start_s < regions$to[r]] <- r
    found[rows] <- region_s
  }
  found
}

extract_dir <- file.path(tempdir(), "check_methylkit_equivalence_extract")
dir.create(extract_dir, showWarnings = FALSE)
extract_file <- file.path(extract_dir, lib$cov_file)
lines_by_group <- removed_by_group <- matrix(0L, n_lib, n_groups)
relative_median <- matrix(NA_real_, n_lib, n_groups)
file_facts <- vector("list", n_lib)
whole05 <- lapply(keys05, function(keys) c(values = 0, different = 0, half = 0))   # part 5, whole tables: tally over the libraries
for (lib_i in seq_len(n_lib)) {
  whole <- read_whole_lines(cov_file[lib_i])
  # Part 5, whole tables: the line of this file for every CpG of the tables of
  # step 05, and the counts of the tables against the counts of those lines.
  on_tested <- which(whole$chr %chin% levels05)
  hit <- match(cpg_key(whole$chr[on_tested], whole$start[on_tested], levels05), all_keys05)
  found <- which(!is.na(hit))
  must(!anyDuplicated(hit[found]), lib$cov_file[lib_i], " holds a CpG of the tables of step 05 on two lines.")
  line_of <- integer(length(all_keys05)); line_of[hit[found]] <- on_tested[found]
  for (mode_name in config$modes) {
    at <- line_of[index05[[mode_name]]]
    must(all(at > 0L), "united_cpg_", mode_name, ".csv.gz holds a CpG that has no line in ", lib$cov_file[lib_i], ".")
    whole05[[mode_name]] <- whole05[[mode_name]] +
      table_against_lines(tables05[[mode_name]], lib_i, whole$nC[at], whole$nT[at], factor05[[mode_name]][lib_i])
  }
  rm(on_tested, hit, found, line_of, at)
  library_median <- stats::median(whole$nC + whole$nT)    # median reads over all lines of the file
  region_index <- region_of(whole$chr, whole$start)
  taken <- whole[region_index > 0L]
  must(nrow(whole) > 1 && region_index[1L] > 0L && !anyDuplicated(taken, by = c("chr", "start")),
       "The extract of ", lib$library_id[lib_i], " lacks the first line of the file or holds a position twice.")
  must(!anyNA(taken) && all(taken$end == taken$start + 2L) && all(taken$nC >= 0L) && all(taken$nT >= 0L) && all(taken$nC + taken$nT >= 1L),
       lib$cov_file[lib_i], " is not a merged CpG coverage file: six columns, end = start + 2, at least one read per line.")
  fwrite(taken, extract_file[lib_i], sep = "\t", col.names = FALSE, quote = FALSE)
  # The extract must hold lines of the original: the same six fields when read
  # back, and a first line that is the first line of the file, byte for byte.
  must(identical(as.list(read_whole_lines(extract_file[lib_i])), as.list(taken)) &&
         identical(readLines(extract_file[lib_i], n = 1L), readLines(cov_file[lib_i], n = 1L)),
       "The extract of ", lib$library_id[lib_i], " does not reproduce the lines of the coverage file.")
  # What each region contributes in this library: lines, median reads relative
  # to the median of the whole coverage file, and lines that corrected mode
  # removes (at least min_coverage reads and at least the cutoff of the extract file).
  group_taken <- region_group[region_index[region_index > 0L]]
  reads <- taken$nC + taken$nT
  cutoff <- upper_coverage_cutoff(reads, config$import_min_coverage, config$upper_coverage_percentile)
  lines_by_group[lib_i, ] <- tabulate(group_taken, nbins = n_groups)
  relative_median[lib_i, ] <- as.numeric(tapply(reads, factor(group_taken, levels = seq_len(n_groups)), stats::median)) / library_median
  removed_by_group[lib_i, ] <- tabulate(group_taken[reads >= config$min_coverage & reads >= cutoff], nbins = n_groups)
  file_facts[[lib_i]] <- data.table(library_id = lib$library_id[lib_i], lines_in_file = nrow(whole), lines_in_extract = nrow(taken))
  log_msg("%-10s %9d lines, %6d in the extract", lib$library_id[lib_i], nrow(whole), nrow(taken))
  rm(whole, taken, region_index, group_taken, reads); gc(verbose = FALSE)
}
file_facts <- rbindlist(file_facts)
must(all(colSums(lines_by_group)[seq_len(n_interval_regions)] > 0), "A region of the extract holds no line in any library.")
# Of the tables of step 05, the later parts need the rows inside the regions only.
n_tested05 <- lapply(tables05, nrow)
tested05_by_region <- lapply(tables05, function(united05) {
  inside <- region_of(united05$chr, united05$start)
  tabulate(inside[inside > 0L], nbins = nrow(regions))
})
tables05 <- lapply(tables05, function(united05) {
  inside <- region_of(united05$chr, united05$start) > 0L
  united05[inside]
})
rm(keys05, all_keys05, index05); gc(verbose = FALSE)

# ---- Facts about one library, from each route ---------------------------------------------
# What methylKit did to one library. `imported` and `kept` are the tables of its
# objects before and after filterByCoverage.
facts_of_methylkit <- function(imported, kept, mode_name, first_chr, first_start) {
  upper <- mode_name == "corrected"
  enough <- imported$coverage >= config$min_coverage
  # filterByCoverage does not return its cutoff. This is the expression in its
  # source, evaluated here on the rows that methRead imported; the rows that
  # filterByCoverage kept must be the rows below this value.
  cutoff <- if (upper) unname(stats::quantile(imported$coverage, config$upper_coverage_percentile / 100)) else NA_real_
  must(nrow(kept) == if (upper) sum(enough & imported$coverage < cutoff) else sum(enough),
       "filterByCoverage did not keep the rows that the rule in its source keeps.")
  data.table(
    rows_imported_methylkit = nrow(imported),
    first_line_imported_methylkit = any(imported$chr == first_chr & imported$start == first_start),
    rows_min_coverage_methylkit = sum(enough),
    upper_cutoff_from_methylkit_rows = cutoff,
    rows_with_cutoff_reads = if (upper) sum(imported$coverage == cutoff) else NA_integer_,
    rows_removed_by_upper_methylkit = sum(enough) - nrow(kept),
    rows_kept_methylkit = nrow(kept),
    first_line_kept_methylkit = any(kept$chr == first_chr & kept$start == first_start),
    median_coverage_methylkit = stats::median(kept$coverage),
    max_coverage_kept_methylkit = max(kept$coverage))
}
# What read_filtered_counts() did to the same library (`helper` is its result),
# and the number of rows it kept that differ from the rows methylKit kept.
facts_of_helper <- function(helper, kept, mode_name, first_chr, first_start) {
  must(helper$record$n_cpg_kept == nrow(helper$counts), "read_filtered_counts reports a number of kept rows that it did not return.")
  position <- key_to_position(helper$counts$pos_key, helper$chr_levels)
  rows_differing <- if (nrow(kept) == nrow(position)) {
    sum(kept$chr != position$chr | kept$start != position$start | kept$coverage != helper$counts$coverage |
          kept$numCs != helper$counts$nC | kept$numTs != helper$counts$nT)
  } else {                                               # not the same number of rows: count the rows that only one route has
    rows_methylkit <- paste(kept$chr, kept$start, kept$numCs, kept$numTs)
    rows_helper <- paste(position$chr, position$start, helper$counts$nC, helper$counts$nT)
    length(setdiff(rows_methylkit, rows_helper)) + length(setdiff(rows_helper, rows_methylkit))
  }
  data.table(
    rows_imported_helper = helper$record$n_cpg_imported,
    rows_min_coverage_helper = helper$record$n_cpg_min_coverage,
    upper_cutoff_helper = if (mode_name == "corrected") helper$record$upper_cutoff else NA_real_,
    rows_removed_by_upper_helper = helper$record$n_removed_by_upper,
    rows_kept_helper = helper$record$n_cpg_kept,
    kept_rows_differing = rows_differing,
    first_line_kept_helper = any(position$chr == first_chr & position$start == first_start),
    median_coverage_helper = helper$record$median_coverage,
    max_coverage_kept_helper = helper$record$max_coverage_kept)
}
# The rows of the result table that compare the two routes library by library.
record_library_facts <- function(part, mode_name, facts) {
  compare <- function(quantity, a, b, rule = "identical", note = "") record_values(part, mode_name, quantity, "library", a, b, rule, note)
  compare("rows imported: methRead against read_filtered_counts", facts$rows_imported_methylkit, facts$rows_imported_helper,
          note = sprintf("the first line of the file is among the rows that methRead imported in %d of %d libraries",
                         sum(facts$first_line_imported_methylkit), nrow(facts)))
  compare("rows with at least min_coverage reads", facts$rows_min_coverage_methylkit, facts$rows_min_coverage_helper)
  if (mode_name == "corrected") {
    compare(paste("upper cutoff: percentile of the rows that methRead imported, which is the expression in the source of filterByCoverage",
                  "evaluated in this script, against upper_coverage_cutoff"),
            facts$upper_cutoff_from_methylkit_rows, facts$upper_cutoff_helper, rule = "tolerance",
            note = sprintf("cutoffs that are whole numbers: %d; that are not: %d; rows with exactly as many reads as their cutoff, which both routes remove: %s",
                           sum(facts$upper_cutoff_from_methylkit_rows %% 1 == 0), sum(facts$upper_cutoff_from_methylkit_rows %% 1 != 0),
                           count_text(sum(facts$rows_with_cutoff_reads))))
    compare("rows removed by the upper filter", facts$rows_removed_by_upper_methylkit, facts$rows_removed_by_upper_helper,
            note = sprintf("%s rows removed in all", count_text(sum(facts$rows_removed_by_upper_methylkit))))
  }
  compare("rows kept: filterByCoverage against read_filtered_counts", facts$rows_kept_methylkit, facts$rows_kept_helper)
  record(part, mode_name, "rows kept, row by row: position, coverage, methylated and unmethylated count", "row",
         sum(pmax(facts$rows_kept_methylkit, facts$rows_kept_helper)), sum(facts$kept_rows_differing))
  compare("first line of the file among the rows kept", facts$first_line_kept_methylkit, facts$first_line_kept_helper,
          note = sprintf("kept in %d of %d libraries", sum(facts$first_line_kept_methylkit), nrow(facts)))
  compare("median coverage of the rows kept", facts$median_coverage_methylkit, facts$median_coverage_helper)
}

# ---- 2. Full libraries --------------------------------------------------------------------
full_set <- which(!duplicated(lib$group) | first_line$reads >= config$min_coverage | records05[[1]]$upper_cutoff %% 1 != 0)
part_full <- sprintf("full libraries (%d of %d)", length(full_set), n_lib)
log_msg("full libraries: %s", paste(lib$library_id[full_set], collapse = ", "))
full_facts <- list()
for (mode_name in config$modes) {
  rows <- list()
  for (lib_i in full_set) {
    # methylKit's own calls, written out for each mode; no helper decides them.
    # Legacy mode leaves methRead at its default header = TRUE.
    read_args <- list(location = cov_file[lib_i], sample.id = sample_ids[lib_i], assembly = "Roslin_v1",
                      pipeline = "bismarkCoverage", mincov = config$import_min_coverage)
    if (mode_name == "corrected") read_args$header <- FALSE
    imported <- suppressMessages(do.call(methRead, read_args))
    gc(verbose = FALSE)
    kept <- getData(if (mode_name == "legacy") filterByCoverage(imported, lo.count = config$min_coverage)
                    else filterByCoverage(imported, lo.count = config$min_coverage, hi.perc = config$upper_coverage_percentile))
    by_methylkit <- facts_of_methylkit(getData(imported), kept, mode_name, first_line$chr[lib_i], first_line$start[lib_i])
    rm(imported); gc(verbose = FALSE)                    # the import is the largest object; it is not needed below
    helper <- read_filtered_counts(cov_file[lib_i], mode_name, config, character(0))
    gc(verbose = FALSE)
    by_helper <- facts_of_helper(helper, kept, mode_name, first_line$chr[lib_i], first_line$start[lib_i])
    rows[[length(rows) + 1L]] <- cbind(
      data.table(scope = "full library", mode = mode_name, library_id = lib$library_id[lib_i],
                 lines_read = file_facts$lines_in_file[lib_i], first_line_reads = first_line$reads[lib_i]),
      by_methylkit, by_helper)
    log_msg("[%s] %-10s full library: %d rows kept by methylKit, %d by the helper, %d rows differ", mode_name, lib$library_id[lib_i],
            by_methylkit$rows_kept_methylkit, by_helper$rows_kept_helper, by_helper$kept_rows_differing)
    rm(kept, helper, by_methylkit, by_helper); gc(verbose = FALSE)
  }
  facts <- rbindlist(rows)
  record_library_facts(part_full, mode_name, facts)
  # methylKit on the full file against what step 05 recorded for the library.
  recorded <- records05[[mode_name]][match(facts$library_id, library_id)]
  pairs <- c(n_cpg_imported = "rows_imported_methylkit", n_cpg_min_coverage = "rows_min_coverage_methylkit",
             n_removed_by_upper = "rows_removed_by_upper_methylkit", n_cpg_kept = "rows_kept_methylkit",
             median_coverage = "median_coverage_methylkit", max_coverage_kept = "max_coverage_kept_methylkit")
  record_values(part_full, mode_name,
                sprintf("methylKit on the full file against library_filtering_%s.csv of step 05: %s", mode_name, paste(names(pairs), collapse = ", ")),
                "value", as.numeric(unlist(lapply(pairs, function(column) facts[[column]]))),
                as.numeric(unlist(lapply(names(pairs), function(column) recorded[[column]]))))
  if (mode_name == "corrected") {
    record_values(part_full, mode_name,
                  "percentile of the rows that methRead imported from the full file against upper_cutoff of library_filtering_corrected.csv of step 05",
                  "value", facts$upper_cutoff_from_methylkit_rows, recorded$upper_cutoff, rule = "tolerance",
                  note = "the table holds 15 significant digits")
  }
  full_facts[[mode_name]] <- facts
}
if (all(c("legacy", "corrected") %in% config$modes)) {
  # I-38 on full files: the import with header = FALSE holds one row more than
  # the default import where the first line has at least import_min_coverage
  # reads, and the same number of rows where it has fewer.
  record_values(part_full, "legacy", "rows that methRead imports with header = FALSE minus rows with its default header = TRUE: 1 or 0, by the reads of the first line",
                "library", full_facts$corrected$rows_imported_methylkit - full_facts$legacy$rows_imported_methylkit,
                as.integer(first_line$reads[full_set] >= config$import_min_coverage),
                note = sprintf("the first line is among the rows imported with the default in %d libraries",
                               sum(full_facts$legacy$first_line_imported_methylkit)))
}

# ---- 3. Both routes on the extract ------------------------------------------------------
# methylKit's own route for one treatment vector. Legacy mode leaves methRead at
# its default header = TRUE and applies no upper filter; corrected mode reads
# every line and passes the percentile under its real name. Both calls are
# written out here, so that this route uses no helper of the route it is
# compared with.
methylkit_route <- function(mode_name, treatment_vector) {
  read_args <- list(location = as.list(extract_file), sample.id = as.list(sample_ids), assembly = "Roslin_v1",
                    treatment = treatment_vector, pipeline = "bismarkCoverage", mincov = config$import_min_coverage)
  if (mode_name == "corrected") read_args$header <- FALSE
  imported <- suppressMessages(do.call(methRead, read_args))
  filtered <- if (mode_name == "legacy") filterByCoverage(imported, lo.count = config$min_coverage)
              else filterByCoverage(imported, lo.count = config$min_coverage, hi.perc = config$upper_coverage_percentile)
  normalized <- normalizeCoverage(filtered)
  list(imported = imported, filtered = filtered, normalized = normalized,
       united = suppressMessages(unite(normalized, destrand = FALSE)))
}

# The helper route, chained as in fit_mode() of 01_pipeline/05_differential_methylation.R.
helper_route <- function(mode_name) {
  kept <- vector("list", n_lib); records <- vector("list", n_lib)
  chr_levels <- character(0); common <- NULL
  for (i in seq_len(n_lib)) {
    r <- read_filtered_counts(extract_file[i], mode_name, config, chr_levels)
    chr_levels <- r$chr_levels
    common <- if (is.null(common)) r$counts$pos_key else common[common %in% r$counts$pos_key]
    kept[[i]] <- r$counts
    records[[i]] <- cbind(library_id = lib$library_id[i], r$record)
  }
  filtering <- rbindlist(records)
  filtering[, scaling_factor := max(median_coverage) / median_coverage]
  pos <- key_to_position(common, chr_levels)
  objects <- lapply(seq_len(n_lib), function(i) {
    k <- kept[[i]][match(common, pos_key)]
    make_methylraw(pos$chr, pos$start, scale_counts_like_methylkit(k$coverage, k$nC, k$nT, filtering$scaling_factor[i]),
                   sample_ids[i])
  })
  raw_list <- new("methylRawList", objects, treatment = treat_ploidy)
  meth_ploidy <- suppressMessages(unite(raw_list, destrand = FALSE))
  meth_heat <- reorganize(meth_ploidy, sample.ids = meth_ploidy@sample.ids, treatment = treat_heat)
  list(kept = kept, chr_levels = chr_levels, filtering = filtering, united = list(ploidy = meth_ploidy, heat = meth_heat))
}

# The covariate of each contrast as step 05 passes it (a factor), and as the
# legacy script passed it: a character vector, "yes" and "no" for heat, "D" and
# "T" for ploidy (08_archive/code/2_WGBS_Methylkit-all-gigas.R, lines 143 to
# 150, where the vectors are typed in by position). Here they come from the
# sample sheet and must be the vectors of that script.
covariates <- list(ploidy = data.frame(heat = lib$treatment), heat = data.frame(ploidy = lib$ploidy))
covariates_legacy <- list(ploidy = data.frame(HEAT = ifelse(lib$treatment == "heat", "yes", "no"), stringsAsFactors = FALSE),
                          heat = data.frame(PLOIDY = ifelse(lib$ploidy == "triploid", "T", "D"), stringsAsFactors = FALSE))
must(identical(covariates_legacy$ploidy$HEAT, rep(c("yes", "no"), c(10, 12))) &&
       identical(covariates_legacy$heat$PLOIDY, rep(c("D", "T", "D", "T"), c(5, 5, 6, 6))),
     "The sample sheet does not give the covariate vectors of the legacy script.")
coding_pair <- "methylKit's route with the legacy coding against the helper route with the coding of step 05"

# methylKit's test with one covariate table. Warnings are counted by message
# and logged.
run_test <- function(united, covariate_table, label) {
  seen <- character(0)
  fit <- withCallingHandlers(
    suppressMessages(calculateDiffMeth(united, covariates = covariate_table, mc.cores = config$cores)),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") })
  tally <- table(seen)
  log_msg("%s: %d CpGs tested; warnings: %s", label, nrow(getData(fit)),
          if (length(seen) == 0) "none" else paste(sprintf("%d x '%s'", as.integer(tally), names(tally)), collapse = "; "))
  fit
}
# A DML selection (getMethylDiff, or dml_table() of the helpers, which calls
# it), with its warnings collected. getMethylDiff warns when no CpG passes (see
# Notes). The selections made, those without a CpG and the warnings are logged
# for each mode.
selections <- list(made = 0L, empty = 0L, warnings = character(0))
quiet_selection <- function(selected) {
  result <- withCallingHandlers(selected, warning = function(w) {
    selections$warnings <<- c(selections$warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  selections$made <<- selections$made + 1L
  selections$empty <<- selections$empty + (nrow(result) == 0L)
  result
}

thresholds <- sort(unique(c(config$dml_difference, config$dml_difference_sweep)))
extract_facts <- list(); united_keys <- list(); imported_keys <- list(); tested_by_region <- list(); misspelled <- NULL
for (mode_name in config$modes) {
  log_msg("[%s] methylKit's route and the helper route on the extract", mode_name)
  native <- methylkit_route(mode_name, treat_ploidy)
  native_united <- list(ploidy = native$united, heat = methylkit_route(mode_name, treat_heat)$united)
  helper <- helper_route(mode_name)

  # -- every library
  facts <- rbindlist(lapply(seq_len(n_lib), function(lib_i) {
    kept <- getData(native$filtered[[lib_i]])
    cbind(data.table(scope = "extract", mode = mode_name, library_id = lib$library_id[lib_i],
                     lines_read = file_facts$lines_in_extract[lib_i], first_line_reads = first_line$reads[lib_i]),
          facts_of_methylkit(getData(native$imported[[lib_i]]), kept, mode_name, first_line$chr[lib_i], first_line$start[lib_i]),
          facts_of_helper(list(counts = helper$kept[[lib_i]], record = helper$filtering[lib_i], chr_levels = helper$chr_levels),
                          kept, mode_name, first_line$chr[lib_i], first_line$start[lib_i]))
  }))
  # normalizeCoverage does not return its factors. This is the expression in its
  # source, evaluated here on the medians of the rows that filterByCoverage
  # kept; the normalized counts compared below show that it was applied.
  facts[, scaling_factor_from_methylkit_medians := max(median_coverage_methylkit) / median_coverage_methylkit]
  facts[, scaling_factor_helper := helper$filtering$scaling_factor]
  record_library_facts(part_extract, mode_name, facts)
  record_values(part_extract, mode_name,
                paste("scaling factor: largest median over own median, which is the expression in the source of normalizeCoverage",
                      "evaluated in this script on the medians of the rows that filterByCoverage kept, against the factor of the helper route"),
                "library", facts$scaling_factor_from_methylkit_medians, facts$scaling_factor_helper, rule = "tolerance",
                note = sprintf("%d distinct factors, from %s to %s", uniqueN(facts$scaling_factor_from_methylkit_medians),
                               format(min(facts$scaling_factor_from_methylkit_medians), digits = 6),
                               format(max(facts$scaling_factor_from_methylkit_medians), digits = 6)))

  # -- normalized counts of every kept row of every library, not only of the united CpGs
  scaled <- rbindlist(lapply(seq_len(n_lib), function(lib_i) {
    by_methylkit <- getData(native$normalized[[lib_i]]); k <- helper$kept[[lib_i]]; f <- helper$filtering$scaling_factor[lib_i]
    n_values <- 3 * max(nrow(by_methylkit), nrow(k))
    if (nrow(by_methylkit) != nrow(k)) return(data.table(n = n_values, different = n_values, half = 0, near = 0))
    by_helper <- scale_counts_like_methylkit(k$coverage, k$nC, k$nT, f)
    unrounded <- c(f * k$coverage, by_helper$coverage * (k$nC / k$coverage), by_helper$coverage * (k$nT / k$coverage))
    data.table(n = n_values,
               different = sum(by_methylkit$coverage != by_helper$coverage) + sum(by_methylkit$numCs != by_helper$numCs) +
                 sum(by_methylkit$numTs != by_helper$numTs),
               half = sum(is_half(unrounded)), near = sum(is_near_half(unrounded)))
  }))
  record(part_extract, mode_name,
         "normalized coverage, numCs and numTs of every kept row: normalizeCoverage against scale_counts_like_methylkit",
         "count", sum(scaled$n), sum(scaled$different),
         note = sprintf("%s values are exactly halfway between two whole numbers before rounding and %s lie within %g of a half",
                        count_text(sum(scaled$half)), count_text(sum(scaled$near)), near_half))

  # -- united tables
  for (contrast_name in names(native_united)) {
    a <- getData(native_united[[contrast_name]]); b <- getData(helper$united[[contrast_name]])
    must(identical(names(a), names(b)) && ncol(a) == 4 + 3 * n_lib, "The united tables of the two routes do not have the same columns.")
    key_a <- cpg_key(a$chr, a$start, extract_levels); key_b <- cpg_key(b$chr, b$start, extract_levels)
    label <- sprintf("united table, treatment vector of the %s contrast", contrast_name)
    record_sets(part_extract, mode_name, paste0(label, ": CpGs"), "CpG", key_a, key_b,
                note = sprintf("%s CpGs by methylKit, %s by the helpers", count_text(length(key_a)), count_text(length(key_b))))
    both <- intersect(key_a, key_b); row_a <- match(both, key_a); row_b <- match(both, key_b)
    same_frame <- identical(key_a, key_b) && identical(a$end, b$end) && identical(a$strand, b$strand) &&
      identical(native_united[[contrast_name]]@sample.ids, helper$united[[contrast_name]]@sample.ids) &&
      isTRUE(all.equal(as.numeric(native_united[[contrast_name]]@treatment), as.numeric(helper$united[[contrast_name]]@treatment)))
    record(part_extract, mode_name, paste0(label, ": order of the rows, end coordinate, strand, sample names, treatment vector"),
           "table", 1, as.integer(!same_frame))
    count_columns <- names(a)[-(1:4)]
    record_values(part_extract, mode_name, paste0(label, ": coverage, numCs and numTs of every library at every CpG"), "count",
                  unlist(a[row_a, count_columns], use.names = FALSE), unlist(b[row_b, count_columns], use.names = FALSE))
  }
  united_now <- getData(helper$united$ploidy)
  united_keys[[mode_name]] <- cpg_key(united_now$chr, united_now$start, extract_levels)
  tested_by_region[[mode_name]] <- tabulate(region_of(united_now$chr, united_now$start), nbins = nrow(regions))
  imported_keys[[mode_name]] <- lapply(seq_len(n_lib), function(lib_i) {
    x <- getData(native$imported[[lib_i]]); cpg_key(x$chr, x$start, extract_levels)
  })
  if (mode_name == "legacy") {
    # The filter as the legacy scripts called it. methylKit knows hi.count and
    # hi.perc; high.count and high.perc are ignored (issue I-05).
    as_called <- filterByCoverage(native$imported, lo.count = config$min_coverage, lo.perc = NULL,
                                  high.count = NULL, high.perc = config$upper_coverage_percentile)
    misspelled <- vapply(seq_len(n_lib), function(lib_i) identical(getData(as_called[[lib_i]]), getData(native$filtered[[lib_i]])), logical(1))
    rm(as_called)
  }

  # -- tests
  dml <- list(); dml_legacy_coding <- list()
  for (contrast_name in names(native_united)) {
    fit_a <- run_test(native_united[[contrast_name]], covariates[[contrast_name]], sprintf("[%s] %s, methylKit's route", mode_name, contrast_name))
    fit_b <- run_test(helper$united[[contrast_name]], covariates[[contrast_name]], sprintf("[%s] %s, helper route", mode_name, contrast_name))
    a <- as.data.table(getData(fit_a)); b <- as.data.table(getData(fit_b))
    key_a <- cpg_key(a$chr, a$start, extract_levels); key_b <- cpg_key(b$chr, b$start, extract_levels)
    both <- intersect(key_a, key_b); row_a <- match(both, key_a); row_b <- match(both, key_b)
    label <- sprintf("%s contrast: ", contrast_name)
    record_values(part_extract, mode_name, paste0(label, "methylation difference (meth.diff)"), "CpG", a$meth.diff[row_a], b$meth.diff[row_b],
                  rule = "tolerance")
    record_values(part_extract, mode_name, paste0(label, "p-value"), "CpG", a$pvalue[row_a], b$pvalue[row_b], rule = "tolerance")
    record_values(part_extract, mode_name, paste0(label, "q-value (SLIM)"), "CpG", a$qvalue[row_a], b$qvalue[row_b], rule = "tolerance")
    # DML: methylKit's getMethylDiff on its own route; the rule as step 05
    # writes it on the helper route, checked against dml_table() as in step 05.
    by_threshold <- lapply(thresholds, function(threshold) {
      called_a <- quiet_selection(getData(getMethylDiff(fit_a, difference = threshold, qvalue = config$dml_qvalue)))
      called_b <- b[qvalue < config$dml_qvalue & abs(meth.diff) > threshold]
      must(nrow(called_b) == nrow(quiet_selection(dml_table(fit_b, threshold, config$dml_qvalue))),
           "The DML rule of step 05 and dml_table() disagree.")
      list(a = cpg_key(called_a$chr, called_a$start, extract_levels), b = cpg_key(called_b$chr, called_b$start, extract_levels))
    })
    names(by_threshold) <- thresholds
    configured <- by_threshold[[as.character(config$dml_difference)]]
    record_sets(part_extract, mode_name,
                sprintf("%sDML at q < %g and a difference above %g points", label, config$dml_qvalue, config$dml_difference),
                "DML", configured$a, configured$b)
    record(part_extract, mode_name,
           sprintf("%sDML at q < %g, at each difference threshold of the sweep (%s points)", label, config$dml_qvalue,
                   paste(range(config$dml_difference_sweep), collapse = " to ")),
           "DML", sum(vapply(by_threshold, function(x) length(union(x$a, x$b)), numeric(1))),
           sum(vapply(by_threshold, function(x) length(setdiff(x$a, x$b)) + length(setdiff(x$b, x$a)), numeric(1))),
           note = paste("DML by threshold:",
                        paste(sprintf("%s: %d", names(by_threshold), vapply(by_threshold, function(x) length(x$a), integer(1))), collapse = "; ")))
    dml[[contrast_name]] <- configured

    # Part 4, the coding of the covariate: methylKit's route once more, with
    # the covariate as the legacy script passed it, against the helper route
    # with the coding of step 05 (fit_b).
    fit_legacy_coding <- run_test(native_united[[contrast_name]], covariates_legacy[[contrast_name]],
                                  sprintf("[%s] %s, methylKit's route, covariate coded as in the legacy script", mode_name, contrast_name))
    coded <- as.data.table(getData(fit_legacy_coding))
    must(identical(cpg_key(coded$chr, coded$start, extract_levels), key_a), "The test with the legacy coding did not return the CpGs of its united table.")
    label_coding <- sprintf("%s contrast, %s: ", contrast_name, coding_pair)
    record_values(part_coding, mode_name, paste0(label_coding, "methylation difference (meth.diff)"), "CpG", coded$meth.diff[row_a], b$meth.diff[row_b],
                  rule = "tolerance", note = "the difference does not depend on the covariate")
    record_values(part_coding, mode_name, paste0(label_coding, "p-value"), "CpG", coded$pvalue[row_a], b$pvalue[row_b], rule = "recorded",
                  note = "one model, with the covariate as an indicator (legacy script) or as the integer codes of a factor (step 05)")
    record_values(part_coding, mode_name, paste0(label_coding, "q-value (SLIM)"), "CpG", coded$qvalue[row_a], b$qvalue[row_b], rule = "recorded")
    called_legacy_coding <- lapply(thresholds, function(threshold) {
      called <- quiet_selection(getData(getMethylDiff(fit_legacy_coding, difference = threshold, qvalue = config$dml_qvalue)))
      cpg_key(called$chr, called$start, extract_levels)
    })
    names(called_legacy_coding) <- thresholds
    record(part_coding, mode_name,
           sprintf("%sDML at q < %g, at each difference threshold of the sweep (%s points)", label_coding, config$dml_qvalue,
                   paste(range(config$dml_difference_sweep), collapse = " to ")),
           "DML", sum(vapply(names(by_threshold), function(x) length(union(called_legacy_coding[[x]], by_threshold[[x]]$b)), numeric(1))),
           sum(vapply(names(by_threshold), function(x) {
             length(setdiff(called_legacy_coding[[x]], by_threshold[[x]]$b)) + length(setdiff(by_threshold[[x]]$b, called_legacy_coding[[x]]))
           }, numeric(1))),
           note = paste("DML by threshold with the legacy coding:",
                        paste(sprintf("%s: %d", names(called_legacy_coding), lengths(called_legacy_coding)), collapse = "; ")))
    dml_legacy_coding[[contrast_name]] <- called_legacy_coding[[as.character(config$dml_difference)]]
    rm(fit_a, fit_b, fit_legacy_coding, a, b, coded)
  }
  record_sets(part_extract, mode_name, "loci in both DML lists", "locus",
              intersect(dml$ploidy$a, dml$heat$a), intersect(dml$ploidy$b, dml$heat$b))
  record_sets(part_coding, mode_name, paste0(coding_pair, ": loci in both DML lists"), "locus",
              intersect(dml_legacy_coding$ploidy, dml_legacy_coding$heat), intersect(dml$ploidy$b, dml$heat$b))
  warned <- table(selections$warnings)
  log_msg("[%s] DML selections: %d made, %d of them without a CpG; warnings: %s", mode_name, selections$made, selections$empty,
          if (length(warned) == 0) "none" else paste(sprintf("%d x '%s'", as.integer(warned), names(warned)), collapse = "; "))
  selections <- list(made = 0L, empty = 0L, warnings = character(0))
  extract_facts[[mode_name]] <- facts
  rm(native, native_united, helper, facts, scaled, united_now, dml, dml_legacy_coding); gc(verbose = FALSE)
}

if ("corrected" %in% config$modes) {
  # The extract table gives, region by region, the lines that the upper filter
  # removes; part 1 counted them with upper_coverage_cutoff. Summed over the
  # regions they must be the rows that filterByCoverage removed in each library.
  record_values(part_extract, "corrected",
                "lines that the upper filter removes, as counted by region for the extract table, against the rows that filterByCoverage removed",
                "library", rowSums(removed_by_group), extract_facts$corrected$rows_removed_by_upper_methylkit)
}

# ---- 4. The legacy call in methylKit's own objects -----------------------------------------
# The coding of the covariate was compared above, next to the tests. Here: the
# two legacy defects.
# I-38: with header = FALSE methRead returns the rows it returns with its default
# header = TRUE, and the first line of the file where that line has at least
# import_min_coverage reads. Nothing else differs.
if (all(c("legacy", "corrected") %in% config$modes)) {
  first_key <- cpg_key(first_line$chr, first_line$start, extract_levels)
  as_expected <- vapply(seq_len(n_lib), function(lib_i) {
    only_complete <- setdiff(imported_keys$corrected[[lib_i]], imported_keys$legacy[[lib_i]])
    only_default <- setdiff(imported_keys$legacy[[lib_i]], imported_keys$corrected[[lib_i]])
    length(only_default) == 0 &&
      identical(only_complete, if (first_line$reads[lib_i] >= config$import_min_coverage) first_key[lib_i] else numeric(0))
  }, logical(1))
  record(part_extract, "legacy",
         "rows that methRead returns with header = FALSE and not with its default header = TRUE: the first line of the file and no other",
         "library", n_lib, sum(!as_expected),
         note = sprintf("the first line has at least %d reads, and so is imported, in %d libraries; it has at least %d reads in %d",
                        config$import_min_coverage, sum(first_line$reads >= config$import_min_coverage),
                        config$min_coverage, sum(first_line$reads >= config$min_coverage)))
  record(part_extract, "legacy",
         "rows kept by filterByCoverage as the legacy scripts called it (high.count, high.perc) against the call with no upper filter",
         "library", n_lib, sum(!misspelled), note = "the misspelled arguments are ignored")
}

# ---- 5. The tables of step 05 ---------------------------------------------------------------
# The whole tables were set against the whole coverage files during the pass of
# part 1 (tally `whole05`). Here: the rows of the tables inside the regions,
# against the lines of the extract files.
for (mode_name in config$modes) {
  united05 <- tables05[[mode_name]]                                    # the rows inside the regions
  recorded <- records05[[mode_name]]
  key05 <- cpg_key(united05$chr, united05$start, extract_levels)
  record_sets(part_tables, mode_name,
              sprintf("CpGs of united_cpg_%s.csv.gz inside the regions against the united CpGs of the extract", mode_name), "CpG",
              key05, united_keys[[mode_name]], rule = if (mode_name == "legacy") "identical" else "recorded",
              note = sprintf("%s CpGs of step 05 inside the regions, %s of them not among the %s united CpGs of the extract%s",
                             count_text(length(key05)), count_text(length(setdiff(key05, united_keys[[mode_name]]))),
                             count_text(length(united_keys[[mode_name]])),
                             if (mode_name == "legacy") "" else "; the cutoffs of the extract are not those of the full files, so the sets may differ"))
  # The CpGs of the regions that pass the rule of the mode in every library. The
  # extract holds every line of the regions, and its first line is the first
  # line of the coverage file. The upper cutoff is the one that step 05 recorded
  # for the whole file, not the cutoff of the extract file.
  inside <- c(values = 0, different = 0, half = 0); passing <- NULL
  for (lib_i in seq_len(n_lib)) {
    lines <- read_cov(extract_file[lib_i])
    key_lines <- cpg_key(lines$chr, lines$start, extract_levels)
    reads <- lines$nC + lines$nT
    passes <- reads >= config$min_coverage
    if (mode_name == "legacy") passes[1L] <- FALSE                     # the first line of the file was not read
    if (mode_name == "corrected") passes <- passes & reads < recorded$upper_cutoff[lib_i]
    passing <- if (is.null(passing)) key_lines[passes] else passing[passing %in% key_lines[passes]]
    at <- match(key05, key_lines)
    must(!anyNA(at), "A CpG that step 05 tested in mode ", mode_name, " has no line in the extract of ", lib$library_id[lib_i], ".")
    inside <- inside + table_against_lines(united05, lib_i, lines$nC[at], lines$nT[at], factor05[[mode_name]][lib_i])
  }
  record_sets(part_tables, mode_name,
              sprintf("CpGs of united_cpg_%s.csv.gz inside the regions against the CpGs of the regions that pass the rule of the mode in every library",
                      mode_name),
              "CpG", key05, passing,
              note = sprintf("the rule: at least %g reads%s; %s CpGs of the regions pass it in every library", config$min_coverage,
                             if (mode_name == "legacy") ", and not the first line of a coverage file"
                             else sprintf(" and fewer reads than upper_cutoff of library_filtering_%s.csv, the cutoff of the whole file", mode_name),
                             count_text(length(passing))))
  factor_note <- sprintf("factor: largest median over own median, from the medians in library_filtering_%s.csv", mode_name)
  record(part_tables, mode_name,
         sprintf("counts of united_cpg_%s.csv.gz at the CpGs inside the regions against the lines of the extract files scaled by scale_counts_like_methylkit",
                 mode_name),
         "count", inside[["values"]], inside[["different"]], note = factor_note)
  record(part_tables, mode_name,
         sprintf("counts of united_cpg_%s.csv.gz at every CpG of the table against the lines of the whole coverage files scaled by scale_counts_like_methylkit",
                 mode_name),
         "count", whole05[[mode_name]][["values"]], whole05[[mode_name]][["different"]],
         note = sprintf("%s CpGs in %d libraries; %s; %s values are exactly halfway between two whole numbers before rounding",
                        count_text(n_tested05[[mode_name]]), n_lib, factor_note, count_text(whole05[[mode_name]][["half"]])))
  rm(united05, recorded, key05, lines, key_lines, reads, passes, passing, at); gc(verbose = FALSE)
}

# ---- 6. Rounding: synthetic counts ----------------------------------------------------------
set.seed(seed)
synthetic <- list(values = 0, different = 0, other_order = 0, half_coverage = 0, half_counts = 0, near = 0, factors = numeric(0))
for (set_i in seq_len(synthetic_sets)) {
  top <- synthetic_top[sample.int(length(synthetic_top), 1L)]
  medians <- c(top, top %/% 3L * 2L, top %/% 5L * 4L, sample(seq.int(top %/% 3L, top), 1L))   # factors 1, 1.5, 1.25 and one arbitrary ratio
  counts <- lapply(medians, function(target) {
    coverage <- rnbinom(synthetic_cpgs, size = 4, mu = target)
    coverage <- pmax(1L, coverage - as.integer(stats::median(coverage)) + target)             # median coverage is now `target`
    methylated <- rbinom(synthetic_cpgs, coverage, rbeta(synthetic_cpgs, 0.4, 0.4))
    data.frame(chr = "synthetic", start = 2L * seq_len(synthetic_cpgs), end = 2L * seq_len(synthetic_cpgs) + 2L, strand = "*",
               coverage = coverage, numCs = methylated, numTs = coverage - methylated, stringsAsFactors = FALSE)
  })
  must(identical(vapply(counts, function(x) stats::median(x$coverage), numeric(1)), as.numeric(medians)),
       "A synthetic library does not have the median coverage it was built for.")
  objects <- lapply(seq_along(counts), function(k) {
    new("methylRaw", counts[[k]], sample.id = paste0("synthetic", k), assembly = "synthetic", context = "CpG", resolution = "base")
  })
  by_methylkit <- normalizeCoverage(new("methylRawList", objects, treatment = c(0, 1, 0, 1)))
  factors <- max(medians) / medians                                                          # as step 05 computes them
  for (k in seq_along(counts)) {
    x <- counts[[k]]; a <- getData(by_methylkit[[k]])
    b <- scale_counts_like_methylkit(x$coverage, x$numCs, x$numTs, factors[k])
    other_c <- round(b$coverage * x$numCs / x$coverage); other_t <- round(b$coverage * x$numTs / x$coverage)   # the order of step 04
    unrounded_counts <- c(b$coverage * (x$numCs / x$coverage), b$coverage * (x$numTs / x$coverage))
    synthetic$values <- synthetic$values + 3 * nrow(x)
    synthetic$different <- synthetic$different + sum(a$coverage != b$coverage) + sum(a$numCs != b$numCs) + sum(a$numTs != b$numTs)
    synthetic$other_order <- synthetic$other_order + sum(a$numCs != other_c) + sum(a$numTs != other_t)
    synthetic$half_coverage <- synthetic$half_coverage + sum(is_half(factors[k] * x$coverage))
    synthetic$half_counts <- synthetic$half_counts + sum(is_half(unrounded_counts))
    synthetic$near <- synthetic$near + sum(is_near_half(c(factors[k] * x$coverage, unrounded_counts)))
  }
  synthetic$factors <- c(synthetic$factors, factors)
}
must(synthetic$half_coverage > 0 && synthetic$half_counts > 0 && synthetic$other_order > 0,
     "The synthetic counts hold no value halfway between two whole numbers, or do not tell the two orders of operations apart.")
record(part_synthetic, NA_character_, "normalized coverage, numCs and numTs: normalizeCoverage against scale_counts_like_methylkit", "count",
       synthetic$values, synthetic$different,
       note = sprintf(paste("exactly halfway between two whole numbers before rounding: %s coverage values and %s counts;",
                            "within %g of a half and not a half: %s values; %d distinct scaling factors"),
                      count_text(synthetic$half_coverage), count_text(synthetic$half_counts), near_half, count_text(synthetic$near),
                      length(unique(synthetic$factors))))
record(part_synthetic, NA_character_, "numCs and numTs: normalizeCoverage against the order of operations of step 04, round(new_cov * nc / cv)",
       "count", synthetic$values / 3 * 2, synthetic$other_order, rule = "recorded",
       note = "a different order of operations is detected by these counts")

# ---- 7. The order of operations of step 04 ---------------------------------------------------
# Rows of every library that step 04 counts as covered at `level` reads (its
# function covered()): at least `level` reads; in corrected mode fewer reads
# than the library's cutoff; in legacy mode not the first line of the file.
rows_of_step04 <- function(mode_name, level) {
  lapply(seq_len(n_lib), function(lib_i) {
    lines <- read_cov(extract_file[lib_i])
    lines[, coverage := nC + nT]
    cutoff <- upper_coverage_cutoff(lines$coverage, config$import_min_coverage, config$upper_coverage_percentile)
    counted <- lines$coverage >= level & (mode_name == "legacy" | lines$coverage < cutoff)
    if (first_line_lost(mode_name)) counted[1L] <- FALSE
    data.table(pos_key = cpg_key(lines$chr[counted], lines$start[counted], extract_levels),
               coverage = lines$coverage[counted], nC = lines$nC[counted])
  })
}
for (mode_name in config$modes) for (level in levels_step04) {
  counted <- rows_of_step04(mode_name, level)
  in_all <- Reduce(function(a, b) a[a %in% b], lapply(counted, function(x) x$pos_key))
  factors <- vapply(counted, function(x) stats::median(x$coverage), numeric(1))
  factors <- max(factors) / factors
  if (level == config$min_coverage) {
    # At this level step 04 and step 05 select by one rule, so the selection made
    # here must give the united CpGs and the scaling factors of the helper route.
    record_sets(part_order(level), mode_name, "CpGs selected by the rule of step 04 against the united CpGs of the helper route", "CpG",
                in_all, united_keys[[mode_name]])
    record_values(part_order(level), mode_name, "scaling factors by the rule of step 04 against the factors of the helper route", "library",
                  factors, extract_facts[[mode_name]]$scaling_factor_helper, rule = "tolerance")
  }
  pm_methylkit <- pm_step04 <- matrix(NA_real_, length(in_all), n_lib)
  differing <- c(numCs = 0, numTs = 0)
  for (lib_i in seq_len(n_lib)) {
    x <- counted[[lib_i]][match(in_all, pos_key)]
    cv <- x$coverage; nc <- x$nC
    by_methylkit <- scale_counts_like_methylkit(cv, nc, cv - nc, factors[lib_i])               # methylKit's order (parts 3 and 6)
    new_cov <- by_methylkit$coverage
    new_c <- round(new_cov * nc / cv); new_t <- round(new_cov * (cv - nc) / cv)                # as written in step 04
    differing <- differing + c(sum(new_c != by_methylkit$numCs), sum(new_t != by_methylkit$numTs))
    pm_methylkit[, lib_i] <- 100 * by_methylkit$numCs / (by_methylkit$numCs + by_methylkit$numTs)
    pm_step04[, lib_i] <- 100 * new_c / (new_c + new_t)
  }
  complete <- stats::complete.cases(pm_methylkit) & stats::complete.cases(pm_step04)
  cells_differing <- sum(xor(is.na(pm_methylkit), is.na(pm_step04))) + sum(pm_methylkit != pm_step04, na.rm = TRUE)
  r_methylkit <- stats::cor(pm_methylkit[complete, , drop = FALSE]); r_step04 <- stats::cor(pm_step04[complete, , drop = FALSE])
  pair <- upper.tri(r_methylkit)
  n_cells <- length(in_all) * n_lib
  record(part_order(level), mode_name, "numCs after scaling: round(new_cov * nc / cv) of step 04 against round(new_cov * (nc / cv)) of methylKit",
         "count", n_cells, differing[["numCs"]], rule = "recorded")
  record(part_order(level), mode_name, "numTs after scaling: the same two orders", "count", n_cells, differing[["numTs"]], rule = "recorded")
  record(part_order(level), mode_name, "percent methylation of a CpG in a library, from the scaled counts", "CpG in a library",
         n_cells, cells_differing, max(abs(pm_methylkit - pm_step04), na.rm = TRUE), rule = "recorded",
         note = sprintf("largest absolute difference in percentage points; %s CpGs", count_text(length(in_all))))
  record(part_order(level), mode_name, "Pearson correlation between two libraries", "pair of libraries", sum(pair),
         sum(r_methylkit[pair] != r_step04[pair]), max(abs(r_methylkit[pair] - r_step04[pair])), rule = "recorded",
         note = sprintf("correlations range from %.4f to %.4f", min(r_methylkit[pair]), max(r_methylkit[pair])))
}
unlink(extract_dir, recursive = TRUE)

# ---- 8. Write -------------------------------------------------------------------------------
# Rows in the order of the argument: the extract, the tables of step 05, the
# full libraries, the coding of the covariate in the legacy script, the
# synthetic counts, the order of operations of step 04.
comparisons <- rbindlist(comparisons)
part_levels <- c(part_extract, part_tables, part_full, part_coding, part_synthetic, vapply(levels_step04, part_order, character(1)))
must(all(comparisons$part %in% part_levels), "A comparison belongs to no part.")
comparisons <- comparisons[order(match(part, part_levels))]          # a stable sort: rows of one part stay in the order they were made
comparisons[, methylkit_version := as.character(utils::packageVersion("methylKit"))]
comparisons[rule == "tolerance", rule := sprintf("tolerance: relative difference at most %g", float_tolerance)]
fwrite(comparisons, out_file())

# Extract first, then the full libraries; inside each, the modes in the order of
# run_config.R and the libraries in the order of the sample sheet.
library_table <- rbind(rbindlist(extract_facts), rbindlist(full_facts), fill = TRUE)
sheet_row <- match(library_table$library_id, lib$library_id)
library_table[, `:=`(ploidy = as.character(lib$ploidy[sheet_row]), treatment = as.character(lib$treatment[sheet_row]))]
setcolorder(library_table, c(
  "scope", "mode", "library_id", "ploidy", "treatment", "lines_read", "first_line_reads",
  "rows_imported_methylkit", "rows_imported_helper", "first_line_imported_methylkit",
  "rows_min_coverage_methylkit", "rows_min_coverage_helper", "upper_cutoff_from_methylkit_rows", "upper_cutoff_helper", "rows_with_cutoff_reads",
  "rows_removed_by_upper_methylkit", "rows_removed_by_upper_helper", "rows_kept_methylkit", "rows_kept_helper", "kept_rows_differing",
  "first_line_kept_methylkit", "first_line_kept_helper", "median_coverage_methylkit", "median_coverage_helper",
  "max_coverage_kept_methylkit", "max_coverage_kept_helper", "scaling_factor_from_methylkit_medians", "scaling_factor_helper"))
fwrite(library_table, out_file("_libraries"))

# One row per region; the positions of first lines outside the other regions form one row.
by_group <- function(x) as.numeric(tapply(x, region_group, sum))
first_region <- region_of(first_line$chr, first_line$start)
extract_table <- data.table(
  region = c(regions$region[seq_len(n_interval_regions)], "positions of first lines outside the regions above"),
  sequence = c(regions$sequence[seq_len(n_interval_regions)], NA_character_),
  from = as.integer(c(regions$from[seq_len(n_interval_regions)], NA)),
  to = as.integer(c(ifelse(is.finite(regions$to[seq_len(n_interval_regions)]), regions$to[seq_len(n_interval_regions)], NA), NA)),   # empty: to the end
  purpose = c(regions$purpose[seq_len(n_interval_regions)], "the first-line rule of legacy mode"),
  positions = c(rep(NA_integer_, n_interval_regions), sum(region_of(regions$sequence, regions$from) > n_interval_regions)),
  libraries_with_lines = colSums(lines_by_group > 0),
  lines_per_library_min = apply(lines_by_group, 2, min), lines_per_library_max = apply(lines_by_group, 2, max),
  files_that_begin_here = as.integer(tabulate(region_group[first_region], nbins = n_groups)),
  median_reads_over_library_median_min = round(apply(relative_median, 2, min, na.rm = TRUE), 2),
  median_reads_over_library_median_max = round(apply(relative_median, 2, max, na.rm = TRUE), 2),
  lines_removed_by_upper_filter = colSums(removed_by_group))
for (mode_name in config$modes) {
  set(extract_table, j = paste0("cpg_united_on_extract_", mode_name), value = by_group(tested_by_region[[mode_name]]))
  set(extract_table, j = paste0("cpg_tested_by_step05_", mode_name), value = by_group(tested05_by_region[[mode_name]]))
}
total <- extract_table[, lapply(.SD, sum), .SDcols = patterns("^cpg_|^files_that|^lines_removed")]
total[, `:=`(region = "all regions", purpose = "the extract", libraries_with_lines = n_lib,
             lines_per_library_min = min(file_facts$lines_in_extract), lines_per_library_max = max(file_facts$lines_in_extract))]
extract_table <- rbind(extract_table, total, fill = TRUE)
for (mode_name in config$modes) {
  set(extract_table, j = paste0("percent_of_all_cpg_tested_by_step05_", mode_name),
      value = round(100 * extract_table[[paste0("cpg_tested_by_step05_", mode_name)]] / n_tested05[[mode_name]], 3))
}
fwrite(extract_table, out_file("_extract"))

# ---- Summary for the log -----------------------------------------------------------------------
options(width = 200)
cat("\nComparisons by part, mode and result\n")
print(comparisons[, .(comparisons = .N, items_compared = sum(n_compared), items_different = sum(n_different)), by = .(part, mode, result)])
for (mode_name in config$modes) {
  row <- extract_table[region == "all regions"]
  cat(sprintf("\n[%s] united CpGs on the extract: %s. CpGs tested by step 05 inside the regions: %s of %s (%.2f%%).\n", mode_name,
              count_text(row[[paste0("cpg_united_on_extract_", mode_name)]]), count_text(row[[paste0("cpg_tested_by_step05_", mode_name)]]),
              count_text(n_tested05[[mode_name]]), row[[paste0("percent_of_all_cpg_tested_by_step05_", mode_name)]]))
}
log_msg("finished in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
failed <- comparisons[result == "DIFFERENT"]
if (nrow(failed) > 0) {
  print(failed[, .(part, mode, quantity, n_compared, n_different, largest_relative_difference)])
  stop(nrow(failed), " comparison(s) that must agree do not; see ", out_file(), call. = FALSE)
}
cat(sprintf(paste("\nAll %d comparisons that must agree do (%d identical, %d within tolerance).",
                  "%d differences are recorded and not judged. %d comparisons had nothing to compare.\n"),
            sum(comparisons$result %in% c("identical", "within tolerance")), sum(comparisons$result == "identical"),
            sum(comparisons$result == "within tolerance"), sum(comparisons$result == "recorded"), sum(comparisons$result == "nothing to compare")))

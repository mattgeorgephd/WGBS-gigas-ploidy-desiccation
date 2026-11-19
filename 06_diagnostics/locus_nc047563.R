# locus_nc047563.R
#
# Purpose   What the locus of issue I-09 is. One region of 2.2 kb on NC_047563.1
#           holds the largest run of DML of each legacy list (ploidy, heat, loci
#           in both lists) and all of the draft's "lncRNA" signal, and every
#           tested CpG in it is removed by the upper coverage filter that the
#           Methods describe (issue I-05). This script documents the region from
#           the raw coverage files and the step 05 tables:
#             a. the block of excess coverage around the region: where it begins
#                and ends, how many CpG sites and tested CpGs it holds, and the
#                coverage of the sequence on either side;
#             b. which libraries carry the excess, how they are spread over
#                ploidy and cohort (Fisher exact tests), and whether the size of
#                the excess differs by ploidy or cohort (Wilcoxon rank-sum tests);
#             c. how many CpGs of the region still differ by more than
#                config$dml_difference points when one library at a time is left
#                out of the read-pooled group means;
#             d. percent methylation of the region in each library, compared
#                between ploidies (Wilcoxon rank-sum tests);
#             e. the genes, exons and transposable elements annotated at the block;
#             f. base composition of the block, how much of its sequence occurs
#                once in the reference, where the other copies lie, and whether
#                coverage reaches the cutoff there as well.
#
# Run       Rscript 06_diagnostics/locus_nc047563.R      (about 2 minutes; 1.5 GB of memory; one core)
#           Needs the 22 coverage files and the reference genome in the cache (step 01).
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files and the reference genome listed in
#             03_input_files/gannet_manifest.csv (cache)
#           03_input_files/genome_features/   gene, exon and transposable element tracks
#           04_output/04_cpg_coverage/tables/library_coverage.csv
#           04_output/05_differential_methylation/tables/
#             united_cpg_legacy.csv.gz, united_cpg_corrected.csv.gz,
#             diffmeth_{ploidy,heat}_legacy.csv.gz, dml_{ploidy,heat,shared}_legacy.csv,
#             library_filtering_legacy.csv
#           04_output/07_genome_features/tables/dml_regions_legacy.csv
# Writes    06_diagnostics/locus_nc047563_summary.csv            one row per fact: section, quantity, value, how defined
#           06_diagnostics/locus_nc047563_per_library.csv        one row per library
#           06_diagnostics/locus_nc047563_leave_one_out.csv      one row per set of libraries left out
#           06_diagnostics/locus_nc047563_coverage_profile.csv   one row per CpG site of the window
#
# Notes     - Positions are the 0-based start of the C, as in the coverage files
#             and in every pipeline table. A CpG occupies bases start + 1 and
#             start + 2 of the sequence in 1-based coordinates.
#           - "DML region" is the run of DML that step 07 reports. "Tested CpGs"
#             are those of legacy mode (at least config$min_coverage reads in all
#             22 libraries). The locus yields DML in legacy mode only, so every
#             DML count here is a legacy count.
#           - Coverage is the raw read count of a CpG in a coverage file.
#             Relative coverage is coverage divided by the median coverage of the
#             library over its CpGs with at least one read (library_coverage.csv).
#             A coverage file lists only CpGs with reads. The median coverage of
#             a library over a set of CpG sites is therefore taken over the sites
#             with reads in that library, as the library median is. In the per-CpG
#             profile a library without reads at a site counts as coverage 0.
#           - Rule for a CpG with excess coverage: at least one library is at or
#             above its upper cutoff there (the percentile of config, the rule by
#             which corrected mode removes a CpG). Rule for the block: CpGs with
#             excess coverage that lie within run_gap_bp of each other form a
#             run, and the block is the run that holds the DML region. The block
#             is where coverage reaches the cutoff, to the nearest CpG. The
#             sequence present in extra copies may reach further: the summary
#             gives the coverage of the CpG site next to each end, and the
#             profile every site.
#           - Rule for a library that carries the excess: at least one tested CpG
#             of the DML region is at or above the library's upper cutoff. This
#             is a threshold on a continuous quantity. Column block_over_flanks
#             of the library table gives the excess of every library, those
#             below the threshold included.
#           - Leave-one-out applies the difference rule only: the model is not
#             refitted and q-values are not recomputed. A CpG must pass the
#             difference rule to be a DML, so each count is the largest number of
#             DML that the region could give without that library. The differences
#             are read-pooled, as methylKit computes meth.diff: methylated reads
#             of all libraries of a group over all their reads, from the
#             normalized counts that step 05 tested. Ploidy differences are
#             triploid minus diploid; heat differences are heat minus control.
#           - Percent methylation of a library is read-pooled over the tested
#             CpGs of the DML region, from raw counts. Each library is taken as
#             one animal; the sample sheet names the animal for the heat cohort
#             only (issue I-36). The Wilcoxon tests are exact; the normal
#             approximation is written next to each, because the two differ at
#             these group sizes (5, 6 or 11 libraries).
#           - Cohort and heat treatment are the same split of the libraries
#             (issue I-01). With 22 libraries the Fisher tests can detect only a
#             strong association; a large p-value does not show independence.
#             The tests are therefore repeated under a stricter rule, and the
#             size of the excess is compared between groups as a continuous
#             quantity.
#           - Part f counts every 30-mer of the block on both strands of every
#             sequence of the reference. For the 30-mers with exactly two copies
#             it reports the places that hold the second copies, and whether any
#             library reaches its upper cutoff there.
#           - No random numbers are used. The 30-mer sample of part f is one
#             30-mer every kmer_sample_step bases, from the first base of the block.
#           - Columns: "region" is the tested CpGs of the DML region, "block" and
#             "flank" are all CpG sites with reads. ploidy_difference and
#             heat_difference of the profile are step 05's meth.diff (legacy).
#             Shares of reads are shares of the normalized reads of a group.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges", "Biostrings"))
setDTthreads(1)                     # run time and memory in the header are for one thread; results do not depend on it

# ---- Settings --------------------------------------------------------------------------
locus_sequence   <- "NC_047563.1"   # the sequence this script is about; the DML region must lie on it
window_margin_bp <- 10000L          # raw coverage is read for the DML region and this much sequence on either side
flank_bp         <- 5000L           # width of each flank of the block
run_gap_bp       <- 1000L           # CpGs with excess coverage closer than this form one run (the gap of step 07's DML regions)
kmer_width       <- 30L             # part f: length of the words (k-mers) of the block that are counted in the reference
kmer_sample_step <- 60L             # part f: the sample takes one word every 60 bases
copy_gap_bp      <- 20000L          # part f: second copies closer than this on one reference sequence are one place
n_places         <- 3L              # part f: number of such places reported

out_file <- function(suffix) here("06_diagnostics", paste0("locus_nc047563_", suffix, ".csv"))
tab05 <- function(f) output_path("05_differential_methylation", "tables", f)

# One row of the summary table. A value is written as text: a number, a label,
# "x of y", or "a to b" for the lowest and highest of several values.
facts <- list()
fact <- function(section, quantity, value, how_defined) {
  stopifnot(length(value) == 1L, !is.na(value))
  facts[[length(facts) + 1L]] <<- data.table(section = section, quantity = quantity, value = as.character(value),
                                              how_defined = how_defined)
  invisible(NULL)
}
whole    <- function(x) formatC(x, format = "d")
decimal  <- function(x, digits) formatC(round(x, digits), format = "f", digits = digits)
# A p-value or q-value as text, with four significant digits. A written value
# that ends in the digit 5 lies on a rounding boundary, and a reader who rounds
# it further can go the wrong way (0.011499 written as 0.0115 becomes 0.012).
# Digits are therefore added until the last digit is not 5; a value that is
# exactly such a number (0.5, 0.25) stays as it is.
p_value  <- function(p) {
  for (digits in 4:15) {
    text <- sprintf("%.*g", digits, p)
    if (!endsWith(sub("e.*$", "", text), "5")) break     # the last digit before any exponent
  }
  text
}
x_of_y   <- function(x, y) paste(whole(x), "of", whole(y))
a_to_b   <- function(x, digits) paste(decimal(min(x), digits), "to", decimal(max(x), digits))
interval <- function(from, to, sequence = locus_sequence) sprintf("%s:%d-%d", sequence, from, to)
# The summary writes positions in two conventions, and each definition says which.
zero_based <- "Both positions are the 0-based start of the C, as in the coverage files."

# ---- 1. Libraries, the DML region, the window ---------------------------------------------
lib <- read_libraries(); n_lib <- nrow(lib)
is_triploid <- lib$ploidy == "triploid"; is_heat <- lib$treatment == "heat"
library_coverage <- fread(output_path("04_cpg_coverage", "tables", "library_coverage.csv"))
filtering <- fread(tab05("library_filtering_legacy.csv"))
stopifnot(identical(library_coverage$library_id, lib$library_id), identical(filtering$library_id, lib$library_id))
library_median <- library_coverage$median_coverage; upper_cutoff <- library_coverage$upper_cutoff

# The DML region: the run of ploidy DML with the most DML (step 07).
dml_runs <- fread(output_path("07_genome_features", "tables", "dml_regions_legacy.csv"))
largest_runs <- dml_runs[, .SD[n_dml == max(n_dml)], by = set]   # the largest run of each DML list
largest_run <- largest_runs[set == "ploidy DML"]
if (nrow(largest_run) != 1L || largest_run$chr != locus_sequence)
  stop("The largest run of legacy ploidy DML is not a single region on ", locus_sequence, "; this script describes that region only.")
same_run <- largest_runs[chr == largest_run$chr & start == largest_run$start & end == largest_run$end, set]
region_from <- largest_run$start; region_to <- largest_run$end - 2L      # start of the first and of the last DML CpG
window_from <- region_from - window_margin_bp; window_to <- region_to + window_margin_bp

dml_list <- lapply(c(ploidy = "ploidy", heat = "heat", shared = "shared"), function(s) fread(tab05(sprintf("dml_%s_legacy.csv", s))))
dml_total <- vapply(dml_list, nrow, integer(1))
dml_window <- lapply(dml_list, function(d) d[chr == locus_sequence & start >= window_from & start <= window_to, start])
dml_region <- lapply(dml_window, function(s) s[s >= region_from & s <= region_to])
stopifnot(length(dml_region$ploidy) == largest_run$n_dml, range(dml_region$ploidy) == c(region_from, region_to))

# ---- 2. Raw coverage of the window, one library at a time -----------------------------------
log_msg("reading %d coverage files for %s", n_lib, interval(window_from, window_to))
manifest <- read_gannet_manifest()
window_rows <- vector("list", n_lib)
cutoff_rows <- vector("list", n_lib)   # every CpG of the genome at or above the library's cutoff; part f looks some of them up
for (lib_i in seq_len(n_lib)) {
  cov_table <- read_cov(fetch_gannet_file(lib$cov_file[lib_i], config, manifest))
  reads <- cov_table$nC + cov_table$nT
  cutoff_here <- upper_coverage_cutoff(reads, config$import_min_coverage, config$upper_coverage_percentile)
  high <- which(reads >= upper_cutoff[lib_i])
  if (median(reads) != library_median[lib_i] || !isTRUE(all.equal(cutoff_here, upper_cutoff[lib_i])) ||
      length(high) != library_coverage$n_cpg_at_or_above_cutoff[lib_i])
    stop("The coverage file of ", lib$library_id[lib_i], " does not give the median, the cutoff and the number of CpGs at the cutoff of ",
         "library_coverage.csv; rerun step 04.")
  keep <- which(cov_table$chr == locus_sequence & cov_table$start >= window_from & cov_table$start <= window_to)
  window_rows[[lib_i]] <- data.table(lib_index = lib_i, start = cov_table$start[keep], nC = cov_table$nC[keep], nT = cov_table$nT[keep])
  cutoff_rows[[lib_i]] <- data.table(lib_index = lib_i, chr = cov_table$chr[high], start = cov_table$start[high])
  rm(cov_table, reads, keep, high); gc(verbose = FALSE)
}
window_rows <- rbindlist(window_rows); cutoff_rows <- rbindlist(cutoff_rows)
stopifnot(nrow(window_rows) > 0L, !anyDuplicated(window_rows, by = c("lib_index", "start")))

# Matrices of raw counts: one row per CpG site with reads in at least one
# library, one column per library (sample sheet order); 0 where a library has no reads.
site_start <- sort(unique(window_rows$start)); n_site <- length(site_start)
cell <- cbind(match(window_rows$start, site_start), window_rows$lib_index)
reads_raw <- matrix(0L, n_site, n_lib, dimnames = list(NULL, lib$library_id)); meth_raw <- reads_raw
reads_raw[cell] <- window_rows$nC + window_rows$nT; meth_raw[cell] <- window_rows$nC
rm(window_rows, cell)
relative <- sweep(reads_raw, 2, library_median, "/")           # coverage over the library median
at_cutoff <- sweep(reads_raw, 2, upper_cutoff, ">=")            # at or above the library's upper cutoff

# ---- 3. Tested CpGs of the window, and checks that tie the tables together -------------------
united <- fread(tab05("united_cpg_legacy.csv.gz"))
n_tested_legacy <- nrow(united)
stopifnot(sum(grepl("^numCs", names(united))) == n_lib)
united_window <- united[chr == locus_sequence & start >= window_from & start <= window_to]
rm(united); invisible(gc(verbose = FALSE))
corrected_window <- fread(tab05("united_cpg_corrected.csv.gz"), select = c("chr", "start"))[
  chr == locus_sequence & start >= window_from & start <= window_to, start]

stopifnot(all(united_window$start %in% site_start), all(corrected_window %in% united_window$start))
tested_legacy <- site_start %in% united_window$start
tested_corrected <- site_start %in% corrected_window
in_region <- site_start >= region_from & site_start <= region_to
region_rows <- which(tested_legacy & in_region)                 # the tested CpGs of the DML region
n_region <- length(region_rows)
stopifnot(all(reads_raw[tested_legacy, ] >= config$min_coverage), all(dml_window$ploidy %in% site_start[tested_legacy]),
          all(dml_window$heat %in% site_start[tested_legacy]), all(dml_window$shared %in% intersect(dml_window$ploidy, dml_window$heat)))

# Column j of the united table must be library j of the sample sheet: the
# normalized counts of every tested CpG of the window are rebuilt from the raw file.
scaling_factor <- max(filtering$median_coverage) / filtering$median_coverage
stopifnot(isTRUE(all.equal(scaling_factor, filtering$scaling_factor)))
united_window <- united_window[match(site_start[tested_legacy], start)]
for (lib_i in seq_len(n_lib)) {
  raw_reads <- reads_raw[tested_legacy, lib_i]; raw_meth <- meth_raw[tested_legacy, lib_i]
  rebuilt <- scale_counts_like_methylkit(raw_reads, raw_meth, raw_reads - raw_meth, scaling_factor[lib_i])
  if (!all(rebuilt$numCs == united_window[[paste0("numCs", lib_i)]]) || !all(rebuilt$numTs == united_window[[paste0("numTs", lib_i)]]))
    stop("The normalized counts of ", lib$library_id[lib_i], " in united_cpg_legacy.csv.gz are not those of its coverage file.")
}
united_region <- united_window[start >= region_from & start <= region_to]
stopifnot(identical(united_region$start, site_start[region_rows]))
norm_meth <- as.matrix(united_region[, paste0("numCs", seq_len(n_lib)), with = FALSE])
norm_reads <- norm_meth + as.matrix(united_region[, paste0("numTs", seq_len(n_lib)), with = FALSE])

# Test results of step 05 at the tested CpGs of the window, in the order of united_window.
diffmeth_window <- lapply(c(ploidy = "ploidy", heat = "heat"), function(s) {
  d <- fread(tab05(sprintf("diffmeth_%s_legacy.csv.gz", s)))[chr == locus_sequence & start >= window_from & start <= window_to]
  stopifnot(setequal(d$start, united_window$start), !anyDuplicated(d$start))
  d[match(united_window$start, start)]
})

sec <- "DML region"
fact(sec, "dml_region", interval(region_from, region_to),
     paste("Run of legacy ploidy DML with the most DML, as step 07 groups them (dml_regions_legacy.csv): its first and its last DML CpG.", zero_based))
fact(sec, "dml_region_length_bp", whole(region_to + 2L - region_from), "From the C of the first DML CpG to the G of the last.")
fact(sec, "lists_whose_largest_run_is_this_region", paste(same_run, collapse = "; "),
     "DML lists of step 07 whose run with the most DML has the same first and last CpG as this region.")
fact(sec, "ploidy_dml_in_region", x_of_y(length(dml_region$ploidy), dml_total[["ploidy"]]),
     "Legacy ploidy DML inside the DML region, of all legacy ploidy DML.")
fact(sec, "heat_dml_in_region", x_of_y(length(dml_region$heat), dml_total[["heat"]]),
     "Legacy heat DML inside the DML region, of all legacy heat DML.")
fact(sec, "loci_in_both_lists_in_region", x_of_y(length(dml_region$shared), dml_total[["shared"]]),
     "Loci in both legacy DML lists inside the DML region, of all such loci.")
fact(sec, "tested_cpgs_in_region_legacy", x_of_y(n_region, n_tested_legacy),
     "Tested CpGs of legacy mode inside the DML region, of all tested CpGs of legacy mode.")
fact(sec, "tested_cpgs_in_region_corrected", whole(sum(tested_corrected & in_region)),
     "Tested CpGs of corrected mode inside the DML region: those that the upper coverage filter leaves.")

# ---- 4. Excess coverage: the libraries that carry it, and the block ------------------------------
# Libraries (rule in the header).
region_median <- apply(reads_raw[region_rows, , drop = FALSE], 2, median)
region_at_cutoff <- colSums(at_cutoff[region_rows, , drop = FALSE])
carries_excess <- region_at_cutoff >= 1L
n_carriers <- sum(carries_excess)
if (n_carriers == 0L || n_carriers == n_lib) stop("All libraries fall on one side of the rule for carrying the excess; this script needs both kinds.")
region_pct_meth <- 100 * colSums(meth_raw[region_rows, , drop = FALSE]) / colSums(reads_raw[region_rows, , drop = FALSE])

# Block (rule in the header).
excess_site <- rowSums(at_cutoff) >= 1L
run_of_site <- rep(NA_integer_, n_site)
run_of_site[excess_site] <- cumsum(c(TRUE, diff(site_start[excess_site]) > run_gap_bp))
block_run <- unique(run_of_site[region_rows])
if (length(block_run) != 1L || is.na(block_run)) stop("The tested CpGs of the DML region do not lie in one run of CpGs with excess coverage.")
block_from <- min(site_start[which(run_of_site == block_run)]); block_to <- max(site_start[which(run_of_site == block_run)])
if (block_from - flank_bp < window_from || block_to + flank_bp > window_to)
  stop("The block and its flanks reach the edge of the window; raise window_margin_bp.")
in_block <- site_start >= block_from & site_start <= block_to
zone <- rep("beyond the flanks", n_site)
zone[site_start >= block_from - flank_bp & site_start < block_from] <- "upstream flank"
zone[site_start > block_to & site_start <= block_to + flank_bp] <- "downstream flank"
zone[in_block] <- "block"
in_flanks <- zone %in% c("upstream flank", "downstream flank")
stopifnot(any(zone == "upstream flank"), any(zone == "downstream flank"))

# Median coverage of each library over the sites of a zone that have reads in it.
median_with_reads <- function(rows) vapply(seq_len(n_lib), function(lib_i) {
  x <- reads_raw[rows, lib_i]; median(x[x > 0])
}, numeric(1))
sites_with_reads <- function(rows) colSums(reads_raw[rows, , drop = FALSE] > 0)
block_median <- median_with_reads(in_block)
up_median <- median_with_reads(zone == "upstream flank"); down_median <- median_with_reads(zone == "downstream flank")
block_over_flanks <- block_median / median_with_reads(in_flanks)
stopifnot(!anyNA(c(block_median, up_median, down_median, block_over_flanks)))   # every library has reads in the block and in each flank

# Per CpG site: relative coverage among the libraries with and without the excess.
across <- function(m, f) apply(m, 1, f)
with_excess_median <- across(relative[, carries_excess, drop = FALSE], median)
without_excess_median <- across(relative[, !carries_excess, drop = FALSE], median)
with_excess_at_cutoff <- rowSums(at_cutoff[, carries_excess, drop = FALSE])
without_excess_at_cutoff <- rowSums(at_cutoff[, !carries_excess, drop = FALSE])
site_before <- max(which(site_start < block_from)); site_after <- min(which(site_start > block_to))
half_at_cutoff <- which(in_block & with_excess_at_cutoff >= n_carriers / 2)
stopifnot(length(half_at_cutoff) > 0L)

sec <- "a. block of excess coverage"
cutoff_text <- paste0("its upper cutoff (", config$upper_coverage_percentile, "th percentile of coverage of the library)")
fact(sec, "window", interval(window_from, window_to),
     paste0("The DML region and ", window_margin_bp, " bp on either side; raw coverage was read for the CpGs whose 0-based start lies between these two positions."))
fact(sec, "window_cpg_sites", whole(n_site), "CpG sites of the window with at least one read in at least one library.")
fact(sec, "runs_of_excess_coverage_in_window", whole(max(run_of_site, na.rm = TRUE)),
     paste0("Runs of CpGs with excess coverage in the window. A CpG has excess coverage when at least one library is at or above ", cutoff_text,
            "; such CpGs within ", run_gap_bp, " bp of each other form a run."))
fact(sec, "block", interval(block_from, block_to),
     paste("The run of CpGs with excess coverage that holds the DML region: its first and its last CpG.", zero_based))
fact(sec, "block_first_and_last_base", interval(block_from + 1L, block_to + 2L),
     paste("The same block as bases of the reference sequence, 1-based: the C of its first CpG and the G of its last.",
           "The positions of genes, 30-mers and second copies in this table are 1-based bases as well."))
fact(sec, "block_length_bp", whole(block_to + 2L - block_from), "From the C of the first CpG of the block to the G of the last.")
fact(sec, "block_cpg_sites", whole(sum(in_block)), "CpG sites inside the block with at least one read in at least one library.")
fact(sec, "block_cpg_sites_with_excess", x_of_y(sum(in_block & excess_site), sum(in_block)),
     "CpG sites of the block at which at least one library is at or above its upper cutoff.")
fact(sec, "block_cpg_sites_with_reads_in_every_library", x_of_y(sum(in_block & rowSums(reads_raw > 0) == n_lib), sum(in_block)),
     paste("CpG sites of the block with at least one read in each of the", n_lib, "libraries."))
fact(sec, "block_tested_cpgs_legacy", whole(sum(in_block & tested_legacy)), "Tested CpGs of legacy mode inside the block.")
fact(sec, "block_tested_cpgs_corrected", whole(sum(in_block & tested_corrected)), "Tested CpGs of corrected mode inside the block.")
fact(sec, "cutoff_reached_outside_block", x_of_y(sum(at_cutoff[!in_block, ]), sum(reads_raw[!in_block, ] > 0)),
     "Pairs of CpG site and library at or above the upper cutoff, of all pairs with reads, in the window outside the block.")
fact(sec, "block_under_stricter_rule", interval(site_start[min(half_at_cutoff)], site_start[max(half_at_cutoff)]),
     paste("For comparison: first and last CpG of the block at which at least half of the libraries with the excess are at or above their upper cutoff.",
           zero_based))
fact(sec, "next_cpg_site_upstream_distance_bp", whole(block_from - site_start[site_before] - 2L),
     "Bases between the first CpG of the block and the nearest CpG site with reads before it.")
fact(sec, "next_cpg_site_upstream_relative_coverage", decimal(with_excess_median[site_before], 2),
     "Coverage over the library median at that site; median of the libraries with the excess.")
fact(sec, "next_cpg_site_downstream_distance_bp", whole(site_start[site_after] - block_to - 2L),
     "Bases between the last CpG of the block and the nearest CpG site with reads after it.")
fact(sec, "next_cpg_site_downstream_relative_coverage", decimal(with_excess_median[site_after], 2),
     "Coverage over the library median at that site; median of the libraries with the excess.")
fact(sec, "block_relative_coverage_with_excess", a_to_b((block_median / library_median)[carries_excess], 2),
     "Median coverage of the CpG sites of the block over the library median; lowest and highest among the libraries with the excess.")
fact(sec, "block_relative_coverage_without_excess", a_to_b((block_median / library_median)[!carries_excess], 2),
     "Same; lowest and highest among the libraries without the excess.")
fact(sec, "flank_cpg_sites_upstream", whole(sum(zone == "upstream flank")),
     paste0("CpG sites with reads in at least one library whose start lies in the ", flank_bp, " bp before the first CpG of the block."))
fact(sec, "flank_cpg_sites_downstream", whole(sum(zone == "downstream flank")),
     paste0("CpG sites with reads in at least one library whose start lies in the ", flank_bp, " bp after the last CpG of the block."))
fact(sec, "flank_relative_coverage_upstream", a_to_b(up_median / library_median, 2),
     paste("Median coverage of the CpG sites of the upstream flank over the library median; lowest and highest among the", n_lib, "libraries."))
fact(sec, "flank_relative_coverage_downstream", a_to_b(down_median / library_median, 2), "Same for the downstream flank.")
fact(sec, "flank_relative_coverage", a_to_b(c(up_median, down_median) / c(library_median, library_median), 2),
     paste("Lowest and highest over both flanks and the", n_lib, "libraries."))
fact(sec, "flank_tested_cpgs_legacy", whole(sum(in_flanks & tested_legacy)), "Tested CpGs of legacy mode in the two flanks.")
fact(sec, "block_over_flanks_with_excess", a_to_b(block_over_flanks[carries_excess], 2),
     paste("Median coverage of the CpG sites of the block over the median coverage of the CpG sites of both flanks, per library;",
           "lowest and highest among the libraries with the excess."))
fact(sec, "block_over_flanks_without_excess", a_to_b(block_over_flanks[!carries_excess], 2),
     "Same; lowest and highest among the libraries without the excess.")
fact(sec, "libraries_with_more_coverage_in_block_than_in_flanks", x_of_y(sum(block_over_flanks > 1), n_lib),
     "Libraries whose median coverage of the CpG sites of the block is above their median coverage of the CpG sites of both flanks.")

sec <- "b. libraries that carry the excess"
fisher_p <- function(in_group, state) fisher.test(table(factor(in_group, c(FALSE, TRUE)), factor(state, c(FALSE, TRUE))))$p.value
with_excess_of <- function(members) x_of_y(sum(carries_excess & members), sum(members))
fact(sec, "libraries_with_excess", x_of_y(n_carriers, n_lib),
     paste0("Libraries with at least one tested CpG of the DML region at or above ", sub("^its", "their", cutoff_text), "."))
fact(sec, "libraries_classified_otherwise_by_median_rule", whole(sum(carries_excess != (region_median >= upper_cutoff))),
     "Libraries that change class under another rule: median coverage of the tested CpGs of the DML region at or above the upper cutoff.")
fact(sec, "region_relative_coverage_with_excess", a_to_b((region_median / library_median)[carries_excess], 2),
     "Median coverage of the tested CpGs of the DML region over the library median; lowest and highest among the libraries with the excess.")
fact(sec, "region_relative_coverage_without_excess", a_to_b((region_median / library_median)[!carries_excess], 2),
     "Same; lowest and highest among the libraries without the excess.")
fact(sec, "region_coverage_over_cutoff_with_excess", a_to_b((region_median / upper_cutoff)[carries_excess], 2),
     "Median coverage of the tested CpGs of the DML region over the upper cutoff of the library; lowest and highest among the libraries with the excess.")
fact(sec, "region_coverage_over_cutoff_without_excess", a_to_b((region_median / upper_cutoff)[!carries_excess], 2),
     "Same; lowest and highest among the libraries without the excess.")
fact(sec, "region_cpgs_at_cutoff_per_library_with_excess", paste(a_to_b(region_at_cutoff[carries_excess], 0), "of", n_region),
     "Tested CpGs of the DML region at or above the upper cutoff, per library; lowest and highest among the libraries with the excess.")
fact(sec, "libraries_without_excess_at_cutoff_in_block",
     x_of_y(sum(colSums(at_cutoff[in_block, !carries_excess, drop = FALSE]) > 0), n_lib - n_carriers),
     "Libraries without the excess that reach their upper cutoff at any CpG site of the block, tested or not.")
fact(sec, "libraries_at_cutoff_per_region_cpg", a_to_b(rowSums(at_cutoff[region_rows, , drop = FALSE]), 0),
     "Libraries at or above their upper cutoff, per tested CpG of the DML region; lowest and highest among those CpGs.")
fact(sec, "with_excess_heat_cohort", with_excess_of(is_heat), "Libraries with the excess among the heat cohort.")
fact(sec, "with_excess_control_cohort", with_excess_of(!is_heat), "Libraries with the excess among the control cohort.")
fact(sec, "with_excess_diploid", with_excess_of(!is_triploid), "Libraries with the excess among diploids.")
fact(sec, "with_excess_triploid", with_excess_of(is_triploid), "Libraries with the excess among triploids.")
for (group_name in levels(lib$group))
  fact(sec, paste0("with_excess_", gsub(" ", "_", group_name)), with_excess_of(lib$group == group_name),
       paste0("Libraries with the excess among the group ", group_name, "."))
fact(sec, "fisher_p_excess_by_ploidy", p_value(fisher_p(is_triploid, carries_excess)),
     paste("Fisher exact test, two-sided, of carrying the excess against ploidy;", n_lib, "libraries."))
fact(sec, "fisher_p_excess_by_cohort", p_value(fisher_p(is_heat, carries_excess)),
     paste("Fisher exact test, two-sided, of carrying the excess against cohort;", n_lib, "libraries."))
# The same counts and tests under a stricter rule, to show how far they depend on the rule.
every_cpg <- region_at_cutoff == n_region
of_every_cpg <- function(members) x_of_y(sum(every_cpg & members), sum(members))
fact(sec, "libraries_with_every_region_cpg_at_cutoff", x_of_y(sum(every_cpg), n_lib),
     "A stricter rule for comparison: libraries in which every tested CpG of the DML region is at or above the upper cutoff.")
fact(sec, "every_region_cpg_at_cutoff_heat_cohort", of_every_cpg(is_heat), "Such libraries among the heat cohort.")
fact(sec, "every_region_cpg_at_cutoff_control_cohort", of_every_cpg(!is_heat), "Such libraries among the control cohort.")
fact(sec, "every_region_cpg_at_cutoff_diploid", of_every_cpg(!is_triploid), "Such libraries among diploids.")
fact(sec, "every_region_cpg_at_cutoff_triploid", of_every_cpg(is_triploid), "Such libraries among triploids.")
fact(sec, "fisher_p_every_region_cpg_at_cutoff_by_ploidy", p_value(fisher_p(is_triploid, every_cpg)),
     "Fisher exact test, two-sided, of that state against ploidy.")
fact(sec, "fisher_p_every_region_cpg_at_cutoff_by_cohort", p_value(fisher_p(is_heat, every_cpg)),
     "Fisher exact test, two-sided, of that state against cohort.")
# The size of the excess as a continuous quantity, by ploidy and by cohort. Two
# measures; the second does not depend on the library median.
size_of_excess <- list(
  region_relative_coverage = list(value = region_median / library_median,
                                  text = "median coverage of the tested CpGs of the DML region over the library median"),
  block_over_flanks = list(value = block_over_flanks,
                           text = "median coverage of the CpG sites of the block over the median coverage of the CpG sites of both flanks"))
two_groups <- list(ploidy = list(in_second = is_triploid, names = c("diploid", "triploid")),
                   cohort = list(in_second = is_heat, names = c("control", "heat")))
for (measure in names(size_of_excess)) for (grouping in names(two_groups)) {
  x <- size_of_excess[[measure]]$value; grp <- two_groups[[grouping]]
  tied <- anyDuplicated(x) > 0                                   # the exact test needs untied values
  fact(sec, paste0(measure, "_median_", grp$names[1]), decimal(median(x[!grp$in_second]), 2),
       paste0("Median among ", grp$names[1], " libraries of the ", size_of_excess[[measure]]$text, "."))
  fact(sec, paste0(measure, "_median_", grp$names[2]), decimal(median(x[grp$in_second]), 2), paste0("Same among ", grp$names[2], " libraries."))
  fact(sec, paste0(measure, "_wilcoxon_p_by_", grouping),
       p_value(wilcox.test(x[!grp$in_second], x[grp$in_second], exact = !tied, correct = TRUE)$p.value),
       paste0("Wilcoxon rank-sum test, two-sided, ", if (tied) "by the normal approximation with continuity correction because values are tied" else "exact",
              ", of that quantity between ", grp$names[1], " and ", grp$names[2], " libraries."))
}

# ---- 5. (c) Leave one library out ------------------------------------------------------------
# Read-pooled difference, group 1 minus group 0, at each tested CpG of the DML
# region, from the libraries in `use`.
pooled_difference <- function(use, in_group1) {
  pooled <- function(members) rowSums(norm_meth[, members, drop = FALSE]) / rowSums(norm_reads[, members, drop = FALSE])
  100 * (pooled(use & in_group1) - pooled(use & !in_group1))
}
every_library <- rep(TRUE, n_lib)
contrast_group1 <- list(ploidy = is_triploid, heat = is_heat)    # triploid minus diploid; heat minus control
diffmeth_region <- lapply(diffmeth_window, function(d) d[start >= region_from & start <= region_to])
for (s in names(contrast_group1)) {
  d <- diffmeth_region[[s]]
  # The pooled difference of all libraries must be step 05's meth.diff, and the two DML rules must give step 05's DML.
  stopifnot(identical(d$start, united_region$start), max(abs(pooled_difference(every_library, contrast_group1[[s]]) - d$meth.diff)) < 1e-8,
            identical(sort(d$start[d$qvalue < config$dml_qvalue & abs(d$meth.diff) > config$dml_difference]), sort(dml_region[[s]])))
}
over_difference <- function(use, s) sum(abs(pooled_difference(use, contrast_group1[[s]])) > config$dml_difference)
share_pct <- function(members) 100 * colSums(norm_reads) / sum(norm_reads[, members])
loo_row <- function(label, left_out) {
  use <- !left_out
  data.table(libraries_left_out = label, n_left_out = sum(left_out),
             diploid_used = sum(use & !is_triploid), triploid_used = sum(use & is_triploid),
             control_used = sum(use & !is_heat), heat_used = sum(use & is_heat), region_cpgs = n_region,
             ploidy_cpgs_over_difference = over_difference(use, "ploidy"),
             ploidy_difference_median = round(median(pooled_difference(use, is_triploid)), 2),
             heat_cpgs_over_difference = over_difference(use, "heat"),
             heat_difference_median = round(median(pooled_difference(use, is_heat)), 2))
}
loo <- rbindlist(c(list(loo_row("none", !every_library)),
                   lapply(seq_len(n_lib), function(lib_i) loo_row(lib$library_id[lib_i], seq_len(n_lib) == lib_i)),
                   list(loo_row("all libraries with the excess", carries_excess))))
one_out <- loo[n_left_out == 1L]
remaining <- one_out$ploidy_cpgs_over_difference + one_out$heat_cpgs_over_difference
if (sum(remaining == min(remaining)) != 1L) stop("No single library has the largest effect in the leave-one-out.")
extreme <- which.min(remaining)                                  # index of the most influential library
without_excess_only <- loo[libraries_left_out == "all libraries with the excess"]

sec <- "c. leave one library out"
difference_text <- paste0("tested CpGs of the DML region whose read-pooled difference (normalized counts) exceeds ", config$dml_difference,
                          " points in absolute value")
for (s in names(contrast_group1)) {
  d <- diffmeth_region[[s]]; counts <- paste0(s, "_cpgs_over_difference")
  fact(sec, paste0(s, "_region_cpgs_failing_qvalue_rule"), x_of_y(sum(d$qvalue >= config$dml_qvalue), n_region),
       paste0("Tested CpGs of the DML region with q at or above ", config$dml_qvalue, " for the ", s, " contrast, all ", n_lib,
              " libraries (diffmeth_", s, "_legacy.csv.gz)."))
  fact(sec, paste0(s, "_region_largest_qvalue"), p_value(max(d$qvalue)), "Largest q-value among those CpGs.")
  fact(sec, paste0(s, "_region_cpgs_over_difference_all_libraries"), x_of_y(loo[n_left_out == 0L][[counts]], n_region),
       paste0("Number of ", difference_text, "; ", s, " contrast, all ", n_lib,
              " libraries. Equal to the number of DML when no CpG fails the q-value rule."))
  fact(sec, paste0(s, "_region_cpgs_over_difference_without_most_influential_library"), whole(one_out[[counts]][extreme]),
       paste("Same count with the most influential library left out. The difference rule only: the model is not refitted,",
             "so this is the largest possible number of DML."))
  fact(sec, paste0(s, "_region_cpgs_over_difference_without_one_other_library"), a_to_b(one_out[[counts]][-extreme], 0),
       paste("Same count with one of the other", n_lib - 1L, "libraries left out; lowest and highest among them."))
  fact(sec, paste0(s, "_region_cpgs_over_difference_libraries_without_excess_only"), whole(without_excess_only[[counts]]),
       sprintf("Same count from the libraries without the excess alone (%d diploid and %d triploid; %d control and %d heat).",
               without_excess_only$diploid_used, without_excess_only$triploid_used, without_excess_only$control_used, without_excess_only$heat_used))
}
fact(sec, "most_influential_library", lib$library_id[extreme],
     "The library whose omission leaves the fewest CpGs over the difference threshold, ploidy and heat contrasts added.")
fact(sec, "most_influential_library_group", as.character(lib$group[extreme]), "Ploidy and cohort of that library.")
fact(sec, "most_influential_library_relative_coverage", decimal(region_median[extreme] / library_median[extreme], 2),
     "Median coverage of the tested CpGs of the DML region over the library median, that library.")
fact(sec, "most_influential_library_pct_methylation", decimal(region_pct_meth[extreme], 1),
     "Read-pooled percent methylation of the tested CpGs of the DML region, raw counts, that library.")
fact(sec, "most_influential_library_pct_of_cohort_normalized_reads", decimal(share_pct(is_heat == is_heat[extreme])[extreme], 1),
     "Its share of the normalized reads of its cohort at the tested CpGs of the DML region.")
fact(sec, "most_influential_library_pct_of_ploidy_normalized_reads", decimal(share_pct(is_triploid == is_triploid[extreme])[extreme], 1),
     "Its share of the normalized reads of its ploidy at the tested CpGs of the DML region.")

# ---- 6. (d) Percent methylation per library ---------------------------------------------------
if (anyDuplicated(region_pct_meth)) stop("Two libraries have the same percent methylation; the exact Wilcoxon test needs untied values.")
sec <- "d. methylation per library"
pct_text <- "read-pooled percent methylation of the tested CpGs of the DML region, raw counts, one value per library"
compare_groups <- function(label, members, grouping, described) {
  in_second <- two_groups[[grouping]]$in_second; group_names <- two_groups[[grouping]]$names
  a <- region_pct_meth[members & !in_second]; b <- region_pct_meth[members & in_second]
  fact(sec, paste0(label, "_group_sizes"), paste(length(a), group_names[1], "and", length(b), group_names[2]),
       paste0("Libraries compared: ", described, "."))
  fact(sec, paste0(label, "_median_", group_names[1]), decimal(median(a), 1),
       paste0("Median of the ", pct_text, "; ", group_names[1], ", ", described, "."))
  fact(sec, paste0(label, "_median_", group_names[2]), decimal(median(b), 1), paste0("Same; ", group_names[2], "."))
  fact(sec, paste0(label, "_range_", group_names[1]), a_to_b(a, 1), paste0("Lowest and highest value; ", group_names[1], ", ", described, "."))
  fact(sec, paste0(label, "_range_", group_names[2]), a_to_b(b, 1), paste0("Same; ", group_names[2], "."))
  fact(sec, paste0(label, "_wilcoxon_p_exact"), p_value(wilcox.test(a, b, exact = TRUE)$p.value),
       paste0("Wilcoxon rank-sum test, two-sided, exact; ", group_names[1], " against ", group_names[2], ", ", described, "."))
  fact(sec, paste0(label, "_wilcoxon_p_normal_approximation"), p_value(wilcox.test(a, b, exact = FALSE, correct = TRUE)$p.value),
       "Same test by the normal approximation with continuity correction.")
}
compare_groups("ploidy_all_libraries", every_library, "ploidy", paste("all", n_lib, "libraries"))
compare_groups("ploidy_heat_cohort", is_heat, "ploidy", "heat cohort")
compare_groups("ploidy_control_cohort", !is_heat, "ploidy", "control cohort")
compare_groups("ploidy_without_most_influential_library", seq_len(n_lib) != extreme, "ploidy", "all libraries but the most influential one of part c")
compare_groups("cohort_all_libraries", every_library, "cohort", paste("all", n_lib, "libraries"))
fact(sec, "library_with_highest_pct_methylation", lib$library_id[which.max(region_pct_meth)],
     paste0("Library with the highest ", pct_text, "."))

# ---- 7. (e) Genes, exons and transposable elements at the block ----------------------------------
log_msg("annotation of the block %s", interval(block_from, block_to))
block_gr <- GRanges(locus_sequence, IRanges(block_from + 1L, block_to + 2L))
block_cpg <- cpg_ranges(locus_sequence, site_start[in_block])
genes <- read_genes()
block_genes <- genes[overlap_pairs(block_gr, GRanges(genes$chr, IRanges(genes$start, genes$end)))$feature]
exon_gr <- read_feature_track("exon", keep_attributes = TRUE)
exon_gr <- exon_gr[overlaps_track(exon_gr, block_gr)]
exon_gene_id <- gff_attribute(S4Vectors::mcols(exon_gr)$attributes, "GeneID:")
sec <- "e. annotation of the block"
fact(sec, "genes_overlapping_block", whole(nrow(block_genes)), "Records of the gene track that overlap the block by at least one base.")
for (gene_i in seq_len(nrow(block_genes))) {
  gene_row <- block_genes[gene_i]; key <- paste0("gene_", gene_row$gene)
  gene_gr <- GRanges(gene_row$chr, IRanges(gene_row$start, gene_row$end))
  gene_exons <- reduce(exon_gr[exon_gene_id == gene_row$gene_id], ignore.strand = TRUE)     # exons of all its transcripts, merged
  exon_part <- GenomicRanges::intersect(gene_exons, block_gr, ignore.strand = TRUE)
  fact(sec, paste0(key, "_biotype"), gene_row$biotype, "Biotype of a gene that overlaps the block (gene track).")
  fact(sec, paste0(key, "_position"), sprintf("%s:%d-%d (%s)", gene_row$chr, gene_row$start, gene_row$end, gene_row$strand),
       "Its first and last base, 1-based, and its strand.")
  fact(sec, paste0(key, "_block_bases_inside_gene"), x_of_y(sum(width(GenomicRanges::intersect(gene_gr, block_gr))), width(block_gr)),
       "Bases of the block inside the gene, of all bases of the block.")
  fact(sec, paste0(key, "_exons_overlapping_block"), whole(length(exon_part)),
       "Exons of the gene that overlap the block, exons of all its transcripts merged (exon track).")
  fact(sec, paste0(key, "_block_bases_in_exons"), whole(sum(width(exon_part))), "Bases of the block inside those exons.")
  if (length(exon_part) > 0L) {
    in_exon <- overlaps_track(block_cpg, gene_exons); in_intron <- overlaps_track(block_cpg, gene_gr) & !in_exon
    fact(sec, paste0(key, "_block_cpg_sites_in_exons"), whole(sum(in_exon)), "CpG sites of the block that overlap an exon of the gene.")
    fact(sec, paste0(key, "_block_cpg_sites_in_introns"), whole(sum(in_intron)), "CpG sites of the block inside the gene and outside its exons.")
    fact(sec, paste0(key, "_relative_coverage_in_exons"), decimal(median(with_excess_median[in_block][in_exon]), 2),
         "Coverage over the library median, median of the libraries with the excess at each CpG site; median of the CpG sites in exons.")
    fact(sec, paste0(key, "_relative_coverage_in_introns"), decimal(median(with_excess_median[in_block][in_intron]), 2),
         "Same; median of the CpG sites in introns.")
  }
}
te_gr <- read_feature_track("TE")
te_here <- te_gr[as.character(seqnames(te_gr)) == locus_sequence]
te_end_before <- end(te_here)[end(te_here) < start(block_gr)]; te_start_after <- start(te_here)[start(te_here) > end(block_gr)]
stopifnot(length(te_end_before) > 0L, length(te_start_after) > 0L)
fact(sec, "transposable_elements_overlapping_block", whole(sum(overlaps_track(te_gr, block_gr))),
     "Records of the transposable element track (RepeatMasker) that overlap the block.")
fact(sec, "block_cpg_sites_in_transposable_elements", x_of_y(sum(overlaps_track(block_cpg, te_gr)), sum(in_block)),
     "CpG sites of the block that overlap a record of that track.")
fact(sec, "nearest_transposable_element_upstream_distance_bp", whole(start(block_gr) - max(te_end_before) - 1L),
     "Bases between the block and the last transposable element that ends before it.")
fact(sec, "nearest_transposable_element_downstream_distance_bp", whole(min(te_start_after) - end(block_gr) - 1L),
     "Bases between the block and the first transposable element that begins after it.")
rm(te_gr, te_here, exon_gr, genes); invisible(gc(verbose = FALSE))

# ---- 8. (f) The block in the reference sequence --------------------------------------------------
log_msg("reading the reference genome")
genome <- readDNAStringSet(fetch_gannet_file("cgigas_uk_roslin_v1_genomic-mito.fa", config, manifest))
names(genome) <- sub(" .*", "", names(genome))
stopifnot(locus_sequence %in% names(genome), !anyDuplicated(names(genome)))
locus_dna <- genome[[locus_sequence]]
if (!all(as.character(Views(locus_dna, start = site_start + 1L, width = 2L)) == "CG"))
  stop("A CpG site of the coverage files is not CG in the reference; the genome in the cache is not the one the reads were aligned to.")
block_dna <- subseq(locus_dna, block_from + 1L, block_to + 2L)
up_dna <- subseq(locus_dna, block_from + 1L - flank_bp, block_from); down_dna <- subseq(locus_dna, block_to + 3L, block_to + 2L + flank_bp)
base_counts <- function(dna) alphabetFrequency(dna, baseOnly = TRUE)
gc_pct <- function(dna) { n <- base_counts(dna); 100 * (n[["C"]] + n[["G"]]) / (n[["A"]] + n[["C"]] + n[["G"]] + n[["T"]]) }
if (base_counts(block_dna)[["other"]] > 0) stop("The block holds bases other than A, C, G, T; the word count assumes none.")

# Every word of kmer_width bases of the block, counted on both strands of every reference sequence.
kmers <- as.character(Views(block_dna, start = seq_len(length(block_dna) - kmer_width + 1L), width = kmer_width))
kmers_rc <- as.character(reverseComplement(DNAStringSet(kmers)))
dictionary <- unique(c(kmers, kmers_rc))
log_msg("counting %d %d-mers of the block in %d reference sequences", length(kmers), kmer_width, length(genome))
hits <- vcountPDict(PDict(DNAStringSet(dictionary)), genome)      # one row per dictionary entry, one column per reference sequence
dimnames(hits) <- list(dictionary, names(genome))
copies_by_sequence <- hits[kmers, , drop = FALSE] + hits[kmers_rc, , drop = FALSE]
own_reverse <- kmers == kmers_rc                                 # a word equal to its reverse complement is counted once
copies_by_sequence[own_reverse, ] <- hits[kmers[own_reverse], , drop = FALSE]
copies <- rowSums(copies_by_sequence)
stopifnot(all(copies >= 1L))                                     # each word is found at least at its own place
sampled <- seq(1L, length(kmers), by = kmer_sample_step)

# Where the second copy lies, for the words with exactly two copies.
two_copies <- which(copies == 2L)
elsewhere <- copies_by_sequence[two_copies, , drop = FALSE]
elsewhere[, locus_sequence] <- elsewhere[, locus_sequence] - 1L   # minus the copy in the block itself
stopifnot(all(rowSums(elsewhere) == 1L))
second_sequence <- colnames(elsewhere)[max.col(elsewhere, ties.method = "first")]
second_copy <- rbindlist(lapply(unique(second_sequence), function(seq_name) {
  patterns <- unique(c(kmers[two_copies[second_sequence == seq_name]], kmers_rc[two_copies[second_sequence == seq_name]]))
  found <- matchPDict(PDict(DNAStringSet(patterns)), genome[[seq_name]])
  data.table(sequence = seq_name, first_base = unlist(startIndex(found)))
}))
# Matches inside the block are the words themselves, or second copies that also lie in the block.
second_copy <- second_copy[!(sequence == locus_sequence & first_base > block_from & first_base + kmer_width - 1L <= block_to + 2L)]
stopifnot(nrow(second_copy) <= length(two_copies))
places <- NULL
if (nrow(second_copy) > 0L) {
  setorder(second_copy, sequence, first_base)
  second_copy[, place := cumsum(c(TRUE, sequence[-1L] != sequence[-.N] | diff(first_base) > copy_gap_bp))]
  places <- second_copy[, .(sequence = sequence[1L], from = min(first_base), to = max(first_base) + kmer_width - 1L, n_kmers = .N), by = place][
    order(-n_kmers, sequence, from)]
}
many_copies <- which(copies >= 10L)                              # classes of the summary: 1, 2, 3 to 9, 10 or more copies

sec <- "f. reference sequence of the block"
kmer_name <- paste0(kmer_width, "-mers"); n_block_kmers <- length(kmers)
fact(sec, "block_bases", whole(length(block_dna)), "Length of the block in the reference: from the C of its first CpG to the G of its last.")
fact(sec, "block_cpg_dinucleotides_in_reference", whole(countPattern("CG", block_dna)), "CG dinucleotides in the reference sequence of the block.")
fact(sec, "block_gc_pct", decimal(gc_pct(block_dna), 1), "G plus C, percent of the bases of the block.")
fact(sec, "flank_gc_pct_upstream", decimal(gc_pct(up_dna), 1), paste0("G plus C, percent of the ", flank_bp, " bases before the block."))
fact(sec, "flank_gc_pct_downstream", decimal(gc_pct(down_dna), 1), paste0("G plus C, percent of the ", flank_bp, " bases after the block."))
fact(sec, "block_kmers", whole(n_block_kmers), paste0(kmer_name, " of the block, one starting at every base."))
fact(sec, "block_kmers_once_in_reference", x_of_y(sum(copies == 1L), n_block_kmers),
     paste0(kmer_name, " of the block that occur once in the reference, both strands of all ", length(genome), " sequences counted."))
fact(sec, "block_kmers_once_in_reference_pct", decimal(100 * mean(copies == 1L), 1), "The same as a percentage.")
fact(sec, "block_kmers_twice_in_reference", x_of_y(sum(copies == 2L), n_block_kmers),
     paste0(kmer_name, " of the block that occur exactly twice in the reference."))
fact(sec, "block_kmers_3_to_9_times_in_reference", x_of_y(sum(copies >= 3L & copies <= 9L), n_block_kmers),
     paste0(kmer_name, " of the block that occur 3 to 9 times in the reference."))
fact(sec, "block_kmers_10_or_more_times_in_reference", x_of_y(length(many_copies), n_block_kmers),
     paste0(kmer_name, " of the block that occur 10 times or more in the reference."))
if (length(many_copies) > 0L) {
  fact(sec, "block_kmers_10_or_more_times_copies", a_to_b(copies[many_copies], 0),
       paste0("Copies in the reference of those ", kmer_name, "; lowest and highest."))
  fact(sec, "block_kmers_10_or_more_times_position", interval(block_from + min(many_copies), block_from + max(many_copies) + kmer_width - 1L),
       paste0("First and last base of the block that those ", kmer_name, " cover, 1-based."))
}
fact(sec, "sampled_kmers_once_in_reference", x_of_y(sum(copies[sampled] == 1L), length(sampled)),
     paste0("The same count for a sample: one of the ", kmer_name, " every ", kmer_sample_step, " bases, from the first base of the block."))
fact(sec, "sampled_kmers_twice_in_reference", x_of_y(sum(copies[sampled] == 2L), length(sampled)),
     paste0(kmer_name, " of the sample that occur exactly twice in the reference."))
fact(sec, "sampled_kmers_3_or_more_times_in_reference", x_of_y(sum(copies[sampled] >= 3L), length(sampled)),
     paste0(kmer_name, " of the sample that occur three times or more."))
fact(sec, "second_copies_inside_block", x_of_y(length(two_copies) - nrow(second_copy), length(two_copies)),
     paste0("Of the ", kmer_name, " with exactly two copies, those whose second copy also lies inside the block."))
fact(sec, "second_copies_elsewhere_on_same_sequence", x_of_y(sum(second_copy$sequence == locus_sequence), length(two_copies)),
     paste0("Of the ", kmer_name, " with exactly two copies, those whose second copy lies outside the block on the sequence of the block."))
for (place_i in seq_len(min(n_places, NROW(places)))) {
  fact(sec, paste0("second_copies_place_", place_i), interval(places$from[place_i], places$to[place_i], places$sequence[place_i]),
       paste0("A place that holds second copies: first and last base, 1-based. Second copies within ", copy_gap_bp,
              " bp of each other on one sequence are one place; places are ranked by the number of ", kmer_name, "."))
  fact(sec, paste0("second_copies_place_", place_i, "_kmers"), x_of_y(places$n_kmers[place_i], length(two_copies)),
       paste0(kmer_name, " with exactly two copies whose second copy lies there."))
  place_dna <- subseq(genome[[places$sequence[place_i]]], places$from[place_i], places$to[place_i])
  place_at_cutoff <- cutoff_rows[chr == places$sequence[place_i] & start + 1L >= places$from[place_i] & start + 2L <= places$to[place_i], .N]
  fact(sec, paste0("second_copies_place_", place_i, "_cutoff_reached"), x_of_y(place_at_cutoff, countPattern("CG", place_dna) * n_lib),
       paste("Pairs of CpG site and library at or above the upper cutoff inside that place, of its CG dinucleotides in the reference times", n_lib,
             "libraries. Compare block_cutoff_reached."))
}
stopifnot(cutoff_rows[chr == locus_sequence & start >= block_from & start <= block_to, .N] == sum(at_cutoff[in_block, ]))
fact(sec, "block_cutoff_reached", x_of_y(sum(at_cutoff[in_block, ]), countPattern("CG", block_dna) * n_lib),
     paste("Pairs of CpG site and library at or above the upper cutoff inside the block, of its CG dinucleotides in the reference times", n_lib, "libraries."))
rm(genome, locus_dna, hits, copies_by_sequence); invisible(gc(verbose = FALSE))

# ---- 9. Output tables ----------------------------------------------------------------------------
summary_table <- rbindlist(facts)
stopifnot(!anyDuplicated(summary_table$quantity))
fwrite(summary_table, out_file("summary"))

per_library <- data.table(
  library_id = lib$library_id, cohort = lib$cohort, ploidy = as.character(lib$ploidy), treatment = as.character(lib$treatment),
  acute_heat_shock = lib$acute_heat_shock,
  library_median_coverage = library_median, upper_cutoff = round(upper_cutoff, 3),
  region_tested_cpgs = n_region, region_median_coverage = region_median,
  region_min_coverage = apply(reads_raw[region_rows, , drop = FALSE], 2, min),
  region_max_coverage = apply(reads_raw[region_rows, , drop = FALSE], 2, max),
  region_coverage_over_library_median = round(region_median / library_median, 2),
  region_coverage_over_cutoff = round(region_median / upper_cutoff, 2),
  region_cpgs_at_or_above_cutoff = region_at_cutoff, carries_excess = carries_excess,
  region_pct_methylation = round(region_pct_meth, 1),
  region_normalized_reads = colSums(norm_reads),
  pct_of_ploidy_normalized_reads = round(ifelse(is_triploid, share_pct(is_triploid), share_pct(!is_triploid)), 1),
  pct_of_cohort_normalized_reads = round(ifelse(is_heat, share_pct(is_heat), share_pct(!is_heat)), 1),
  block_cpg_sites_with_reads = sites_with_reads(in_block), block_median_coverage = block_median,
  block_coverage_over_library_median = round(block_median / library_median, 2),
  block_cpg_sites_at_or_above_cutoff = colSums(at_cutoff[in_block, , drop = FALSE]),
  upstream_flank_cpg_sites_with_reads = sites_with_reads(zone == "upstream flank"), upstream_flank_median_coverage = up_median,
  upstream_flank_coverage_over_library_median = round(up_median / library_median, 2),
  downstream_flank_cpg_sites_with_reads = sites_with_reads(zone == "downstream flank"), downstream_flank_median_coverage = down_median,
  downstream_flank_coverage_over_library_median = round(down_median / library_median, 2),
  block_over_flanks = round(block_over_flanks, 2))
fwrite(per_library, out_file("per_library"))

loo_groups <- data.table(libraries_left_out = lib$library_id, ploidy = as.character(lib$ploidy), treatment = as.character(lib$treatment),
                         carries_excess = carries_excess, pct_of_ploidy_normalized_reads = per_library$pct_of_ploidy_normalized_reads,
                         pct_of_cohort_normalized_reads = per_library$pct_of_cohort_normalized_reads)
loo_table <- merge(loo, loo_groups, by = "libraries_left_out", all.x = TRUE, sort = FALSE)
setcolorder(loo_table, c("libraries_left_out", "n_left_out", names(loo_groups)[-1L]))
fwrite(loo_table, out_file("leave_one_out"))

tested_value <- function(x, digits) { out <- rep(NA_real_, n_site); out[tested_legacy] <- round(x, digits); out }
profile <- data.table(
  sequence = locus_sequence, start = site_start, zone = zone, in_dml_region = in_region,
  tested_legacy = tested_legacy, tested_corrected = tested_corrected,
  ploidy_dml = site_start %in% dml_window$ploidy, heat_dml = site_start %in% dml_window$heat,
  in_both_dml_lists = site_start %in% dml_window$shared,
  ploidy_difference = tested_value(diffmeth_window$ploidy$meth.diff, 1), heat_difference = tested_value(diffmeth_window$heat$meth.diff, 1),
  libraries_with_reads = rowSums(reads_raw > 0),
  libraries_at_cutoff_with_excess = with_excess_at_cutoff, libraries_at_cutoff_without_excess = without_excess_at_cutoff,
  relative_coverage_with_excess_median = round(with_excess_median, 2),
  relative_coverage_with_excess_min = round(across(relative[, carries_excess, drop = FALSE], min), 2),
  relative_coverage_with_excess_max = round(across(relative[, carries_excess, drop = FALSE], max), 2),
  relative_coverage_without_excess_median = round(without_excess_median, 2),
  relative_coverage_without_excess_min = round(across(relative[, !carries_excess, drop = FALSE], min), 2),
  relative_coverage_without_excess_max = round(across(relative[, !carries_excess, drop = FALSE], max), 2))
fwrite(profile, out_file("coverage_profile"))

# ---- Summary for the log ---------------------------------------------------------------------------
options(width = 200)
cat("\nLibraries (coverage of the tested CpGs of the DML region):\n")
print(per_library[, .(library_id, ploidy, treatment, region_median_coverage, over_median = region_coverage_over_library_median,
                      over_cutoff = region_coverage_over_cutoff, at_cutoff = region_cpgs_at_or_above_cutoff, carries_excess,
                      pct_meth = region_pct_methylation, block_over_flanks)])
cat("\nLeave one library out (CpGs of the region over the difference threshold):\n")
print(loo_table[, .(libraries_left_out, ploidy_cpgs_over_difference, ploidy_difference_median, heat_cpgs_over_difference, heat_difference_median)])
cat("\nFacts:\n")
print(as.data.frame(summary_table[, .(section, quantity, value)]), right = FALSE, row.names = FALSE)
log_msg("wrote %s", paste(basename(vapply(c("summary", "per_library", "leave_one_out", "coverage_profile"), out_file, "")), collapse = ", "))

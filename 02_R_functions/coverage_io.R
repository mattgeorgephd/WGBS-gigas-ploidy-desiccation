# coverage_io.R
# Reading Bismark merged CpG coverage files.
#
# File format (coverage2cytosine --merge_CpG --zero_based), tab separated, no header:
#   sequence, start (0-based position of the C), end (start + 2), percent
#   methylation, methylated read count, unmethylated read count.
# One row per CpG dinucleotide, reads from both strands combined.

# Read one coverage file as a data.table with columns chr, start, nC, nT.
#
# drop_first_line = TRUE reproduces what the legacy analysis read. The legacy
# scripts called methylKit::methRead() without `header = FALSE`; its default is
# header = TRUE, so the first line of every (headerless) coverage file was taken
# as column names and that CpG was lost (issue I-38). Use first_line_lost(mode)
# to decide.
read_cov <- function(path, drop_first_line = FALSE) {
  d <- data.table::fread(path, header = FALSE, sep = "\t", select = c(1, 2, 5, 6),
                         col.names = c("chr", "start", "nC", "nT"), showProgress = FALSE)
  if (drop_first_line) d <- d[-1L]
  d
}

# Whether a mode reads the coverage files the way the legacy scripts did.
first_line_lost <- function(mode) {
  if (!mode %in% c("legacy", "corrected")) stop("Unknown mode: ", mode)
  mode == "legacy"
}

# A CpG position as one number, so that positions can be matched and intersected
# quickly: index of the sequence in `chr_levels` times 1e9, plus the start
# coordinate. Exact in double precision (sequences are shorter than 1e9 bp).
cpg_key <- function(chr, start, chr_levels) {
  i <- match(chr, chr_levels)
  if (anyNA(i)) stop("Sequence name not in chr_levels")
  i * 1e9 + start
}

key_to_position <- function(key, chr_levels) {
  i <- key %/% 1e9
  data.table::data.table(chr = chr_levels[i], start = as.integer(key - i * 1e9))
}

# Upper coverage cutoff of one library as methylKit::filterByCoverage(hi.perc = p)
# computes it: the p-th percentile of coverage among rows with at least
# `import_min_coverage` reads. CpGs with coverage at or above the cutoff are removed.
upper_coverage_cutoff <- function(coverage, import_min_coverage, percentile) {
  unname(stats::quantile(coverage[coverage >= import_min_coverage], percentile / 100))
}

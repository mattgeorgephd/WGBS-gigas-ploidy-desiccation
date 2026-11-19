# fastq_read_names.R
#
# Purpose   Which sequencer run, flowcell and lane the reads of each of the 22
#           libraries come from, taken from the read names at the start of each
#           trimmed FASTQ file (read 1). Methods 2.2 of the manuscript says that
#           the libraries were pooled into one lane of a single NovaSeq flowcell
#           (issues I-01 and I-28).
#
# Run       Rscript 06_diagnostics/fastq_read_names.R     (about one minute; needs
#           the network, because the FASTQ files are not pipeline inputs and are
#           read from gannet as a stream; nothing is saved to the cache)
#
# Reads     03_input_files/metadata/wgbs_libraries.csv
#           the first reads of <library>_R1.fastp-trim.<date>.fq.gz in the two
#           folders named below
# Writes    06_diagnostics/fastq_read_names.csv   one row per library and
#           combination of instrument, run, flowcell and lane found
#
# Notes     - An Illumina read name starts with instrument:run:flowcell:lane.
#           - Only the first 10,000 reads of each file are read. A library that
#             was also sequenced on a further lane whose reads stand later in the
#             file would not show it; the table therefore proves that two
#             flowcells were used, and it does not prove that each project was
#             confined to one lane.
#           - Read 2 files are not read; a pair shares its name.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))

n_reads <- 10000L
fastq_folder <- c(
  zr3534 = "https://gannet.fish.washington.edu/Atumefaciens/20201202_cgig_fastp-10bp-5-3-prime_ronit-ploidy-wgbs",
  zr3644 = "https://gannet.fish.washington.edu/Atumefaciens/20201206_cgig_fastp-10bp-5-3-prime_ploidy-pH-wgbs"
)
fastq_date <- c(zr3534 = "20201202", zr3644 = "20201206")

lib <- read_libraries()
stopifnot(all(lib$zymo_project %in% names(fastq_folder)))

# Names of the first n reads of a gzipped FASTQ file, read as a stream.
first_read_names <- function(u, n) {
  con <- gzcon(url(u, open = "rb"))
  on.exit(close(con))
  x <- readLines(con, n = 4L * n)
  stopifnot(length(x) == 4L * n)
  x[seq(1L, length(x), by = 4L)]
}

old <- options(timeout = config$download_timeout_s)
found <- rbindlist(lapply(seq_len(nrow(lib)), function(i) {
  project <- as.character(lib$zymo_project[i])
  file <- sprintf("%s_R1.fastp-trim.%s.fq.gz", lib$library_id[i], fastq_date[[project]])
  name <- first_read_names(file.path(fastq_folder[[project]], file), n_reads)
  stopifnot(all(startsWith(name, "@")))
  part <- tstrsplit(sub("^@", "", name), ":", fixed = TRUE, keep = 1:4)
  tab <- data.table(instrument = part[[1]], run = part[[2]], flowcell = part[[3]], lane = part[[4]])[, .(n_reads = .N), by = .(instrument, run, flowcell, lane)]
  cbind(data.table(library_id = lib$library_id[i], zymo_project = project, cohort = lib$cohort[i],
                   treatment = lib$treatment[i], ploidy = lib$ploidy[i], fastq_file = file, reads_examined = n_reads), tab)
}))
options(old)

fwrite(found, here("06_diagnostics", "fastq_read_names.csv"))

by_treatment <- found[, .(libraries = uniqueN(library_id)), by = .(treatment, instrument, run, flowcell, lane)]
print(by_treatment)
cat(sprintf("%d libraries, %d combinations of run, flowcell and lane; wrote 06_diagnostics/fastq_read_names.csv\n",
            uniqueN(found$library_id), nrow(unique(found[, .(run, flowcell, lane)]))))

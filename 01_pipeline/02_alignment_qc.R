# 02_alignment_qc.R
#
# Purpose   One table of alignment statistics per library, from the Bismark
#           reports of the June 2021 alignment (the run that produced the
#           coverage files), and a comparison of the two sequencing batches.
#           Percent CpG methylation per library, used in Figure 1B, comes from
#           this table.
#
# Run       Rscript 01_pipeline/02_alignment_qc.R     (needs no other step; under a minute)
#
# Inputs    03_input_files/metadata/wgbs_libraries.csv
#           03_input_files/bismark_reports/061021-big_bismark_summary_report.txt
#           03_input_files/bismark_reports/<library>_R1_bismark_bt2_PE_report.txt
# Outputs   04_output/02_alignment_qc/tables/bismark_by_library.csv
#           04_output/02_alignment_qc/tables/bismark_by_cohort.csv
#           04_output/02_alignment_qc/figures/alignment_qc.(png|pdf)
#           04_output/02_alignment_qc/session_info.txt
#
# Manuscript items: Methods 2.4 (alignment); input to Figure 1B (step 03).
# Notes     Two kinds of methylation counts exist. The summary report counts
#           calls after deduplication (the reads the coverage files are built
#           from); the alignment report counts calls before deduplication. Both
#           are kept; pct_mCpG is after deduplication. See issue I-13.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "ggplot2", "patchwork"))

step <- "02_alignment_qc"
t0 <- start_step(step, "alignment statistics per library")
require_run_kind("paper", "this step does not depend on the analysis and writes the paper's outputs")
lib <- read_libraries()

# ---- Summary report (after deduplication) ------------------------------------
smry <- fread(input_path("bismark_reports", "061021-big_bismark_summary_report.txt"), sep = "\t")
setnames(smry, c("file", "read_pairs", "aligned_pairs", "unaligned_pairs", "ambiguous_pairs",
                 "no_genomic_sequence", "duplicate_pairs", "unique_pairs", "total_C",
                 "meth_CpG", "unmeth_CpG", "meth_CHG", "unmeth_CHG", "meth_CHH", "unmeth_CHH"))
smry[, library_id := sub("_R1_bismark_bt2_pe\\.bam$", "", file)]
smry <- smry[library_id %in% lib$library_id]
stopifnot(nrow(smry) == nrow(lib))

# ---- Per-library alignment reports (before deduplication) --------------------
read_pe_report <- function(library_id) {
  x <- readLines(input_path("bismark_reports", paste0(library_id, "_R1_bismark_bt2_PE_report.txt")))
  num <- function(pattern) as.numeric(sub(".*:\\s*", "", grep(pattern, x, value = TRUE)[1]))
  data.table(
    library_id       = library_id,
    bismark_version  = sub(".*version: (v[0-9.]+).*", "\\1", x[1]),
    bowtie2_options  = sub(".*specified options: ", "", x[2]),
    non_directional  = any(grepl("--non_directional", x)),
    pe_read_pairs    = num("^Sequence pairs analysed in total"),
    meth_CpG_predup   = num("^Total methylated C's in CpG context"),
    unmeth_CpG_predup = num("^Total unmethylated C's in CpG context")
  )
}
pe <- rbindlist(lapply(lib$library_id, read_pe_report))

# ---- Combine -------------------------------------------------------------------
qc <- merge(as.data.table(lib[, c("library_id", "legacy_sample_index", "zymo_project", "cohort", "ploidy",
                                  "treatment", "acute_heat_shock", "animal_id")]),
            smry[, !"file"], by = "library_id")
qc <- merge(qc, pe, by = "library_id")
stopifnot(all(qc$read_pairs == qc$pe_read_pairs))     # both reports describe the same run
qc[, `:=`(
  mapping_efficiency_pct = 100 * aligned_pairs / read_pairs,
  duplication_pct        = 100 * duplicate_pairs / aligned_pairs,
  pct_mCpG               = 100 * meth_CpG / (meth_CpG + unmeth_CpG),
  pct_mCHG               = 100 * meth_CHG / (meth_CHG + unmeth_CHG),
  pct_mCHH               = 100 * meth_CHH / (meth_CHH + unmeth_CHH),
  pct_mCpG_predup        = 100 * meth_CpG_predup / (meth_CpG_predup + unmeth_CpG_predup)
)]
qc[, pe_read_pairs := NULL]
setorder(qc, legacy_sample_index)
write_table(qc, output_path(step, "tables", "bismark_by_library.csv"))

# ---- Batch comparison: the heat and control libraries are different cohorts
#      and different sequencing projects (issue I-01) -------------------------
vars <- c("read_pairs", "mapping_efficiency_pct", "duplication_pct", "unique_pairs",
          "pct_mCpG", "pct_mCHG", "pct_mCHH")
long <- qc[, c(list(library_id = library_id, cohort = cohort, zymo_project = zymo_project),
                   lapply(.SD, as.numeric)), .SDcols = vars]
by_cohort <- melt(long, id.vars = c("library_id", "cohort", "zymo_project"), measure.vars = vars)[
  , .(n = .N, mean = mean(value), min = min(value), max = max(value)),
  by = .(variable, cohort, zymo_project)]
write_table(by_cohort, output_path(step, "tables", "bismark_by_cohort.csv"))

stopifnot(length(unique(qc$bismark_version)) == 1, all(qc$non_directional))
log_msg("%d libraries, Bismark %s, bowtie2 options: %s", nrow(qc), qc$bismark_version[1], qc$bowtie2_options[1])
print(dcast(by_cohort, variable ~ cohort, value.var = "mean"), digits = 4)

# ---- Figure ----------------------------------------------------------------------
panel <- function(y, ylab) {
  ggplot(qc, aes(x = treatment, y = .data[[y]], colour = ploidy, shape = acute_heat_shock)) +
    geom_point(size = 2.6, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5, seed = 1)) +
    scale_colour_manual(values = config$colours$ploidy) +
    scale_shape_manual(values = c(no = 16, yes = 17), name = "acute heat shock") +
    labs(x = NULL, y = ylab) + theme_wgbs()
}
fig <- (panel("read_pairs", "Read pairs") + scale_y_continuous(labels = function(v) paste0(v / 1e6, " M"))) +
  panel("duplication_pct", "Duplicate alignments (%)") +
  panel("pct_mCpG", "CpG methylation (%), deduplicated") +
  panel("pct_mCHH", "CHH methylation (%)") +
  plot_layout(ncol = 2, guides = "collect") + plot_annotation(tag_levels = "A")
save_figure(fig, output_path(step, "figures", "alignment_qc"), width = 9, height = 7)

finish_step(step, t0)

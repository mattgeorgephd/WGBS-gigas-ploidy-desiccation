# check_legacy_outputs.R
#
# Purpose   Compare what the pipeline regenerates, in legacy mode, with the
#           outputs of the legacy analysis kept in 08_archive/. Every comparison
#           becomes one row of 06_diagnostics/legacy_check_results.csv with the
#           legacy value, the regenerated value, and a status.
#
# Run       Rscript 06_diagnostics/check_legacy_outputs.R
#           Checks whose pipeline outputs do not exist yet are reported as
#           "not run" rather than failing.
#
# Reads     04_output/<step>/tables/*      (pipeline outputs)
#           08_archive/...                 (legacy outputs; never modified)
# Writes    06_diagnostics/legacy_check_results.csv

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "readxl"))

results <- list()
record <- function(check, item, legacy, regenerated, tolerance = 0, note = "") {
  status <- if (is.na(regenerated) || is.na(legacy)) "not run" else
    if (isTRUE(abs(as.numeric(legacy) - as.numeric(regenerated)) <= tolerance)) "match" else "differs"
  results[[length(results) + 1]] <<- data.table(check = check, item = item, legacy = as.character(legacy),
                                                regenerated = as.character(regenerated), status = status, note = note)
}
have <- function(step, file) file.exists(here("04_output", step, "tables", file))
out  <- function(step, file) here("04_output", step, "tables", file)

# ---- Step 03: ELISA values and model table ------------------------------------------
if (have("03_global_methylation", "elisa_percent_5mC.csv")) {
  new <- fread(out("03_global_methylation", "elisa_percent_5mC.csv"))
  old <- as.data.table(read_excel(archive_path("data", "ELISA-percent-mCPG-summary.xlsx"), sheet = "summary"))
  m <- merge(new[, .(sample = paste0(animal_id, "-C"), percent_5mC)], old[, .(sample, percent)], by = "sample")
  record("ELISA percent 5-mC", "animals matched", nrow(old), nrow(m))
  record("ELISA percent 5-mC", "largest per-animal difference (percentage points)", 0,
         signif(max(abs(m$percent_5mC - m$percent)), 3), tolerance = 0.001,
         note = "legacy sheet used constants rounded to 0.224 and 0.6105")
  aov_new <- fread(out("03_global_methylation", "anova_by_scale.csv"))[scale == "arcsine" & grepl("^ELISA", assay)]
  aov_old <- read.table(archive_path("stats", "AOV_5-mC.csv"), header = TRUE, check.names = FALSE)
  rownames(aov_old) <- trimws(rownames(aov_old))
  for (tm in c("ploidy", "treatment", "ploidy:treatment")) {
    record("ELISA two-way ANOVA (arcsine scale)", paste("F,", tm), round(aov_old[tm, "F value"], 2),
           round(aov_new[term == tm, F_sequential], 2), tolerance = 0.005,
           note = "regenerated from raw absorbance with unrounded constants")
  }
  mc_new <- fread(out("03_global_methylation", "anova_by_scale.csv"))[scale == "arcsine" & grepl("before", assay)]
  mc_old <- read.table(archive_path("stats", "AOV_%mCpG.csv"), header = TRUE, check.names = FALSE)
  rownames(mc_old) <- trimws(rownames(mc_old))
  record("Percent mCpG two-way ANOVA (arcsine scale)", "F, ploidy:treatment, counts before deduplication",
         round(mc_old["ploidy:treatment", "F value"], 2),
         round(mc_new[term == "ploidy:treatment", F_sequential], 2), tolerance = 0.005,
         note = "legacy sheet mixed two alignment runs (issue I-13); exact agreement is not expected")
} else record("ELISA percent 5-mC", "step 03 outputs", NA, NA)

# ---- Step 04: CpG counts by coverage ---------------------------------------------------
if (have("04_cpg_coverage", "cpg_counts_by_coverage.csv")) {
  new <- fread(out("04_cpg_coverage", "cpg_counts_by_coverage.csv"))[mode == "legacy"]
  old <- as.data.table(read_excel(archive_path("analyses", "coverage_summary.xlsx"), sheet = "total_cpg"))
  old[, sample_set := fifelse(ploidy == "both", "all", ploidy)]
  m <- merge(old, new, by.x = c("sample_set", "coverage"), by.y = c("sample_set", "min_coverage"))
  for (i in seq_len(nrow(m))) {
    record("CpGs covered in all libraries of a set", sprintf("%s, %dx", m$sample_set[i], m$coverage[i]),
           m$cpg[i], m$cpg_in_all_libraries[i])
  }
} else record("CpGs covered in all libraries of a set", "step 04 outputs", NA, NA)

# ---- Step 05: DML lists ---------------------------------------------------------------
if (have("05_differential_methylation", "dml_ploidy_legacy.csv")) {
  compare_dml <- function(contrast, legacy_file, label) {
    new <- fread(out("05_differential_methylation", sprintf("dml_%s_legacy.csv", contrast)))
    old <- fread(archive_path("tables", legacy_file)); setnames(old, 1, "legacy_row")
    record(label, "number of DML", nrow(old), nrow(new))
    record(label, "hyper", sum(old$meth.diff > 0), sum(new$meth.diff > 0))
    record(label, "hypo", sum(old$meth.diff < 0), sum(new$meth.diff < 0))
    m <- merge(old, new, by = c("chr", "start"), suffixes = c("_old", "_new"))
    record(label, "loci present in both lists", nrow(old), nrow(m))
    if (nrow(m) > 0) {
      record(label, "largest difference in meth.diff (percentage points)", 0,
             signif(max(abs(m$meth.diff_old - m$meth.diff_new)), 3), tolerance = 1e-6)
      record(label, "largest relative difference in q-value", 0,
             signif(max(abs(m$qvalue_old - m$qvalue_new) / m$qvalue_old), 3), tolerance = 1e-4,
             note = "SLIM q-values agree to five significant digits across methylKit versions")
      if ("united_row" %in% names(m)) record(label, "loci at the same row of the united table", nrow(m),
                                             sum(m$legacy_row == m$united_row))
    }
  }
  compare_dml("ploidy", "DML-getMethylDiff-PLOIDY-Cov10-20.csv", "Ploidy DML (20%, q < 0.01)")
  compare_dml("heat", "DML-getMethylDiff-HEAT-Cov10-20.csv", "Heat DML (20%, q < 0.01)")
  sh_new <- fread(out("05_differential_methylation", "dml_shared_legacy.csv"))
  sh_old <- fread(archive_path("tables", "DML-getMethylDiff-INTERACTION-Cov10-20.csv"))
  record("Loci in both lists (legacy 'interaction')", "number", nrow(sh_old), nrow(sh_new))
  record("Loci in both lists (legacy 'interaction')", "present in both versions", nrow(sh_old),
         nrow(merge(sh_old, sh_new, by = c("chr", "start"))))
  record("Loci in both lists (legacy 'interaction')", "both hyper (manuscript value)", 28, sum(sh_new$class == "both_hyper"))
  record("Loci in both lists (legacy 'interaction')", "both hypo (manuscript value)", 49, sum(sh_new$class == "both_hypo"))
  record("Loci in both lists (legacy 'interaction')", "opposite sign (manuscript value)", 123, sum(sh_new$class == "opposite"))

  sw_new <- fread(out("05_differential_methylation", "dml_counts_by_threshold_legacy.csv"))
  sw_old <- as.data.table(read_excel(archive_path("analyses", "DML_difference_threshold_comparison.xlsx"), sheet = "plot"))
  sw_old[, set := c(desiccation = "heat", ploidy = "ploidy", interaction = "shared")[trt]]
  m <- merge(sw_old, sw_new, by = c("set", "difference"))
  for (i in seq_len(nrow(m))) {
    record("DML count by difference threshold", sprintf("%s, %d%%", m$set[i], m$difference[i]), m$DML[i], m$n[i])
    record("DML count by difference threshold", sprintf("%s, %d%%, hyper", m$set[i], m$difference[i]), m$hyper[i], m$n_hyper[i],
           note = if (m$set[i] == "heat" && m$difference[i] == 50)
             "legacy workbook has hyper and hypo swapped at this threshold (1 and 3); the legacy script comments say 3 and 1" else "")
  }
} else record("DML lists", "step 05 outputs", NA, NA)

res <- rbindlist(results)
fwrite(res, here("06_diagnostics", "legacy_check_results.csv"))
cat("\nLegacy reproduction checks:\n"); print(res[, .N, by = .(check, status)])
if (any(res$status == "differs")) { cat("\nDifferences:\n"); print(res[status == "differs"]) }

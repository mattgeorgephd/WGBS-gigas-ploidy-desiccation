# compare_analyses.R
#
# One table comparing every analysis that has been run: the paper's
# ("legacy") and each optional analysis under 04_output/not_adopted/<id>/,
# at the paper's 20-point and, where run, the 50-point difference threshold.
#
# Run:     Rscript 01_pipeline/optional/compare_analyses.R   (from the repository root)
# Needs:   04_output/not_adopted/comparisons/analyses_included.csv (written by every run)
#          04_output/not_adopted/comparisons/05_differential_methylation/tables/dml_counts_by_mode.csv
# Writes:  04_output/not_adopted/comparisons/analysis_comparison.csv
#          04_output/not_adopted/comparisons/analysis_comparison.md
#
# Reads only tables the pipeline wrote (and rewrites analyses_included.csv from the
# analysis.csv stamps, as the steps do); computes nothing new from the data. The DML rule is
# SLIM q < 0.01 and an absolute difference strictly above the threshold (methylKit
# getMethylDiff). Venkataraman et al. 2022 (BMC Genomics 23:556,
# doi:10.1186/s12864-022-08781-5) used 5x in each sample, the 99.9th percentile filter,
# MN overdispersion with the chi-squared test, C-to-T variant loci removed, and 50 points.

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))

# Rebuild analyses_included.csv from the stamps on disk first: steps 04, 05, 08 and 10
# write it, step 09 does not, so after a run that ended with step 09 the file can lag
# its last stamp (found by the final error check of 2026-11-21).
write_analyses_included()

cmp_dir <- here("04_output", "not_adopted", "comparisons")
inc <- fread(file.path(cmp_dir, "analyses_included.csv"), colClasses = "character")
cnt <- fread(file.path(cmp_dir, "05_differential_methylation", "tables", "dml_counts_by_mode.csv"))
setnames(cnt, "mode", "id")

if (!"dml_difference" %in% names(inc)) inc[, dml_difference := NA_character_]
inc[is.na(dml_difference) | dml_difference == "", dml_difference := "20"]

x <- merge(inc[, .(id, filter, min_coverage, overdispersion, test, ct_excluded, dml_difference, fit_id)],
           cnt, by = "id", all = TRUE)
missing_counts <- x[is.na(cpg_analysed), id]
missing_stamp <- x[is.na(filter), id]
if (length(missing_stamp)) stop("DML counts without a stamp in analyses_included.csv: ", paste(missing_stamp, collapse = ", "))
x <- x[!is.na(cpg_analysed)]   # e.g. an analysis whose step 05 has not run

x[, `:=`(
  upper_filter       = ifelse(filter == "corrected", "yes", "no (as the draft)"),
  overdispersion_lab = ifelse(overdispersion == "MN", "MN, chi-squared", "none (as the draft)"),
  ct_variants        = ifelse(toupper(ct_excluded) == "TRUE", "excluded", "kept"),
  difference         = as.integer(dml_difference)
)]
x[, venkataraman_2022_settings := upper_filter == "yes" & min_coverage == "5" & overdispersion == "MN" &
                                   ct_variants == "excluded" & difference == 50L]

paper <- x[id == "legacy"]
if (nrow(paper) != 1) stop("the paper's row (legacy) is missing")
x[, `:=`(
  ploidy_per_1000_cpg   = round(1000 * ploidy_dml / cpg_analysed, 2),
  heat_per_1000_cpg     = round(1000 * heat_dml / cpg_analysed, 2),
  shared_opposite_pct   = ifelse(shared > 0, round(100 * shared_opposite / shared, 1), NA_real_),
  ploidy_dml_vs_paper   = round(ploidy_dml / paper$ploidy_dml, 3),
  heat_dml_vs_paper     = round(heat_dml / paper$heat_dml, 3)
)]

# Fixed order: by threshold, then coverage, then the options in token order.
x[, ord := (difference == 50) * 1000 + (min_coverage == "5") * 100 + (upper_filter == "yes") * 10 +
           (overdispersion == "MN") * 2 + (ct_variants == "excluded") * 1]
setorder(x, ord, id)

out <- x[, .(analysis = id, adopted = ifelse(id == "legacy", "yes (the paper)", "no"),
             min_coverage, upper_filter, overdispersion = overdispersion_lab, ct_variants,
             difference_threshold = difference, qvalue_threshold = 0.01, fit_id,
             cpg_tested = cpg_analysed,
             ploidy_dml, ploidy_hyper, ploidy_hypo, heat_dml, heat_hyper, heat_hypo,
             shared, shared_both_hyper, shared_both_hypo, shared_opposite, shared_opposite_pct,
             ploidy_per_1000_cpg, heat_per_1000_cpg, ploidy_dml_vs_paper, heat_dml_vs_paper,
             venkataraman_2022_settings)]
fwrite(out, file.path(cmp_dir, "analysis_comparison.csv"))

fmt <- function(v) ifelse(is.na(v), "", format(v, big.mark = ",", scientific = FALSE, trim = TRUE))
md <- c(
  "# Comparison of all analyses",
  "",
  "Written by `01_pipeline/optional/compare_analyses.R` from `analyses_included.csv` and `05_differential_methylation/tables/dml_counts_by_mode.csv` in this folder; do not edit by hand. Only `legacy` is the paper's analysis; every other row is an optional analysis, not adopted. A DML has SLIM q < 0.01 and an absolute difference of more than the threshold (percentage points). \"Upper filter\" is the 99.9th percentile coverage filter that the draft describes but did not apply (I-05). Shared: loci in both the ploidy and the heat list; opposite: different directions in the two lists. The row marked V2022 has the settings that the Methods of Venkataraman et al. 2022 state for coverage, the test, variant removal and the threshold (the variant rule counts reads with their full CIGAR, I-08).",
  "",
  "| Analysis | Coverage | Upper filter | Overdispersion | C-to-T loci | Threshold | CpGs tested | Ploidy DML (hyper/hypo) | Heat DML (hyper/hypo) | Shared (both hyper/both hypo/opposite) | Opposite % |",
  "|---|---|---|---|---|---|---|---|---|---|---|",
  out[, sprintf("| %s%s | %sx | %s | %s | %s | %d | %s | %s (%s/%s) | %s (%s/%s) | %s (%s/%s/%s) | %s |",
                analysis, ifelse(venkataraman_2022_settings, " (V2022)", ifelse(analysis == "legacy", " (paper)", "")),
                min_coverage, sub(" \\(as the draft\\)", "", upper_filter), sub(" \\(as the draft\\)", "", overdispersion),
                ct_variants, difference_threshold, fmt(cpg_tested),
                fmt(ploidy_dml), fmt(ploidy_hyper), fmt(ploidy_hypo), fmt(heat_dml), fmt(heat_hyper), fmt(heat_hypo),
                fmt(shared), fmt(shared_both_hyper), fmt(shared_both_hypo), fmt(shared_opposite), fmt(shared_opposite_pct))],
  "",
  if (length(missing_counts)) paste0("Stamped but without DML counts (step 05 not run): ", paste(missing_counts, collapse = ", "), ".") else NULL
)
writeLines(md, file.path(cmp_dir, "analysis_comparison.md"))
cat(sprintf("%d analyses written to %s\n", nrow(out), file.path(cmp_dir, "analysis_comparison.csv")))

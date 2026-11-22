# compare_analyses.R
#
# One table comparing every analysis that has been run: the paper's
# ("legacy") and each optional analysis under 04_output/not_adopted/<id>/,
# at the paper's 20-point and, where run, the 50-point difference threshold, and
# at the paper's q-value cutoff of 0.01 and, where run, 0.05 (option q05), with the
# paper's contrasts and, where run, the control-only ploidy contrast (option ctrl) or the
# tested interaction (option inter, whose interaction list is counted in four columns).
#
# Run:     Rscript 01_pipeline/optional/compare_analyses.R   (from the repository root)
# Needs:   04_output/not_adopted/comparisons/analyses_included.csv (written by every run)
#          04_output/not_adopted/comparisons/05_differential_methylation/tables/dml_counts_by_mode.csv
# Writes:  04_output/not_adopted/comparisons/analysis_comparison.csv
#          04_output/not_adopted/comparisons/analysis_comparison.md
#
# Reads only tables the pipeline wrote (and rewrites analyses_included.csv from the
# analysis.csv stamps, as the steps do); computes nothing new from the data. The DML rule is
# SLIM q below the cutoff of the analysis (0.01, or 0.05 with option q05) and an absolute
# difference strictly above the threshold (methylKit getMethylDiff). Both come from the
# stamps. Venkataraman et al. 2022 (BMC Genomics 23:556,
# doi:10.1186/s12864-022-08781-5) used 5x in each sample, the 99.9th percentile filter,
# MN overdispersion with the chi-squared test, C-to-T variant loci removed, and 50 points at
# q < 0.01.

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
if (!"dml_qvalue" %in% names(inc)) inc[, dml_qvalue := NA_character_]
inc[is.na(dml_qvalue) | dml_qvalue == "", dml_qvalue := "0.01"]
if (!"contrast_variant" %in% names(inc)) inc[, contrast_variant := NA_character_]
inc[is.na(contrast_variant) | contrast_variant == "", contrast_variant := "none"]
# The interaction counts are in dml_counts_by_mode.csv of the analyses with option inter only
# (empty for the others); without any such analysis the columns are absent.
for (col in c("interaction_dml", "interaction_hyper", "interaction_hypo", "interaction_in_shared")) {
  if (!col %in% names(cnt)) cnt[, (col) := NA_integer_]
}

x <- merge(inc[, .(id, filter, min_coverage, overdispersion, test, contrast_variant, ct_excluded, dml_difference, dml_qvalue, fit_id)],
           cnt, by = "id", all = TRUE)
missing_counts <- x[is.na(cpg_analysed), id]
missing_stamp <- x[is.na(filter), id]
if (length(missing_stamp)) stop("DML counts without a stamp in analyses_included.csv: ", paste(missing_stamp, collapse = ", "))
x <- x[!is.na(cpg_analysed)]   # e.g. an analysis whose step 05 has not run

x[, `:=`(
  upper_filter       = ifelse(filter == "corrected", "yes", "no (as the draft)"),
  overdispersion_lab = ifelse(overdispersion == "MN", "MN, chi-squared", "none (as the draft)"),
  ct_variants        = ifelse(toupper(ct_excluded) == "TRUE", "excluded", "kept"),
  contrasts          = ifelse(contrast_variant == "ctrl", "ploidy in controls only (ctrl)",
                              ifelse(contrast_variant == "inter", "paper's and interaction (inter)", "paper's")),
  difference         = as.integer(dml_difference),
  qvalue             = as.numeric(dml_qvalue)
)]
x[, venkataraman_2022_settings := upper_filter == "yes" & min_coverage == "5" & overdispersion == "MN" &
                                   ct_variants == "excluded" & difference == 50L & qvalue == 0.01 & contrast_variant == "none"]

paper <- x[id == "legacy"]
if (nrow(paper) != 1) stop("the paper's row (legacy) is missing")
x[, `:=`(
  ploidy_per_1000_cpg   = round(1000 * ploidy_dml / cpg_analysed, 2),
  heat_per_1000_cpg     = round(1000 * heat_dml / cpg_analysed, 2),
  shared_opposite_pct   = ifelse(shared > 0, round(100 * shared_opposite / shared, 1), NA_real_),
  ploidy_dml_vs_paper   = round(ploidy_dml / paper$ploidy_dml, 3),
  heat_dml_vs_paper     = round(heat_dml / paper$heat_dml, 3)
)]

# Fixed order: by the contrasts (the paper's, then ctrl, then inter), then q-value cutoff, then
# threshold, then coverage, then the options in token order. The paper's row is first.
x[, ord := (contrast_variant == "ctrl") * 100000 + (contrast_variant == "inter") * 200000 +
           (qvalue == 0.05) * 10000 + (difference == 50) * 1000 + (min_coverage == "5") * 100 + (upper_filter == "yes") * 10 +
           (overdispersion == "MN") * 2 + (ct_variants == "excluded") * 1]
setorder(x, ord, id)

out <- x[, .(analysis = id, adopted = ifelse(id == "legacy", "yes (the paper)", "no"),
             min_coverage, upper_filter, overdispersion = overdispersion_lab, contrasts, ct_variants,
             difference_threshold = difference, qvalue_threshold = qvalue, fit_id,
             cpg_tested = cpg_analysed,
             ploidy_dml, ploidy_hyper, ploidy_hypo, heat_dml, heat_hyper, heat_hypo,
             shared, shared_both_hyper, shared_both_hypo, shared_opposite, shared_opposite_pct,
             interaction_dml, interaction_hyper, interaction_hypo, interaction_in_shared,
             ploidy_per_1000_cpg, heat_per_1000_cpg, ploidy_dml_vs_paper, heat_dml_vs_paper,
             venkataraman_2022_settings)]
fwrite(out, file.path(cmp_dir, "analysis_comparison.csv"))

fmt <- function(v) ifelse(is.na(v), "", format(v, big.mark = ",", scientific = FALSE, trim = TRUE))
md <- c(
  "# Comparison of all analyses",
  "",
  "Written by `01_pipeline/optional/compare_analyses.R` from `analyses_included.csv` and `05_differential_methylation/tables/dml_counts_by_mode.csv` in this folder; do not edit by hand. Only `legacy` is the paper's analysis; every other row is an optional analysis, not adopted. A DML has a SLIM q-value below the cutoff of its row (column Q: 0.01 in the paper, 0.05 with option `q05`) and an absolute difference of more than the threshold of its row (percentage points). \"Upper filter\" is the 99.9th percentile coverage filter that the draft describes but did not apply (I-05). Shared: loci in both the ploidy and the heat list; opposite: different directions in the two lists. Contrasts: the paper's (ploidy with heat as covariate and heat with ploidy as covariate, all 22 libraries), `ctrl` (the ploidy contrast on the 12 control libraries alone; I-07), or `inter` (the paper's, and a tested ploidy by heat interaction whose DML list is in the column Interaction: its DML (hyper/hypo) and how many of them are among the shared loci; I-02). The row marked V2022 has the settings that the Methods of Venkataraman et al. 2022 state for coverage, the test, variant removal and the threshold, at q < 0.01 (the variant rule counts reads with their full CIGAR, I-08).",
  "",
  "| Analysis | Coverage | Upper filter | Overdispersion | Contrasts | C-to-T loci | Threshold | Q | CpGs tested | Ploidy DML (hyper/hypo) | Heat DML (hyper/hypo) | Shared (both hyper/both hypo/opposite) | Opposite % | Interaction DML (hyper/hypo; in shared) |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
  out[, sprintf("| %s%s | %sx | %s | %s | %s | %s | %d | %g | %s | %s (%s/%s) | %s (%s/%s) | %s (%s/%s/%s) | %s | %s |",
                analysis, ifelse(venkataraman_2022_settings, " (V2022)", ifelse(analysis == "legacy", " (paper)", "")),
                min_coverage, sub(" \\(as the draft\\)", "", upper_filter), sub(" \\(as the draft\\)", "", overdispersion),
                contrasts, ct_variants, difference_threshold, qvalue_threshold, fmt(cpg_tested),
                fmt(ploidy_dml), fmt(ploidy_hyper), fmt(ploidy_hypo), fmt(heat_dml), fmt(heat_hyper), fmt(heat_hypo),
                fmt(shared), fmt(shared_both_hyper), fmt(shared_both_hypo), fmt(shared_opposite), fmt(shared_opposite_pct),
                ifelse(is.na(interaction_dml), "",
                       sprintf("%s (%s/%s; %s)", fmt(interaction_dml), fmt(interaction_hyper), fmt(interaction_hypo), fmt(interaction_in_shared))))],
  "",
  if (length(missing_counts)) paste0("Stamped but without DML counts (step 05 not run): ", paste(missing_counts, collapse = ", "), ".") else NULL
)
writeLines(md, file.path(cmp_dir, "analysis_comparison.md"))
cat(sprintf("%d analyses written to %s\n", nrow(out), file.path(cmp_dir, "analysis_comparison.csv")))

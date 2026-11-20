# table1_group_means.R
#
# Purpose   What the "methylation difference" values of the manuscript's Table 1
#           are made of. For each of the 20 loci of Table 1 this script computes
#           the read-pooled percent methylation of the four groups (diploid
#           control, diploid heat, triploid control, triploid heat) and from them
#             - the two differences the legacy analysis reported (triploid minus
#               diploid over all 22 libraries; heat cohort minus control cohort
#               over all 22 libraries), to show which of them Table 1 prints;
#             - the ploidy difference within each cohort (issue I-07: values
#               described as "under control conditions" were computed across
#               both cohorts);
#             - the cohort difference within each ploidy (issue I-02: the
#               "interaction" values are the pooled heat contrast and describe
#               neither ploidy).
#           It also writes the per-library percent methylation at each locus, so
#           that the number of animals behind a group difference can be seen.
#           Everything here is descriptive: groups hold 5 or 6 libraries and no
#           test is made.
#
# Run       Rscript 06_diagnostics/table1_group_means.R      (under a minute)
#
# Reads     04_output/09_go_enrichment_legacy/tables/table_1_legacy.csv
#           04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz
#           03_input_files/metadata/wgbs_libraries.csv
# Writes    06_diagnostics/table1_group_means.csv        one row per row of Table 1
#           06_diagnostics/table1_library_values.csv     one row per locus and library
#
# Notes     - "Read-pooled" is the definition methylKit uses for meth.diff: the
#             methylated reads of all libraries of a group divided by all their
#             reads, after coverage normalization. The counts are the normalized
#             counts that step 05 tested.
#           - Legacy mode only, because Table 1 comes from the legacy analysis.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))

lib <- read_libraries()
table1 <- fread(diagnostic_input("09_go_enrichment_legacy", "tables", "table_1_legacy.csv"), colClasses = c(gene_id = "character"))
united <- fread(diagnostic_input("05_differential_methylation", "tables", "united_cpg_legacy.csv.gz"))
stopifnot(sum(grepl("^numCs", names(united))) == nrow(lib))

# Percent methylation of a set of libraries at one CpG, reads pooled.
pooled <- function(row, idx) {
  cs <- sum(unlist(row[, paste0("numCs", idx), with = FALSE]))
  ts <- sum(unlist(row[, paste0("numTs", idx), with = FALSE]))
  100 * cs / (cs + ts)
}
groups <- levels(lib$group)                      # diploid control, diploid heat, triploid control, triploid heat
idx <- lapply(setNames(groups, groups), function(g) which(lib$group == g))

per_locus <- list(); per_library <- list()
for (i in seq_len(nrow(table1))) {
  row <- united[chr == table1$chromosome[i] & start == table1$start[i]]
  stopifnot(nrow(row) == 1)
  g <- vapply(idx, function(x) pooled(row, x), numeric(1))
  ploidy_all <- pooled(row, which(lib$ploidy == "triploid")) - pooled(row, which(lib$ploidy == "diploid"))
  heat_all   <- pooled(row, which(lib$treatment == "heat")) - pooled(row, which(lib$treatment == "control"))
  reported <- table1$methylation_difference[i]
  per_locus[[i]] <- data.table(
    table1_row = i, dml_set = table1$dml_set[i], gene_id = table1$gene_id[i], gene_name = table1$gene_name[i],
    sequence = table1$chromosome[i], start = table1$start[i], table1_value = reported,
    pct_diploid_control = g[["diploid control"]], pct_diploid_heat = g[["diploid heat"]],
    pct_triploid_control = g[["triploid control"]], pct_triploid_heat = g[["triploid heat"]],
    triploid_minus_diploid_all_libraries = ploidy_all, heat_minus_control_all_libraries = heat_all,
    table1_value_is = fifelse(abs(reported - ploidy_all) < 0.01, "triploid minus diploid, all libraries",
                      fifelse(abs(reported - heat_all) < 0.01, "heat cohort minus control cohort, all libraries", "neither")),
    triploid_minus_diploid_control_cohort = g[["triploid control"]] - g[["diploid control"]],
    triploid_minus_diploid_heat_cohort = g[["triploid heat"]] - g[["diploid heat"]],
    heat_minus_control_diploids = g[["diploid heat"]] - g[["diploid control"]],
    heat_minus_control_triploids = g[["triploid heat"]] - g[["triploid control"]])
  cs <- unlist(row[, paste0("numCs", seq_len(nrow(lib))), with = FALSE])
  ts <- unlist(row[, paste0("numTs", seq_len(nrow(lib))), with = FALSE])
  per_library[[i]] <- data.table(table1_row = i, dml_set = table1$dml_set[i], gene_name = table1$gene_name[i],
                                 sequence = table1$chromosome[i], start = table1$start[i],
                                 library_id = lib$library_id, group = as.character(lib$group),
                                 normalized_reads = cs + ts, pct_methylation = 100 * cs / (cs + ts))
}
per_locus <- rbindlist(per_locus); per_library <- rbindlist(per_library)
num <- names(per_locus)[vapply(per_locus, is.double, logical(1))]
per_locus[, (num) := lapply(.SD, round, 2), .SDcols = num]
per_library[, pct_methylation := round(pct_methylation, 1)]
fwrite(per_locus, here("06_diagnostics", "table1_group_means.csv"))
fwrite(per_library, here("06_diagnostics", "table1_library_values.csv"))

# ---- Summary for the log ---------------------------------------------------------------
cat("What each Table 1 value is:\n"); print(per_locus[, .N, by = .(dml_set, table1_value_is)])
ploidy_rows <- per_locus[dml_set == "ploidy-DML"]
cat(sprintf("\nPloidy rows whose difference exceeds %g points in the control cohort alone: %d of %d\n",
            config$dml_difference, sum(abs(ploidy_rows$triploid_minus_diploid_control_cohort) > config$dml_difference), nrow(ploidy_rows)))
cat("\nInteraction rows, heat cohort minus control cohort within each ploidy:\n")
print(per_locus[dml_set == "interaction-DML", .(gene_name, table1_value, heat_minus_control_diploids, heat_minus_control_triploids)])

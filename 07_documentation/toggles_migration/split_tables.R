# split_tables.R
#
# One-time migration of the committed outputs to the settings layout of
# 2026-11-20 (07_documentation/agent_workflows/scope_inputs/settings_build_spec.md,
# section 7.2, step 3). Run once from the repository root, after the renames
# (git mv) of step 2, with the BASE commit as argument:
#   Rscript 07_documentation/toggles_migration/split_tables.R 66c7121
#
# 1. The 8 tables that listed the legacy and the corrected mode side by side
#    are split: the paper's copy keeps the header and the legacy lines, the copy
#    of the analysis 'corrected' (04_output/not_adopted/corrected/) the header and
#    the corrected lines, and the full table goes to 04_output/not_adopted/comparisons/.
#    Lines are copied as text, so every line keeps its bytes.
# 2. Step 11: the five gene_level_ columns leave named_terms_results_3_4.csv and
#    the gene-level rows leave semantic_summary.csv and semantic_terms_dropped.csv;
#    they go to the _gene_level files of 04_output/not_adopted/legacy/11_semantic_similarity/.
# 3. Stamps (analysis.csv) of the folders not_adopted/legacy/ and not_adopted/corrected/.
# Reads the BASE versions with git show; writes only the files named above.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
base <- args[1]
suppressPackageStartupMessages(library(data.table))
git_lines <- function(path) {
  x <- system2("git", c("show", paste0(base, ":", path)), stdout = TRUE)
  stopifnot(length(x) >= 1)
  x
}
put <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, "wb"); writeLines(lines, con, sep = "\n"); close(con)
}
spanning <- c("04_cpg_coverage/tables/cpg_counts_by_coverage.csv", "04_cpg_coverage/tables/cpg_counts_by_feature_track.csv",
              "04_cpg_coverage/tables/methylation_distribution_by_library.csv", "04_cpg_coverage/tables/methylation_distribution_by_group.csv",
              "04_cpg_coverage/tables/library_correlation_summary.csv", "04_cpg_coverage/tables/pca_variance.csv",
              "05_differential_methylation/tables/dml_counts_by_mode.csv", "08_gene_annotation/tables/annotation_summary.csv")
for (f in spanning) {
  l <- git_lines(file.path("04_output", f))
  mode <- fread(text = paste0(paste(l, collapse = "\n"), "\n"), colClasses = "character")$mode
  stopifnot(length(mode) == length(l) - 1, all(mode %in% c("legacy", "corrected")))
  put(c(l[1], l[-1][mode == "legacy"]), file.path("04_output", f))
  put(c(l[1], l[-1][mode == "corrected"]), file.path("04_output/not_adopted/corrected", f))
  put(l, file.path("04_output/not_adopted/comparisons", f))
  cat(sprintf("%-60s %5d lines: %5d legacy, %5d corrected\n", f, length(l) - 1, sum(mode == "legacy"), sum(mode == "corrected")))
}

# Step 11
s11 <- "04_output/11_semantic_similarity/tables"
g11 <- "04_output/not_adopted/legacy/11_semantic_similarity/tables"
gene_cols <- c("gene_level_genes_annotated", "gene_level_genes_with_dml", "gene_level_p_wallenius", "gene_level_q_wallenius", "gene_level_result")
named <- fread(text = paste(git_lines(file.path(s11, "named_terms_results_3_4.csv")), collapse = "\n"), colClasses = "character", na.strings = NULL)
dir.create(g11, recursive = TRUE, showWarnings = FALSE)
fwrite(named[, setdiff(names(named), gene_cols), with = FALSE], file.path(s11, "named_terms_results_3_4.csv"), eol = "\n", na = "")
fwrite(named[, c("set", "phrase", "go_id", gene_cols), with = FALSE], file.path(g11, "named_terms_results_3_4_gene_level.csv"), eol = "\n", na = "")
for (f in c("semantic_summary.csv", "semantic_terms_dropped.csv")) {
  l <- git_lines(file.path(s11, f))
  tl <- fread(text = paste0(paste(l, collapse = "\n"), "\n"), colClasses = "character")$term_list
  stopifnot(length(tl) == length(l) - 1)
  gene <- startsWith(tl, "gene level")
  put(c(l[1], l[-1][!gene]), file.path(s11, f))
  put(c(l[1], l[-1][gene]), file.path(g11, sub("\\.csv$", "_gene_level.csv", f)))
  cat(sprintf("%-60s %5d lines: %5d paper, %5d gene level\n", f, length(l) - 1, sum(!gene), sum(gene)))
}

# Step 10: the tables that listed both modes are in comparisons/ (moved with git mv);
# each analysis folder gets its own copy, as step 10 now writes one per analysis.
s10 <- "04_output/10_go_enrichment_gene_level/tables"
for (f in c("go_gene_level_summary.csv", "dml_probability_by_cpg_count.csv", "go_enriched_terms.csv", "go_enriched_term_genes.csv",
            "go_terms_by_analysis.csv")) {
  l <- git_lines(file.path(s10, f))
  d <- fread(text = paste0(paste(l, collapse = "\n"), "\n"), colClasses = "character")
  stopifnot(nrow(d) == length(l) - 1)
  mode <- if ("mode" %in% names(d)) d$mode else ifelse(d$analysis == "Same, DML called with the coverage filter", "corrected", "legacy")
  for (m in c("legacy", "corrected")) put(c(l[1], l[-1][mode == m]), file.path("04_output/not_adopted", m, "10_go_enrichment_gene_level/tables", f))
  cat(sprintf("%-60s %5d lines: %5d legacy, %5d corrected\n", f, length(l) - 1, sum(mode == "legacy"), sum(mode == "corrected")))
}
top <- fread(text = paste(git_lines(file.path(s10, "legacy_top10_at_gene_level.csv")), collapse = "\n"))
corrected_cols <- grep("_eligible_corrected_dml$", names(top), value = TRUE)
fwrite(top[, setdiff(names(top), corrected_cols), with = FALSE],
       "04_output/not_adopted/legacy/10_go_enrichment_gene_level/tables/legacy_top10_at_gene_level.csv", eol = "\n")
corr <- top[, setdiff(names(top), c(grep("_(tested|eligible)$", names(top), value = TRUE), "gene_level_status")), with = FALSE]
corr[, gene_level_status := fifelse(top$gene_level_status == "not in the current ontology", "not in the current ontology",
                            fifelse(is.na(genes_annotated_eligible_corrected_dml), "no eligible gene carries it",
                            fifelse(genes_annotated_eligible_corrected_dml < 5, "fewer annotated genes than go_min_genes",
                            fifelse(q_wallenius_eligible_corrected_dml < 0.05, "enriched", "not enriched"))))]
fwrite(corr, "04_output/not_adopted/corrected/10_go_enrichment_gene_level/tables/legacy_top10_at_gene_level.csv", eol = "\n")
cat("legacy_top10_at_gene_level.csv: one copy per analysis folder\n")

# Stamps
stamp <- function(id, filter, steps) {
  data.table(id = id, filter = filter, min_coverage = "10", overdispersion = "none", test = "not passed", ct_excluded = "FALSE",
             fit_id = id, fit_fingerprint_md5 = "", ct_variant_list_md5 = "", git_commit = base,
             written_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), steps_done = steps)
}
fwrite(stamp("corrected", "corrected", "04;05;06;07;08;10"), "04_output/not_adopted/corrected/analysis.csv", eol = "\n")
fwrite(stamp("legacy", "legacy", "10;11"), "04_output/not_adopted/legacy/analysis.csv", eol = "\n")
cat("stamps written\n")

# compare_downstream.R
#
# How the downstream results of the paper change across the 192 analyses of the
# comparison table (author's request of 2026-11-25): steps 07 (genome features), 08
# (genes with DML) and 09 (the paper's GO procedure, reproduced), which every analysis
# has run. Not adopted: a description of the optional analyses. Reads only tables
# the pipeline wrote; computes nothing new from the data.
#
#   A. Genome features (step 07). Every odds ratio of feature_association_tests_<id>.csv
#      side by side, and per background, set and track: the paper's odds ratio and how
#      many analyses with at least `min_dml` DML in the set find the same direction at
#      q < 0.05, the opposite direction at q < 0.05, or neither.
#   B. Genes with DML (step 08). Per analysis and set: genes with a DML, how many of
#      them are genes of the paper's list of that set, and the share of the paper's genes
#      kept. Per set: the genes that hold a DML in the most analyses.
#   C. The loci of the draft's Table 1 (step 09, table_1_legacy.csv): in how many analyses
#      each locus is a DML of its set (by sequence and start).
#   D. The paper's GO procedure (step 09; issues I-03, I-21, I-39, I-40, I-43, I-48): per
#      analysis, set and ontology, how many of the paper's ten top-ranked terms are among
#      the analysis's ten top-ranked terms. The procedure is the draft's and is not a valid
#      enrichment test; this measures only how stable its output is across settings.
#
# The shared set ("interaction-DML" of the draft) is the overlap of the ploidy and heat
# lists, as in the draft; with option ctrl the ploidy list is the control-only list. The
# tested interaction list (option inter) is annotated by step 08 and tested by step 10
# (2026-11-25) but is not part of these tables.
#
# Run:     Rscript 01_pipeline/optional/compare_downstream.R   (from the repository root,
#          after compare_analyses.R)
# Reads:   04_output/not_adopted/comparisons/analysis_comparison.csv
#          per analysis (the paper under 04_output/<step>/, the others under
#          04_output/not_adopted/<id>/<step>/): 07 feature_association_tests_<id>.csv,
#          08 genes_with_dml_<set>_<id>.csv, 05 dml_<set>_<id>.csv,
#          09 legacy_fisher_<set>_<ontology>.csv and (the paper) table_1_legacy.csv
# Writes:  04_output/not_adopted/comparisons/downstream/tables/
#            feature_enrichment_by_analysis.csv, feature_enrichment_summary.csv,
#            genes_by_analysis.csv, recurring_genes.csv,
#            table1_loci_by_analysis.csv, table1_loci_summary.csv,
#            legacy_go_top10_by_analysis.csv, legacy_go_top10_summary.csv

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))

min_dml <- 30          # part A: sets with fewer DML are listed but not counted in the summary
sets <- c("ploidy", "heat", "shared")
ontologies <- c("BP", "CC", "MF")
out_dir <- here("04_output", "not_adopted", "comparisons", "downstream", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x, name) write_table(x, file.path(out_dir, name))

cmp <- fread(here("04_output", "not_adopted", "comparisons", "analysis_comparison.csv"))
if (nrow(cmp) != 192) stop("Expected the 192 analyses of analysis_comparison.csv; found ", nrow(cmp), call. = FALSE)
settings <- cmp[, .(analysis, contrasts, min_coverage, upper_filter, overdispersion, ct_variants, difference_threshold, qvalue_threshold,
                    ploidy_dml, heat_dml, shared)]
tab <- function(id, step, file) {
  f <- if (id == "legacy") here("04_output", step, "tables", file) else here("04_output", "not_adopted", id, step, "tables", file)
  if (!file.exists(f)) stop("Missing ", f, call. = FALSE)
  f
}

# ---- A. Genome features ----------------------------------------------------------------------
fe <- rbindlist(lapply(cmp$analysis, function(id) {
  x <- fread(tab(id, "07_genome_features", sprintf("feature_association_tests_%s.csv", id)))
  x[, analysis := id][, mode := NULL]
}))
fe[, set := fcase(set == "ploidy DML", "ploidy", set == "heat DML", "heat", set == "shared", "shared")]
if (anyNA(fe$set)) stop("Unexpected set label in feature_association_tests.", call. = FALSE)
out(merge(settings, fe, by = "analysis"), "feature_enrichment_by_analysis.csv")
paper_fe <- fe[analysis == "legacy", .(background, set, track, paper_odds_ratio = odds_ratio, paper_q = q_BH)]
fes <- merge(fe[analysis != "legacy" & dml_total >= min_dml], paper_fe, by = c("background", "set", "track"))
fe_sum <- fes[, .(analyses_with_enough_dml = .N,
                  same_direction_q05 = sum(q_BH < 0.05 & sign(log(odds_ratio)) == sign(log(paper_odds_ratio)), na.rm = TRUE),
                  opposite_direction_q05 = sum(q_BH < 0.05 & sign(log(odds_ratio)) == -sign(log(paper_odds_ratio)), na.rm = TRUE),
                  median_odds_ratio = stats::median(odds_ratio, na.rm = TRUE),
                  min_odds_ratio = min(odds_ratio, na.rm = TRUE), max_odds_ratio = max(odds_ratio, na.rm = TRUE)),
              by = .(background, set, track, paper_odds_ratio, paper_q)]
fe_sum[, neither := analyses_with_enough_dml - same_direction_q05 - opposite_direction_q05]
fe_sum[, min_dml_in_set := min_dml]
setorder(fe_sum, background, set, track)
out(fe_sum, "feature_enrichment_summary.csv")

# ---- B. Genes with DML ------------------------------------------------------------------------
genes <- rbindlist(lapply(cmp$analysis, function(id) rbindlist(lapply(sets, function(s) {
  g <- fread(tab(id, "08_gene_annotation", sprintf("genes_with_dml_%s_%s.csv", s, id)), select = c("gene_id", "gene", "biotype", "product", "n_dml"))
  if (nrow(g)) g[, `:=`(analysis = id, set = s)] else data.table(gene_id = integer(), gene = character(), biotype = character(),
                                                                      product = character(), n_dml = integer(), analysis = character(), set = character())
}))))
paper_genes <- genes[analysis == "legacy", .(gene_id, set)]
gb <- genes[, .(n_genes = .N, n_protein_coding = sum(biotype == "protein_coding")), by = .(analysis, set)]
gb <- merge(CJ(analysis = cmp$analysis, set = sets), gb, by = c("analysis", "set"), all.x = TRUE)
gb[is.na(n_genes), `:=`(n_genes = 0L, n_protein_coding = 0L)]
ov <- merge(genes[, .(analysis, set, gene_id)], paper_genes[, .(set, gene_id, in_paper = TRUE)], by = c("set", "gene_id"), all.x = TRUE)
ov <- ov[, .(n_genes_in_paper_list = sum(!is.na(in_paper))), by = .(analysis, set)]
gb <- merge(gb, ov, by = c("analysis", "set"), all.x = TRUE)
gb[is.na(n_genes_in_paper_list), n_genes_in_paper_list := 0L]
gb <- merge(gb, paper_genes[, .(paper_genes = .N), by = set], by = "set")
gb[, `:=`(share_of_paper_genes_kept = n_genes_in_paper_list / paper_genes,
          share_of_genes_in_paper_list = fifelse(n_genes > 0, n_genes_in_paper_list / n_genes, NA_real_))]
out(merge(settings, gb, by = "analysis")[order(set, -n_genes)], "genes_by_analysis.csv")
rec <- genes[, .(analyses = uniqueN(analysis), analyses_corrected_test = uniqueN(analysis[grepl("mnchisq", analysis)]),
                 gene = gene[1], biotype = biotype[1], product = product[1]), by = .(set, gene_id)]
rec <- merge(rec, paper_genes[, .(set, gene_id, in_paper_list = TRUE)], by = c("set", "gene_id"), all.x = TRUE)
rec[is.na(in_paper_list), in_paper_list := FALSE]
rec <- rec[order(set, -analyses, -analyses_corrected_test)][, head(.SD, 50), by = set]
out(rec, "recurring_genes.csv")

# ---- C. The loci of Table 1 -------------------------------------------------------------------
t1 <- fread(here("04_output", "09_go_enrichment_legacy", "tables", "table_1_legacy.csv"))
t1[, set := fcase(dml_set == "ploidy-DML", "ploidy", dml_set == "heat-DML", "heat", dml_set == "interaction-DML", "shared")]
t1 <- unique(t1[, .(set, gene_name, chromosome, start)])
dml_keys <- rbindlist(lapply(cmp$analysis, function(id) rbindlist(lapply(sets, function(s) {
  d <- fread(tab(id, "05_differential_methylation", sprintf("dml_%s_%s.csv", s, id)), select = c("chr", "start"))
  if (!nrow(d)) return(NULL)
  d[, .(analysis = id, set = s, chromosome = chr, start)]
}))))
tl <- CJ(analysis = cmp$analysis, row = seq_len(nrow(t1)))
tl <- cbind(tl[, .(analysis)], t1[tl$row])
tl <- merge(tl, unique(dml_keys)[, in_list := TRUE], by = c("analysis", "set", "chromosome", "start"), all.x = TRUE)
tl[is.na(in_list), in_list := FALSE]
if (!all(tl[analysis == "legacy", in_list])) stop("A Table 1 locus is not a DML of the paper's list of its set.", call. = FALSE)
out(merge(settings[, .(analysis, contrasts, min_coverage, upper_filter, overdispersion, ct_variants, difference_threshold, qvalue_threshold)], tl, by = "analysis"),
    "table1_loci_by_analysis.csv")
tls <- tl[analysis != "legacy", .(analyses_keeping = sum(in_list), analyses_corrected_test_keeping = sum(in_list & grepl("mnchisq", analysis)),
                                  analyses_corrected_test = sum(grepl("mnchisq", analysis)), analyses = .N),
          by = .(set, gene_name, chromosome, start)]
out(tls[order(set, -analyses_keeping)], "table1_loci_summary.csv")

# ---- D. The paper's GO procedure ----------------------------------------------------------------
top10 <- function(id, s, o) {
  x <- fread(tab(id, "09_go_enrichment_legacy", sprintf("legacy_fisher_%s_%s.csv", s, o)))
  if (!nrow(x) || !"top10" %in% names(x)) return(character())
  x[top10 == TRUE, go_id]
}
paper_top <- setNames(lapply(sets, function(s) setNames(lapply(ontologies, function(o) top10("legacy", s, o)), ontologies)), sets)
go <- rbindlist(lapply(cmp$analysis[cmp$analysis != "legacy"], function(id) rbindlist(lapply(sets, function(s) rbindlist(lapply(ontologies, function(o) {
  t <- top10(id, s, o)
  data.table(analysis = id, set = s, ontology = o, n_top_terms = length(t), paper_top_terms = length(paper_top[[s]][[o]]),
             paper_top_terms_kept = length(intersect(t, paper_top[[s]][[o]])))
}))))))
go <- merge(settings[, .(analysis, contrasts, min_coverage, upper_filter, overdispersion, ct_variants, difference_threshold, qvalue_threshold,
                         n_dml = NA_integer_)], go, by = "analysis")
go[, n_dml := fcase(set == "ploidy", cmp$ploidy_dml[match(analysis, cmp$analysis)], set == "heat", cmp$heat_dml[match(analysis, cmp$analysis)],
                    set == "shared", cmp$shared[match(analysis, cmp$analysis)])]
out(go, "legacy_go_top10_by_analysis.csv")
gos <- go[, .(analyses = .N, with_any_top_terms = sum(n_top_terms > 0),
              kept_median = stats::median(as.double(paper_top_terms_kept)), kept_max = max(paper_top_terms_kept),
              kept_at_least_5 = sum(paper_top_terms_kept >= 5), kept_none = sum(paper_top_terms_kept == 0)), by = .(set, ontology, overdispersion)]
out(gos[order(set, ontology, overdispersion)], "legacy_go_top10_summary.csv")

cat("written:", paste(list.files(out_dir), collapse = ", "), "\n")

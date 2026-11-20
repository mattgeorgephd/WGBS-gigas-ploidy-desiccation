# 11_semantic_similarity.R
#
# Purpose   Semantic similarity of GO terms: the computation behind Figures 6, 7
#           and 8 of the draft, redone in R so that it can be traced, and the same
#           procedure applied to the terms the gene-level analysis calls enriched.
#
#           The draft made the three figures on the REVIGO web server (Supek et
#           al. 2011, PLoS ONE 6:e21800, doi:10.1371/journal.pone.0021800). A web
#           server cannot be a pipeline step. This step uses rrvgo (Sayols 2023,
#           MicroPublication Biology, doi:10.17912/micropub.biology.000811), which
#           follows the same approach in three parts:
#             1. Similarity. rrvgo::calculateSimMatrix gives the pairwise similarity
#                of the terms with GOSemSim (Yu et al. 2010, Bioinformatics 26:976,
#                doi:10.1093/bioinformatics/btq064). The measure is Rel (Schlicker
#                et al. 2006, BMC Bioinformatics 7:302, doi:10.1186/1471-2105-7-302),
#                the SimRel of the Methods. It needs the information content of
#                each term, taken from the annotation of one organism
#                (config$semantic_similarity$organism_db).
#             2. Clusters. rrvgo::reduceSimMatrix clusters 1 - similarity by
#                complete linkage, cuts the tree at config$semantic_similarity$threshold,
#                and names each cluster after its highest-scoring term, the
#                representative (ties: see Notes).
#             3. Layout. Classical multidimensional scaling of 1 - similarity
#                (stats::cmdscale; Torgerson 1952, Psychometrika 17:401; Gower 1966,
#                Biometrika 53:325). It has no random start, unlike UMAP. Each axis
#                is oriented so that its most extreme term has a positive coordinate.
#
#           Three kinds of term list are clustered:
#             legacy      terms with Fisher p < config$go_pvalue in the legacy
#                         enrichment (step 09), per DML set and ontology. This is
#                         the list the Methods and the figure captions describe.
#                         Score: -log10 p.
#             as drawn    every BP term annotated to a transcript of a gene that
#                         holds a DML, as the legacy reader returned the annotation
#                         table. This is the list the draft's figures were drawn
#                         from (see Notes). It has no p-values; the score is the
#                         number of genes of the organism database annotated to the
#                         term or its children (rrvgo scores = "size"), so the
#                         broadest term names a cluster.
#             gene level  terms with q < config$go_fdr in step 10. At present there
#                         are none; the step then writes a one-row table that says
#                         so and draws nothing.
#
# Run       Rscript 01_pipeline/11_semantic_similarity.R
#           (after steps 08 and 09; about 5 minutes and 1.6 GB of memory. The
#           gene-level part runs instead, in an optional run with --gene-level-go,
#           after step 10; about a minute)
#
# Inputs    04_output/09_go_enrichment_legacy/tables/legacy_fisher_<set>_<ontology>.csv, legacy_reader_loss.csv
#           04_output/08_gene_annotation/tables/dml_go_annotation_<set>_legacy.tsv.gz, annotation_summary.csv
#           gene-level part: go_enriched_terms.csv and go_gene_level_<set>_<id>.csv of step 10 of the
#             analysis (04_output/not_adopted/<id>/10_go_enrichment_gene_level/), and the paper's
#             named_terms_results_3_4.csv
#           03_input_files/annotation/legacy/revigo_exports/code_Revigo.R, code_OLD_Revigo.R, Revigo-PLOIDY-BP.R
#             (plot scripts exported by REVIGO, data embedded; read as text, never run)
#           03_input_files/annotation/legacy/fisher_other_analysis/ploidy-<ontology>-FisherTestResults.csv
#             (Fisher tables of another analysis; comparison only)
#           GO.db and the organism database of the R installation
# Outputs   04_output/11_semantic_similarity/tables/
#             semantic_clusters_legacy_<set>_<ontology>.csv    one row per term: p, score, cluster, representative, coordinates
#             semantic_clusters_as_drawn_<set>_BP.csv          the same for the lists the draft's figures were drawn from
#             similarity_<list>_<set>_<ontology>.csv           the similarity matrices (.csv.gz when larger than gzip_above_mb)
#             semantic_summary.csv                             per list: terms in, terms dropped and why, clusters, checks
#             semantic_terms_dropped.csv                       the dropped terms, with the reason
#             term_lists_compared.csv                          size and overlap of the legacy and the as-drawn list
#             revigo_exports.csv                               each REVIGO export against each candidate term list
#             named_terms_results_3_4.csv                      every GO term named in Results 3.4, traced
#           04_output/11_semantic_similarity/figures/
#             figure_6_legacy, figure_7_legacy, figure_8_legacy   panel A: the list as drawn; panel B: the list the Methods describe
#             figure_revigo_export_<name>                          each distinct REVIGO export redrawn from its embedded data
#           04_output/11_semantic_similarity/session_info.txt
#           Gene-level part (NOT ADOPTED; only in an optional run with --gene-level-go, after
#           step 10), in 04_output/not_adopted/<id>/11_semantic_similarity/ (<id> = legacy for
#           the paper's analysis):
#             tables/semantic_clusters_gene_level.csv          the clusters of the gene-level terms of step 10;
#                                                              one row with a note when there are none
#             tables/named_terms_results_3_4_gene_level.csv    the terms named in Results 3.4 (set, phrase, go_id)
#                                                              at gene level: five gene_level_ columns
#             tables/semantic_summary_gene_level.csv, semantic_terms_dropped_gene_level.csv
#             figures/figure_semantic_gene_level_<set>_<ontology>_<id>   only when terms are enriched
#
# Manuscript items: Figures 6, 7 and 8; Results 3.4 (the three paragraphs on semantic
#           clustering); Methods 2.6, last sentence.
# Notes     - Run kinds. In the paper run the step redraws figures of the draft
#             from the legacy enrichment of the paper's DML (sections 1 to 5 and
#             7). In an optional run with --gene-level-go it does the gene-level
#             part only (sections 6 and 7) for the analysis of the run, writing
#             under 04_output/not_adopted/; in any other optional run it ends at
#             once with status 0. When the summary of step 08 holds no row of the
#             paper for a DML set, the step stops with an error, as step 09 does.
#           - Legacy files. The REVIGO exports and the three Fisher tables are
#             copies of files in 08_archive/ (code/Revigo.R, code/OLD/Revigo.R,
#             figures/revigo_figure/Revigo-PLOIDY-BP.R, tables/ploidy-<ontology>-
#             FisherTestResults.csv), kept under 03_input_files/ so that no step
#             reads the archive. The two files named code_ are the same export.
#             revigo_exports.csv names each file by its path in the repository.
#           - Settings. Every setting named in this header and in the comments
#             below (threshold, n_labels, order_check_n, gzip_above_mb,
#             cluster_palette, label_rings and the others) is an element of the
#             list config$semantic_similarity in run_config.R, where each one
#             is explained. The code reads the list as cfg. The checks of
#             run_config.R stop the step when n_labels is larger than the
#             number of colours in cluster_palette.
#           - Which terms the draft plotted. The REVIGO export behind Figure 6
#             (Revigo-PLOIDY-BP.R) holds no p-values, and its terms are not the
#             terms below the p-value cutoff: every one of them is a BP term
#             annotated to a transcript of a gene with a ploidy DML
#             (revigo_exports.csv). For Figures 7 and 8 no export survives; the
#             terms that Results 3.4 names for them point the same way
#             (named_terms_results_3_4.csv). Each figure is therefore regenerated
#             twice: as drawn (panel A) and as described (panel B).
#           - Labels are chosen by rule, not by hand (issue I-22).
#               Which clusters are labelled. In each panel, the
#               config$semantic_similarity$n_labels clusters whose representative
#               has the highest score. Ties go to the larger cluster, then to the
#               lower GO ID. Clusters are numbered in that order, so the labelled
#               clusters are 1 to n_labels. The label is the name of the
#               representative.
#               Which term is the representative. The term with the highest score.
#               When several terms of a cluster share that score: the one least
#               similar to the other terms of the cluster (the highest uniqueness
#               within the cluster, column uniqueness_in_cluster). When that is
#               shared too: the higher GO ID. This is the order rrvgo applies. The
#               step applies the rule itself and stops if rrvgo names a cluster
#               after a term the rule excludes.
#               Where a label is drawn. place_labels() tries a fixed set of
#               positions around the circle of the cluster and takes the nearest
#               one where the label touches no other label and no labelled circle;
#               a position over an unlabelled circle counts as farther away.
#               ggrepel then draws the labels and the leader lines, with a fixed
#               seed, and moves a label only if it still overlaps something.
#           - What the threshold means. rrvgo cuts the complete-linkage tree of
#             1 - similarity at height `threshold`, so any two terms of a cluster
#             are at least 1 - threshold similar. REVIGO's "allowed similarity"
#             works the other way: a term is removed when it is more similar than
#             that value to a term that stays. semantic_summary.csv gives the
#             number of clusters under both readings. The step stops if rrvgo's
#             clusters are not the cut described here.
#           - The clusters depend on the order of the terms. GOSemSim rounds
#             similarities to three decimals, so equal distances are common, and
#             stats::hclust breaks such ties by the position of the terms. The step
#             sorts the terms by GO ID, which fixes the result; the order of the
#             input rows does not matter. Another order of the same terms can give
#             a slightly different partition. The step measures how different:
#             each list is clustered again in order_check_n random orders, and
#             semantic_summary.csv gives the number of orders with another
#             partition, the lowest agreement with the reported partition
#             (adjusted Rand index; Hubert and Arabie 1985, Journal of
#             Classification 2:193, doi:10.1007/BF01908075), the range of the
#             number of clusters, and the fewest labelled representatives that
#             stay labelled.
#           - Terms are dropped when the GO.db of the R installation does not hold
#             them (obsolete, or merged into another term; the legacy lists come
#             from the 2021 ontology), when it files them under another ontology,
#             when GOSemSim has no entry for them, or when no gene of the organism
#             database is annotated to them or to a child, which leaves the
#             information content undefined. The organism database is human, so
#             terms used only for other taxa are lost. REVIGO carries a merged term
#             under its current ID; this step drops it and counts it.
#           - Colour. Each labelled cluster has its own colour: its terms, the
#             lines from them to the representative, its circle, the leader line
#             of its label and the tint of the label. All other clusters are grey;
#             a list has dozens of clusters, and that many colours cannot be told
#             apart. The colours (cluster_palette) were chosen so that every pair
#             of them stays apart under normal colour vision and under simulated
#             protanopia and deuteranopia (Machado et al. 2009, IEEE Transactions
#             on Visualization and Computer Graphics 15:1291,
#             doi:10.1109/TVCG.2009.113), and so that each has a contrast of at
#             least 3 to 1 against white. Colours are assigned within a panel:
#             the same colour in two panels does not mark the same cluster.
#           - A similarity matrix is written as plain CSV. One larger than
#             gzip_above_mb is written gzipped, and its file name ends in .csv.gz.
#           - The legacy enrichment is unsound (issues I-03, I-39, I-40) and step 10
#             calls no term enriched at gene level at its configured settings
#             (I-41). The legacy figures are regenerated for traceability, not as
#             results.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))

step <- "11_semantic_similarity"
t0 <- start_step(step, "semantic similarity of GO terms (Figures 6 to 8)")
# The paper run does the paper's sections; an optional run with --gene-level-go the gene-level part only (see Notes).
paper_part <- identical(config$run_kind, "paper")
gene_part  <- identical(config$run_kind, "optional") && isTRUE(config$addons$gene_level_go)
if (!paper_part) require_run_kind("optional_gene_level", "step 11 in an optional run does the gene-level part of step 10's terms")
an <- config$analysis
if (gene_part) check_stamp(an, needs = c("04", "05", "08", "10"))

# Rcpp is in the list because GOSemSim calls compiled code of Rcpp without loading it.
load_packages(c("data.table", "GO.db", "Rcpp", "GOSemSim", "rrvgo", config$semantic_similarity$organism_db, "ggplot2", "ggrepel", "patchwork"))
if (as.character(packageVersion("GO.db")) == config$legacy_go_db$version) {
  stop("GO.db ", config$legacy_go_db$version, " (the legacy release) is loaded. Run this step in a fresh R process.")
}

# ---- Settings ------------------------------------------------------------------------
# All settings of this step are in run_config.R, list config$semantic_similarity:
# the similarity measure and the threshold, the facts about the draft and the
# legacy files, the seeds, and the sizes and costs that lay out the figures.
cfg  <- config$semantic_similarity
sets <- c(ploidy = "ploidy-DML", heat = "heat-DML", shared = "interaction-DML")   # names used in the draft
onts <- c(BP = "P", CC = "C", MF = "F")                                           # ontology and its code in the annotation tables
tab  <- function(name) if (paper_part) output_path(step, "tables", name) else analysis_path(an, step, "tables", name)
fig  <- function(name) if (paper_part) output_path(step, "figures", name) else analysis_path(an, step, "figures", name)
go_db_version <- as.character(packageVersion("GO.db"))
log_msg("GO.db %s (ontology of %s), %s %s, rrvgo %s, GOSemSim %s; measure %s, threshold %g", go_db_version, go_source_date(),
        cfg$organism_db, as.character(packageVersion(cfg$organism_db)), as.character(packageVersion("rrvgo")),
        as.character(packageVersion("GOSemSim")), cfg$measure, cfg$threshold)

# Every input file, checked before any work starts.
repo_relative      <- function(path) substring(path, nchar(here()) + 2L)          # a path as it reads from the repository root
legacy_fisher_file <- function(s_i, o_i) paper_path("09_go_enrichment_legacy", "tables", sprintf("legacy_fisher_%s_%s.csv", s_i, o_i))
reader_loss_file   <- paper_path("09_go_enrichment_legacy", "tables", "legacy_reader_loss.csv")
annotation_file    <- function(s_i) paper_path("08_gene_annotation", "tables", sprintf("dml_go_annotation_%s_legacy.tsv.gz", s_i))
gene_level_file    <- function(s_i) analysis_input(an, "10_go_enrichment_gene_level", "tables", sprintf("go_gene_level_%s_%s.csv", s_i, an$id))
enriched_file      <- analysis_input(an, "10_go_enrichment_gene_level", "tables", "go_enriched_terms.csv")
named_terms_file   <- paper_path(step, "tables", "named_terms_results_3_4.csv")
export_file        <- function(e_i) input_path(cfg$revigo_exports[[e_i]])
other_fisher_file  <- function(o_i) input_path(sprintf(cfg$other_fisher, o_i))
input_files <- if (paper_part) {
  c(unlist(lapply(names(sets), function(s_i) vapply(names(onts), function(o_i) legacy_fisher_file(s_i, o_i), ""))), reader_loss_file,
    vapply(names(sets), annotation_file, ""), vapply(names(cfg$revigo_exports), export_file, ""), vapply(names(onts), other_fisher_file, ""))
} else c(vapply(names(sets), gene_level_file, ""), enriched_file, named_terms_file)
if (!all(file.exists(input_files))) {
  stop("Input file(s) not found:\n  ", paste(repo_relative(input_files[!file.exists(input_files)]), collapse = "\n  "),
       if (paper_part) "\nSteps 08 and 09 write the files under 04_output; run them first." else
         "\nThe paper run (step 11) and step 10 of the analysis write them; run those first.", call. = FALSE)
}
if (paper_part) {
  # Annotation tables of step 08 without a row of the paper per DML set stop the step here, as they stop step 09.
  invisible(legacy_annotation_summary(names(sets)))
}
if (gene_part) {
  enriched <- fread(enriched_file)               # the gene-level terms (section 6); read here so that a faulty table stops the step at once
  need <- c("mode", "set", "ontology", "go_id", "term", "p_wallenius", "q_wallenius")
  if (!all(need %in% names(enriched))) stop("go_enriched_terms.csv lacks column(s): ", paste(setdiff(need, names(enriched)), collapse = ", "))
  if (anyDuplicated(enriched[, .(mode, set, ontology, go_id)])) stop("A GO ID occurs twice for one mode, set and ontology in go_enriched_terms.csv")
  # The gene-level figures and matrices exist only when terms are enriched; remove those of an earlier run.
  unlink(list.files(analysis_dir(an, step, "figures"), pattern = "^figure_semantic_gene_level_", full.names = TRUE))
  unlink(list.files(analysis_dir(an, step, "tables"), pattern = "^similarity_gene_level_", full.names = TRUE))
}

# =============================================================================
# 1. Term lists
# =============================================================================
# The ontology of the R installation: terms, and secondary IDs (IDs merged into another term).
go_terms <- as.data.table(suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db),
                                                                 columns = c("TERM", "ONTOLOGY"), keytype = "GOID")))
setnames(go_terms, c("go_id", "term", "ontology"))
secondary_ids <- unique(as.data.table(AnnotationDbi::toTable(GO.db::GOSYNONYM))[, .(secondary_id = synonym, primary_id = go_id)])
obsolete_ids <- AnnotationDbi::keys(GO.db::GOOBSOLETE)
primary_id <- function(ids) {            # the current ID of a term given under a secondary ID; other IDs unchanged
  hit <- secondary_ids$primary_id[match(ids, secondary_ids$secondary_id)]
  out <- as.character(ids); out[!is.na(hit)] <- hit[!is.na(hit)]
  out
}

if (paper_part) {   # sections 1 to 5 (their work, not their functions) belong to the paper run
# Legacy enrichment: every tested term of the nine Fisher tables.
legacy <- rbindlist(lapply(names(sets), function(s_i) rbindlist(lapply(names(onts), function(o_i) {
  x <- fread(legacy_fisher_file(s_i, o_i))
  need <- c("set", "ontology", "rank", "go_id", "term", "annotated", "significant", "p_fisher", "genes_significant", "below_p_cutoff")
  if (!all(need %in% names(x))) stop("legacy_fisher_", s_i, "_", o_i, ".csv lacks column(s): ", paste(setdiff(need, names(x)), collapse = ", "))
  if (!identical(x$below_p_cutoff, x$p_fisher < config$go_pvalue)) {
    stop("legacy_fisher_", s_i, "_", o_i, ".csv was written with another p-value cutoff than config$go_pvalue. Rerun step 09.")
  }
  x[, ..need]
}))))
if (anyDuplicated(legacy[, .(set, ontology, go_id)])) stop("A GO ID occurs twice in a legacy Fisher table")
legacy_names <- unique(legacy[, .(go_id, ontology, term_legacy_release = term)])   # names in the release the legacy analysis ran on
if (anyDuplicated(legacy_names$go_id)) stop("A GO ID has two names or two ontologies in the legacy Fisher tables")

# The lists the draft plotted: GO IDs of the DML annotation tables, by ontology code, read as the legacy script
# read them (issue I-40). That reader ran some lines together; such a row keeps its own GO ID but carries the
# ontology code of a later line, so a list can hold a few terms of another ontology. Term names are taken from
# a correct read of the same file.
annotated_lists <- rbindlist(lapply(names(sets), function(s_i) {
  a <- suppressWarnings(read_dml_annotation(annotation_file(s_i), "legacy"))
  a <- a[GOcat %chin% onts & grepl("^GO:[0-9]{7}$", GOID), .(annotation_rows = .N), by = .(ontology = names(onts)[match(GOcat, onts)], go_id = GOID)]
  true_names <- unique(read_dml_annotation(annotation_file(s_i), "exact")[, .(go_id = GOID, term_annotation = gsub('^"|"$', "", GOterm))])
  a[, term_annotation := true_names$term_annotation[match(go_id, true_names$go_id)]]
  cbind(data.table(set = rep(s_i, nrow(a))), a)
}))
if (anyDuplicated(annotated_lists[, .(set, ontology, go_id)])) stop("A GO ID occurs twice in an annotated list")
reader_loss <- fread(reader_loss_file)[grepl("legacy", reader)]
for (s_i in names(sets)) for (o_i in names(onts)) {
  n_here <- annotated_lists[set == s_i & ontology == o_i, .N]; n_step09 <- reader_loss[set == s_i][[paste0("go_ids_", o_i)]]
  if (length(n_step09) != 1 || n_here != n_step09) {
    stop(sprintf("Annotated %s terms of the %s set: %d here, %s in legacy_reader_loss.csv of step 09", o_i, s_i, n_here,
                 paste(n_step09, collapse = "/")))
  }
}

lists_compared <- rbindlist(lapply(names(sets), function(s_i) rbindlist(lapply(names(onts), function(o_i) {
  p_list <- legacy[set == s_i & ontology == o_i & below_p_cutoff == TRUE, go_id]; a_list <- annotated_lists[set == s_i & ontology == o_i, go_id]
  data.table(set = s_i, ontology = o_i, terms_tested_legacy = legacy[set == s_i & ontology == o_i, .N], terms_legacy_p_below_cutoff = length(p_list),
             terms_annotated_legacy_reader = length(a_list), terms_in_both = length(intersect(p_list, a_list)),
             annotated_terms_tested_legacy = sum(a_list %chin% legacy[set == s_i & ontology == o_i, go_id]))
}))))
lists_compared[, `:=`(percent_of_annotated_with_p_below_cutoff = round(100 * terms_in_both / terms_annotated_legacy_reader, 1),
                      draft_figure = fifelse(ontology == cfg$figure_ontology, paste("Figure", cfg$figure_numbers[set]), ""))]
write_table(lists_compared, tab("term_lists_compared.csv"))
}   # end of the paper part of section 1

# =============================================================================
# 2. Similarity, clusters, layout
# =============================================================================
semantic_cache <- new.env()
semantic_data <- function(ont_code) {      # information content of the terms of one ontology, computed once
  if (is.null(semantic_cache[[ont_code]])) {
    semantic_cache[[ont_code]] <- suppressMessages(GOSemSim::godata(cfg$organism_db, ont = ont_code, keytype = "ENTREZID"))
  }
  semantic_cache[[ont_code]]
}
# Genes of the organism database annotated to a term or its children, counted as rrvgo counts them.
term_gene_count <- function(ids) {
  x <- suppressMessages(AnnotationDbi::select(getExportedValue(cfg$organism_db, cfg$organism_db), keytype = "GOALL", keys = ids,
                                              columns = "ENTREZID"))
  vapply(ids, function(g) uniqueN(x$ENTREZID[x$GOALL == g]), integer(1))
}
drop_reasons <- c(obsolete = "obsolete in GO.db", merged = "merged into another term in GO.db (secondary ID)", unknown_id = "not in GO.db",
                  other_ontology = "filed under another ontology in GO.db", not_in_gosemsim = "no entry in the term table of GOSemSim",
                  no_gene = "no gene of the organism database annotated to it or to a child")
cluster_columns <- c("go_id", "term", "p_value", "score", "cluster", "cluster_size", "is_representative", "parent_go_id", "parent_term",
                     "labelled", "term_size", "uniqueness", "uniqueness_in_cluster", "mds_1", "mds_2")

# Name and number the clusters of one partition by the rule in the Notes.
#   ids         GO IDs of the terms (row names of `sim`, in any order)
#   membership  the cluster of each term
#   scores      the score of each term
# Returns `clusters` (one row per cluster: representative, its score, number of terms, and the cluster number, 1 being the
# cluster whose representative scores highest), `apart` (uniqueness of each term within its cluster, as rrvgo computes it:
# 1 minus the mean similarity to the terms of the cluster, itself counted as 0) and `eligible` (the terms the rule allows as
# representative before the GO ID decides).
name_clusters <- function(ids, membership, scores, sim) {
  d <- data.table(term_id = ids, at = match(ids, rownames(sim)), member_of = membership, term_score = scores)
  d[, apart := { m <- sim[at, at, drop = FALSE]; diag(m) <- 0; 1 - rowMeans(m) }, by = member_of]
  d[, eligible := term_score == max(term_score), by = member_of]
  d[eligible == TRUE, eligible := apart >= max(apart) - cfg$tie_tolerance, by = member_of]
  named <- d[eligible == TRUE, .(representative = max(term_id), top_score = term_score[1]), by = member_of]
  named <- merge(named, d[, .(terms = .N), by = member_of], by = "member_of")
  if (nrow(named) != uniqueN(membership)) stop("A cluster has no representative or more than one")
  setorder(named, -top_score, -terms, representative)
  named[, number := seq_len(.N)]
  list(clusters = named, apart = d$apart, eligible = d[eligible == TRUE, term_id])
}
# Agreement of two partitions of the same terms: adjusted Rand index (Hubert and Arabie 1985). 1 means the same partition.
adjusted_rand <- function(part_a, part_b) {
  both <- table(part_a, part_b); pairs_in <- function(v) sum(choose(v, 2))
  expected <- pairs_in(rowSums(both)) * pairs_in(colSums(both)) / choose(sum(both), 2)
  highest <- (pairs_in(rowSums(both)) + pairs_in(colSums(both))) / 2
  if (highest == expected) 1 else (pairs_in(both) - expected) / (highest - expected)
}
# The same terms clustered in random orders (see Notes). `reference` is the reported cluster of each term, `scores` its score,
# both in the row order of `sim`; `labelled_reference` holds the representatives of the labelled clusters.
order_check <- function(sim, scores, reference, labelled_reference) {
  ids <- rownames(sim)
  set.seed(cfg$order_check_seed)
  one_order <- function(k) {
    shuffled <- sample(length(ids))
    part <- stats::cutree(stats::hclust(stats::as.dist(1 - sim[shuffled, shuffled])), h = cfg$threshold)[ids]
    named <- name_clusters(ids, part, scores, sim)$clusters
    c(adjusted_rand(part, reference), nrow(named), sum(named[number <= cfg$n_labels, representative] %chin% labelled_reference))
  }
  found <- vapply(seq_len(cfg$order_check_n), one_order, numeric(3))      # rows: index, clusters, labelled representatives kept
  list(tried = cfg$order_check_n, other = sum(found[1, ] < 1 - 1e-12), rand = min(found[1, ]), fewest = min(found[2, ]), most = max(found[2, ]),
       kept = min(found[3, ]))
}

# Cluster one list of GO terms.
#   ids       GO IDs of the list
#   p_values  p-values named by GO ID, or NULL for a list without p-values
#   ont_code  "BP", "CC" or "MF"
# Returns the table of placed terms (cluster_columns), the dropped terms with the reason, the similarity matrix,
# and a one-row summary.
cluster_term_list <- function(ids, p_values, ont_code) {
  ids <- sort(unique(ids))
  ic <- semantic_data(ont_code)@IC
  filed <- go_terms$ontology[match(ids, go_terms$go_id)]
  absent <- fifelse(ids %chin% secondary_ids$secondary_id, drop_reasons[["merged"]],
            fifelse(ids %chin% obsolete_ids, drop_reasons[["obsolete"]], drop_reasons[["unknown_id"]]))
  reason <- fifelse(is.na(filed), absent,
            fifelse(filed != ont_code, drop_reasons[["other_ontology"]],
            fifelse(!ids %chin% names(ic), drop_reasons[["not_in_gosemsim"]],
            fifelse(!is.finite(ic[ids]), drop_reasons[["no_gene"]], "kept"))))
  kept <- ids[reason == "kept"]; n_kept <- length(kept)
  dropped <- data.table(go_id = ids[reason != "kept"], reason = reason[reason != "kept"])
  score_is <- if (is.null(p_values)) "genes annotated in the organism database" else "-log10 p"
  if (!is.null(p_values) && (anyNA(p_values[kept]) || any(p_values[kept] < 0))) stop("A term of the list has no valid p-value")
  score_of <- function(g, sizes) if (is.null(p_values)) as.numeric(sizes) else -log10(pmax(unname(p_values[g]), .Machine$double.xmin))
  checks <- list(clusters = 0L, singletons = 0L, largest = 0L, min_within = NA_real_, complement = 0L, fit = NA_real_)
  orders <- list(tried = 0L, other = NA_integer_, rand = NA_real_, fewest = NA_integer_, most = NA_integer_, kept = NA_integer_)
  sim <- NULL

  if (n_kept == 0) {
    out <- data.table(go_id = character(0), term = character(0), p_value = numeric(0), score = numeric(0), cluster = integer(0),
                      cluster_size = integer(0), is_representative = logical(0), parent_go_id = character(0), parent_term = character(0),
                      labelled = logical(0), term_size = integer(0), uniqueness = numeric(0), uniqueness_in_cluster = numeric(0),
                      mds_1 = numeric(0), mds_2 = numeric(0))
  } else if (n_kept == 1) {                       # nothing to compare: one cluster of one term
    size_1 <- term_gene_count(kept)
    out <- data.table(go_id = kept, term = go_terms$term[match(kept, go_terms$go_id)], score = score_of(kept, size_1), cluster = 1L,
                      cluster_size = 1L, parent_go_id = kept, term_size = size_1, uniqueness = 1, uniqueness_in_cluster = 1, mds_1 = 0, mds_2 = 0)
    checks <- list(clusters = 1L, singletons = 1L, largest = 1L, min_within = NA_real_, complement = 1L, fit = NA_real_)
  } else {
    sim <- withCallingHandlers(
      suppressMessages(rrvgo::calculateSimMatrix(kept, orgdb = cfg$organism_db, ont = ont_code, method = cfg$measure,
                                                 semdata = semantic_data(ont_code))),
      warning = function(w) stop("rrvgo::calculateSimMatrix warned: ", conditionMessage(w), call. = FALSE))
    if (!is.matrix(sim) || !setequal(rownames(sim), kept) || nrow(sim) != n_kept) {
      stop("rrvgo::calculateSimMatrix kept other terms than this script expected (", ont_code, "). Check the rrvgo and GOSemSim versions.")
    }
    sim <- sim[kept, kept]
    if (anyNA(sim) || !isSymmetric(sim)) stop("The similarity matrix has missing values or is not symmetric")
    scores <- if (is.null(p_values)) "size" else setNames(score_of(kept, NULL), kept)
    red <- as.data.table(suppressMessages(rrvgo::reduceSimMatrix(sim, scores, threshold = cfg$threshold, orgdb = cfg$organism_db)))
    if (nrow(red) != n_kept || !setequal(red$go, kept)) stop("rrvgo::reduceSimMatrix returned other terms than it was given")

    # rrvgo's clusters must be the complete-linkage cut that the header describes.
    tree <- stats::hclust(stats::as.dist(1 - sim))
    own_cut <- stats::cutree(tree, h = cfg$threshold)
    if (uniqueN(own_cut) != uniqueN(red$cluster) || uniqueN(data.table(a = own_cut[red$go], b = red$cluster)) != uniqueN(red$cluster)) {
      stop("rrvgo::reduceSimMatrix did not return the complete-linkage clusters cut at the threshold. Check the rrvgo version.")
    }
    members <- split(red$go, red$cluster)
    within <- vapply(members, function(g) if (length(g) > 1) min(sim[g, g][upper.tri(diag(length(g)))]) else NA_real_, numeric(1))

    # The representative and the number of each cluster, by the rule in the Notes. rrvgo's own choice must be a term the rule allows.
    named <- name_clusters(red$go, red$cluster, red$score, sim)
    if (!all(unique(red$parent) %chin% named$eligible)) {
      stop("rrvgo::reduceSimMatrix named a cluster after a term that the rule in the header excludes. Check the rrvgo version.")
    }
    in_named <- match(red$cluster, named$clusters$member_of)

    # Layout: classical multidimensional scaling; each axis oriented by its most extreme term.
    mds <- stats::cmdscale(stats::as.dist(1 - sim), k = min(2L, n_kept - 1L), eig = TRUE)
    pts <- cbind(mds$points, 0, 0)[, 1:2, drop = FALSE]; rownames(pts) <- kept
    for (j in 1:2) if (any(pts[, j] != 0) && pts[which.max(abs(pts[, j])), j] < 0) pts[, j] <- -pts[, j]

    out <- red[, .(go_id = go, term, score, term_size = as.integer(size), uniqueness = termUniqueness)]
    out[, `:=`(cluster = named$clusters$number[in_named], cluster_size = named$clusters$terms[in_named],
               parent_go_id = named$clusters$representative[in_named], uniqueness_in_cluster = named$apart,
               mds_1 = pts[go_id, 1], mds_2 = pts[go_id, 2])]
    checks <- list(clusters = length(members), singletons = sum(lengths(members) == 1), largest = max(lengths(members)),
                   min_within = if (all(is.na(within))) NA_real_ else min(within, na.rm = TRUE),
                   complement = uniqueN(stats::cutree(tree, h = 1 - cfg$threshold)), fit = mds$GOF[1])
    orders <- order_check(sim, unname(setNames(red$score, red$go)[kept]), own_cut[kept],
                          named$clusters[number <= cfg$n_labels, representative])
  }

  if (n_kept > 0) {
    out[, `:=`(parent_term = term[match(parent_go_id, go_id)], is_representative = go_id == parent_go_id, labelled = cluster <= cfg$n_labels,
               p_value = if (is.null(p_values)) NA_real_ else unname(p_values[go_id]),
               score = round(score, 4), uniqueness = round(uniqueness, 4), uniqueness_in_cluster = round(uniqueness_in_cluster, 4),
               mds_1 = round(mds_1, 5), mds_2 = round(mds_2, 5))]
    if (anyNA(out[, .(term, cluster, parent_term)])) stop("A placed term lacks a name, a cluster, or a representative")
    if (out[is_representative == TRUE, .N] != uniqueN(out$cluster)) stop("A cluster has no representative or more than one")
    out <- out[order(cluster, -is_representative, -score, go_id), ..cluster_columns]
  }
  n_drop <- function(key) sum(dropped$reason == drop_reasons[[key]])
  summary_row <- data.table(
    terms_in = length(ids), dropped_obsolete = n_drop("obsolete"), dropped_merged = n_drop("merged"), dropped_unknown_id = n_drop("unknown_id"),
    dropped_other_ontology = n_drop("other_ontology"), dropped_not_in_gosemsim = n_drop("not_in_gosemsim"),
    dropped_no_gene_in_organism_db = n_drop("no_gene"), terms_placed = n_kept,
    score = score_is, clusters = checks$clusters, clusters_of_one_term = checks$singletons, largest_cluster = checks$largest,
    lowest_similarity_within_a_cluster = checks$min_within, clusters_if_cut_at_similarity_of_threshold = checks$complement,
    mds_goodness_of_fit = round(checks$fit, 3), orders_tried = orders$tried, orders_with_another_partition = orders$other,
    lowest_adjusted_rand_index = round(orders$rand, 3), fewest_clusters_in_another_order = orders$fewest,
    most_clusters_in_another_order = orders$most, fewest_labelled_representatives_kept = orders$kept)
  if (summary_row[, terms_in != dropped_obsolete + dropped_merged + dropped_unknown_id + dropped_other_ontology + dropped_not_in_gosemsim +
                    dropped_no_gene_in_organism_db + terms_placed]) {
    stop("Terms in do not equal terms dropped plus terms placed")
  }
  list(table = out, dropped = dropped, sim = sim, summary = summary_row)
}

# Write one similarity matrix as a table, gzipped only when the plain file is larger than gzip_above_mb. NULL: nothing to write.
write_similarity <- function(sim, name) {
  plain <- tab(paste0(name, ".csv")); zipped <- paste0(plain, ".gz")
  unlink(c(plain, zipped))
  if (is.null(sim)) return(invisible(NULL))
  as_table <- data.table(go_id = rownames(sim), as.data.table(sim))
  write_table(as_table, plain)
  if (file.size(plain) > cfg$gzip_above_mb * 1e6) { write_table(as_table, zipped); unlink(plain) }
  invisible(NULL)
}

results <- list(); summary_rows <- list(); dropped_rows <- list()
run_list <- function(list_name, list_label, s_i, o_i, ids, p_values, term_names) {
  res <- cluster_term_list(ids, p_values, o_i)
  key <- paste(list_name, s_i, o_i, sep = "_")
  write_table(cbind(data.table(set = rep(s_i, nrow(res$table)), ontology = rep(o_i, nrow(res$table))), res$table),
              tab(sprintf("semantic_clusters_%s.csv", key)))
  write_similarity(res$sim, paste0("similarity_", key))
  results[[key]] <<- res
  summary_rows[[key]] <<- cbind(data.table(term_list = list_label, mode = "legacy", set = s_i, ontology = o_i), res$summary)
  dropped_rows[[key]] <<- cbind(data.table(term_list = rep(list_label, nrow(res$dropped)), mode = rep("legacy", nrow(res$dropped)),
                                           set = rep(s_i, nrow(res$dropped)), ontology = rep(o_i, nrow(res$dropped))),
                                res$dropped[, .(go_id, term = unname(term_names[go_id]), reason)])
  log_msg("%-22s %-6s %s: %4d terms, %3d dropped, %4d placed, %3d clusters", list_label, s_i, o_i, res$summary$terms_in,
          res$summary$terms_in - res$summary$terms_placed, res$summary$terms_placed, res$summary$clusters)
  invisible(key)
}
label_legacy <- sprintf("legacy, p < %g", config$go_pvalue)
label_drawn  <- "as drawn in the draft"
if (paper_part) {
for (s_i in names(sets)) for (o_i in names(onts)) {
  x <- legacy[set == s_i & ontology == o_i & below_p_cutoff == TRUE]
  run_list("legacy", label_legacy, s_i, o_i, x$go_id, setNames(x$p_fisher, x$go_id), setNames(x$term, x$go_id))
}
for (s_i in names(sets)) {
  x <- annotated_lists[set == s_i & ontology == cfg$figure_ontology]
  run_list("as_drawn", label_drawn, s_i, cfg$figure_ontology, x$go_id, NULL, setNames(x$term_annotation, x$go_id))
}
}   # end of the paper part of section 2

# =============================================================================
# 3. Figures 6, 7 and 8
# =============================================================================
wrap_text   <- function(x, width) paste(strwrap(x, width), collapse = "\n")
wrap_label  <- function(x) vapply(x, wrap_text, character(1), width = cfg$label_width, USE.NAMES = FALSE)
size_breaks <- function(size_max) { b <- pretty(c(1, size_max), n = 4); sort(unique(c(1, b[b > 1 & b <= size_max & b == round(b)]))) }
n_of        <- function(n, noun) sprintf("%s %s%s", format(n, big.mark = ","), noun, if (n == 1) "" else "s")       # "1 term", "2 terms"
circle_size <- function(n_terms, size_max) cfg$max_circle * sqrt(n_terms / size_max)       # area in proportion to the number of terms
tint_of <- function(colour) {                                                           # the colour mixed with white
  mixed <- grDevices::col2rgb(colour) * cfg$tint_share + 255 * (1 - cfg$tint_share)
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ], maxColorValue = 255)
}
# Radius in inches of a ggplot2 point of shape 21: R draws the circle with a radius of 0.375 font sizes.
circle_radius_in <- function(size, stroke) 0.375 * (size * ggplot2::.pt + stroke * ggplot2::.stroke / 2) / 72
# Width and height in inches of the box ggrepel draws around a label, estimated from the text (no graphics device is needed).
label_box_in <- function(text, font_mm) {
  parts <- strsplit(text, "\n", fixed = TRUE); font_pt <- font_mm * ggplot2::.pt; pad_pt <- cfg$label_padding * 14.4
  data.table(w = (vapply(parts, function(v) max(nchar(v)), numeric(1)) * cfg$label_box_em[["char"]] * font_pt + 2 * pad_pt) / 72,
             h = ((cfg$label_box_em[["first_line"]] + cfg$label_box_em[["next_line"]] * (lengths(parts) - 1)) * font_pt + 2 * pad_pt) / 72)
}

# Positions of the labels of one panel, chosen by rule. All lengths are in inches from the lower left corner of the panel.
#   lab      one row per label: id, u, v (centre of its circle), r (radius of the circle), w, h (size of its box)
#   circles  every circle of the panel: id, u, v, r, hard (TRUE: a label may not touch it; FALSE: it costs charge_circle)
#   panel    width and height of the panel
# Labels are placed one after the other. Each takes the cheapest allowed position among label_rings x label_angles positions
# around its circle. A position is allowed when the box is inside the panel, touches no label already placed and no hard
# circle, and when its leader line and the leader lines already drawn cross neither each other nor a label box. The cost is the
# distance from the circle, plus a charge for each other circle under the box or the leader. When no position is allowed, the
# leader rule is dropped, then the circle rule, then the label rule.
# The result depends on which label is placed first, so six orders are tried: the most crowded places first, the least crowded
# first, and from the left, the right, the bottom and the top. The order that drops the fewest rules is kept; among equals,
# the one with the lowest total cost, then the earlier in this list.
# Returns the middle of each box and the number of rules that had to be dropped for it.
place_labels <- function(lab, circles, panel) {
  n_lab <- nrow(lab)
  if (n_lab == 0) return(data.table(label_u = numeric(0), label_v = numeric(0), rules_dropped = integer(0)))
  fixed <- circles[hard == TRUE]; loose <- circles[hard == FALSE]
  crowd <- vapply(seq_len(n_lab), function(i) sum((fixed$u - lab$u[i])^2 + (fixed$v - lab$v[i])^2 < cfg$crowd_radius^2), numeric(1))
  middle <- c(mean(fixed$u), mean(fixed$v))
  cand <- CJ(ring = cfg$label_rings, angle = cfg$label_angles, sorted = FALSE); n_cand <- nrow(cand)
  ux <- cos(cand$angle); uy <- sin(cand$angle)
  touches <- function(x1, y1, x2, y2, cu, cv, cr) {               # which boxes touch the circle
    dx <- pmax(x1 - cu, 0, cu - x2); dy <- pmax(y1 - cv, 0, cv - y2); dx^2 + dy^2 < cr^2
  }
  turn_of <- function(pu, pv, qu, qv, ru, rv) sign((qu - pu) * (rv - pv) - (qv - pv) * (ru - pu))   # side of p-q on which r lies
  # Every position of every label, and what does not depend on the labels placed before it.
  positions <- lapply(seq_len(n_lab), function(i) {
    p <- lab[i]
    reach <- p$r + cfg$label_gap + cand$ring + p$w / 2 * abs(ux) + p$h / 2 * abs(uy)          # the box clears the label's own circle
    bu <- p$u + reach * ux; bv <- p$v + reach * uy                                         # middle of the box
    x1 <- bu - p$w / 2; x2 <- bu + p$w / 2; y1 <- bv - p$h / 2; y2 <- bv + p$h / 2
    lu <- p$u + outer(bu - p$u, cfg$leader_steps); lv <- p$v + outer(bv - p$v, cfg$leader_steps)     # leader: from the circle to the middle of the box
    free <- !(lu > x1 & lu < x2 & lv > y1 & lv < y2)                                         # the part of the leader not under its own box
    inside <- x1 >= cfg$label_gap & x2 <= panel[1] - cfg$label_gap & y1 >= cfg$label_gap & y2 <= panel[2] - cfg$label_gap
    hard_hit <- rep(FALSE, n_cand); leader_hit <- numeric(n_cand); loose_hit <- numeric(n_cand)
    for (k in which(fixed$id != p$id)) {
      hard_hit <- hard_hit | touches(x1, y1, x2, y2, fixed$u[k], fixed$v[k], fixed$r[k] + cfg$label_gap)
      leader_hit <- leader_hit + (rowSums(free & (lu - fixed$u[k])^2 + (lv - fixed$v[k])^2 < fixed$r[k]^2) > 0)
    }
    for (k in seq_len(nrow(loose))) loose_hit <- loose_hit + touches(x1, y1, x2, y2, loose$u[k], loose$v[k], loose$r[k])
    outward <- atan2(p$v - middle[2], p$u - middle[1])
    turn <- abs(atan2(sin(cand$angle - outward), cos(cand$angle - outward))) / pi          # 0: away from the middle of the labelled circles
    list(pu = p$u, pv = p$v, bu = bu, bv = bv, x1 = x1, y1 = y1, x2 = x2, y2 = y2, lu = lu, lv = lv, free = free, inside = inside,
         hard_hit = hard_hit, cost = cand$ring + cfg$charge_circle * loose_hit + cfg$charge_leader * leader_hit + cfg$charge_turn * turn)
  })
  place_in_order <- function(turns) {
    at_u <- rep(NA_real_, n_lab); at_v <- rep(NA_real_, n_lab); dropped <- rep(NA_integer_, n_lab); paid <- 0
    placed <- list()
    for (i in turns) {
      o <- positions[[i]]
      label_hit <- rep(FALSE, n_cand); leader_bad <- numeric(n_cand)
      for (q in placed) {
        label_hit <- label_hit | (o$x1 < q$x2 + cfg$label_gap & o$x2 > q$x1 - cfg$label_gap & o$y1 < q$y2 + cfg$label_gap & o$y2 > q$y1 - cfg$label_gap)
        # two leader lines cross when the ends of each lie on both sides of the other
        crossing <- turn_of(o$pu, o$pv, o$bu, o$bv, q$pu, q$pv) * turn_of(o$pu, o$pv, o$bu, o$bv, q$bu, q$bv) < 0 &
                    turn_of(q$pu, q$pv, q$bu, q$bv, o$pu, o$pv) * turn_of(q$pu, q$pv, q$bu, q$bv, o$bu, o$bv) < 0
        under_label <- rowSums(o$lu > q$x1 & o$lu < q$x2 & o$lv > q$y1 & o$lv < q$y2) > 0
        over_leader <- rowSums(outer(o$x1, q$lu, "<") & outer(o$x2, q$lu, ">") & outer(o$y1, q$lv, "<") & outer(o$y2, q$lv, ">")) > 0
        leader_bad <- leader_bad + crossing + under_label + over_leader
      }
      cost <- o$cost + cfg$charge_leader * leader_bad
      allowed <- list(o$inside & !label_hit & !o$hard_hit & leader_bad == 0, o$inside & !label_hit & !o$hard_hit, o$inside & !label_hit, o$inside,
                      rep(TRUE, n_cand))
      level <- which(vapply(allowed, any, logical(1)))[1]
      best <- which(allowed[[level]])[which.min(cost[allowed[[level]]])]
      at_u[i] <- o$bu[best]; at_v[i] <- o$bv[best]; dropped[i] <- level - 1L; paid <- paid + cost[best]
      placed[[length(placed) + 1]] <- list(x1 = o$x1[best], y1 = o$y1[best], x2 = o$x2[best], y2 = o$y2[best], pu = o$pu, pv = o$pv, bu = o$bu[best],
                                           bv = o$bv[best], lu = o$lu[best, o$free[best, ]], lv = o$lv[best, o$free[best, ]])
    }
    list(spot = data.table(label_u = at_u, label_v = at_v, rules_dropped = dropped), dropped = sum(dropped), paid = paid)
  }
  first <- seq_len(n_lab)
  tries <- lapply(list(order(-crowd, first), order(crowd, first), order(lab$u, first), order(-lab$u, first), order(lab$v, first),
                       order(-lab$v, first)), place_in_order)
  tries[[order(vapply(tries, `[[`, numeric(1), "dropped"), vapply(tries, `[[`, numeric(1), "paid"))[1]]]$spot
}

# Say in the log when labels could not be placed within every rule.
placement_note <- function(what, spot) {
  if (any(spot$rules_dropped > 0)) {
    log_msg("%s: %d of %d labels placed with the leader rule dropped, %d with a further rule dropped", what, sum(spot$rules_dropped == 1),
            nrow(spot), sum(spot$rules_dropped > 1))
  }
}

# One panel: every term a dot, every cluster a circle at its representative, each labelled cluster in its own colour.
semantic_panel <- function(terms, title, subtitle, size_max) {
  subtitle <- wrap_text(subtitle, cfg$subtitle_width)
  if (nrow(terms) < 3) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "Fewer than three terms placed.\nNo layout is drawn.", size = 3.5) +
             labs(title = title, subtitle = subtitle) + theme_void(base_size = 11) +
             theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5)))
  }
  terms <- copy(terms)[, hue := fifelse(labelled, cfg$cluster_palette[pmin(cluster, length(cfg$cluster_palette))], cfg$grey_other)]
  reps <- terms[is_representative == TRUE][order(-cluster_size, cluster)]                  # large circles first: small ones stay in view
  reps[, `:=`(circle = circle_size(cluster_size, size_max), label_text = fifelse(labelled, wrap_label(term), ""), tint = tint_of(hue))]
  rays <- merge(terms[labelled == TRUE & is_representative == FALSE], reps[, .(cluster, rep_1 = mds_1, rep_2 = mds_2)], by = "cluster")
  # A square panel with one scale on both axes (a map of distances), whatever the shape of the cloud of terms.
  half <- max(diff(range(terms$mds_1)), diff(range(terms$mds_2))) / 2
  half <- (if (half > 0) half else 1) * (1 + 2 * cfg$panel_margin)
  lim_1 <- mean(range(terms$mds_1)) + c(-half, half); lim_2 <- mean(range(terms$mds_2)) + c(-half, half)
  # Label positions by rule; ggrepel starts from them.
  inches <- function(x, lim) (x - lim[1]) / (2 * half) * cfg$panel_in
  circles <- reps[, .(id = .I, u = inches(mds_1, lim_1), v = inches(mds_2, lim_2), r = circle_radius_in(circle, cfg$circle_stroke), hard = labelled)]
  lab <- circles[hard == TRUE]
  spot <- place_labels(cbind(lab[, .(id, u, v, r)], label_box_in(reps$label_text[lab$id], cfg$label_mm[["cluster"]])), circles, c(cfg$panel_in, cfg$panel_in))
  placement_note(title, spot)
  callouts <- reps[lab$id]
  callouts[, `:=`(nudge_1 = (spot$label_u - lab$u) / cfg$panel_in * 2 * half, nudge_2 = (spot$label_v - lab$v) / cfg$panel_in * 2 * half)]
  breaks <- size_breaks(size_max)
  ggplot(terms, aes(mds_1, mds_2)) +
    geom_point(data = terms[labelled == FALSE], colour = cfg$grey_other, size = 0.7) +
    geom_segment(data = rays, aes(xend = rep_1, yend = rep_2, colour = hue), linewidth = 0.25, alpha = 0.5) +
    geom_point(data = reps[labelled == FALSE], aes(size = circle), shape = 21, fill = cfg$grey_other, colour = "white", stroke = 0.5) +
    geom_point(data = terms[labelled == TRUE], aes(colour = hue), size = 0.9) +
    geom_point(data = reps[labelled == TRUE], aes(size = circle, fill = hue), shape = 21, colour = "white", stroke = cfg$circle_stroke) +
    ggrepel::geom_label_repel(data = callouts, aes(label = label_text, point.size = circle, segment.colour = hue, fill = tint),
                              nudge_x = callouts$nudge_1, nudge_y = callouts$nudge_2, size = cfg$label_mm[["cluster"]], lineheight = 0.9,
                              colour = "black", label.size = 0.2,
                              label.padding = cfg$label_padding, box.padding = 0.1, point.padding = 0.1, min.segment.length = 0, segment.size = 0.5,
                              max.overlaps = Inf, max.iter = 20000, max.time = Inf, seed = cfg$label_seed, show.legend = FALSE) +
    scale_colour_identity() + scale_fill_identity() +
    scale_size_identity(guide = "legend", limits = c(0, cfg$max_circle), breaks = circle_size(breaks, size_max), labels = breaks,
                        name = "Terms in the cluster") +
    coord_equal(xlim = lim_1, ylim = lim_2, expand = FALSE) +
    guides(size = guide_legend(override.aes = list(fill = "grey45"))) +
    labs(x = "MDS axis 1", y = "MDS axis 2", title = title, subtitle = subtitle) +
    theme_wgbs(base_size = 11) + theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5))
}
panel_subtitle <- function(res, representative_rule) {
  sprintf("%s in the list, %s placed, %s. Each cluster is named after %s.", n_of(res$summary$terms_in, "term"),
          format(res$summary$terms_placed, big.mark = ","), n_of(res$summary$clusters, "cluster"), representative_rule)
}
layout_fit <- function(res) {
  fit <- res$summary$mds_goodness_of_fit
  if (is.na(fit)) "" else sprintf("Goodness of fit of the layout in two dimensions: %.2f.", fit)
}
# The caption of a figure, wrapped to its width. `score_note` says which score ranks the clusters, `panel_note` what holds across
# panels, and `closing_note` ends the caption.
figure_caption <- function(score_note, panel_note, closing_note, figure_width) {
  lines <- c(
    "Each dot is a GO term. Each circle is a cluster of terms, drawn at its representative term; its area shows the number of terms.",
    paste0(sprintf("Colour: the clusters whose representative scores highest%s, at most %d, one colour each. ", score_note, cfg$n_labels),
           "Each is labelled with the name of its representative, and lines join its terms to the representative. ", panel_note,
           "Grey: all other clusters. ",
           "When terms share the highest score of a cluster, the one least similar to the others names it, then the higher GO ID."),
    sprintf("Layout: classical multidimensional scaling of 1 - similarity (%s measure, %s). Clusters: rrvgo, threshold %g.%s", cfg$measure,
            cfg$organism_db, cfg$threshold, closing_note))
  paste(vapply(lines, wrap_text, character(1), width = floor(cfg$caption_per_in * figure_width)), collapse = "\n")
}
caption_theme <- theme(plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.1))
reproduction_note <- " The legacy enrichment is unsound (issues I-03, I-39, I-40); this figure is a reproduction, not a result."
if (paper_part) {   # sections 3 (figures), 4 and 5
for (s_i in names(sets)) {
  drawn <- results[[paste("as_drawn", s_i, cfg$figure_ontology, sep = "_")]]; described <- results[[paste("legacy", s_i, cfg$figure_ontology, sep = "_")]]
  both <- lists_compared[set == s_i & ontology == cfg$figure_ontology]
  size_max <- max(1L, drawn$table$cluster_size, described$table$cluster_size)
  log_msg("drawing figure_%d_legacy (%s)", cfg$figure_numbers[[s_i]], sets[[s_i]])
  panel_a <- semantic_panel(drawn$table, sprintf("A  As drawn: all %s terms annotated to genes with a DML of this set", cfg$figure_ontology),
                            paste(panel_subtitle(drawn, "its term with the most annotated genes"),
                                  sprintf("%d of the %s terms of the list have p < %g in the legacy test.", both$terms_in_both,
                                          format(both$terms_annotated_legacy_reader, big.mark = ","), config$go_pvalue), layout_fit(drawn)), size_max)
  panel_b <- semantic_panel(described$table, sprintf("B  As described: %s terms with p < %g in the legacy test", cfg$figure_ontology, config$go_pvalue),
                            paste(panel_subtitle(described, "its term with the smallest p-value"),
                                  sprintf("%d of the %d terms of the list are also in the list of panel A.", both$terms_in_both,
                                          both$terms_legacy_p_below_cutoff), layout_fit(described)), size_max)
  figure <- (panel_a | panel_b) + plot_layout(guides = "collect") +
    plot_annotation(title = sprintf("Figure %d of the draft, regenerated in R: semantic similarity of GO terms, %s", cfg$figure_numbers[[s_i]],
                                    sets[[s_i]]),
                    caption = figure_caption(" (panel A: most annotated genes; panel B: smallest p-value)",
                                             "The same colour in the two panels does not mark the same cluster. ", reproduction_note,
                                             cfg$figure_size[1]),
                    theme = theme(plot.title = element_text(size = 12, face = "bold")) + caption_theme) &
    theme(legend.position = "bottom")
  save_figure(figure, fig(sprintf("figure_%d_legacy", cfg$figure_numbers[[s_i]])), width = cfg$figure_size[1], height = cfg$figure_size[2])
}

# =============================================================================
# 4. The REVIGO exports that survive
# =============================================================================
# An export is an R script with one c("GO:...", "name", numbers) record per term.
read_revigo_export <- function(path) {
  txt <- readLines(path, warn = FALSE)
  header <- grep("^revigo\\.names <- c\\(", txt, value = TRUE)
  rows <- grep('c\\("GO:[0-9]{7}",', txt, value = TRUE)
  hit <- regmatches(rows, regexec('c\\("(GO:[0-9]{7})","(.*)",((?:-?[0-9.eE+-]+,)*-?[0-9.eE+-]+)\\)', rows, perl = TRUE))
  if (length(header) != 1 || length(rows) == 0 || any(lengths(hit) != 4)) stop("Unexpected format in REVIGO export ", path)
  columns <- gsub('"', "", regmatches(header, gregexpr('"[^"]+"', header))[[1]])
  numbers <- do.call(rbind, strsplit(vapply(hit, `[`, "", 4), ",", fixed = TRUE))
  if (ncol(numbers) != length(columns) - 2) stop("Records and column names do not agree in REVIGO export ", path)
  d <- data.table(vapply(hit, `[`, "", 2), vapply(hit, `[`, "", 3), as.data.table(apply(numbers, 2, as.numeric)))
  setnames(d, columns)
  if (anyNA(d) || anyDuplicated(d$term_ID)) stop("Missing values or repeated terms in REVIGO export ", path)
  rule <- regmatches(txt, regexec("one\\.data\\$dispensability < ([0-9.]+)", txt)); rule <- rule[lengths(rule) == 2]
  if (length(rule) != 1) stop("The label rule of REVIGO export ", path, " was not found")
  list(data = d, label_below_dispensability = as.numeric(rule[[1]][2]), md5 = unname(tools::md5sum(path)))
}
# A Fisher table of the other analysis: no quoting, commas inside term names, p as printed. A record that is
# cut short (the BP file holds runs of NUL bytes) does not fit the pattern and is not read. The lines are matched
# as bytes: 19 lines of the BP file hold bytes that are not valid UTF-8, and under a UTF-8 locale a match on
# characters skips them, so that 2,430 records were read there where the C locale reads 2,449.
read_other_fisher <- function(path) {
  txt <- grep("^GO:[0-9]{7},", readLines(path, warn = FALSE), value = TRUE, useBytes = TRUE)
  hit <- regmatches(txt, regexec("^(GO:[0-9]{7}),.*,[0-9]+,[0-9]+,[0-9.]+,(< ?)?([0-9.eE+-]+),?$", txt, useBytes = TRUE))
  good <- lengths(hit) == 4
  d <- data.table(go_id = vapply(hit[good], `[`, "", 2), p = as.numeric(vapply(hit[good], `[`, "", 4)))
  if (nrow(d) == 0 || anyNA(d$p) || anyDuplicated(d$go_id)) stop("Unexpected format in ", path)
  list(data = d, lines_not_read = sum(!good))
}

exports <- lapply(setNames(names(cfg$revigo_exports), names(cfg$revigo_exports)), function(e_i) read_revigo_export(export_file(e_i)))
candidates <- list()
add_candidate <- function(list_name, source_file, ids, p = NULL) {
  candidates[[length(candidates) + 1]] <<- list(list_name = list_name, source_file = repo_relative(source_file), ids = ids, p = p)
}
# Candidate term lists: the nine legacy lists, the nine annotated lists, and the lists of the other analysis.
for (s_i in names(sets)) for (o_i in names(onts)) {
  x <- legacy[set == s_i & ontology == o_i]
  add_candidate(sprintf("legacy test, p < %g: %s, %s", config$go_pvalue, sets[[s_i]], o_i), legacy_fisher_file(s_i, o_i),
                x[below_p_cutoff == TRUE, go_id], setNames(x$p_fisher, x$go_id))
}
for (s_i in names(sets)) for (o_i in names(onts)) {
  add_candidate(sprintf("all terms annotated to genes with a DML (legacy reader): %s, %s", sets[[s_i]], o_i), annotation_file(s_i),
                annotated_lists[set == s_i & ontology == o_i, go_id])
}
for (o_i in names(onts)) {
  x <- read_other_fisher(other_fisher_file(o_i))
  if (x$lines_not_read > 0) log_msg("%s: %d damaged line(s) not read", repo_relative(other_fisher_file(o_i)), x$lines_not_read)
  add_candidate(sprintf("another analysis, not the manuscript's, p < %g: ploidy, %s", config$go_pvalue, o_i), other_fisher_file(o_i),
                x$data[p < config$go_pvalue, go_id], setNames(x$data$p, x$data$go_id))
}
# REVIGO gives a term under its current ID, so a list term can be in an export under the ID it was merged into.
export_path <- function(e_i) repo_relative(export_file(e_i))
first_with_md5 <- vapply(exports, function(e) names(exports)[match(e$md5, vapply(exports, `[[`, "", "md5"))], "")
export_rows <- rbindlist(lapply(names(exports), function(e_i) {
  d <- exports[[e_i]]$data; ids <- d$term_ID; has_p <- "value" %in% names(d)
  rbindlist(lapply(candidates, function(cand) {
    same_id <- ids %chin% cand$ids
    by_secondary <- !same_id & ids %chin% primary_id(cand$ids)
    same_p <- if (is.null(cand$p) || !has_p) NA_integer_ else sum(abs(d$value - log10(cand$p[ids])) < cfg$value_tolerance, na.rm = TRUE)
    data.table(export = export_path(e_i), same_content_as = if (first_with_md5[[e_i]] == e_i) "" else export_path(first_with_md5[[e_i]]),
               export_terms = length(ids), export_has_p_values = has_p, term_list = cand$list_name, term_list_source = cand$source_file,
               list_terms = length(cand$ids), export_terms_in_list = sum(same_id), export_terms_in_list_by_secondary_id = sum(by_secondary),
               list_terms_in_export = sum(primary_id(cand$ids) %chin% ids | cand$ids %chin% ids), export_values_equal_to_log10_p = same_p)
  }))
}))
export_rows[, percent_of_export_in_list := round(100 * (export_terms_in_list + export_terms_in_list_by_secondary_id) / export_terms, 1)]
export_rows[, best_match := seq_len(.N) == order(-(export_terms_in_list + export_terms_in_list_by_secondary_id), list_terms)[1], by = export]
write_table(export_rows, tab("revigo_exports.csv"))
print(export_rows[best_match == TRUE, .(export, export_terms, term_list, list_terms, export_terms_in_list, export_terms_in_list_by_secondary_id,
                                        list_terms_in_export)])

# Redraw each distinct export from its own data, with the label rule of its own plot code. Label positions: place_labels().
for (e_i in names(exports)) {
  if (first_with_md5[[e_i]] != e_i) {
    log_msg("%s has the same content as %s; not drawn twice", export_path(e_i), export_path(first_with_md5[[e_i]]))
    next
  }
  d <- copy(exports[[e_i]]$data)[order(-log_size)]; label_cut <- exports[[e_i]]$label_below_dispensability
  has_p <- "value" %in% names(d)
  spread <- diff(range(d$log_size))
  bubble_size <- function(v) if (spread > 0) cfg$export_bubble[1] + diff(cfg$export_bubble) * sqrt((v - min(d$log_size)) / spread) else mean(cfg$export_bubble)
  d[, `:=`(shade = if (has_p) value else uniqueness, label_text = fifelse(dispensability < label_cut, wrap_label(description), ""),
           bubble = bubble_size(log_size))]
  lim_x <- range(d$plot_X) + c(-1, 1) * cfg$export_margin[1] * diff(range(d$plot_X))
  lim_y <- range(d$plot_Y) + c(-1, 1) * cfg$export_margin[2] * diff(range(d$plot_Y))
  circles <- d[, .(id = .I, u = (plot_X - lim_x[1]) / diff(lim_x) * cfg$export_panel_in[1], v = (plot_Y - lim_y[1]) / diff(lim_y) * cfg$export_panel_in[2],
                   r = circle_radius_in(bubble, cfg$export_stroke), hard = label_text != "")]
  lab <- circles[hard == TRUE]
  spot <- place_labels(cbind(lab[, .(id, u, v, r)], label_box_in(d$label_text[lab$id], cfg$label_mm[["export"]])), circles, cfg$export_panel_in)
  placement_note(export_path(e_i), spot)
  callouts <- d[lab$id]
  callouts[, `:=`(nudge_1 = (spot$label_u - lab$u) / cfg$export_panel_in[1] * diff(lim_x),
                  nudge_2 = (spot$label_v - lab$v) / cfg$export_panel_in[2] * diff(lim_y))]
  size_keys <- pretty(d$log_size, n = 3); size_keys <- size_keys[size_keys >= min(d$log_size) & size_keys <= max(d$log_size)]
  best <- export_rows[export == export_path(e_i) & best_match == TRUE]
  n_secondary <- best$export_terms_in_list_by_secondary_id
  p_export <- ggplot(d, aes(plot_X, plot_Y)) +
    geom_point(aes(size = bubble, fill = shade), shape = 21, colour = "white", stroke = cfg$export_stroke, alpha = 0.9) +
    ggrepel::geom_label_repel(data = callouts, aes(label = label_text, point.size = bubble), nudge_x = callouts$nudge_1, nudge_y = callouts$nudge_2,
                              size = cfg$label_mm[["export"]], lineheight = 0.9, colour = "black", fill = scales::alpha("white", 0.9), label.size = 0.2,
                              label.padding = cfg$label_padding, box.padding = 0.1, point.padding = 0.1, min.segment.length = 0,
                              segment.colour = "grey25", segment.size = 0.3, max.overlaps = Inf, max.iter = 20000, max.time = Inf,
                              seed = cfg$label_seed) +
    scale_fill_gradient(low = cfg$revigo_colours[["low"]], high = cfg$revigo_colours[["high"]],
                        name = if (has_p) "Value in the export\n(log10 p)" else "Uniqueness\n(REVIGO)") +
    scale_size_identity(guide = "legend", limits = c(0, cfg$export_bubble[2]), breaks = bubble_size(size_keys), labels = size_keys,
                        name = "Log size\n(REVIGO)") +
    coord_cartesian(xlim = lim_x, ylim = lim_y, expand = FALSE) +
    guides(size = guide_legend(order = 1, override.aes = list(fill = "grey45")), fill = guide_colourbar(order = 2)) +
    labs(x = "Semantic space X (REVIGO)", y = "Semantic space Y (REVIGO)",
         title = sprintf("REVIGO export %s, redrawn from its embedded data", export_path(e_i)),
         subtitle = paste0(sprintf("%d terms. %s. Labels: terms with dispensability below %g, the rule of the export's own plot code.\n", nrow(d),
                                   if (has_p) "Colour: the value column of the export" else
                                     "The export holds no p-values; colour is REVIGO's uniqueness", label_cut),
                           sprintf("Closest term list: %s (%s terms). It holds %d of the export's terms%s.", best$term_list,
                                   format(best$list_terms, big.mark = ","), best$export_terms_in_list + n_secondary,
                                   if (n_secondary > 0) sprintf(", %d of them under a secondary ID", n_secondary) else ""))) +
    theme_wgbs(base_size = 11) + theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5))
  save_figure(p_export, fig(paste0("figure_revigo_export_", e_i)), width = cfg$export_size[1], height = cfg$export_size[2])
}

# =============================================================================
# 5. The GO terms named in Results 3.4
# =============================================================================
# Quoted from the manuscript draft (Results 3.4, paragraphs on ploidy-DML, heat-DML and interaction-DML). `read_as` is the GO
# term name each phrase is taken to mean; where it differs from the phrase the match is "close" and the reason is given.
# Phrases that name a theme and not a term (for example "genomic maintenance and replication") are not listed.
named <- rbindlist(list(
  data.table(set = "ploidy", phrase = c("DNA replication", "homologous chromosome pairing at meiosis", "chromosome segregation",
                                        "ATP generation from ADP", "proton transmembrane transport", "calcium ion transport")),
  data.table(set = "heat", phrase = c("the MAPK cascade", "the MAPK cascade", "ER overload response", "signal transduction",
                                      "physiological responses to osmotic and oxidative stress",
                                      "physiological responses to osmotic and oxidative stress",
                                      "mRNA methylation", "mRNA export from the nucleus in response to heat stress", "rRNA catabolic processes"),
             read_as = c("MAPK cascade", "stress-activated MAPK cascade", NA, NA, "response to osmotic stress", "response to oxidative stress", NA,
                         "mRNA export from nucleus in response to heat stress", "rRNA catabolic process"),
             reading_note = c("article dropped", "the label on the archived image of Figure 7 (08_archive/figures/revigo_figure/HEAT-BP.png)", NA, NA,
                              "one phrase for two terms", "one phrase for two terms", NA, "article dropped", "plural")),
  data.table(set = "shared", phrase = c("tRNA modification", "RNA 3'-end processing", "post-translational protein modification",
                                        "calcium-mediated signaling", "ER overload response", "response to osmotic stress",
                                        "positive regulation of gene expression", "gene silencing by RNA"))),
  fill = TRUE)
named[is.na(read_as), read_as := phrase]
named[is.na(reading_note), reading_note := ""]
named[, row_in_text := .I]

# Match on term names: first the GO.db of the R installation, then the release the legacy analysis ran on.
name_current <- go_terms[, .(go_id, name_key = tolower(term))]; name_legacy <- legacy_names[, .(go_id, name_key = tolower(term_legacy_release))]
if (anyDuplicated(name_current$name_key) || anyDuplicated(name_legacy$name_key)) stop("Two GO terms share a name")
named[, id_current := name_current$go_id[match(tolower(read_as), name_current$name_key)]]
named[, id_legacy := name_legacy$go_id[match(tolower(read_as), name_legacy$name_key)]]
if (named[!is.na(id_current) & !is.na(id_legacy) & id_current != id_legacy, .N] > 0) stop("A name points to different GO IDs in the two releases")
named[, go_id := fifelse(!is.na(id_current), id_current, id_legacy)]
named[, `:=`(match_type = fifelse(is.na(go_id), "not found", fifelse(tolower(phrase) == tolower(read_as), "exact", "close")),
             name_matched_in = fifelse(!is.na(id_current), paste("GO.db", go_db_version),
                               fifelse(!is.na(id_legacy), paste("GO.db", config$legacy_go_db$version, "(release of the legacy analysis)"), "")),
             term_current = go_terms$term[match(go_id, go_terms$go_id)], ontology = legacy_names$ontology[match(go_id, legacy_names$go_id)])]
named[is.na(ontology), ontology := go_terms$ontology[match(go_id, go_terms$go_id)]]
named[, term_legacy_release := legacy_names$term_legacy_release[match(go_id, legacy_names$go_id)]]

# The legacy test, the list the figure was drawn from, and the gene-level test, each for the DML set of the paragraph.
n_named <- nrow(named)
named <- merge(named, legacy[, .(set, go_id, legacy_rank = rank, legacy_p = p_fisher, in_legacy_p_list = below_p_cutoff,
                                 transcripts_annotated = annotated, transcripts_of_interest = significant,
                                 genes_behind_transcripts_of_interest = genes_significant)],
               by = c("set", "go_id"), all.x = TRUE, sort = FALSE)
named[is.na(in_legacy_p_list), in_legacy_p_list := FALSE]
drawn_ids <- function(s_i) annotated_lists[set == s_i & ontology == cfg$figure_ontology, go_id]
named[, in_list_as_drawn := unname(mapply(function(s_i, g) {
  if (is.na(g)) return("no")
  ids <- drawn_ids(s_i); via <- ids[primary_id(ids) == g & ids != g]
  if (g %chin% ids) "yes" else if (length(via) > 0) paste0("under a secondary ID (", paste(via, collapse = ", "), ")") else "no"
}, set, go_id))]
drawn_clusters <- rbindlist(lapply(names(sets), function(s_i) {       # where the term sits in panel A of the regenerated figure
  x <- results[[paste("as_drawn", s_i, cfg$figure_ontology, sep = "_")]]$table
  cbind(data.table(set = rep(s_i, nrow(x))), x[, .(go_id, as_drawn_cluster = cluster, as_drawn_cluster_named_after = parent_term)])
}))
named <- merge(named, drawn_clusters, by = c("set", "go_id"), all.x = TRUE, sort = FALSE)
if (nrow(named) != n_named) stop("A join changed the number of named terms")
named[, `:=`(manuscript_set = unname(sets[set]), draft_figure = paste("Figure", cfg$figure_numbers[set]))]
named_out <- named[order(row_in_text),
                   .(set, manuscript_set, draft_figure, phrase, read_as, match_type, reading_note, name_matched_in, go_id, term_current,
                     term_legacy_release, ontology, in_legacy_p_list, legacy_rank, legacy_p, transcripts_annotated, transcripts_of_interest,
                     genes_behind_transcripts_of_interest, in_list_as_drawn, as_drawn_cluster, as_drawn_cluster_named_after)]
write_table(named_out, tab("named_terms_results_3_4.csv"))
print(named_out[, .(named_terms = .N, matched = sum(match_type != "not found"), in_legacy_p_list = sum(in_legacy_p_list),
                    in_list_as_drawn = sum(in_list_as_drawn != "no")), by = set])
}   # end of the paper part (sections 3 to 5)

# =============================================================================
# 6. The same procedure on the gene-level terms (optional run, --gene-level-go)
# =============================================================================
if (gene_part) {
# The terms named in Results 3.4 (the paper's table), tested at gene level in the analysis of the run.
named <- fread(named_terms_file, select = c("set", "phrase", "go_id"))
gene_level <- rbindlist(lapply(names(sets), function(s_i) fread(gene_level_file(s_i))[universe == "eligible"]))
if (anyDuplicated(gene_level[, .(set, go_id)])) stop("A GO ID occurs twice for one set in the gene-level tables")
hit <- match(paste(named$set, named$go_id), paste(gene_level$set, gene_level$go_id))
named[, `:=`(gene_level_genes_annotated = gene_level$genes_annotated[hit], gene_level_genes_with_dml = gene_level$genes_with_dml[hit],
             gene_level_p_wallenius = gene_level$p_wallenius[hit], gene_level_q_wallenius = gene_level$q_wallenius[hit])]
tested_term <- gene_level$tested[hit]
named[, gene_level_result := fifelse(is.na(go_id) | go_id == "", "no GO term matched",
                             fifelse(!go_id %chin% go_terms$go_id, paste("not tested: the term is not in GO.db", go_db_version),
                             fifelse(is.na(gene_level_genes_annotated), "not tested: no eligible gene carries the term",
                             fifelse(!tested_term, sprintf("not tested: fewer than %d eligible genes carry the term", config$go_min_genes),
                             fifelse(gene_level_q_wallenius < config$go_fdr, "enriched", "not enriched")))))]
write_table(named, tab("named_terms_results_3_4_gene_level.csv"))
print(named[, .(named_terms = .N, enriched_at_gene_level = sum(gene_level_result == "enriched")), by = set])

gene_level_columns <- c("mode", "set", "ontology", cluster_columns, "q_wallenius", "note")
label_gene <- sprintf("gene level, q < %g", config$go_fdr)
if (nrow(enriched) == 0) {
  # Step 10 calls no term at its settings; that is not the same as no term being enriched (issues.md, I-41).
  no_terms <- sprintf(paste("Step 10 calls no GO term at gene level at its configured settings (q < %g, eligible universe, Wallenius test)",
                            "in any DML set, ontology or mode; issues.md, I-41, gives the settings and tests under which terms pass.",
                            "Nothing was clustered and no figure was drawn."), config$go_fdr)
  gene_table <- as.data.table(c(setNames(rep(list(NA), length(gene_level_columns) - 1), head(gene_level_columns, -1)), list(note = no_terms)))
  summary_rows[["gene_level"]] <- cbind(data.table(term_list = label_gene, mode = "all", set = "all", ontology = "all"),
                                        cluster_term_list(character(0), numeric(0), cfg$figure_ontology)$summary)
  log_msg("%s", no_terms)
} else {
  groups <- unique(enriched[, .(mode, set, ontology)])
  gene_table <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    m_i <- groups$mode[i]; s_i <- groups$set[i]; o_i <- groups$ontology[i]
    x <- enriched[mode == m_i & set == s_i & ontology == o_i]
    res <- cluster_term_list(x$go_id, setNames(x$p_wallenius, x$go_id), o_i)
    key <- paste(s_i, o_i, m_i, sep = "_")
    write_similarity(res$sim, paste0("similarity_gene_level_", key))
    summary_rows[[paste0("gene_level_", key)]] <<- cbind(data.table(term_list = label_gene, mode = m_i, set = s_i, ontology = o_i), res$summary)
    n_out <- nrow(res$dropped)
    dropped_rows[[paste0("gene_level_", key)]] <<- cbind(data.table(term_list = rep(label_gene, n_out), mode = rep(m_i, n_out), set = rep(s_i, n_out),
                                                                    ontology = rep(o_i, n_out)),
                                                         res$dropped[, .(go_id, term = x$term[match(go_id, x$go_id)], reason)])
    if (nrow(res$table) >= 3) {
      log_msg("drawing figure_semantic_gene_level_%s", key)
      p_gene <- semantic_panel(res$table, sprintf("%s terms enriched at gene level: %s set, %s mode", o_i, s_i, m_i),
                               paste(panel_subtitle(res, "its term with the smallest p-value"), layout_fit(res)), max(res$table$cluster_size)) +
        labs(caption = figure_caption(" (smallest p-value)", "", "", cfg$gene_figure_size[1])) + caption_theme + theme(legend.position = "bottom")
      save_figure(p_gene, fig(sprintf("figure_semantic_gene_level_%s", key)), width = cfg$gene_figure_size[1], height = cfg$gene_figure_size[2])
    }
    log_msg("%-22s %-6s %s (%s): %d terms, %d placed, %d clusters", label_gene, s_i, o_i, m_i, res$summary$terms_in, res$summary$terms_placed,
            res$summary$clusters)
    out <- cbind(data.table(mode = rep(m_i, nrow(res$table)), set = rep(s_i, nrow(res$table)), ontology = rep(o_i, nrow(res$table))), res$table)
    out[, `:=`(q_wallenius = x$q_wallenius[match(go_id, x$go_id)], note = "")]
    out[, ..gene_level_columns]
  }))
  if (nrow(gene_table) == 0) {
    gene_table <- as.data.table(c(setNames(rep(list(NA), length(gene_level_columns) - 1), head(gene_level_columns, -1)),
                                  list(note = "Terms are enriched at gene level, but none could be placed: see semantic_terms_dropped.csv.")))
  }
}
write_table(gene_table, tab("semantic_clusters_gene_level.csv"))
}   # end of section 6

# =============================================================================
# 7. Summary
# =============================================================================
semantic_summary <- rbindlist(summary_rows)
semantic_summary[, `:=`(measure = cfg$measure, rrvgo_threshold = cfg$threshold,
                        organism_db = paste(cfg$organism_db, as.character(packageVersion(cfg$organism_db))), go_db = go_db_version,
                        go_source_date = go_source_date(), rrvgo = as.character(packageVersion("rrvgo")),
                        gosemsim = as.character(packageVersion("GOSemSim")))]
write_table(semantic_summary, tab(if (paper_part) "semantic_summary.csv" else "semantic_summary_gene_level.csv"))
dropped_terms <- rbindlist(dropped_rows)
if (ncol(dropped_terms) == 0) {                    # no term dropped from any list (the gene-level part without terms)
  dropped_terms <- data.table(term_list = character(0), mode = character(0), set = character(0), ontology = character(0),
                              go_id = character(0), term = character(0), reason = character(0))
}
setorder(dropped_terms, term_list, mode, set, ontology, reason, go_id)
write_table(dropped_terms, tab(if (paper_part) "semantic_terms_dropped.csv" else "semantic_terms_dropped_gene_level.csv"))
if (nrow(dropped_terms) != semantic_summary[, sum(terms_in - terms_placed)]) stop("The dropped terms do not add up to the summary")

print(semantic_summary[, .(term_list, set, ontology, terms_in, obsolete = dropped_obsolete, merged = dropped_merged, unknown = dropped_unknown_id,
                           other_ontology = dropped_other_ontology, not_in_gosemsim = dropped_not_in_gosemsim,
                           no_gene = dropped_no_gene_in_organism_db, placed = terms_placed, clusters,
                           clusters_other_reading = clusters_if_cut_at_similarity_of_threshold)])
print(semantic_summary[orders_tried > 0, .(term_list, set, ontology, orders_tried, another_partition = orders_with_another_partition,
                                           lowest_rand = lowest_adjusted_rand_index, clusters_from = fewest_clusters_in_another_order,
                                           clusters_to = most_clusters_in_another_order, labelled_kept = fewest_labelled_representatives_kept)])
if (gene_part) write_stamp(an, step)
finish_step(step, t0)

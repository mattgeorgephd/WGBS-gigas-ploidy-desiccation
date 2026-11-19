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
#                and names each cluster after its highest-scoring term.
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
# Inputs    04_output/09_go_enrichment_legacy/tables/legacy_fisher_<set>_<ontology>.csv, legacy_reader_loss.csv
#           04_output/08_gene_annotation/tables/dml_go_annotation_<set>_legacy.tsv.gz
#           04_output/10_go_enrichment_gene_level/tables/go_enriched_terms.csv, go_gene_level_<set>_legacy.csv
#           08_archive/code/Revigo.R, 08_archive/code/OLD/Revigo.R,
#             08_archive/figures/revigo_figure/Revigo-PLOIDY-BP.R   (plot scripts exported by REVIGO, data embedded)
#           08_archive/tables/ploidy-<ontology>-FisherTestResults.csv (Fisher tables of an earlier analysis; comparison only)
#           GO.db and the organism database of the R installation
# Outputs   04_output/11_semantic_similarity/tables/
#             semantic_clusters_legacy_<set>_<ontology>.csv    one row per term: p, score, cluster, representative, coordinates
#             semantic_clusters_as_drawn_<set>_BP.csv          the same for the lists the draft's figures were drawn from
#             semantic_clusters_gene_level.csv                 the same for the gene-level terms; one row with a note when there are none
#             similarity_<list>_<set>_<ontology>.csv.gz        the similarity matrices
#             semantic_summary.csv                             per list: terms in, terms dropped and why, clusters, checks
#             semantic_terms_dropped.csv                       the dropped terms, with the reason
#             term_lists_compared.csv                          size and overlap of the legacy and the as-drawn list
#             revigo_exports.csv                               each REVIGO export against each candidate term list
#             named_terms_results_3_4.csv                      every GO term named in Results 3.4, traced
#           04_output/11_semantic_similarity/figures/
#             figure_6_legacy, figure_7_legacy, figure_8_legacy   panel A: the list as drawn; panel B: the list the Methods describe
#             figure_revigo_export_<name>                          each distinct REVIGO export redrawn from its embedded data
#             figure_semantic_gene_level_<set>_<ontology>_<mode>   only when the gene-level analysis has enriched terms
#
# Manuscript items: Figures 6, 7 and 8; Results 3.4 (the three paragraphs on semantic
#           clustering); Methods 2.6, last sentence.
# Notes     - Which terms the draft plotted. The REVIGO export behind Figure 6
#             (Revigo-PLOIDY-BP.R) holds no p-values, and its terms are not the
#             terms below the p-value cutoff: every one of them is a BP term
#             annotated to a transcript of a gene with a ploidy DML
#             (revigo_exports.csv). For Figures 7 and 8 no export survives; the
#             terms that Results 3.4 names for them point the same way
#             (named_terms_results_3_4.csv). Each figure is therefore regenerated
#             twice: as drawn (panel A) and as described (panel B).
#           - Labels are chosen by rule, not by hand (issue I-22). In each panel the
#             config$semantic_similarity$n_labels clusters whose representative has
#             the highest score are labelled with the name of the representative.
#             Ties go to the larger cluster, then to the lower GO ID. Clusters are
#             numbered in that order, so the labelled clusters are 1 to n_labels.
#           - What the threshold means. rrvgo cuts the complete-linkage tree of
#             1 - similarity at height `threshold`, so any two terms of a cluster
#             are at least 1 - threshold similar. REVIGO's "allowed similarity"
#             works the other way: a term is removed when it is more similar than
#             that value to a term that stays. semantic_summary.csv gives the
#             number of clusters under both readings. The step stops if rrvgo's
#             clusters are not the cut described here.
#           - Terms are dropped when the GO.db of the R installation does not hold
#             them (obsolete, or merged into another term; the legacy lists come
#             from the 2021 ontology), when it files them under another ontology,
#             when GOSemSim has no entry for them, or when no gene of the organism
#             database is annotated to them or to a child, which leaves the
#             information content undefined. The organism database is human, so
#             terms used only for other taxa are lost. REVIGO carries a merged term
#             under its current ID; this step drops it and counts it.
#           - Colour marks the labelled clusters and nothing else. One hue per
#             cluster cannot be told apart when a list has dozens of clusters.
#             Cluster membership is in the tables; for the labelled clusters it is
#             drawn as lines from each term to the representative.
#           - Results do not depend on the order of the input rows: terms are
#             sorted by GO ID before clustering.
#           - The legacy enrichment is unsound (issues I-03, I-39, I-40) and no term
#             is enriched at gene level (I-41). The legacy figures are regenerated
#             for traceability, not as results.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GO.db", "GOSemSim", "rrvgo", config$semantic_similarity$organism_db, "ggplot2", "ggrepel", "patchwork"))
if (as.character(packageVersion("GO.db")) == config$legacy_go_db$version) {
  stop("GO.db ", config$legacy_go_db$version, " (the legacy release) is loaded. Run this step in a fresh R process.")
}

step <- "11_semantic_similarity"
t0 <- start_step(step, "semantic similarity of GO terms (Figures 6 to 8)")

# ---- Settings of this step that are not in run_config.R -----------------------------
# Facts about the draft and about the archive, and presentation choices.
figure_numbers  <- c(ploidy = 6L, heat = 7L, shared = 8L)   # figure of the draft for each DML set (manuscript)
figure_ontology <- "BP"                                     # the three figures show biological process only (their captions)
label_seed      <- 11L                                      # seed of the label placement (ggrepel)
label_width     <- 28L                                      # characters per line in a label
subtitle_width  <- 78L                                      # characters per line in a panel subtitle
value_tolerance <- 0.025                                    # a REVIGO value "equals" log10 p within this (p rounded to two digits)
revigo_exports  <- c(code_revigo      = "code/Revigo.R",    # plot scripts exported by REVIGO, relative to 08_archive/
                     code_old_revigo  = "code/OLD/Revigo.R",
                     revigo_ploidy_bp = "figures/revigo_figure/Revigo-PLOIDY-BP.R")
earlier_fisher  <- "tables/ploidy-%s-FisherTestResults.csv" # Fisher tables of an earlier analysis, relative to 08_archive/

cfg  <- config$semantic_similarity
sets <- c(ploidy = "ploidy-DML", heat = "heat-DML", shared = "interaction-DML")   # names used in the draft
onts <- c(BP = "P", CC = "C", MF = "F")                                           # ontology and its code in the annotation tables
tab  <- function(name) output_path(step, "tables", name)
fig  <- function(name) output_path(step, "figures", name)
go_db_version <- as.character(packageVersion("GO.db"))
log_msg("GO.db %s (ontology of %s), %s %s, rrvgo %s, GOSemSim %s; measure %s, threshold %g", go_db_version, go_source_date(),
        cfg$organism_db, as.character(packageVersion(cfg$organism_db)), as.character(packageVersion("rrvgo")),
        as.character(packageVersion("GOSemSim")), cfg$measure, cfg$threshold)
# The gene-level figures and matrices exist only when terms are enriched; remove those of an earlier run.
unlink(list.files(output_dir(step, "figures"), pattern = "^figure_semantic_gene_level_", full.names = TRUE))
unlink(list.files(output_dir(step, "tables"), pattern = "^similarity_gene_level_", full.names = TRUE))

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
  ifelse(is.na(hit), ids, hit)
}

# Legacy enrichment: every tested term of the nine Fisher tables.
legacy <- rbindlist(lapply(names(sets), function(s_i) rbindlist(lapply(names(onts), function(o_i) {
  x <- fread(output_path("09_go_enrichment_legacy", "tables", sprintf("legacy_fisher_%s_%s.csv", s_i, o_i)))
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
annotated <- rbindlist(lapply(names(sets), function(s_i) {
  f <- output_path("08_gene_annotation", "tables", sprintf("dml_go_annotation_%s_legacy.tsv.gz", s_i))
  a <- suppressWarnings(read_dml_annotation(f, "legacy"))
  a <- a[GOcat %chin% onts & grepl("^GO:[0-9]{7}$", GOID), .(annotation_rows = .N), by = .(ontology = names(onts)[match(GOcat, onts)], go_id = GOID)]
  true_names <- unique(read_dml_annotation(f, "exact")[, .(go_id = GOID, term_annotation = gsub('^"|"$', "", GOterm))])
  a[, term_annotation := true_names$term_annotation[match(go_id, true_names$go_id)]]
  cbind(data.table(set = rep(s_i, nrow(a))), a)
}))
if (anyDuplicated(annotated[, .(set, ontology, go_id)])) stop("A GO ID occurs twice in an annotated list")
reader_loss <- fread(output_path("09_go_enrichment_legacy", "tables", "legacy_reader_loss.csv"))[grepl("legacy", reader)]
for (s_i in names(sets)) for (o_i in names(onts)) {
  n_here <- annotated[set == s_i & ontology == o_i, .N]; n_step09 <- reader_loss[set == s_i][[paste0("go_ids_", o_i)]]
  if (length(n_step09) != 1 || n_here != n_step09) {
    stop(sprintf("Annotated %s terms of the %s set: %d here, %s in legacy_reader_loss.csv of step 09", o_i, s_i, n_here, paste(n_step09, collapse = "/")))
  }
}

lists_compared <- rbindlist(lapply(names(sets), function(s_i) rbindlist(lapply(names(onts), function(o_i) {
  p_list <- legacy[set == s_i & ontology == o_i & below_p_cutoff == TRUE, go_id]; a_list <- annotated[set == s_i & ontology == o_i, go_id]
  data.table(set = s_i, ontology = o_i, terms_tested_legacy = legacy[set == s_i & ontology == o_i, .N], terms_legacy_p_below_cutoff = length(p_list),
             terms_annotated_legacy_reader = length(a_list), terms_in_both = length(intersect(p_list, a_list)),
             annotated_terms_tested_legacy = sum(a_list %chin% legacy[set == s_i & ontology == o_i, go_id]))
}))))
lists_compared[, `:=`(percent_of_annotated_with_p_below_cutoff = round(100 * terms_in_both / terms_annotated_legacy_reader, 1),
                      draft_figure = fifelse(ontology == figure_ontology, paste("Figure", figure_numbers[set]), ""))]
write_table(lists_compared, tab("term_lists_compared.csv"))

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
  x <- suppressMessages(AnnotationDbi::select(getExportedValue(cfg$organism_db, cfg$organism_db), keytype = "GOALL", keys = ids, columns = "ENTREZID"))
  vapply(ids, function(g) uniqueN(x$ENTREZID[x$GOALL == g]), integer(1))
}
drop_reasons <- c(obsolete = "obsolete in GO.db", merged = "merged into another term in GO.db (secondary ID)", unknown_id = "not in GO.db",
                  other_ontology = "filed under another ontology in GO.db", not_in_gosemsim = "no entry in the term table of GOSemSim",
                  no_gene = "no gene of the organism database annotated to it or to a child")
cluster_columns <- c("go_id", "term", "p_value", "score", "cluster", "cluster_size", "is_representative", "parent_go_id", "parent_term",
                     "labelled", "term_size", "uniqueness", "mds_1", "mds_2")

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
  sim <- NULL

  if (n_kept == 0) {
    out <- data.table(go_id = character(0), term = character(0), p_value = numeric(0), score = numeric(0), cluster = integer(0),
                      cluster_size = integer(0), is_representative = logical(0), parent_go_id = character(0), parent_term = character(0),
                      labelled = logical(0), term_size = integer(0), uniqueness = numeric(0), mds_1 = numeric(0), mds_2 = numeric(0))
  } else if (n_kept == 1) {                       # nothing to compare: one cluster of one term
    size_1 <- term_gene_count(kept)
    out <- data.table(go_id = kept, term = go_terms$term[match(kept, go_terms$go_id)], score = score_of(kept, size_1), cluster_rrvgo = 1L,
                      parent_go_id = kept, term_size = size_1, uniqueness = 1, mds_1 = 0, mds_2 = 0)
    checks <- list(clusters = 1L, singletons = 1L, largest = 1L, min_within = NA_real_, complement = 1L, fit = NA_real_)
  } else {
    sim <- withCallingHandlers(
      suppressMessages(rrvgo::calculateSimMatrix(kept, orgdb = cfg$organism_db, ont = ont_code, method = cfg$measure, semdata = semantic_data(ont_code))),
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

    # Layout: classical multidimensional scaling; each axis oriented by its most extreme term.
    mds <- stats::cmdscale(stats::as.dist(1 - sim), k = min(2L, n_kept - 1L), eig = TRUE)
    pts <- cbind(mds$points, 0, 0)[, 1:2, drop = FALSE]; rownames(pts) <- kept
    for (j in 1:2) if (any(pts[, j] != 0) && pts[which.max(abs(pts[, j])), j] < 0) pts[, j] <- -pts[, j]

    out <- red[, .(go_id = go, term, score, cluster_rrvgo = cluster, parent_go_id = parent, term_size = as.integer(size), uniqueness = termUniqueness)]
    out[, `:=`(mds_1 = pts[go_id, 1], mds_2 = pts[go_id, 2])]
    checks <- list(clusters = length(members), singletons = sum(lengths(members) == 1), largest = max(lengths(members)),
                   min_within = if (all(is.na(within))) NA_real_ else min(within, na.rm = TRUE),
                   complement = uniqueN(stats::cutree(tree, h = 1 - cfg$threshold)), fit = mds$GOF[1])
  }

  if (n_kept > 0) {
    # Number the clusters by the score of the representative; ties: larger cluster, then lower GO ID.
    out[, cluster_size := .N, by = cluster_rrvgo]
    reps <- out[go_id == parent_go_id][order(-score, -cluster_size, go_id)]
    if (nrow(reps) != uniqueN(out$cluster_rrvgo)) stop("A cluster has no representative or more than one")
    reps[, cluster := seq_len(.N)]
    out[, cluster := reps$cluster[match(cluster_rrvgo, reps$cluster_rrvgo)]]
    out[, `:=`(parent_term = term[match(parent_go_id, go_id)], is_representative = go_id == parent_go_id, labelled = cluster <= cfg$n_labels,
               p_value = if (is.null(p_values)) NA_real_ else unname(p_values[go_id]),
               score = round(score, 4), uniqueness = round(uniqueness, 4), mds_1 = round(mds_1, 5), mds_2 = round(mds_2, 5))]
    if (anyNA(out[, .(term, cluster, parent_term)])) stop("A placed term lacks a name, a cluster, or a representative")
    out <- out[order(cluster, -is_representative, -score, go_id), ..cluster_columns]
  }
  n_drop <- function(key) sum(dropped$reason == drop_reasons[[key]])
  summary_row <- data.table(
    terms_in = length(ids), dropped_obsolete = n_drop("obsolete"), dropped_merged = n_drop("merged"), dropped_unknown_id = n_drop("unknown_id"),
    dropped_other_ontology = n_drop("other_ontology"), dropped_not_in_gosemsim = n_drop("not_in_gosemsim"),
    dropped_no_gene_in_organism_db = n_drop("no_gene"), terms_placed = n_kept,
    score = score_is, clusters = checks$clusters, clusters_of_one_term = checks$singletons, largest_cluster = checks$largest,
    lowest_similarity_within_a_cluster = checks$min_within, clusters_if_cut_at_similarity_of_threshold = checks$complement,
    mds_goodness_of_fit = round(checks$fit, 3))
  if (summary_row[, terms_in != dropped_obsolete + dropped_merged + dropped_unknown_id + dropped_other_ontology + dropped_not_in_gosemsim +
                    dropped_no_gene_in_organism_db + terms_placed]) {
    stop("Terms in do not equal terms dropped plus terms placed")
  }
  list(table = out, dropped = dropped, sim = sim, summary = summary_row)
}

write_similarity <- function(sim, path) {
  write_table(data.table(go_id = rownames(sim), as.data.table(sim)), path)
}

results <- list(); summary_rows <- list(); dropped_rows <- list()
run_list <- function(list_name, list_label, s_i, o_i, ids, p_values, term_names) {
  res <- cluster_term_list(ids, p_values, o_i)
  key <- paste(list_name, s_i, o_i, sep = "_")
  write_table(cbind(data.table(set = rep(s_i, nrow(res$table)), ontology = rep(o_i, nrow(res$table))), res$table),
              tab(sprintf("semantic_clusters_%s.csv", key)))
  if (is.null(res$sim)) unlink(tab(sprintf("similarity_%s.csv.gz", key))) else write_similarity(res$sim, tab(sprintf("similarity_%s.csv.gz", key)))
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
for (s_i in names(sets)) for (o_i in names(onts)) {
  x <- legacy[set == s_i & ontology == o_i & below_p_cutoff == TRUE]
  run_list("legacy", label_legacy, s_i, o_i, x$go_id, setNames(x$p_fisher, x$go_id), setNames(x$term, x$go_id))
}
for (s_i in names(sets)) {
  x <- annotated[set == s_i & ontology == figure_ontology]
  run_list("as_drawn", label_drawn, s_i, figure_ontology, x$go_id, NULL, setNames(x$term_annotation, x$go_id))
}

# =============================================================================
# 3. Figures 6, 7 and 8
# =============================================================================
grey_other <- "grey72"
wrap_label <- function(x) vapply(x, function(v) paste(strwrap(v, label_width), collapse = "\n"), character(1), USE.NAMES = FALSE)
size_breaks <- function(size_max) { b <- pretty(c(1, size_max), n = 4); sort(unique(c(1, b[b > 1 & b <= size_max]))) }

# One panel: every term a dot, every cluster a circle at its representative, the labelled clusters in colour.
semantic_panel <- function(terms, accent, title, subtitle, size_max) {
  subtitle <- paste(strwrap(subtitle, subtitle_width), collapse = "\n")
  if (nrow(terms) < 3) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "Fewer than three terms placed.\nNo layout is drawn.", size = 3.5) +
             labs(title = title, subtitle = subtitle) + theme_void(base_size = 11) +
             theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5)))
  }
  reps <- terms[is_representative == TRUE]
  rays <- merge(terms[labelled == TRUE & is_representative == FALSE], reps[, .(cluster, rep_1 = mds_1, rep_2 = mds_2)], by = "cluster")
  reps[, label_text := fifelse(labelled, wrap_label(term), "")]
  ggplot(terms, aes(mds_1, mds_2)) +
    geom_segment(data = rays, aes(xend = rep_1, yend = rep_2), colour = accent, linewidth = 0.2, alpha = 0.35) +
    geom_point(data = terms[labelled == FALSE], colour = grey_other, size = 0.7) +
    geom_point(data = terms[labelled == TRUE], colour = accent, size = 0.8) +
    geom_point(data = reps[labelled == FALSE], aes(size = cluster_size), shape = 21, fill = grey_other, colour = "white", stroke = 0.5) +
    geom_point(data = reps[labelled == TRUE], aes(size = cluster_size), shape = 21, fill = accent, colour = "white", stroke = 0.6) +
    ggrepel::geom_label_repel(data = reps, aes(label = label_text), size = 2.5, lineheight = 0.9, colour = "black", fill = scales::alpha("white", 0.9),
                              label.size = 0.2, label.padding = 0.15, box.padding = 0.5, point.padding = 0.3, min.segment.length = 0,
                              segment.colour = "grey25", segment.size = 0.3, force = 3, max.overlaps = Inf, max.iter = 20000, max.time = 30,
                              seed = label_seed) +
    scale_size_area(max_size = 9, limits = c(0, size_max), breaks = size_breaks(size_max), name = "Terms in the cluster") +
    scale_x_continuous(expand = expansion(mult = 0.1)) + scale_y_continuous(expand = expansion(mult = 0.1)) +
    coord_equal() +                                                     # a map of distances: one scale on both axes
    guides(size = guide_legend(override.aes = list(fill = "grey45"))) +
    labs(x = "MDS axis 1", y = "MDS axis 2", title = title, subtitle = subtitle) +
    theme_wgbs(base_size = 11) + theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5))
}
panel_subtitle <- function(res, representative_rule) {
  sprintf("%s terms in the list, %s placed, %d clusters. Each cluster is named after %s.", format(res$summary$terms_in, big.mark = ","),
          format(res$summary$terms_placed, big.mark = ","), res$summary$clusters, representative_rule)
}
layout_fit <- function(res) {
  if (is.na(res$summary$mds_goodness_of_fit)) "" else sprintf("Goodness of fit of the layout in two dimensions: %.2f.", res$summary$mds_goodness_of_fit)
}
figure_caption <- function(score_note, closing_note) paste0(
  "Each dot is a GO term. Each circle is a cluster of terms, drawn at its representative term; its area shows the number of terms.\n",
  sprintf("Colour: the %d clusters whose representative scores highest%s. ", cfg$n_labels, score_note),
  "They are labelled, and lines join their terms to the representative. Grey: all other clusters.\n",
  sprintf("Layout: classical multidimensional scaling of 1 - similarity (%s measure, %s). Clusters: rrvgo, threshold %g.%s", cfg$measure, cfg$organism_db,
          cfg$threshold, closing_note))
for (s_i in names(sets)) {
  drawn <- results[[paste("as_drawn", s_i, figure_ontology, sep = "_")]]; described <- results[[paste("legacy", s_i, figure_ontology, sep = "_")]]
  both <- lists_compared[set == s_i & ontology == figure_ontology]
  size_max <- max(1L, drawn$table$cluster_size, described$table$cluster_size)
  accent <- config$colours$dml_set[[s_i]]
  panel_a <- semantic_panel(drawn$table, accent, sprintf("A  As drawn: all %s terms annotated to genes with a DML of this set", figure_ontology),
                            paste(panel_subtitle(drawn, "its term with the most annotated genes"),
                                  sprintf("%d of the %s terms of the list have p < %g in the legacy test.", both$terms_in_both,
                                          format(both$terms_annotated_legacy_reader, big.mark = ","), config$go_pvalue), layout_fit(drawn)), size_max)
  panel_b <- semantic_panel(described$table, accent, sprintf("B  As described: %s terms with p < %g in the legacy test", figure_ontology, config$go_pvalue),
                            paste(panel_subtitle(described, "its term with the smallest p-value"),
                                  sprintf("%d of the %d terms of the list are also in the list of panel A.", both$terms_in_both, both$terms_legacy_p_below_cutoff),
                                  layout_fit(described)), size_max)
  figure <- (panel_a | panel_b) + plot_layout(guides = "collect") +
    plot_annotation(title = sprintf("Figure %d of the draft, regenerated in R: semantic similarity of GO terms, %s", figure_numbers[[s_i]], sets[[s_i]]),
                    caption = figure_caption(" (panel A: most annotated genes; panel B: smallest p-value)",
                                             " The legacy enrichment is unsound (issues I-03, I-39, I-40); this figure is a reproduction, not a result."),
                    theme = theme(plot.title = element_text(size = 12, face = "bold"), plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.1))) &
    theme(legend.position = "bottom")
  save_figure(figure, fig(sprintf("figure_%d_legacy", figure_numbers[[s_i]])), width = 15, height = 8.6)
}

# =============================================================================
# 4. The REVIGO exports that survive
# =============================================================================
# An export is an R script with one c("GO:...", "name", numbers) record per term.
read_revigo_export <- function(path) {
  if (!file.exists(path)) stop("REVIGO export not found: ", path)
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
# A Fisher table of the earlier analysis: no quoting, commas inside term names, p as printed. A record that is
# cut short (the BP file holds runs of NUL bytes) does not fit the pattern and is not read.
read_earlier_fisher <- function(path) {
  if (!file.exists(path)) stop("Archived table not found: ", path)
  txt <- grep("^GO:[0-9]{7},", readLines(path, warn = FALSE), value = TRUE)
  hit <- regmatches(txt, regexec("^(GO:[0-9]{7}),.*,[0-9]+,[0-9]+,[0-9.]+,(< ?)?([0-9.eE+-]+),?$", txt))
  good <- lengths(hit) == 4
  d <- data.table(go_id = vapply(hit[good], `[`, "", 2), p = as.numeric(vapply(hit[good], `[`, "", 4)))
  if (nrow(d) == 0 || anyNA(d$p) || anyDuplicated(d$go_id)) stop("Unexpected format in ", path)
  list(data = d, lines_not_read = sum(!good))
}

exports <- lapply(revigo_exports, function(f) read_revigo_export(archive_path(f)))
candidates <- list()
add_candidate <- function(list_name, source_file, ids, p = NULL) {
  candidates[[length(candidates) + 1]] <<- list(list_name = list_name, source_file = source_file, ids = ids, p = p)
}
# Candidate term lists: the nine legacy lists, the nine annotated lists, and the lists of the earlier analysis.
for (s_i in names(sets)) for (o_i in names(onts)) {
  x <- legacy[set == s_i & ontology == o_i]
  add_candidate(sprintf("legacy test, p < %g: %s, %s", config$go_pvalue, sets[[s_i]], o_i),
                sprintf("04_output/09_go_enrichment_legacy/tables/legacy_fisher_%s_%s.csv", s_i, o_i),
                x[below_p_cutoff == TRUE, go_id], setNames(x$p_fisher, x$go_id))
}
for (s_i in names(sets)) for (o_i in names(onts)) {
  add_candidate(sprintf("all terms annotated to genes with a DML (legacy reader): %s, %s", sets[[s_i]], o_i),
                sprintf("04_output/08_gene_annotation/tables/dml_go_annotation_%s_legacy.tsv.gz", s_i), annotated[set == s_i & ontology == o_i, go_id])
}
for (o_i in names(onts)) {
  x <- read_earlier_fisher(archive_path(sprintf(earlier_fisher, o_i)))
  if (x$lines_not_read > 0) log_msg("%s: %d damaged line(s) not read", sprintf(earlier_fisher, o_i), x$lines_not_read)
  add_candidate(sprintf("earlier analysis, not in the manuscript, p < %g: ploidy, %s", config$go_pvalue, o_i),
                file.path("08_archive", sprintf(earlier_fisher, o_i)), x$data[p < config$go_pvalue, go_id], setNames(x$data$p, x$data$go_id))
}
# REVIGO gives a term under its current ID, so a list term can be in an export under the ID it was merged into.
export_path <- function(e_i) file.path("08_archive", revigo_exports[[e_i]])
first_with_md5 <- vapply(exports, function(e) names(exports)[match(e$md5, vapply(exports, `[[`, "", "md5"))], "")
export_rows <- rbindlist(lapply(names(exports), function(e_i) {
  d <- exports[[e_i]]$data; ids <- d$term_ID; has_p <- "value" %in% names(d)
  rbindlist(lapply(candidates, function(cand) {
    same_id <- ids %chin% cand$ids
    by_secondary <- !same_id & ids %chin% primary_id(cand$ids)
    same_p <- if (is.null(cand$p) || !has_p) NA_integer_ else sum(abs(d$value - log10(cand$p[ids])) < value_tolerance, na.rm = TRUE)
    data.table(export = export_path(e_i), same_content_as = if (first_with_md5[[e_i]] == e_i) "" else export_path(first_with_md5[[e_i]]),
               export_terms = length(ids), export_has_p_values = has_p, term_list = cand$list_name, term_list_source = cand$source_file,
               list_terms = length(cand$ids), export_terms_in_list = sum(same_id), export_terms_in_list_by_secondary_id = sum(by_secondary),
               list_terms_in_export = sum(primary_id(cand$ids) %chin% ids | cand$ids %chin% ids), export_values_equal_to_log10_p = same_p)
  }))
}))
export_rows[, percent_of_export_in_list := round(100 * (export_terms_in_list + export_terms_in_list_by_secondary_id) / export_terms, 1)]
export_rows[, best_match := seq_len(.N) == order(-(export_terms_in_list + export_terms_in_list_by_secondary_id), list_terms)[1], by = export]
write_table(export_rows, tab("revigo_exports.csv"))
print(export_rows[best_match == TRUE, .(export, export_terms, term_list, list_terms, export_terms_in_list, export_terms_in_list_by_secondary_id, list_terms_in_export)])

# Redraw each distinct export from its own data, with the label rule of its own plot code.
revigo_low <- "#b30000"; revigo_high <- "#f7d038"          # red for low values, yellow for high, as in the archived REVIGO screenshots
for (e_i in names(exports)) {
  if (first_with_md5[[e_i]] != e_i) {
    log_msg("%s has the same content as %s; not drawn twice", revigo_exports[[e_i]], revigo_exports[[first_with_md5[[e_i]]]])
    next
  }
  d <- copy(exports[[e_i]]$data); label_cut <- exports[[e_i]]$label_below_dispensability
  has_p <- "value" %in% names(d)
  d[, `:=`(shade = if (has_p) value else uniqueness, label_text = fifelse(dispensability < label_cut, wrap_label(description), ""))]
  best <- export_rows[export == export_path(e_i) & best_match == TRUE]
  n_secondary <- best$export_terms_in_list_by_secondary_id
  p_export <- ggplot(d[order(-log_size)], aes(plot_X, plot_Y)) +
    geom_point(aes(size = log_size, fill = shade), shape = 21, colour = "white", stroke = 0.4, alpha = 0.9) +
    ggrepel::geom_label_repel(aes(label = label_text), size = 2.4, lineheight = 0.9, colour = "black", fill = scales::alpha("white", 0.9), label.size = 0.2,
                              label.padding = 0.15, box.padding = 0.45, point.padding = 0.2, min.segment.length = 0, segment.colour = "grey25",
                              segment.size = 0.3, force = 2, max.overlaps = Inf, max.iter = 20000, max.time = 30, seed = label_seed) +
    scale_fill_gradient(low = revigo_low, high = revigo_high, name = if (has_p) "Value in the export\n(log10 p)" else "Uniqueness\n(REVIGO)") +
    scale_size(range = c(1.2, 9), name = "Log size\n(REVIGO)") +
    scale_x_continuous(expand = expansion(mult = 0.08)) + scale_y_continuous(expand = expansion(mult = 0.06)) +
    guides(size = guide_legend(order = 1, override.aes = list(fill = "grey45")), fill = guide_colourbar(order = 2)) +
    labs(x = "Semantic space X (REVIGO)", y = "Semantic space Y (REVIGO)",
         title = sprintf("REVIGO export %s, redrawn from its embedded data", export_path(e_i)),
         subtitle = paste0(sprintf("%d terms. %s. Labels: terms with dispensability below %g, the rule of the export's own plot code.\n", nrow(d),
                                   if (has_p) "Colour: the value column of the export" else "The export holds no p-values; colour is REVIGO's uniqueness", label_cut),
                           sprintf("Closest term list: %s (%s terms). It holds %d of the export's terms%s.", best$term_list, format(best$list_terms, big.mark = ","),
                                   best$export_terms_in_list + n_secondary,
                                   if (n_secondary > 0) sprintf(", %d of them under a secondary ID", n_secondary) else ""))) +
    theme_wgbs(base_size = 11) + theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.5))
  save_figure(p_export, fig(paste0("figure_revigo_export_", e_i)), width = 10.5, height = 8)
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
                                      "physiological responses to osmotic and oxidative stress", "physiological responses to osmotic and oxidative stress",
                                      "mRNA methylation", "mRNA export from the nucleus in response to heat stress", "rRNA catabolic processes"),
             read_as = c("MAPK cascade", "stress-activated MAPK cascade", NA, NA, "response to osmotic stress", "response to oxidative stress", NA,
                         "mRNA export from nucleus in response to heat stress", "rRNA catabolic process"),
             reading_note = c("article dropped", "the label on the archived image of Figure 7 (08_archive/figures/revigo_figure/HEAT-BP.png)", NA, NA,
                              "one phrase for two terms", "one phrase for two terms", NA, "article dropped", "plural")),
  data.table(set = "shared", phrase = c("tRNA modification", "RNA 3'-end processing", "post-translational protein modification", "calcium-mediated signaling",
                                        "ER overload response", "response to osmotic stress", "positive regulation of gene expression", "gene silencing by RNA"))),
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
named <- merge(named, legacy[, .(set, go_id, legacy_rank = rank, legacy_p = p_fisher, in_legacy_p_list = below_p_cutoff, transcripts_annotated = annotated,
                                 transcripts_of_interest = significant, genes_behind_transcripts_of_interest = genes_significant)],
               by = c("set", "go_id"), all.x = TRUE, sort = FALSE)
named[is.na(in_legacy_p_list), in_legacy_p_list := FALSE]
drawn_ids <- function(s_i) annotated[set == s_i & ontology == figure_ontology, go_id]
named[, in_list_as_drawn := unname(mapply(function(s_i, g) {
  if (is.na(g)) return("no")
  ids <- drawn_ids(s_i); via <- ids[primary_id(ids) == g & ids != g]
  if (g %chin% ids) "yes" else if (length(via) > 0) paste0("under a secondary ID (", paste(via, collapse = ", "), ")") else "no"
}, set, go_id))]
drawn_clusters <- rbindlist(lapply(names(sets), function(s_i) {       # where the term sits in panel A of the regenerated figure
  x <- results[[paste("as_drawn", s_i, figure_ontology, sep = "_")]]$table
  cbind(data.table(set = rep(s_i, nrow(x))), x[, .(go_id, as_drawn_cluster = cluster, as_drawn_cluster_named_after = parent_term)])
}))
named <- merge(named, drawn_clusters, by = c("set", "go_id"), all.x = TRUE, sort = FALSE)
gene_level <- rbindlist(lapply(names(sets), function(s_i) {
  fread(output_path("10_go_enrichment_gene_level", "tables", sprintf("go_gene_level_%s_legacy.csv", s_i)))[universe == "eligible"]
}))
if (anyDuplicated(gene_level[, .(set, go_id)])) stop("A GO ID occurs twice for one set in the gene-level tables")
named <- merge(named, gene_level[, .(set, go_id, gene_level_genes_annotated = genes_annotated, gene_level_genes_with_dml = genes_with_dml,
                                     gene_level_p_wallenius = p_wallenius, gene_level_q_wallenius = q_wallenius, gene_level_tested = tested)],
               by = c("set", "go_id"), all.x = TRUE, sort = FALSE)
if (nrow(named) != n_named) stop("A join changed the number of named terms")
named[, gene_level_result := fifelse(is.na(go_id), "no GO term matched",
                             fifelse(!go_id %chin% go_terms$go_id, paste("not tested: the term is not in GO.db", go_db_version),
                             fifelse(is.na(gene_level_genes_annotated), "not tested: no eligible gene carries the term",
                             fifelse(!gene_level_tested, sprintf("not tested: fewer than %d eligible genes carry the term", config$go_min_genes),
                             fifelse(gene_level_q_wallenius < config$go_fdr, "enriched", "not enriched")))))]
named[, `:=`(manuscript_set = unname(sets[set]), draft_figure = paste("Figure", figure_numbers[set]))]
named_out <- named[order(row_in_text),
                   .(set, manuscript_set, draft_figure, phrase, read_as, match_type, reading_note, name_matched_in, go_id, term_current,
                     term_legacy_release, ontology, in_legacy_p_list, legacy_rank, legacy_p, transcripts_annotated, transcripts_of_interest,
                     genes_behind_transcripts_of_interest, in_list_as_drawn, as_drawn_cluster, as_drawn_cluster_named_after,
                     gene_level_genes_annotated, gene_level_genes_with_dml, gene_level_p_wallenius, gene_level_q_wallenius, gene_level_result)]
write_table(named_out, tab("named_terms_results_3_4.csv"))
print(named_out[, .(named_terms = .N, matched = sum(match_type != "not found"), in_legacy_p_list = sum(in_legacy_p_list),
                    in_list_as_drawn = sum(in_list_as_drawn != "no"), enriched_at_gene_level = sum(gene_level_result == "enriched")), by = set])

# =============================================================================
# 6. The same procedure on the gene-level terms
# =============================================================================
enriched <- fread(output_path("10_go_enrichment_gene_level", "tables", "go_enriched_terms.csv"))
need <- c("mode", "set", "ontology", "go_id", "term", "p_wallenius", "q_wallenius")
if (!all(need %in% names(enriched))) stop("go_enriched_terms.csv lacks column(s): ", paste(setdiff(need, names(enriched)), collapse = ", "))
gene_level_columns <- c("mode", "set", "ontology", cluster_columns, "q_wallenius", "note")
label_gene <- sprintf("gene level, q < %g", config$go_fdr)
if (nrow(enriched) == 0) {
  no_terms <- sprintf(paste("No GO term is enriched at gene level (q < %g, eligible universe, Wallenius test) in any DML set, ontology or mode.",
                            "Nothing was clustered and no figure was drawn."), config$go_fdr)
  gene_table <- as.data.table(c(setNames(rep(list(NA), length(gene_level_columns) - 1), head(gene_level_columns, -1)), list(note = no_terms)))
  summary_rows[["gene_level"]] <- cbind(data.table(term_list = label_gene, mode = "all", set = "all", ontology = "all"),
                                        cluster_term_list(character(0), numeric(0), figure_ontology)$summary)
  log_msg("%s", no_terms)
} else {
  groups <- unique(enriched[, .(mode, set, ontology)])
  gene_table <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    m_i <- groups$mode[i]; s_i <- groups$set[i]; o_i <- groups$ontology[i]
    x <- enriched[mode == m_i & set == s_i & ontology == o_i]
    res <- cluster_term_list(x$go_id, setNames(x$p_wallenius, x$go_id), o_i)
    key <- paste(s_i, o_i, m_i, sep = "_")
    if (!is.null(res$sim)) write_similarity(res$sim, tab(sprintf("similarity_gene_level_%s.csv.gz", key)))
    summary_rows[[paste0("gene_level_", key)]] <<- cbind(data.table(term_list = label_gene, mode = m_i, set = s_i, ontology = o_i), res$summary)
    dropped_rows[[paste0("gene_level_", key)]] <<- cbind(data.table(term_list = rep(label_gene, nrow(res$dropped)), mode = rep(m_i, nrow(res$dropped)),
                                                                    set = rep(s_i, nrow(res$dropped)), ontology = rep(o_i, nrow(res$dropped))),
                                                         res$dropped[, .(go_id, term = x$term[match(go_id, x$go_id)], reason)])
    if (nrow(res$table) >= 3) {
      p_gene <- semantic_panel(res$table, config$colours$dml_set[[s_i]], sprintf("%s terms enriched at gene level: %s set, %s mode", o_i, s_i, m_i),
                               paste(panel_subtitle(res, "its term with the smallest p-value"), layout_fit(res)), max(res$table$cluster_size)) +
        labs(caption = figure_caption(" (smallest p-value)", "")) +
        theme(legend.position = "bottom", plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.1))
      save_figure(p_gene, fig(sprintf("figure_semantic_gene_level_%s", key)), width = 9, height = 8.6)
    }
    log_msg("%-22s %-6s %s (%s): %d terms, %d placed, %d clusters", label_gene, s_i, o_i, m_i, res$summary$terms_in, res$summary$terms_placed, res$summary$clusters)
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

# =============================================================================
# 7. Summary
# =============================================================================
semantic_summary <- rbindlist(summary_rows)
semantic_summary[, `:=`(measure = cfg$measure, rrvgo_threshold = cfg$threshold, organism_db = paste(cfg$organism_db, as.character(packageVersion(cfg$organism_db))),
                        go_db = go_db_version, go_source_date = go_source_date(), rrvgo = as.character(packageVersion("rrvgo")),
                        gosemsim = as.character(packageVersion("GOSemSim")))]
write_table(semantic_summary, tab("semantic_summary.csv"))
dropped_terms <- rbindlist(dropped_rows)
setorder(dropped_terms, term_list, mode, set, ontology, reason, go_id)
write_table(dropped_terms, tab("semantic_terms_dropped.csv"))
if (nrow(dropped_terms) != semantic_summary[, sum(terms_in - terms_placed)]) stop("The dropped terms do not add up to the summary")

print(semantic_summary[, .(term_list, set, ontology, terms_in, obsolete = dropped_obsolete, merged = dropped_merged, unknown = dropped_unknown_id,
                           other_ontology = dropped_other_ontology, not_in_gosemsim = dropped_not_in_gosemsim, no_gene = dropped_no_gene_in_organism_db,
                           placed = terms_placed, clusters, clusters_other_reading = clusters_if_cut_at_similarity_of_threshold)])
finish_step(step, t0)

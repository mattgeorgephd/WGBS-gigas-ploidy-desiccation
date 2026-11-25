# 10_go_enrichment_gene_level.R
#
# Purpose   GO enrichment of the genes that hold DML, done at the level of genes
#           and against the genes that hold an eligible CpG. This is the analysis
#           the audit recommends in place of step 09; it changes five things:
#             1. Unit. One gene counts once, whatever its number of transcripts
#                (issue I-03). A gene carries the union of its mRNAs' GO terms.
#             2. Universe. Two are used and reported side by side:
#                  tested    protein-coding genes with GO terms and at least one
#                            CpG that was tested for differential methylation
#                  eligible  those of them with at least one eligible CpG
#                            (mean methylation from half the threshold of the
#                            analysis to 100 minus that: 10% to 90% at the
#                            paper's 20 points, 25% to 75% at the 50 points of
#                            option d50; a screen and not a bound; see
#                            dml_eligible() and step 08). Gene-body methylation in this
#                            genome sits in a minority of genes, mostly broadly
#                            expressed ones, so against the tested universe any
#                            DML list looks enriched for housekeeping functions.
#                            The eligible universe is the primary one.
#             3. Every line of the annotation is read (issue I-40) and the
#                universe is complete (issue I-39).
#             4. Multiple testing. Benjamini-Hochberg across the terms of an
#                ontology with at least config$go_min_genes annotated genes; a
#                term is called enriched at q < config$go_fdr (issue I-21).
#             5. Bias from the number of CpGs. A gene with more eligible CpGs is
#                more likely to hold a DML. Next to the classic Fisher test the
#                step reports the Wallenius test of goseq (Young et al. 2010,
#                Genome Biology 11:R14), with the gene's CpG count as the bias
#                variable, as recommended for methylation data (Geeleher et al.
#                2013, Bioinformatics 29:1851).
#
# Run       NOT ADOPTED in the paper: an optional step, run by
#             Rscript run_pipeline.R --gene-level-go [with the switches of an analysis]
#           or alone with WGBS_GENE_LEVEL_GO=TRUE (and WGBS_ANALYSIS=<id>); in any other
#           configuration it ends at once with status 0. Needs steps 08 and 09 of the
#           paper (the draft's terms) and step 08 of the analysis. About 3 minutes per
#           analysis, about 1.6 GB of memory.
#
# Inputs    04_output/08_gene_annotation/tables/gene_go_map.csv.gz (the paper run)
#           gene_table_<id>.csv.gz of step 08 of the analysis (paper folder for the paper's
#             analysis, 04_output/not_adopted/<id>/08_gene_annotation/ for another)
#           04_output/09_go_enrichment_legacy/tables/legacy_fisher_<set>_<ontology>.csv (the paper run)
# Outputs   04_output/not_adopted/<id>/10_go_enrichment_gene_level/tables/ (the paper's
#           analysis: <id> = legacy)
#             go_gene_level_<set>_<id>.csv          every term, both universes, both tests
#             go_gene_level_summary.csv             terms tested and enriched (rows of the analysis)
#             go_enriched_terms.csv                 terms with q < go_fdr (eligible universe, Wallenius)
#             go_enriched_term_genes.csv            the DML genes behind those terms
#             legacy_top10_at_gene_level.csv        each "enriched" term of the draft, retested (for the paper's
#                                                   analysis in both universes; for another, eligible universe,
#                                                   columns ending in _eligible_<id>_dml)
#             dml_probability_by_cpg_count.csv      the bias the Wallenius test adjusts for
#             go_terms_by_analysis.csv              terms called enriched by the analysis (for the paper's
#                                                   analysis also by the draft's procedure)
#           .../figures/figure_go_bias_<id>, figure_go_enriched_<id> (when any term is enriched)
#           04_output/not_adopted/comparisons/10_go_enrichment_gene_level/: the last six tables of
#             every analysis on disk side by side, the paper's first (legacy_top10_at_gene_level.csv
#             side by side in columns), and figures/figure_go_terms_by_analysis
#           .../session_info.txt
#
# Manuscript items: replacement candidates for Results 3.4 and 3.5, Figures 5 to 9, Table 1.
# Notes     - The GO annotation is the frozen one of the legacy analysis (blastx of
#             2021 against Swiss-Prot, GO terms from a 2017 UniProt table). The
#             ontology is the GO.db of the R installation, recorded in
#             go_gene_level_summary.csv; GO IDs it no longer holds are dropped.
#           - Known defect, left as it is by the scope decision of 2026-11-19: the
#             830 GO IDs of the annotation that are not primary IDs of the GO.db
#             in use are dropped without a message, 395 of them secondary IDs of
#             current terms, and with those replaced by their primary IDs the
#             step calls one term (organelle; heat set, CC) where it now calls
#             none (issues.md, I-41).
#           - With option inter (not adopted) the tested interaction list is a fourth
#             set, "interaction" (step 08 annotates it). It is tested in a pass of its
#             own, whose eligible universe adds the genes of its DML; the three sets of
#             the paper and their universe are those of the analysis without inter. The
#             draft's terms (legacy_top10_at_gene_level.csv) concern the three sets only.
#             The interaction rests mostly on five libraries with a raised overall
#             methylation level (issues I-02, I-49): its genes describe those libraries
#             as much as a response to heat.
#           - The heat contrast is confounded with cohort and sequencing batch
#             (issue I-01). Enrichment of "heat" DML genes describes genes that
#             differ between the two cohorts, whatever the cause.
#           - A DML is assigned to a gene when it lies in the gene body.
#           - A set none of whose DML genes carries a term of the ontology (a set
#             without DML, or whose DML lie in genes without GO terms; at 50
#             points, option d50, the ploidy set and the loci in both lists) is
#             not tested in that ontology and universe: the log says so, its
#             table go_gene_level_<set>_<id>.csv has no row for them, and its
#             summary row has dml_genes 0, no term tested and no top term; in
#             legacy_top10_at_gene_level.csv every term of it has the status
#             "set not tested: no DML gene carries a term of the ontology".

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "topGO", "goseq", "ggplot2", "patchwork"))
if (as.character(packageVersion("GO.db")) == config$legacy_go_db$version) {
  stop("GO.db ", config$legacy_go_db$version, " (the legacy release) is loaded. Run this step in a fresh R process.")
}

step <- "10_go_enrichment_gene_level"
t0 <- start_step(step, "gene-level GO enrichment against tested and eligible genes (NOT ADOPTED)")
require_run_kind("optional_gene_level", "step 10 (gene-level GO enrichment, not adopted in the paper) runs only when asked for")
an <- config$analysis
check_stamp(an, needs = c("04", "05", "08"))
tab <- function(name) analysis_path(an, step, "tables", name)
ann_in <- function(name) analysis_input(an, "08_gene_annotation", "tables", name)
set_labels <- c(ploidy = "ploidy DML", heat = "heat DML", shared = "in both lists", interaction = "tested interaction DML")
paper_sets <- c("ploidy", "heat", "shared")
# Option inter (not adopted; author's request of 2026-11-25): the tested interaction list,
# annotated by step 08 as a fourth set, is tested in a pass of its own (below), so that
# the three sets of the paper keep exactly the universe they have without inter.
inter_set <- identical(an$contrast_variant, "inter")
sets <- set_labels[c(paper_sets, if (inter_set) "interaction")]
onts <- c("BP", "CC", "MF")
universes <- c(tested = "n_cpg_tested", eligible = "n_cpg_eligible")
log_msg("GO.db %s (ontology of %s), topGO %s, goseq %s", as.character(packageVersion("GO.db")), go_source_date(),
        as.character(packageVersion("topGO")), as.character(packageVersion("goseq")))

gene_go <- fread(paper_path("08_gene_annotation", "tables", "gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
gene2go_all <- split(gene_go$go_id, gene_go$gene_id)

results <- list(); summary_rows <- list(); bias_rows <- list(); term_genes <- list(); not_tested <- list()
for (mode in an$id) {
  gene_tab <- fread(ann_in(sprintf("gene_table_%s.csv.gz", mode)), colClasses = c(gene_id = "character"))
  gene_tab <- gene_tab[biotype == "protein_coding"]
  any_dml_paper <- gene_tab$n_dml_ploidy + gene_tab$n_dml_heat + gene_tab$n_dml_shared > 0
  # Passes: the paper's three sets; with inter, the interaction set on its own, its
  # eligible universe adding the genes of its own DML (many interaction DML lie at
  # CpGs whose mean is outside the eligible range).
  passes <- list(list(sets = paper_sets, any_dml = any_dml_paper))
  if (inter_set) {
    if (!"n_dml_interaction" %in% names(gene_tab)) stop("gene_table_", mode, ".csv.gz has no interaction set: run step 08 of the analysis again.", call. = FALSE)
    passes[[2]] <- list(sets = "interaction", any_dml = any_dml_paper | gene_tab$n_dml_interaction > 0)
  }
  for (pass in passes) {
  pass_sets <- pass$sets; any_dml <- pass$any_dml

  for (u in names(universes)) {
    # Genes of the universe. A gene that holds a DML belongs to the eligible
    # universe by definition; a few DML can have a mean just outside the eligible
    # range (2 of the paper's DML lie outside 10% to 90%), and their genes are
    # kept with a count of one.
    n_cpg <- gene_tab[[universes[[u]]]]
    if (u == "eligible") n_cpg <- pmax(n_cpg, as.integer(any_dml))
    uni <- cbind(gene_tab[n_cpg > 0], bias = n_cpg[n_cpg > 0])
    annotated_ids <- uni[n_go_terms > 0, gene_id]

    # Probability that a gene holds a DML, by number of CpGs (all protein-coding genes of the universe).
    uni[, cpg_bin := cut(bias, c(0, 1, 2, 5, 10, 20, 50, 100, Inf), labels = c("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"))]
    for (s in pass_sets) {
      bias_rows[[length(bias_rows) + 1]] <- uni[, .(mode = mode, universe = u, set = s, genes = .N,
        genes_with_dml = sum(get(paste0("n_dml_", s)) > 0), percent_with_dml = 100 * mean(get(paste0("n_dml_", s)) > 0)),
        by = cpg_bin][order(cpg_bin)]
    }

    for (ont in onts) {
      td <- go_gene_level_object(ont, annotated_ids, gene2go_all)
      for (s in pass_sets) {
        dml_genes <- uni[get(paste0("n_dml_", s)) > 0, gene_id]
        if (!any(topGO::genes(td) %in% dml_genes)) {
          log_msg("[%s] %-8s %s %-6s: no gene that holds a DML of the set carries a term of the ontology; not tested", mode, u, ont, s)
          not_tested[[length(not_tested) + 1]] <- data.table(mode = mode, set = s, ontology = ont, universe = u)
          summary_rows[[length(summary_rows) + 1]] <- data.table(
            mode = mode, set = s, ontology = ont, universe = u, universe_genes = numGenes(td), dml_genes = 0L, terms_tested = 0L,
            fisher_p_below_cutoff = 0L, fisher_q_below_fdr = 0L, wallenius_p_below_cutoff = 0L, wallenius_q_below_fdr = 0L,
            smallest_q_wallenius = NA_real_, top_term = NA_character_)
          next
        }
        fit <- go_gene_level_test(td, uni, dml_genes, config)
        res <- cbind(data.table(mode = mode, set = s, ontology = ont, universe = u), fit$table)
        results[[length(results) + 1]] <- res

        summary_rows[[length(summary_rows) + 1]] <- data.table(
          mode = mode, set = s, ontology = ont, universe = u,
          universe_genes = numGenes(fit$td), dml_genes = numSigGenes(fit$td), terms_tested = sum(res$tested),
          fisher_p_below_cutoff = res[tested == TRUE, sum(p_fisher < config$go_pvalue)],
          fisher_q_below_fdr = res[tested == TRUE, sum(q_fisher < config$go_fdr)],
          wallenius_p_below_cutoff = res[tested == TRUE, sum(p_wallenius < config$go_pvalue)],
          wallenius_q_below_fdr = res[tested == TRUE, sum(q_wallenius < config$go_fdr)],
          smallest_q_wallenius = res[tested == TRUE, min(q_wallenius)],
          top_term = res[tested == TRUE][1, term])

        hit <- res[tested == TRUE & q_wallenius < config$go_fdr, go_id]
        if (u == "eligible" && length(hit) > 0) {
          term_genes[[length(term_genes) + 1]] <- rbindlist(lapply(hit, function(g) {
            data.table(mode = mode, set = s, ontology = ont, go_id = g, gene_id = intersect(fit$members[[g]], dml_genes))
          }))
        }
      }
      log_msg("[%s] %-8s %s: %5d genes; terms with q < %g (Wallenius): %s", mode, u, ont, length(annotated_ids), config$go_fdr,
              paste(vapply(tail(summary_rows, length(pass_sets)), function(r) sprintf("%s %d", r$set, r$wallenius_q_below_fdr), ""), collapse = ", "))
    }
  }
  }
}
# With no set tested at all, the table of results keeps its columns.
results <- if (length(results)) rbindlist(results) else
  data.table(mode = character(), set = character(), ontology = character(), universe = character(), go_id = character(),
             term = character(), genes_annotated = integer(), genes_with_dml = integer(), expected = numeric(), p_fisher = numeric(),
             p_wallenius = numeric(), tested = logical(), q_fisher = numeric(), q_wallenius = numeric(), enrichment = numeric())
go_summary <- rbindlist(summary_rows)
go_summary[, `:=`(go_db = as.character(packageVersion("GO.db")), go_source_date = go_source_date())]
for (mode_i in an$id) for (s in names(sets)) {
  write_table(results[mode == mode_i & set == s], tab(sprintf("go_gene_level_%s_%s.csv", s, mode_i)))
}
write_table(go_summary, tab("go_gene_level_summary.csv"))
write_table(rbindlist(bias_rows), tab("dml_probability_by_cpg_count.csv"))

enriched <- results[universe == "eligible" & tested == TRUE & q_wallenius < config$go_fdr]
setorder(enriched, mode, set, ontology, q_wallenius)
write_table(enriched, tab("go_enriched_terms.csv"))
term_genes <- if (length(term_genes) > 0) rbindlist(term_genes) else
  data.table(mode = character(0), set = character(0), ontology = character(0), go_id = character(0), gene_id = character(0))
gene_names <- unique(rbindlist(lapply(an$id, function(m) fread(ann_in(sprintf("gene_table_%s.csv.gz", m)),
                                                                       colClasses = c(gene_id = "character"))[, .(gene_id, gene, product)])))
term_genes <- merge(term_genes, gene_names, by = "gene_id", all.x = TRUE, sort = FALSE)
term_genes <- merge(term_genes, unique(results[, .(go_id, term)]), by = "go_id", all.x = TRUE, sort = FALSE)
setcolorder(term_genes, c("mode", "set", "ontology", "go_id", "term", "gene_id", "gene", "product"))
setorder(term_genes, mode, set, ontology, go_id, gene_id)
write_table(term_genes, tab("go_enriched_term_genes.csv"))

# =============================================================================
# The "enriched" terms of the draft, retested at gene level
# =============================================================================
legacy_dir <- function(name) paper_path("09_go_enrichment_legacy", "tables", name)     # the draft's terms
legacy_top <- rbindlist(lapply(paper_sets, function(s) rbindlist(lapply(onts, function(o) {
  fread(legacy_dir(sprintf("legacy_fisher_%s_%s.csv", s, o)))[top10 == TRUE]
}))))
known_go <- AnnotationDbi::keys(GO.db::GO.db)
pick <- function(m, u, suffix) {
  x <- results[mode == m & universe == u, .(set, ontology, go_id, genes_annotated, genes_with_dml, p_fisher, q_fisher, p_wallenius, q_wallenius)]
  setnames(x, 4:9, paste0(names(x)[4:9], "_", suffix)); x
}
top_check <- legacy_top[, .(set, ontology, legacy_rank = rank, go_id, term, legacy_transcripts_annotated = annotated,
                            legacy_transcripts_of_interest = significant, legacy_genes_behind_them = genes_significant,
                            legacy_p = p_fisher)]
# The paper's analysis: both universes, suffixes _tested and _eligible. Another
# analysis: the eligible universe, suffix _eligible_<id>_dml, so that the
# comparison table can set the analyses side by side in columns.
suffix <- if (an$paper) "eligible" else sprintf("eligible_%s_dml", an$id)
if (an$paper) top_check <- merge(top_check, pick(an$id, "tested", "tested"), by = c("set", "ontology", "go_id"), all.x = TRUE, sort = FALSE)
top_check <- merge(top_check, pick(an$id, "eligible", suffix), by = c("set", "ontology", "go_id"), all.x = TRUE, sort = FALSE)
# A set that was not tested in an ontology (see the Notes) has its own status, for every term.
not_tested <- if (length(not_tested)) rbindlist(not_tested) else data.table(mode = character(), set = character(), ontology = character(), universe = character())
set_untested <- paste(top_check$set, top_check$ontology) %in% not_tested[mode == an$id & universe == "eligible", paste(set, ontology)]
status_of <- function(annotated, q) fifelse(!top_check$go_id %in% known_go, "not in the current ontology",
                                    fifelse(set_untested, "set not tested: no DML gene carries a term of the ontology",
                                    fifelse(is.na(annotated), "no eligible gene carries it",
                                    fifelse(annotated < config$go_min_genes, "fewer annotated genes than go_min_genes",
                                    fifelse(q < config$go_fdr, "enriched", "not enriched")))))
top_check[, gene_level_status := status_of(get(paste0("genes_annotated_", suffix)), get(paste0("q_wallenius_", suffix)))]
setorder(top_check, set, ontology, legacy_rank)
write_table(top_check, tab("legacy_top10_at_gene_level.csv"))
print(top_check[, .N, by = .(set, gene_level_status)][order(set, gene_level_status)])

# =============================================================================
# Figures
# =============================================================================
bias <- rbindlist(bias_rows)
bias[, `:=`(set_label = factor(sets[set], levels = sets), universe = factor(universe, levels = names(universes)))]
# The three colours of run_config.R, and one for the tested interaction (option inter).
set_cols <- setNames(c(unname(config$colours$dml_set), "mediumpurple3"), set_labels)
# The eligible range of the analysis (dml_eligible(): half the threshold to 100 minus that).
eligible_label <- sprintf("Eligible CpGs in the gene (mean %g%% to %g%%)", an$dml_difference / 2, 100 - an$dml_difference / 2)
for (mode_i in an$id) {
  pb <- ggplot(bias[mode == mode_i], aes(cpg_bin, percent_with_dml, colour = set_label, group = set_label)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    facet_wrap(~universe, labeller = as_labeller(c(tested = "CpGs tested in the gene", eligible = eligible_label))) +
    scale_colour_manual(values = set_cols, name = NULL) +
    labs(x = "Number of CpGs in the gene body", y = "Genes that hold a DML (%)") +
    theme_wgbs() + theme(legend.position = "top")
  save_figure(pb, analysis_path(an, step, "figures", sprintf("figure_go_bias_%s", mode_i)), width = 9, height = 4.5)
}

# How many terms each analysis calls enriched.
legacy_counts <- rbindlist(lapply(paper_sets, function(s) rbindlist(lapply(onts, function(o) {
  x <- fread(legacy_dir(sprintf("legacy_fisher_%s_%s.csv", s, o)))
  data.table(set = s, ontology = o, analysis = "Draft: transcripts, p < 0.01, no correction", n_terms = sum(x$below_p_cutoff))
}))))
gl <- function(m, u, col, label) go_summary[mode == m & universe == u, .(set, ontology, analysis = label, n_terms = get(col))]
by_analysis <- if (an$paper) rbind(
  legacy_counts,
  gl("legacy", "tested", "fisher_q_below_fdr", "Genes, tested universe, Fisher, FDR"),
  gl("legacy", "eligible", "fisher_q_below_fdr", "Genes, eligible universe, Fisher, FDR"),
  gl("legacy", "eligible", "wallenius_q_below_fdr", "Genes, eligible universe, Wallenius, FDR")) else
  gl(an$id, "eligible", "wallenius_q_below_fdr",
     if (an$id == "corrected") "Same, DML called with the coverage filter" else sprintf("Same, DML of the optional analysis %s", an$id))
by_analysis[, set_label := factor(sets[set], levels = sets)]
write_table(by_analysis, tab("go_terms_by_analysis.csv"))

# The six tables of the analyses side by side (04_output/not_adopted/comparisons/),
# rebuilt from every analysis folder whose stamp lists step 10, the paper's first.
# The stamp of this analysis therefore lists step 10 from here on.
write_stamp(an, step)
ids <- analysis_order(c("legacy", analysis_folders()))
rebuild_comparison(step, "go_gene_level_summary.csv")
rebuild_comparison(step, "dml_probability_by_cpg_count.csv")
rebuild_comparison(step, "go_enriched_terms.csv", function(b) order(b$mode, b$set, b$ontology, as.numeric(b$q_wallenius), method = "radix"))
rebuild_comparison(step, "go_enriched_term_genes.csv", function(b) order(b$mode, b$set, b$ontology, b$go_id, b$gene_id, method = "radix"))
rebuild_comparison(step, "go_terms_by_analysis.csv")
top_files <- vapply(ids, function(id) analysis_input(id, step, "tables", "legacy_top10_at_gene_level.csv"), "")
top_files <- top_files[file.exists(top_files) & vapply(ids, stamp_lists_step, logical(1), step = step)]
if (length(top_files) > 0) {
  tops <- lapply(top_files, fread)
  side <- tops[[1]]                                 # the paper's, when step 10 has run for it
  status <- setNames(side$gene_level_status, paste(side$set, side$ontology, side$go_id)); side[, gene_level_status := NULL]
  for (k in seq_along(tops)[-1]) {
    cols <- grep("_eligible_.*_dml$", names(tops[[k]]), value = TRUE)
    side <- merge(side, tops[[k]][, c("set", "ontology", "go_id", cols), with = FALSE], by = c("set", "ontology", "go_id"), all.x = TRUE, sort = FALSE)
  }
  side[, gene_level_status := unname(status[paste(set, ontology, go_id)])]
  setorder(side, set, ontology, legacy_rank)
  write_table(side, not_adopted_path("comparisons", step, "tables", "legacy_top10_at_gene_level.csv"))
}
by_analysis <- fread(not_adopted_path("comparisons", step, "tables", "go_terms_by_analysis.csv"))
by_analysis[, `:=`(analysis = factor(analysis, levels = rev(unique(analysis))), set_label = factor(set_labels[set], levels = set_labels))]
pa <- ggplot(by_analysis, aes(n_terms, analysis, fill = set_label)) +
  geom_col(position = position_dodge2(width = 0.85, reverse = TRUE), width = 0.8) +
  geom_text(aes(label = n_terms), position = position_dodge2(width = 0.85, reverse = TRUE), hjust = -0.2, size = 2.8) +
  facet_wrap(~ontology, scales = "free_x") +
  scale_fill_manual(values = set_cols, name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "GO terms called enriched", y = NULL) + theme_wgbs() +
  theme(legend.position = "top", panel.spacing.x = unit(1.4, "lines"))
dir.create(not_adopted_path("comparisons", step, "figures"), recursive = TRUE, showWarnings = FALSE)
save_figure(pa, not_adopted_path("comparisons", step, "figures", "figure_go_terms_by_analysis"), width = 11, height = 4.8)

# The enriched terms themselves (eligible universe, Wallenius), strongest ten per set and ontology.
for (mode_i in an$id) {
  e <- enriched[mode == mode_i]
  if (nrow(e) == 0) { log_msg("[%s] no term is enriched at q < %g in the eligible universe", mode_i, config$go_fdr); next }
  e <- e[, head(.SD[order(q_wallenius)], 10), by = .(set, ontology)]
  e[, `:=`(label = paste0(ifelse(nchar(term) > 55, paste0(substr(term, 1, 52), "..."), term), "  (", go_id, ")"),
           set_label = factor(sets[set], levels = sets))]
  e[, label := factor(label, levels = unique(label[order(-q_wallenius)]))]
  pe <- ggplot(e, aes(enrichment, label, size = genes_with_dml, colour = -log10(q_wallenius))) + geom_point() +
    facet_grid(ontology ~ set_label, scales = "free_y", space = "free_y") +
    scale_colour_viridis_c(name = expression(-log[10] ~ q)) + scale_size_area(max_size = 5, name = "Genes with a DML") +
    labs(x = "Observed / expected genes with a DML", y = NULL) + theme_wgbs() +
    theme(axis.text.y = element_text(size = 7), strip.text.y = element_text(angle = 0))
  save_figure(pe, analysis_path(an, step, "figures", sprintf("figure_go_enriched_%s", mode_i)), width = 14,
              height = max(4, 0.19 * nrow(e) + 2))
}

print(dcast(go_summary, mode + set + ontology ~ universe, value.var = c("terms_tested", "wallenius_q_below_fdr")))
write_analyses_included()
finish_step(step, t0)

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
#                            (mean methylation 10% to 90%, a screen and not a
#                            bound; see dml_eligible()). Gene-body methylation in this
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
# Inputs    04_output/08_gene_annotation/tables/gene_table_<mode>.csv.gz, gene_go_map.csv.gz
#           04_output/09_go_enrichment_legacy/tables/legacy_fisher_<set>_<ontology>.csv
# Outputs   04_output/10_go_enrichment_gene_level/tables/
#             go_gene_level_<set>_<mode>.csv        every term, both universes, both tests
#             go_gene_level_summary.csv             terms tested and enriched, per analysis
#             go_enriched_terms.csv                 terms with q < go_fdr (eligible universe, Wallenius)
#             go_enriched_term_genes.csv            the DML genes behind those terms
#             legacy_top10_at_gene_level.csv        each "enriched" term of the draft, retested
#             dml_probability_by_cpg_count.csv      the bias the Wallenius test adjusts for
#           04_output/10_go_enrichment_gene_level/figures/
#             figure_go_bias_<mode>, figure_go_terms_by_analysis, figure_go_enriched_<mode> (when any term is enriched)
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
#           - The heat contrast is confounded with cohort and sequencing batch
#             (issue I-01). Enrichment of "heat" DML genes describes genes that
#             differ between the two cohorts, whatever the cause.
#           - A DML is assigned to a gene when it lies in the gene body.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "topGO", "goseq", "ggplot2", "patchwork"))
if (as.character(packageVersion("GO.db")) == config$legacy_go_db$version) {
  stop("GO.db ", config$legacy_go_db$version, " (the legacy release) is loaded. Run this step in a fresh R process.")
}

step <- "10_go_enrichment_gene_level"
t0 <- start_step(step, "gene-level GO enrichment against tested and eligible genes")
tab <- function(name) output_path(step, "tables", name)
ann_in <- function(name) output_path("08_gene_annotation", "tables", name)
sets <- c(ploidy = "ploidy DML", heat = "heat DML", shared = "in both lists")
onts <- c("BP", "CC", "MF")
universes <- c(tested = "n_cpg_tested", eligible = "n_cpg_eligible")
log_msg("GO.db %s (ontology of %s), topGO %s, goseq %s", as.character(packageVersion("GO.db")), go_source_date(),
        as.character(packageVersion("topGO")), as.character(packageVersion("goseq")))

gene_go <- fread(ann_in("gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
gene2go_all <- split(gene_go$go_id, gene_go$gene_id)

results <- list(); summary_rows <- list(); bias_rows <- list(); term_genes <- list()
for (mode in config$modes) {
  gene_tab <- fread(ann_in(sprintf("gene_table_%s.csv.gz", mode)), colClasses = c(gene_id = "character"))
  gene_tab <- gene_tab[biotype == "protein_coding"]
  any_dml <- gene_tab$n_dml_ploidy + gene_tab$n_dml_heat + gene_tab$n_dml_shared > 0

  for (u in names(universes)) {
    # Genes of the universe. A gene that holds a DML belongs to the eligible
    # universe by definition; one or two DML have a mean just outside the 10% to
    # 90% range, and their genes are kept with a count of one.
    n_cpg <- gene_tab[[universes[[u]]]]
    if (u == "eligible") n_cpg <- pmax(n_cpg, as.integer(any_dml))
    uni <- cbind(gene_tab[n_cpg > 0], bias = n_cpg[n_cpg > 0])
    annotated_ids <- uni[n_go_terms > 0, gene_id]

    # Probability that a gene holds a DML, by number of CpGs (all protein-coding genes of the universe).
    uni[, cpg_bin := cut(bias, c(0, 1, 2, 5, 10, 20, 50, 100, Inf), labels = c("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"))]
    for (s in names(sets)) {
      bias_rows[[length(bias_rows) + 1]] <- uni[, .(mode = mode, universe = u, set = s, genes = .N,
        genes_with_dml = sum(get(paste0("n_dml_", s)) > 0), percent_with_dml = 100 * mean(get(paste0("n_dml_", s)) > 0)),
        by = cpg_bin][order(cpg_bin)]
    }

    for (ont in onts) {
      td <- go_gene_level_object(ont, annotated_ids, gene2go_all)
      for (s in names(sets)) {
        dml_genes <- uni[get(paste0("n_dml_", s)) > 0, gene_id]
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
              paste(vapply(tail(summary_rows, 3), function(r) sprintf("%s %d", r$set, r$wallenius_q_below_fdr), ""), collapse = ", "))
    }
  }
}
results <- rbindlist(results); go_summary <- rbindlist(summary_rows)
go_summary[, `:=`(go_db = as.character(packageVersion("GO.db")), go_source_date = go_source_date())]
for (mode_i in config$modes) for (s in names(sets)) {
  write_table(results[mode == mode_i & set == s], tab(sprintf("go_gene_level_%s_%s.csv", s, mode_i)))
}
write_table(go_summary, tab("go_gene_level_summary.csv"))
write_table(rbindlist(bias_rows), tab("dml_probability_by_cpg_count.csv"))

enriched <- results[universe == "eligible" & tested == TRUE & q_wallenius < config$go_fdr]
setorder(enriched, mode, set, ontology, q_wallenius)
write_table(enriched, tab("go_enriched_terms.csv"))
term_genes <- if (length(term_genes) > 0) rbindlist(term_genes) else
  data.table(mode = character(0), set = character(0), ontology = character(0), go_id = character(0), gene_id = character(0))
gene_names <- unique(rbindlist(lapply(config$modes, function(m) fread(ann_in(sprintf("gene_table_%s.csv.gz", m)),
                                                                       colClasses = c(gene_id = "character"))[, .(gene_id, gene, product)])))
term_genes <- merge(term_genes, gene_names, by = "gene_id", all.x = TRUE, sort = FALSE)
term_genes <- merge(term_genes, unique(results[, .(go_id, term)]), by = "go_id", all.x = TRUE, sort = FALSE)
setcolorder(term_genes, c("mode", "set", "ontology", "go_id", "term", "gene_id", "gene", "product"))
setorder(term_genes, mode, set, ontology, go_id, gene_id)
write_table(term_genes, tab("go_enriched_term_genes.csv"))

# =============================================================================
# The "enriched" terms of the draft, retested at gene level
# =============================================================================
legacy_dir <- function(name) output_path("09_go_enrichment_legacy", "tables", name)
legacy_top <- rbindlist(lapply(names(sets), function(s) rbindlist(lapply(onts, function(o) {
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
top_check <- merge(top_check, pick("legacy", "tested", "tested"), by = c("set", "ontology", "go_id"), all.x = TRUE, sort = FALSE)
top_check <- merge(top_check, pick("legacy", "eligible", "eligible"), by = c("set", "ontology", "go_id"), all.x = TRUE, sort = FALSE)
if ("corrected" %in% config$modes) {
  top_check <- merge(top_check, pick("corrected", "eligible", "eligible_corrected_dml"), by = c("set", "ontology", "go_id"),
                     all.x = TRUE, sort = FALSE)
}
top_check[, gene_level_status := fifelse(!go_id %in% known_go, "not in the current ontology",
                                 fifelse(is.na(genes_annotated_eligible), "no eligible gene carries it",
                                 fifelse(genes_annotated_eligible < config$go_min_genes, "fewer annotated genes than go_min_genes",
                                 fifelse(q_wallenius_eligible < config$go_fdr, "enriched", "not enriched"))))]
setorder(top_check, set, ontology, legacy_rank)
write_table(top_check, tab("legacy_top10_at_gene_level.csv"))
print(top_check[, .N, by = .(set, gene_level_status)][order(set, gene_level_status)])

# =============================================================================
# Figures
# =============================================================================
bias <- rbindlist(bias_rows)
bias[, `:=`(set_label = factor(sets[set], levels = sets), universe = factor(universe, levels = names(universes)))]
set_cols <- setNames(unname(config$colours$dml_set), sets)
for (mode_i in config$modes) {
  pb <- ggplot(bias[mode == mode_i], aes(cpg_bin, percent_with_dml, colour = set_label, group = set_label)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    facet_wrap(~universe, labeller = as_labeller(c(tested = "CpGs tested in the gene", eligible = "Eligible CpGs in the gene (mean 10% to 90%)"))) +
    scale_colour_manual(values = set_cols, name = NULL) +
    labs(x = "Number of CpGs in the gene body", y = "Genes that hold a DML (%)") +
    theme_wgbs() + theme(legend.position = "top")
  save_figure(pb, output_path(step, "figures", sprintf("figure_go_bias_%s", mode_i)), width = 9, height = 4.5)
}

# How many terms each analysis calls enriched.
legacy_counts <- rbindlist(lapply(names(sets), function(s) rbindlist(lapply(onts, function(o) {
  x <- fread(legacy_dir(sprintf("legacy_fisher_%s_%s.csv", s, o)))
  data.table(set = s, ontology = o, analysis = "Draft: transcripts, p < 0.01, no correction", n_terms = sum(x$below_p_cutoff))
}))))
gl <- function(m, u, col, label) go_summary[mode == m & universe == u, .(set, ontology, analysis = label, n_terms = get(col))]
by_analysis <- rbind(
  legacy_counts,
  gl("legacy", "tested", "fisher_q_below_fdr", "Genes, tested universe, Fisher, FDR"),
  gl("legacy", "eligible", "fisher_q_below_fdr", "Genes, eligible universe, Fisher, FDR"),
  gl("legacy", "eligible", "wallenius_q_below_fdr", "Genes, eligible universe, Wallenius, FDR"),
  if ("corrected" %in% config$modes) gl("corrected", "eligible", "wallenius_q_below_fdr",
                                        "Same, DML called with the coverage filter"))
by_analysis[, `:=`(analysis = factor(analysis, levels = rev(unique(analysis))), set_label = factor(sets[set], levels = sets))]
write_table(by_analysis, tab("go_terms_by_analysis.csv"))
pa <- ggplot(by_analysis, aes(n_terms, analysis, fill = set_label)) +
  geom_col(position = position_dodge2(width = 0.85, reverse = TRUE), width = 0.8) +
  geom_text(aes(label = n_terms), position = position_dodge2(width = 0.85, reverse = TRUE), hjust = -0.2, size = 2.8) +
  facet_wrap(~ontology, scales = "free_x") +
  scale_fill_manual(values = set_cols, name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "GO terms called enriched", y = NULL) + theme_wgbs() +
  theme(legend.position = "top", panel.spacing.x = unit(1.4, "lines"))
save_figure(pa, output_path(step, "figures", "figure_go_terms_by_analysis"), width = 11, height = 4.8)

# The enriched terms themselves (eligible universe, Wallenius), strongest ten per set and ontology.
for (mode_i in config$modes) {
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
  save_figure(pe, output_path(step, "figures", sprintf("figure_go_enriched_%s", mode_i)), width = 14,
              height = max(4, 0.19 * nrow(e) + 2))
}

print(dcast(go_summary, mode + set + ontology ~ universe, value.var = c("terms_tested", "wallenius_q_below_fdr")))
finish_step(step, t0)

# 09_go_enrichment_legacy.R
#
# Purpose   Repeat the GO enrichment exactly as it was run for the manuscript
#           draft, so that every GO number in the draft can be traced, and record
#           what that procedure did. It is a reproduction, not a recommendation:
#           step 10 holds the gene-level analysis.
#
#           The legacy procedure (5_WGBS_Functional-enrichement-analysis-BOTH.R):
#             1. "Genes" are transcripts. The universe is the transcript-to-GO map
#                geneid2go-union1x-BOTH.tab, from which 15,666 annotated mRNAs
#                that are adjacent when sorted by accession are missing
#                (XM_011455141.3 to XM_034461018.1; issue I-39).
#             2. The transcripts of interest are those named in the DML
#                annotation table, read with read.table() defaults. Apostrophes in
#                GO term names made that reader skip most of the ploidy table and
#                a quarter of the heat table (issue I-40).
#             3. topGO classic Fisher test per ontology, no correction for the
#                number of terms. The ten top-ranked terms per ontology were
#                treated as "enriched" for Figure 5, Figure 9 and Table 1. The
#                legacy script filtered p < 0.01 only for its simplifyEnrichment
#                PDFs, which are not in the draft (issue I-21); Figures 6 to 8
#                were drawn from unfiltered annotation (issue I-48).
#             4. The ontology was GO.db 3.13.0. Under any other release the
#                tables change; this step loads that release from a side library.
#             5. topGO received one name per line of the map, and names repeat:
#                48,263 lines for 27,469 transcripts. It counts every line in the
#                totals of each Fisher test and each name once per term, so the
#                Expected column and the p-values mix two units (issue I-43).
#                This step keeps the repeats in order to match the archive.
#
# Inputs    04_output/08_gene_annotation/tables/
#             dml_go_annotation_{ploidy,heat,shared}_legacy.tsv.gz, transcript_table.csv.gz,
#             gene_table_legacy.csv.gz, annotation_summary.csv
#           03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz
#           GO.db 3.13.0 (Bioconductor archive; installed into the cache on first use)
# Outputs   04_output/09_go_enrichment_legacy/tables/
#             legacy_fisher_<set>_<ontology>.csv          every term: counts, p, rank, genes behind it
#             legacy_top10_term_transcripts_<set>_<ontology>.csv
#             legacy_dml_in_top10_genes_<set>_<ontology>.csv
#             legacy_reader_loss.csv, legacy_universe_by_chromosome.csv, legacy_funnel.csv
#             figure_5_values_legacy.csv, figure_9_values_legacy.csv, table_1_legacy.csv
#             legacy_run_info.csv
#           04_output/09_go_enrichment_legacy/figures/figure_5_legacy.(png|pdf), figure_9_legacy.(png|pdf)
#
# Manuscript items: Results 3.4 and 3.5; Figure 5; Figure 9; Table 1.
# Notes     - Runs on the legacy-mode DML only.
#           - "shared" is the set the draft calls interaction-DML: loci in both
#             the ploidy and the heat list (issue I-02). Its methylation
#             difference is the heat one, as in the legacy table.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
go_db_version <- use_legacy_go_db(config)            # must come before topGO is loaded
load_packages(c("data.table", "topGO", "ggplot2", "patchwork"))

step <- "09_go_enrichment_legacy"
t0 <- start_step(step, "GO enrichment as run for the manuscript draft")
tab <- function(name) output_path(step, "tables", name)
ann_in <- function(name) output_path("08_gene_annotation", "tables", name)
sets <- c(ploidy = "ploidy-DML", heat = "heat-DML", shared = "interaction-DML")
onts <- c(BP = "P", CC = "C", MF = "F")
log_msg("GO.db %s (ontology of %s), topGO %s", go_db_version, go_source_date(), as.character(packageVersion("topGO")))

tx_table <- fread(ann_in("transcript_table.csv.gz"), colClasses = c(gene_id = "character"))
gene_table <- fread(ann_in("gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))
tx_gene <- setNames(tx_table$gene_id, tx_table$transcript_id)

# =============================================================================
# 1. The DML annotation tables, as written and as the legacy script read them
# =============================================================================
ann_file <- function(s) ann_in(sprintf("dml_go_annotation_%s_legacy.tsv.gz", s))
ann_exact  <- lapply(setNames(names(sets), names(sets)), function(s) read_dml_annotation(ann_file(s), "exact"))
ann_legacy <- lapply(setNames(names(sets), names(sets)), function(s) suppressWarnings(read_dml_annotation(ann_file(s), "legacy")))
count_side <- function(d) {
  data.table(rows = nrow(d), transcripts = uniqueN(d$transcript), genes = uniqueN(d$geneID),
             dml = uniqueN(paste(d$chr, d$start)),
             rows_BP = sum(d$GOcat == "P"), rows_CC = sum(d$GOcat == "C"), rows_MF = sum(d$GOcat == "F"),
             go_ids_BP = uniqueN(d$GOID[d$GOcat == "P"]), go_ids_CC = uniqueN(d$GOID[d$GOcat == "C"]),
             go_ids_MF = uniqueN(d$GOID[d$GOcat == "F"]))
}
reader_loss <- rbindlist(lapply(names(sets), function(s) rbind(
  cbind(set = s, reader = "every line (correct)", count_side(ann_exact[[s]])),
  cbind(set = s, reader = "read.table defaults (legacy)", count_side(ann_legacy[[s]])))))
write_table(reader_loss, tab("legacy_reader_loss.csv"))
for (s in names(sets)) {
  log_msg("%-6s annotation table: %6d lines, %4d transcripts; the legacy reader kept %6d rows, %4d transcripts",
          s, nrow(ann_exact[[s]]), uniqueN(ann_exact[[s]]$transcript), nrow(ann_legacy[[s]]), uniqueN(ann_legacy[[s]]$transcript))
}

# =============================================================================
# 2. The legacy universe and what is missing from it
# =============================================================================
geneID2GO <- topGO::readMappings(input_path("annotation", "geneid2go-union1x-BOTH.tab.gz"))
geneNames <- names(geneID2GO)                  # one name per line of the map; names repeat, as in the legacy run
universe <- unique(geneNames)

annotated <- tx_table[in_legacy_go_map == TRUE]          # mRNAs the legacy route could annotate
annotated[, sequence := fifelse(chr %in% config$linkage_groups, chr, "unplaced scaffolds")]
by_chr <- annotated[, .(mRNAs_with_GO_terms = .N, in_legacy_universe = sum(in_legacy_universe)), by = sequence][order(sequence)]
by_chr <- rbind(by_chr, data.table(sequence = "all", mRNAs_with_GO_terms = nrow(annotated),
                                   in_legacy_universe = sum(annotated$in_legacy_universe)))
by_chr[, percent_in_universe := round(100 * in_legacy_universe / mRNAs_with_GO_terms, 1)]
# Nearly all missing transcripts are adjacent when the annotated mRNAs are sorted by accession.
setorder(annotated, transcript_id)
runs <- rle(annotated$in_legacy_universe)
longest <- which.max(ifelse(runs$values, 0L, runs$lengths))
block_end <- sum(runs$lengths[seq_len(longest)]); block_start <- block_end - runs$lengths[longest] + 1L
by_chr[, `:=`(largest_missing_block_n = runs$lengths[longest],
              largest_missing_block_first = annotated$transcript_id[block_start],
              largest_missing_block_last = annotated$transcript_id[block_end])]
write_table(by_chr, tab("legacy_universe_by_chromosome.csv"))
log_msg("legacy universe: %d transcripts; %d of the %d mRNAs with GO terms (%.0f%%); largest missing block %d transcripts, %s to %s",
        length(universe), sum(annotated$in_legacy_universe), nrow(annotated),
        100 * mean(annotated$in_legacy_universe), runs$lengths[longest],
        annotated$transcript_id[block_start], annotated$transcript_id[block_end])

# From genes with a DML to transcripts that entered the test.
ann_summary <- fread(ann_in("annotation_summary.csv"))[mode == "legacy"]
funnel <- rbindlist(lapply(names(sets), function(s) {
  e <- unique(ann_exact[[s]][, .(transcript, geneID)]); l <- unique(ann_legacy[[s]][, .(transcript, geneID)])
  tested <- l[transcript %chin% universe]
  data.table(set = s,
             stage = c("1 genes holding a DML", "2 with GO terms by the legacy route", "3 kept by the legacy reader",
                       "4 in the legacy universe: tested as genes of interest"),
             genes = c(ann_summary[set == s, n_genes_with_dml], uniqueN(e$geneID), uniqueN(l$geneID), uniqueN(tested$geneID)),
             transcripts = c(NA, uniqueN(e$transcript), uniqueN(l$transcript), uniqueN(tested$transcript)))
}))
funnel[, percent_of_genes_with_GO_terms := round(100 * genes / genes[2], 1), by = set]
funnel[stage == "1 genes holding a DML", percent_of_genes_with_GO_terms := NA]
write_table(funnel, tab("legacy_funnel.csv"))

# =============================================================================
# 3. topGO, classic Fisher test, as in the legacy script
# =============================================================================
fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
fig5 <- list(); fig9 <- list(); table1 <- list()
for (ont in names(onts)) {
  td <- NULL
  for (s in names(sets)) {
    gene_list <- factor(as.integer(geneNames %in% ann_legacy[[s]]$transcript)); names(gene_list) <- geneNames
    td <- if (is.null(td)) {
      suppressMessages(new("topGOdata", ontology = ont, allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = geneID2GO))
    } else updateGenes(td, gene_list)
    rf <- suppressMessages(getSigGroups(td, fisher)); p <- score(rf)
    gt <- GenTable(td, classic = rf, ranksOf = "classic", orderBy = "classic", topNodes = length(p), numChar = 1000)

    # Distinct genes behind the transcript counts of each term.
    in_term <- genesInTerm(td)[gt$GO.ID]; sig_tx <- sigGenes(td)
    res <- data.table(
      set = s, ontology = ont, rank = seq_len(nrow(gt)), go_id = gt$GO.ID, term = gt$Term,
      annotated = gt$Annotated, significant = gt$Significant, expected = gt$Expected, p_fisher = unname(p[gt$GO.ID]),
      genes_annotated = vapply(in_term, function(v) uniqueN(tx_gene[v]), integer(1)),
      genes_significant = vapply(in_term, function(v) uniqueN(tx_gene[v[v %chin% sig_tx]]), integer(1)))
    res[, `:=`(top10 = rank <= config$go_legacy_topn, below_p_cutoff = p_fisher < config$go_pvalue)]
    write_table(res, tab(sprintf("legacy_fisher_%s_%s.csv", s, ont)))

    # ---- the ten top-ranked terms, and what the legacy script derived from them ----
    top <- res[top10 == TRUE]
    members <- genesInTerm(td, whichGO = top$go_id)                 # universe transcripts carrying each term
    term_tx <- data.table(GOID = rep(names(members), lengths(members)), transcript = unlist(members, use.names = FALSE))
    term_tx <- merge(top[, .(GOID = go_id, topGOterm = term, p.value = p_fisher)], term_tx, by = "GOID", sort = FALSE)
    ann_ont <- ann_legacy[[s]][GOcat == onts[[ont]]]
    annot <- merge(term_tx, ann_ont, by = c("GOID", "transcript"), all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
    annot <- unique(annot[, !"GOcat"])
    write_table(annot, tab(sprintf("legacy_top10_term_transcripts_%s_%s.csv", s, ont)))

    # DML "within genes with enriched GO terms": every DML of every transcript on
    # the list above, provided the transcript has a row of this ontology.
    in_genes <- unique(merge(unique(annot[, .(transcript)]), ann_legacy[[s]], by = "transcript")[
      GOcat == onts[[ont]], .(geneID, chr, start, end, meth.diff)])
    in_genes <- merge(in_genes[, geneID := as.character(geneID)], gene_table[, .(geneID = gene_id, gene, product)],
                      by = "geneID", all.x = TRUE, sort = FALSE)
    setorder(in_genes, geneID, chr, start)
    write_table(in_genes, tab(sprintf("legacy_dml_in_top10_genes_%s_%s.csv", s, ont)))

    fig5[[length(fig5) + 1]] <- data.table(
      set = s, ontology = ont,
      published_label = c("Unique transcripts", "GO terms", "Unique genes", "Gene products"),
      value = c(uniqueN(annot$transcript), nrow(ann_ont), uniqueN(na.omit(annot$geneID)), nrow(annot)),
      what_it_counts = c("transcripts of the universe that carry one of the ten top-ranked terms, with or without a DML",
                         "rows of the DML x transcript x GO term table for this ontology, as the legacy reader returned it",
                         "genes that hold a DML and are annotated directly with one of the ten top-ranked terms",
                         paste("rows of the table that pairs each of the ten top-ranked terms with the universe transcripts that",
                               "carry it, a transcript being repeated for each of its DML and each GOslim bin of the term")))
    fig9[[length(fig9) + 1]] <- data.table(set = s, ontology = ont, total = nrow(in_genes),
                                           hyper = sum(in_genes$meth.diff >= 0), hypo = sum(in_genes$meth.diff <= 0),
                                           genes = uniqueN(in_genes$geneID))
    if (ont == "BP") table1[[s]] <- cbind(dml_set = sets[[s]], in_genes)
    log_msg("%-6s %s: %5d terms; %4d with p < %g; transcripts of interest %d; genes behind the ten top-ranked terms: %s",
            s, ont, nrow(res), sum(res$below_p_cutoff), config$go_pvalue, uniqueN(sig_tx),
            paste(top$genes_significant, collapse = " "))
  }
}
fig5 <- rbindlist(fig5); fig9 <- rbindlist(fig9); table1 <- rbindlist(table1)
write_table(fig5, tab("figure_5_values_legacy.csv"))
write_table(fig9, tab("figure_9_values_legacy.csv"))
table1[, meth.diff := round(meth.diff, 2)]
setnames(table1, c("geneID", "product", "chr", "meth.diff"), c("gene_id", "gene_name", "chromosome", "methylation_difference"))
write_table(table1[, .(dml_set, gene_id, gene_name, chromosome, start, end, methylation_difference)], tab("table_1_legacy.csv"))
write_table(data.table(item = c("GO.db", "GO source date", "topGO", "universe transcripts", "universe map lines"),
                       value = c(go_db_version, go_source_date(), as.character(packageVersion("topGO")),
                                 length(universe), length(geneNames))), tab("legacy_run_info.csv"))

# =============================================================================
# 4. Figures 5 and 9 with the values of the draft
# =============================================================================
set_cols <- setNames(unname(config$colours$dml_set), sets)
fig5[, `:=`(set_label = factor(sets[set], levels = sets),
            published_label = factor(published_label, levels = c("Unique transcripts", "GO terms", "Unique genes", "Gene products")))]
p5 <- ggplot(fig5, aes(ontology, value, fill = set_label)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = format(value, big.mark = ",", trim = TRUE)), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.8) +
  facet_wrap(~published_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = set_cols, name = NULL) +
  scale_y_continuous(labels = function(v) format(v, big.mark = ","), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Ontology", y = "Number", caption = paste(
    "Values as in the draft. The panel titles are the draft's; figure_5_values_legacy.csv states what each number counts.")) +
  theme_wgbs() + theme(legend.position = "top", plot.caption = element_text(size = 7, hjust = 0))
save_figure(p5, output_path(step, "figures", "figure_5_legacy"), width = 9, height = 7.5)

long9 <- melt(fig9, id.vars = c("set", "ontology"), measure.vars = c("hyper", "hypo"), variable.name = "direction", value.name = "n")
long9[, set_label := factor(sets[set], levels = sets)]
p9 <- ggplot(long9, aes(direction, n, fill = set_label)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  facet_wrap(~ontology, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = set_cols, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "DML in genes that carry one of\nthe ten top-ranked GO terms") +
  theme_wgbs() + theme(legend.position = "top")
save_figure(p9, output_path(step, "figures", "figure_9_legacy"), width = 9, height = 4)

print(dcast(fig9, set ~ ontology, value.var = "total"))
finish_step(step, t0)

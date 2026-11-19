# 08_gene_annotation.R
#
# Purpose   Link the DML to genes, transcripts, and GO terms. Replaces the
#           bedtools, cut, tr, and join commands of the legacy notebook
#           (4_WGBS_GOterm-Annotation-BOTH.ipynb) with R.
#             1. Annotation tables that do not depend on the DML: transcripts and
#                their genes, the transcript-to-GO map with repaired transcript
#                names, and a gene-to-GO map.
#             2. For each mode and DML set: every (DML, gene) overlap, a per-gene
#                summary, and a gene table with the number of tested CpGs and of
#                DML per gene (the input of the gene-level enrichment in step 09).
#             3. For the legacy mode: the ten-column annotation tables of the
#                legacy notebook (one line per DML, transcript, and GO term), in
#                the legacy line order (the input of the legacy enrichment in
#                step 09).
#
# Inputs    04_output/05_differential_methylation/tables/
#             united_cpg_<mode>.csv.gz, dml_{ploidy,heat,shared}_<mode>.csv
#           03_input_files/genome_features/ gene and mRNA tracks
#           03_input_files/annotation/ geneid2go.tab.gz, geneid2go-union1x-BOTH.tab.gz,
#             20210605-cgigas-roslin-mito-blastx.outfmt6.gz, legacy_go_goslim_terms.tsv.gz
# Outputs   04_output/08_gene_annotation/tables/
#             transcript_table.csv.gz           transcript, gene, product, GO flags
#             gene_go_map.csv.gz                gene_id, go_id (union over the gene's mRNAs)
#             go_map_id_repair.csv              what the legacy transcript names stand for
#             gene_table_<mode>.csv.gz          one row per gene: tested CpGs, DML per set
#             dml_gene_overlap_<set>_<mode>.csv one row per (DML, gene)
#             genes_with_dml_<set>_<mode>.csv   one row per gene that holds a DML
#             dml_go_annotation_<set>_legacy.tsv.gz   legacy ten-column table, no header
#             annotation_summary.csv            counts per mode and set
#
# Manuscript items: Results 3.4 and 3.5 (inputs); supplementary lists of genes with DML.
# Notes     - A DML is assigned to every gene whose body it lies in, as in the
#             legacy analysis. Genes overlap, so a DML can be assigned twice.
#           - Legacy transcript route: gene -> all of its mRNAs -> GO terms of
#             each mRNA. Three things were lost on the way and are reproduced in
#             the legacy tables only: mRNAs of the 9 genes whose name holds a
#             hyphen, the 123 curated (NM_) mRNAs with GO terms, and every GO ID
#             missing from the 2017 GOslim table. Issue I-03.
#           - Gene route (both modes): gene -> union of the GO terms of its mRNAs,
#             with the NM_ names repaired.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges"))

step <- "08_gene_annotation"
t0 <- start_step(step, "DML to genes, transcripts, and GO terms")
tab <- function(name) output_path(step, "tables", name)
dml_file <- function(f) output_path("05_differential_methylation", "tables", f)
dml_sets <- c("ploidy", "heat", "shared")

# =============================================================================
# 1. Annotation tables
# =============================================================================
genes <- read_genes(); tx <- read_transcripts()
stopifnot(!anyDuplicated(genes$gene_id), !anyDuplicated(tx$transcript_id), all(tx$gene_id %in% genes$gene_id))
gene_gr <- GRanges(genes$chr, IRanges(genes$start, genes$end))

go_full <- read_go_map("full")                       # map_id, go_id
go_universe_ids <- unique(read_go_map("legacy_universe")$map_id)
repair <- merge(unique(go_full[, .(map_id)]), true_transcript_ids(), by = "map_id", all.x = TRUE)
repair[, on_mrna_track := transcript_id %in% tx$transcript_id]
repair[, matched_by_legacy := map_id %in% tx$transcript_id]   # the legacy join used the map name as is
write_table(repair, tab("go_map_id_repair.csv"))
log_msg("GO map: %d transcript names; %d are mRNAs on the track; %d of those were renamed and lost in the legacy join",
        nrow(repair), sum(repair$on_mrna_track), sum(repair$on_mrna_track & !repair$matched_by_legacy))

# Transcript-level GO, true names, mRNAs only.
tx_go <- merge(go_full, repair[on_mrna_track == TRUE, .(map_id, transcript_id)], by = "map_id")[, .(transcript_id, go_id)]
gene_go <- unique(merge(tx_go, tx[, .(transcript_id, gene_id)], by = "transcript_id")[, .(gene_id, go_id)])
setorder(gene_go, gene_id, go_id)
write_table(gene_go, tab("gene_go_map.csv.gz"))

tx[, `:=`(n_go_terms = tx_go[, .N, by = transcript_id][.SD, on = "transcript_id", N],
          in_legacy_go_map = transcript_id %in% go_full$map_id,
          in_legacy_universe = transcript_id %in% go_universe_ids)]
tx[is.na(n_go_terms), n_go_terms := 0L]
write_table(tx, tab("transcript_table.csv.gz"))

gene_base <- merge(genes, gene_products(tx), by = "gene_id", all.x = TRUE, sort = FALSE)
gene_base[is.na(n_transcripts), n_transcripts := 0L]
gene_base <- merge(gene_base, gene_go[, .(n_go_terms = .N), by = gene_id], by = "gene_id", all.x = TRUE, sort = FALSE)
gene_base[is.na(n_go_terms), n_go_terms := 0L]
gene_base <- merge(gene_base, tx[, .(in_legacy_universe = any(in_legacy_universe)), by = gene_id], by = "gene_id", all.x = TRUE, sort = FALSE)
gene_base[is.na(in_legacy_universe), in_legacy_universe := FALSE]
gene_base <- gene_base[match(genes$gene_id, gene_id)]          # keep the order of gene_gr
log_msg("%d genes (%d protein coding); %d with GO terms; %d of those have a transcript in the legacy universe",
        nrow(gene_base), sum(gene_base$biotype == "protein_coding"), sum(gene_base$n_go_terms > 0),
        sum(gene_base$n_go_terms > 0 & gene_base$in_legacy_universe))

# =============================================================================
# 2. DML to genes, per mode and set
# =============================================================================
legacy_terms <- read_legacy_go_terms()
summary_rows <- list()
for (mode in config$modes) {
  united <- fread(dml_file(sprintf("united_cpg_%s.csv.gz", mode)), select = c("chr", "start"))
  gene_tab <- copy(gene_base)
  gene_tab[, n_cpg_tested := suppressWarnings(countOverlaps(gene_gr, cpg_ranges(united$chr, united$start), ignore.strand = TRUE))]

  for (set in dml_sets) {
    dml <- fread(dml_file(sprintf("dml_%s_%s.csv", set, mode)))
    # The loci in both lists carry two differences. The legacy "interaction" table
    # kept the heat one, so meth.diff is heat minus control for this set; the
    # ploidy difference is carried along in its own column.
    if (set == "shared") dml[, `:=`(meth.diff = meth.diff_heat, direction = direction_heat)]
    pr <- overlap_pairs(cpg_ranges(dml$chr, dml$start), gene_gr)
    pairs <- cbind(dml[pr$cpg, .(chr, start, end, meth.diff, direction)],
                   gene_tab[pr$feature, .(gene_id, gene, biotype, product, n_go_terms)])
    if (set == "shared") pairs <- cbind(pairs, dml[pr$cpg, .(meth.diff_ploidy, direction_ploidy, class)])
    setorder(pairs, chr, start, gene_id)
    write_table(pairs, tab(sprintf("dml_gene_overlap_%s_%s.csv", set, mode)))

    per_gene <- pairs[, .(n_dml = .N, n_hyper = sum(meth.diff > 0), n_hypo = sum(meth.diff < 0),
                          min_meth_diff = min(meth.diff), max_meth_diff = max(meth.diff),
                          dml_positions = paste(start, collapse = ";")), by = gene_id]
    gene_tab[, (paste0("n_dml_", set)) := per_gene[.SD, on = "gene_id", n_dml]]
    set(gene_tab, which(is.na(gene_tab[[paste0("n_dml_", set)]])), paste0("n_dml_", set), 0L)
    with_dml <- merge(gene_tab[, .(gene_id, gene, biotype, product, chr, gene_start = start, gene_end = end, strand,
                                   n_transcripts, n_go_terms, in_legacy_universe, n_cpg_tested)],
                      per_gene, by = "gene_id")
    setorder(with_dml, -n_dml, gene_id)
    write_table(with_dml, tab(sprintf("genes_with_dml_%s_%s.csv", set, mode)))

    key <- paste(pairs$chr, pairs$start)
    pc <- pairs$biotype == "protein_coding"
    row <- data.table(
      mode = mode, set = set, n_dml = nrow(dml),
      n_dml_in_gene = uniqueN(key), n_dml_in_protein_coding_gene = uniqueN(key[pc]),
      n_dml_in_lncRNA_gene = uniqueN(key[pairs$biotype == "lncRNA"]),
      n_dml_in_two_or_more_genes = sum(table(key) > 1),
      n_dml_gene_pairs = nrow(pairs),
      n_genes_with_dml = uniqueN(pairs$gene_id), n_protein_coding_genes_with_dml = uniqueN(pairs$gene_id[pc]),
      n_genes_with_dml_and_go = uniqueN(pairs$gene_id[pairs$n_go_terms > 0]),
      n_genes_with_dml_and_go_in_legacy_universe = with_dml[n_go_terms > 0 & in_legacy_universe == TRUE, .N])

    # ---- legacy ten-column annotation table (legacy mode only) ---------------
    if (mode == "legacy") {
      via_tx <- unique(merge(pairs[, .(gene_id, chr, start, end, meth.diff)],
                             tx[legacy_gene_id_lost == FALSE, .(gene_id, transcript_id)],
                             by = "gene_id", allow.cartesian = TRUE))
      ann <- merge(via_tx, go_full, by.x = "transcript_id", by.y = "map_id", allow.cartesian = TRUE)
      ann <- merge(ann, legacy_terms, by.x = "go_id", by.y = "GOID", allow.cartesian = TRUE)
      lines <- ann[, paste(transcript_id, gene_id, chr, start, end, format(meth.diff, digits = 15, trim = TRUE),
                           go_id, GOterm, GOSlim, GOcat, sep = "\t")]
      lines <- sort(unique(lines), method = "radix")          # byte order, as the legacy `sort`
      con <- gzfile(tab(sprintf("dml_go_annotation_%s_legacy.tsv.gz", set)), "w"); writeLines(lines, con); close(con)
      row <- cbind(row, data.table(
        legacy_dml_transcript_rows = nrow(via_tx), legacy_transcripts = uniqueN(via_tx$transcript_id),
        legacy_annotation_rows = length(lines), legacy_transcripts_with_go = uniqueN(ann$transcript_id),
        legacy_genes_with_go = uniqueN(ann$gene_id),
        legacy_go_ids_BP = ann[GOcat == "P", uniqueN(go_id)], legacy_go_ids_CC = ann[GOcat == "C", uniqueN(go_id)],
        legacy_go_ids_MF = ann[GOcat == "F", uniqueN(go_id)]))
    }
    summary_rows[[length(summary_rows) + 1]] <- row
    log_msg("[%s] %-6s %4d DML; %4d in a gene; %4d genes, %4d with GO terms", mode, set, nrow(dml),
            row$n_dml_in_gene, row$n_genes_with_dml, row$n_genes_with_dml_and_go)
  }
  setorder(gene_tab, chr, start)
  write_table(gene_tab, tab(sprintf("gene_table_%s.csv.gz", mode)))
  log_msg("[%s] genes with at least one tested CpG: %d of %d; with GO terms as well: %d", mode,
          sum(gene_tab$n_cpg_tested > 0), nrow(gene_tab), sum(gene_tab$n_cpg_tested > 0 & gene_tab$n_go_terms > 0))
}
write_table(rbindlist(summary_rows, fill = TRUE), tab("annotation_summary.csv"))

finish_step(step, t0)

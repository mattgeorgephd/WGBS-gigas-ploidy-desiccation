# check_legacy_outputs.R
#
# Purpose   Compare what the pipeline regenerates, in legacy mode, with the
#           outputs of the legacy analysis kept in 08_archive/. Every comparison
#           becomes one row of 06_diagnostics/legacy_check_results.csv with the
#           legacy value, the regenerated value, and a status.
#
# Run       Rscript 06_diagnostics/check_legacy_outputs.R
#           Checks whose pipeline outputs do not exist yet are reported as
#           "not run" rather than failing.
#
# Reads     04_output/<step>/tables/*            (pipeline outputs)
#           08_archive/...                       (legacy outputs; never modified)
#           03_input_files/annotation/legacy/    (legacy tables that exist only on gannet)
# Writes    06_diagnostics/legacy_check_results.csv

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "readxl"))

results <- list()
record <- function(check, item, legacy, regenerated, tolerance = 0, note = "") {
  status <- if (is.na(regenerated) || is.na(legacy)) "not run" else
    if (isTRUE(abs(as.numeric(legacy) - as.numeric(regenerated)) <= tolerance)) "match" else "differs"
  results[[length(results) + 1]] <<- data.table(check = check, item = item, legacy = as.character(legacy),
                                                regenerated = as.character(regenerated), status = status, note = note)
}
have <- function(step, file) file.exists(here("04_output", step, "tables", file))
out  <- function(step, file) here("04_output", step, "tables", file)

# ---- Step 03: ELISA values and model table ------------------------------------------
if (have("03_global_methylation", "elisa_percent_5mC.csv")) {
  new <- fread(out("03_global_methylation", "elisa_percent_5mC.csv"))
  old <- as.data.table(read_excel(archive_path("data", "ELISA-percent-mCPG-summary.xlsx"), sheet = "summary"))
  m <- merge(new[, .(sample = paste0(animal_id, "-C"), percent_5mC)], old[, .(sample, percent)], by = "sample")
  record("ELISA percent 5-mC", "animals matched", nrow(old), nrow(m))
  record("ELISA percent 5-mC", "largest per-animal difference (percentage points)", 0,
         signif(max(abs(m$percent_5mC - m$percent)), 3), tolerance = 0.001,
         note = "legacy sheet used constants rounded to 0.224 and 0.6105")
  aov_new <- fread(out("03_global_methylation", "anova_by_scale.csv"))[scale == "arcsine" & grepl("^ELISA", assay)]
  aov_old <- read.table(archive_path("stats", "AOV_5-mC.csv"), header = TRUE, check.names = FALSE)
  rownames(aov_old) <- trimws(rownames(aov_old))
  for (tm in c("ploidy", "treatment", "ploidy:treatment")) {
    record("ELISA two-way ANOVA (arcsine scale)", paste("F,", tm), round(aov_old[tm, "F value"], 2),
           round(aov_new[term == tm, F_sequential], 2), tolerance = 0.005,
           note = "regenerated from raw absorbance with unrounded constants")
  }
  mc_new <- fread(out("03_global_methylation", "anova_by_scale.csv"))[scale == "arcsine" & grepl("before", assay)]
  mc_old <- read.table(archive_path("stats", "AOV_%mCpG.csv"), header = TRUE, check.names = FALSE)
  rownames(mc_old) <- trimws(rownames(mc_old))
  record("Percent mCpG two-way ANOVA (arcsine scale)", "F, ploidy:treatment, counts before deduplication",
         round(mc_old["ploidy:treatment", "F value"], 2),
         round(mc_new[term == "ploidy:treatment", F_sequential], 2), tolerance = 0.005,
         note = "legacy sheet mixed two alignment runs (issue I-13); exact agreement is not expected")
} else record("ELISA percent 5-mC", "step 03 outputs", NA, NA)

# ---- Step 04: CpG counts by coverage ---------------------------------------------------
if (have("04_cpg_coverage", "cpg_counts_by_coverage.csv")) {
  new <- fread(out("04_cpg_coverage", "cpg_counts_by_coverage.csv"))[mode == "legacy"]
  old <- as.data.table(read_excel(archive_path("analyses", "coverage_summary.xlsx"), sheet = "total_cpg"))
  old[, sample_set := fifelse(ploidy == "both", "all", ploidy)]
  m <- merge(old, new, by.x = c("sample_set", "coverage"), by.y = c("sample_set", "min_coverage"))
  for (i in seq_len(nrow(m))) {
    record("CpGs covered in all libraries of a set", sprintf("%s, %dx", m$sample_set[i], m$coverage[i]),
           m$cpg[i], m$cpg_in_all_libraries[i])
  }
} else record("CpGs covered in all libraries of a set", "step 04 outputs", NA, NA)

# ---- Step 05: DML lists ---------------------------------------------------------------
if (have("05_differential_methylation", "dml_ploidy_legacy.csv")) {
  compare_dml <- function(contrast, legacy_file, label) {
    new <- fread(out("05_differential_methylation", sprintf("dml_%s_legacy.csv", contrast)))
    old <- fread(archive_path("tables", legacy_file)); setnames(old, 1, "legacy_row")
    record(label, "number of DML", nrow(old), nrow(new))
    record(label, "hyper", sum(old$meth.diff > 0), sum(new$meth.diff > 0))
    record(label, "hypo", sum(old$meth.diff < 0), sum(new$meth.diff < 0))
    m <- merge(old, new, by = c("chr", "start"), suffixes = c("_old", "_new"))
    record(label, "loci present in both lists", nrow(old), nrow(m))
    if (nrow(m) > 0) {
      record(label, "largest difference in meth.diff (percentage points)", 0,
             signif(max(abs(m$meth.diff_old - m$meth.diff_new)), 3), tolerance = 1e-6)
      record(label, "largest relative difference in q-value", 0,
             signif(max(abs(m$qvalue_old - m$qvalue_new) / m$qvalue_old), 3), tolerance = 1e-4,
             note = "SLIM q-values agree to five significant digits across methylKit versions")
      if ("united_row" %in% names(m)) record(label, "loci at the same row of the united table", nrow(m),
                                             sum(m$legacy_row == m$united_row))
    }
  }
  compare_dml("ploidy", "DML-getMethylDiff-PLOIDY-Cov10-20.csv", "Ploidy DML (20%, q < 0.01)")
  compare_dml("heat", "DML-getMethylDiff-HEAT-Cov10-20.csv", "Heat DML (20%, q < 0.01)")
  sh_new <- fread(out("05_differential_methylation", "dml_shared_legacy.csv"))
  sh_old <- fread(archive_path("tables", "DML-getMethylDiff-INTERACTION-Cov10-20.csv"))
  record("Loci in both lists (legacy 'interaction')", "number", nrow(sh_old), nrow(sh_new))
  record("Loci in both lists (legacy 'interaction')", "present in both versions", nrow(sh_old),
         nrow(merge(sh_old, sh_new, by = c("chr", "start"))))
  record("Loci in both lists (legacy 'interaction')", "both hyper (manuscript value)", 28, sum(sh_new$class == "both_hyper"))
  record("Loci in both lists (legacy 'interaction')", "both hypo (manuscript value)", 49, sum(sh_new$class == "both_hypo"))
  record("Loci in both lists (legacy 'interaction')", "opposite sign (manuscript value)", 123, sum(sh_new$class == "opposite"))

  sw_new <- fread(out("05_differential_methylation", "dml_counts_by_threshold_legacy.csv"))
  sw_old <- as.data.table(read_excel(archive_path("analyses", "DML_difference_threshold_comparison.xlsx"), sheet = "plot"))
  sw_old[, set := c(desiccation = "heat", ploidy = "ploidy", interaction = "shared")[trt]]
  m <- merge(sw_old, sw_new, by = c("set", "difference"))
  for (i in seq_len(nrow(m))) {
    record("DML count by difference threshold", sprintf("%s, %d%%", m$set[i], m$difference[i]), m$DML[i], m$n[i])
    record("DML count by difference threshold", sprintf("%s, %d%%, hyper", m$set[i], m$difference[i]), m$hyper[i], m$n_hyper[i],
           note = if (m$set[i] == "heat" && m$difference[i] == 50)
             "legacy workbook has hyper and hypo swapped at this threshold (1 and 3); the legacy script comments say 3 and 1" else "")
  }
} else record("DML lists", "step 05 outputs", NA, NA)

# ---- Step 07: feature overlap counts of the DML sets ------------------------------------
if (have("07_genome_features", "feature_overlap_by_track_legacy.csv")) {
  new <- fread(out("07_genome_features", "feature_overlap_by_track_legacy.csv"))
  old <- as.data.table(read_excel(archive_path("analyses", "DML_locations-all-datasets.xlsx"), sheet = "data",
                                  col_names = FALSE, .name_repair = "minimal"))
  # Row 2 of the sheet names the columns; columns 11 to 16 are the 22-library DML sets.
  cols <- c(`ploidy DML` = 11, `heat DML` = 12, shared = 13, `shared, both hypo` = 14, `shared, both hyper` = 15, `shared, opposite` = 16)
  track_of <- c(exon_UTR = "exon_UTR", CDS = "CDS", intron = "intron", upstream_flank = "upstream_flank",
                downstrea_flank = "downstream_flank", lncRNA = "lncRNA", TE = "TE", intergenic = "intergenic",
                gene = "gene", exon = "exon", mRNA = "mRNA", nonCDS = "nonCDS", flanks = "flanks")
  for (i in which(old[[1]] %in% names(track_of))) for (st in names(cols)) {
    record("DML overlapping a genome feature track", sprintf("%s, %s", st, track_of[[old[[1]][i]]]),
           as.numeric(old[[cols[[st]]]][i]), new[set == st & track == track_of[[old[[1]][i]]], n_overlap])
  }
} else record("DML overlapping a genome feature track", "step 07 outputs", NA, NA)

# ---- Step 08: annotation tables against the legacy intermediates ------------------------
legacy_set <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")
if (have("08_gene_annotation", "transcript_table.csv.gz")) {
  tx <- fread(out("08_gene_annotation", "transcript_table.csv.gz"), colClasses = c(gene_id = "character"))
  leg <- fread(input_path("annotation", "legacy", "cgigas_uk_roslin_v1_mRNA.transcriptID.geneID.gz"), header = FALSE,
               col.names = c("gene_legacy", "transcript_id"))
  m <- merge(tx, leg, by = "transcript_id")
  record("Transcript to gene map", "transcripts", nrow(leg), nrow(m))
  record("Transcript to gene map", "same gene ID", sum(!tx$legacy_gene_id_lost), sum(m$gene_legacy == m$gene_id),
         note = "the others are mRNAs of genes with a hyphen in their name, which the legacy parsing assigned to the word GeneID")
  for (s in names(legacy_set)) {
    label <- sprintf("DML annotation table, %s", s)
    new <- read_dml_annotation(out("08_gene_annotation", sprintf("dml_go_annotation_%s_legacy.tsv.gz", s)))
    old_f <- input_path("annotation", "legacy", sprintf("DML-Cov10-20-BOTH-%s.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz", legacy_set[[s]]))
    old <- read_dml_annotation(old_f)
    record(label, "lines", nrow(old), nrow(new))
    if (nrow(old) == nrow(new)) {
      text_cols <- c("transcript", "geneID", "chr", "start", "end", "GOID", "GOterm", "GOSlim", "GOcat")
      record(label, "text columns identical, line for line (of 9)", 9,
             sum(vapply(text_cols, function(k) identical(as.character(new[[k]]), as.character(old[[k]])), logical(1))))
      record(label, "largest difference in meth.diff", 0, signif(max(abs(new$meth.diff - old$meth.diff)), 3), tolerance = 1e-7,
             note = "the legacy files hold 10 significant digits for ploidy and shared")
    }
    old_tx <- fread(input_path("annotation", "legacy", sprintf("DML-Cov10-20-BOTH-%s.geneOverlap.transcriptIDs.gz", legacy_set[[s]])),
                    header = FALSE, col.names = c("transcript", "geneID", "chr", "start", "end", "meth.diff"))
    ann_sum <- fread(out("08_gene_annotation", "annotation_summary.csv"))[mode == "legacy" & set == s]
    record(label, "DML x transcript rows before the GO join", nrow(old_tx), ann_sum$legacy_dml_transcript_rows)
    # What the legacy script's reader returned; the legacy value is the row count noted in that script.
    noted <- c(ploidy = 10979, heat = 81109, shared = 6941)
    record(label, "rows returned by the legacy reader", noted[[s]],
           nrow(suppressWarnings(read_dml_annotation(out("08_gene_annotation", sprintf("dml_go_annotation_%s_legacy.tsv.gz", s)), "legacy"))),
           note = "legacy value from comments in 5_WGBS_Functional-enrichement-analysis-BOTH.R")
  }
} else record("DML annotation tables", "step 08 outputs", NA, NA)

# ---- Step 09: topGO tables, Figure 5, Figure 9, gene lists ------------------------------
# The legacy tables were written without quoting, and term names hold commas, so
# they are parsed from both ends of each line.
read_legacy_fisher <- function(f) {
  l <- readLines(f)[-1]
  m <- regmatches(l, regexec("^(GO:[0-9]+),(.*),([0-9]+),([0-9]+),([0-9.eE+-]+),([<> 0-9.eE+-]+)$", l))
  stopifnot(all(lengths(m) == 7))
  d <- as.data.table(do.call(rbind, m)[, -1, drop = FALSE]); setnames(d, c("go_id", "term", "annotated", "significant", "expected", "classic"))
  d[, `:=`(annotated = as.integer(annotated), significant = as.integer(significant), rank = .I)]
  d
}
# The values typed into the legacy plotting script (first data frame: Figure 5; second: Figure 9).
legacy_plot_tables <- function() {
  l <- readLines(archive_path("code", "gene-plots.Rmd"))
  starts <- grep("^data <- data.frame\\(", l)
  lapply(starts, function(i) {
    j <- i + which(trimws(l[(i + 1):length(l)]) == ")")[1]
    e <- new.env(); eval(parse(text = l[i:j]), e); as.data.table(e$data)[treatment == "control + heat"]
  })
}
if (have("09_go_enrichment_legacy", "figure_5_values_legacy.csv")) {
  plot_tabs <- legacy_plot_tables()
  set_of <- c(ploidy = "ploidy", heat = "heat", interaction = "shared")
  for (s in names(legacy_set)) for (o in c("BP", "CC", "MF")) {
    label <- sprintf("topGO table, %s %s", s, o)
    new <- fread(out("09_go_enrichment_legacy", sprintf("legacy_fisher_%s_%s.csv", s, o)))
    old <- read_legacy_fisher(archive_path("GOterms", paste0("BOTH-", legacy_set[[s]]), sprintf("%s-%s-FisherTestResults.csv", legacy_set[[s]], o)))
    m <- merge(old, new, by = "go_id", suffixes = c("_old", "_new"))
    record(label, "terms", nrow(old), nrow(new))
    record(label, "terms in both", nrow(old), nrow(m))
    record(label, "terms with the same annotated and significant counts", nrow(old),
           sum(m$annotated_old == m$annotated_new & m$significant_old == m$significant_new))
    record(label, "ten top-ranked terms in common", 10, length(intersect(old$go_id[1:10], new[top10 == TRUE, go_id])))
    p_old <- suppressWarnings(as.numeric(old$classic))
    record(label, "terms with p < 0.01", sum(p_old < 0.01 | grepl("<", old$classic)), sum(new$below_p_cutoff))
    # Gene list behind Figure 9 and Table 1
    g_new <- fread(out("09_go_enrichment_legacy", sprintf("legacy_dml_in_top10_genes_%s_%s.csv", s, o)))
    g_old <- fread(archive_path("GOterms", paste0("BOTH-", legacy_set[[s]]), sprintf("GENELIST-%s-EnrichedGO-DMLOnly.csv", o)))
    record(sprintf("DML in genes with a top-ten term, %s %s", s, o), "rows", nrow(g_old), nrow(g_new))
    record(sprintf("DML in genes with a top-ten term, %s %s", s, o), "rows in both (gene, position)", nrow(g_old),
           nrow(merge(g_old[, .(geneID, chr, start)], g_new[, .(geneID, chr, start)], by = c("geneID", "chr", "start"))))
  }
  f5 <- fread(out("09_go_enrichment_legacy", "figure_5_values_legacy.csv"))
  o5 <- melt(plot_tabs[[1]], id.vars = c("treatment", "factor", "variable"), variable.name = "ontology", value.name = "legacy")
  for (i in seq_len(nrow(o5))) {
    record("Figure 5 value", sprintf("%s, %s, %s", o5$variable[i], o5$factor[i], o5$ontology[i]), o5$legacy[i],
           f5[published_label == o5$variable[i] & set == set_of[[o5$factor[i]]] & ontology == o5$ontology[i], value],
           note = "legacy value typed in gene-plots.Rmd")
  }
  f9 <- fread(out("09_go_enrichment_legacy", "figure_9_values_legacy.csv"))
  o9 <- plot_tabs[[2]]
  for (i in seq_len(nrow(o9))) for (k in c("total", "hyper", "hypo")) {
    record("Figure 9 value", sprintf("%s, %s, %s", k, o9$factor[i], o9$type[i]), o9[[k]][i],
           f9[set == set_of[[o9$factor[i]]] & ontology == o9$type[i]][[k]], note = "legacy value typed in gene-plots.Rmd")
  }
} else record("topGO tables", "step 09 outputs", NA, NA)

res <- rbindlist(results)
fwrite(res, here("06_diagnostics", "legacy_check_results.csv"))
cat("\nLegacy reproduction checks:\n"); print(res[, .N, by = .(check, status)])
if (any(res$status == "differs")) { cat("\nDifferences:\n"); print(res[status == "differs"]) }

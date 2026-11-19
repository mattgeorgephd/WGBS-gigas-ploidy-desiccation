# Checks of what each number of Figure 5, Figure 9 and Table 1 counts, using the literal port (A_port.rds),
# the archived legacy outputs (08_archive/GOterms), the typed values (gene-plots.Rmd) and step 09's tables.
.libPaths(c("/home/user/gannet_cache/legacy_r_lib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(GO.db) })
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scratch, "A_port.rds"))
out <- P$out; L <- lapply(P$DMLGOterms, as.data.table); E <- lapply(P$DMLGOterms_exact, as.data.table)
geneID2GO <- P$geneID2GO
sets <- c(PLOIDY = "ploidy", HEAT = "heat", COMBO = "shared"); onts <- c(BP = "P", CC = "C", MF = "F")
hdr <- function(x) cat("\n==================== ", x, " ====================\n", sep = "")

# ---------------------------------------------------------------------------------------------
hdr("1. Typed values of gene-plots.Rmd (rows 'control + heat') vs literal port vs step 09")
BP <- c(34086, 514086, 50696, 14503, 4379, 39688, 3332, 950, 17007, 412, 19517, 128, 93, 55, 29, 28, 8, 1, 1, 3, 1, 2736, 68452, 1236, 109458, 199, 1441, 165)
CC <- c(21772, 375681, 39736, 9588, 4145, 26084, 2337, 2705, 31297, 12821, 39913, 601, 727, 610, 21, 1313, 348, 331, 12, 17, 10, 3407, 259092, 20617, 39913, 1420, 1698, 1145)
MF <- c(14724, 235261, 23665, 5706, 2443, 15213, 1262, 1358, 22189, 1548, 13304, 178, 1841, 207, 52, 71, 27, 23, 6, 49, 4, 2288, 61978, 3861, 37810, 392, 5801, 318)
treatment <- c("control", "heat", "heat", "heat", "control + heat", "control + heat", "control + heat")
fac <- c("ploidy", "ploidy", "heat", "interaction", "ploidy", "heat", "interaction")
variable <- rep(c("GO terms", "Unique transcripts", "Unique genes", "Gene products"), each = 7)
typed5 <- data.table(treatment = rep(treatment, 4), factor = rep(fac, 4), variable, BP, CC, MF)[treatment == "control + heat"]
typed5 <- melt(typed5, id.vars = c("treatment", "factor", "variable"), variable.name = "ont", value.name = "typed")
typed5[, set := c(ploidy = "PLOIDY", heat = "HEAT", interaction = "COMBO")[factor]]
cnt <- fread(file.path(scratch, "A_counts.csv"))
port5 <- melt(cnt[, .(set, ont, `GO terms` = GO_terms_panel, `Unique transcripts` = unique_transcripts_panel,
                      `Unique genes` = unique_genes_panel, `Gene products` = gene_products_panel)],
              id.vars = c("set", "ont"), variable.name = "variable", value.name = "port")
s9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables/figure_5_values_legacy.csv"))
s9[, set := c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")[set]]
m5 <- merge(merge(typed5, port5, by = c("set", "ont", "variable")), s9[, .(set, ont = ontology, variable = published_label, step09 = value)],
            by = c("set", "ont", "variable"))
cat("Figure 5: values compared:", nrow(m5), "; typed == port:", sum(m5$typed == m5$port), "; typed == step09:", sum(m5$typed == m5$step09), "\n")

total <- c(59, 10766, 43, 568, 10, 7, 3, 147, 16605, 1054, 953, 24, 54, 16, 106, 12993, 125, 385, 15, 134, 8)
hyper <- c(52, 1597, 29, 277, 5, 2, 2, 124, 2560, 633, 428, 12, 30, 8, 87, 1827, 66, 170, 5, 87, 5)
hypo <- c(7, 9169, 14, 291, 5, 5, 1, 23, 14045, 421, 525, 12, 24, 8, 19, 11166, 59, 215, 10, 47, 3)
typed9 <- data.table(treatment = rep(treatment, 3), factor = rep(fac, 3), ont = rep(c("BP", "CC", "MF"), each = 7), total, hyper, hypo)[treatment == "control + heat"]
typed9[, set := c(ploidy = "PLOIDY", heat = "HEAT", interaction = "COMBO")[factor]]
s9b <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables/figure_9_values_legacy.csv"))
s9b[, set := c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")[set]]
m9 <- merge(merge(typed9, cnt[, .(set, ont, fig9_total, fig9_hyper, fig9_hypo, fig9_strict_hyper, fig9_strict_hypo)], by = c("set", "ont")),
            s9b[, .(set, ont = ontology, s9_total = total, s9_hyper = hyper, s9_hypo = hypo, s9_genes = genes)], by = c("set", "ont"))
cat("Figure 9: 27 typed values (total, hyper, hypo x 9): typed == port:",
    sum(m9$total == m9$fig9_total) + sum(m9$hyper == m9$fig9_hyper) + sum(m9$hypo == m9$fig9_hypo),
    "; typed == step09:", sum(m9$total == m9$s9_total) + sum(m9$hyper == m9$s9_hyper) + sum(m9$hypo == m9$s9_hypo), "\n")
cat("  any meth.diff exactly 0 (would be counted in both hyper and hypo):",
    sum(m9$fig9_hyper != m9$fig9_strict_hyper) + sum(m9$fig9_hypo != m9$fig9_strict_hypo), "\n")
print(m9[, .(set, ont, total, hyper, hypo, s9_genes)])
cat("Shares of the sum over the three ontologies (as Results 3.5 uses them):\n")
print(m9[, .(ont, total, share_percent = round(100 * total / sum(total), 2), sum_total = sum(total)), by = set])

# ---------------------------------------------------------------------------------------------
hdr("2. Top ten: port vs archived Fisher tables vs step 09; ties at the boundary; p < 0.01")
arch_dir <- c(PLOIDY = "BOTH-PLOIDY", HEAT = "BOTH-HEAT", COMBO = "BOTH-COMBO")
for (s in names(sets)) for (o in names(onts)) {
  k <- paste(s, o); r <- out[[k]]
  a <- readLines(file.path(repo, "08_archive/GOterms", arch_dir[[s]], sprintf("%s-%s-FisherTestResults.csv", s, o)))
  a_go <- sub(",.*$", "", a[-1]); a_p <- sub("^.*,", "", a[-1]); a_pn <- as.numeric(sub("< ", "", a_p))
  # annotated, significant, expected are the three fields before the last
  parts <- strsplit(a[-1], ",", fixed = TRUE)
  a_tail <- t(vapply(parts, function(v) tail(v, 4), character(4)))
  f9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_fisher_%s_%s.csv", sets[[s]], o)))
  pv <- r$pval[r$allRes$GO.ID]
  cat(sprintf("%-6s %s: terms port %d / archive %d / step09 %d | same order as archive (all rows): %s | top10 set port==archive: %s, step09==archive: %s | counts equal (Annotated, Significant, Expected): %s | p<0.01 port %d, archive (printed p) %d, step09 %d | p rank10 %.4g rank11 %.4g tie: %s\n",
              s, o, nrow(r$allRes), length(a_go), nrow(f9),
              identical(r$allRes$GO.ID, a_go), setequal(r$top10, a_go[1:10]), setequal(f9[rank <= 10, go_id], a_go[1:10]),
              all(as.character(r$allRes$Annotated) == a_tail[, 1] & as.character(r$allRes$Significant) == a_tail[, 2] &
                    as.numeric(r$allRes$Expected) == as.numeric(a_tail[, 3])),
              sum(pv < 0.01), sum(a_pn < 0.01), sum(f9$p_fisher < 0.01), pv[10], pv[11], isTRUE(all.equal(pv[[10]], pv[[11]]))))
  if (isTRUE(all.equal(pv[[10]], pv[[11]]))) {
    tied <- names(pv)[abs(pv - pv[[10]]) < 1e-12 * pv[[10]] | pv == pv[[10]]]
    g <- genesInTermList <- NULL
    cat("       terms tied with rank 10:", paste(tied, collapse = " "), "| ranks", paste(match(tied, names(pv)), collapse = " "), "\n")
  }
}

# ---------------------------------------------------------------------------------------------
hdr("3. Content: port vs archived legacy outputs vs step 09")
for (s in names(sets)) for (o in names(onts)) {
  k <- paste(s, o); r <- out[[k]]
  # archived GENELIST (DML only)
  g <- fread(file.path(repo, "08_archive/GOterms", arch_dir[[s]], sprintf("GENELIST-%s-EnrichedGO-DMLOnly.csv", o)))
  keyf <- function(d) sort(paste(d$geneID, d$chr, d$start, d$end, sprintf("%.6f", d$meth.diff)))
  s9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_dml_in_top10_genes_%s_%s.csv", sets[[s]], o)))
  # archived withTranscript: written with quote = FALSE, so parse by pattern
  w <- readLines(file.path(repo, "08_archive/GOterms", arch_dir[[s]], sprintf("all-%s-EnrichedGO-DML-withTranscript.csv", o)))[-1]
  w_go <- sub(",.*$", "", w)
  w_tx <- regmatches(w, regexpr("(?<=,)[XN][MR]_[0-9]+\\.[0-9]+(?=,)", w, perl = TRUE))
  w_gene <- sub("^.*,[XN][MR]_[0-9]+\\.[0-9]+,([^,]*),.*$", "\\1", w)
  a_pairs <- sort(paste(w_go, w_tx, w_gene)); p_pairs <- sort(paste(r$annot$GOID, r$annot$transcript, r$annot$geneID))
  t9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_top10_term_transcripts_%s_%s.csv", sets[[s]], o)))
  s_pairs <- sort(paste(t9$GOID, t9$transcript, t9$geneID))
  cat(sprintf("%-6s %s: DML list rows port %d archive %d step09 %d; identical port==archive %s, step09==archive %s | term-transcript rows port %d archive %d step09 %d; identical (GOID,transcript,geneID) port==archive %s step09==archive %s\n",
              s, o, nrow(r$methdiff), nrow(g), nrow(s9), identical(keyf(r$methdiff), keyf(g)), identical(keyf(s9), keyf(g)),
              nrow(r$annot), length(w), nrow(t9), identical(a_pairs, p_pairs), identical(a_pairs, s_pairs)))
}

# ---------------------------------------------------------------------------------------------
hdr("4. What a row of the DML annotation table is; the 'GO terms' panel")
uni <- unique(names(geneID2GO))
for (s in names(sets)) {
  l <- L[[s]]; e <- E[[s]]
  cat(sprintf("%-6s lines in file %d; rows legacy reader %d (%.1f%%)\n", s, nrow(e), nrow(l), 100 * nrow(l) / nrow(e)))
  for (o in names(onts)) {
    lo <- l[GOcat == onts[[o]]]; eo <- e[GOcat == onts[[o]]]
    cat(sprintf("   %s rows: legacy %d (panel 'GO terms'), exact %d | distinct GO IDs legacy %d, exact %d | distinct (DML,transcript,GOID) legacy %d -> rows per combination max %d | distinct DML legacy %d | rows of transcripts NOT in the legacy universe: %d (%d transcripts)\n",
                o, nrow(lo), nrow(eo), uniqueN(lo$GOID), uniqueN(eo$GOID),
                uniqueN(lo[, .(chr, start, transcript, GOID)]), max(lo[, .N, by = .(chr, start, transcript, GOID)]$N),
                uniqueN(lo[, .(chr, start)]), sum(!lo$transcript %chin% uni), uniqueN(lo$transcript[!lo$transcript %chin% uni])))
  }
  cat(sprintf("   GOcat values in legacy-read table: %s\n", paste(names(table(l$GOcat, useNA = "ifany")), table(l$GOcat, useNA = "ifany"), collapse = "; ")))
}

# ---------------------------------------------------------------------------------------------
hdr("5. 'Unique transcripts', 'Unique genes', 'Gene products': decomposition")
txg <- fread("/home/user/gannet_cache/legacy_annot/analyses_GOterm-annotation_cgigas_uk_roslin_v1_mRNA.transcriptID.geneID",
             header = FALSE, col.names = c("gene", "transcript"), colClasses = "character")
tx2gene <- setNames(txg$gene, txg$transcript)
anc <- list(BP = as.list(GOBPANCESTOR), CC = as.list(GOCCANCESTOR), MF = as.list(GOMFANCESTOR))
first_line <- geneID2GO[!duplicated(names(geneID2GO))]     # topGO uses the first line of a repeated name
res5 <- list()
for (s in names(sets)) for (o in names(onts)) {
  k <- paste(s, o); r <- out[[k]]; l <- L[[s]]; e <- E[[s]]
  lo <- l[GOcat == onts[[o]]]; eo <- e[GOcat == onts[[o]]]
  pairs <- data.table(GOID = rep(names(r$inTerm), lengths(r$inTerm)), transcript = unlist(r$inTerm, use.names = FALSE))
  stopifnot(nrow(pairs) == nrow(r$topGO))
  pairs[, sig := transcript %chin% l$transcript]                       # "of interest" in the legacy test
  pairs[, direct_in_map := mapply(function(g, t) g %in% first_line[[t]], GOID, transcript)]
  pairs[, in_legacy_rows := paste(GOID, transcript) %chin% paste(lo$GOID, lo$transcript)]
  pairs[, in_exact_rows := paste(GOID, transcript) %chin% paste(eo$GOID, eo$transcript)]
  pairs[, in_exact_any := transcript %chin% e$transcript]              # gene holds a DML, every line read
  pairs[, gene := tx2gene[transcript]]
  # independent closure: is the term the transcript's own term or an ancestor of one of them?
  chk <- mapply(function(g, t) { d <- first_line[[t]]; g %in% d || any(vapply(d, function(x) g %in% anc[[o]][[x]], logical(1))) },
                pairs$GOID, pairs$transcript)
  ut <- unique(pairs$transcript)
  # genes
  genes_direct_legacy <- uniqueN(lo[paste(GOID, transcript) %chin% paste(pairs$GOID, pairs$transcript), geneID])
  genes_direct_exact  <- uniqueN(eo[paste(GOID, transcript) %chin% paste(pairs$GOID, pairs$transcript), geneID])
  genes_any_sig <- uniqueN(l[transcript %chin% ut, geneID])
  genes_any_sig_ontrows <- uniqueN(lo[transcript %chin% ut, geneID])
  genes_any_exact <- uniqueN(e[transcript %chin% ut, geneID])
  # gene products: NA rows + matched rows
  an <- as.data.table(r$annot)
  n_na <- sum(is.na(an$geneID)); n_match <- sum(!is.na(an$geneID))
  mult <- an[!is.na(geneID), .(rows = .N, dml = uniqueN(paste(chr, start)), slim = uniqueN(GOSlim), genes = uniqueN(geneID)), by = .(GOID, transcript)]
  res5[[k]] <- data.table(
    set = s, ont = o, pairs = nrow(pairs), closure_ok = all(chk),
    unique_tx = length(ut), tx_sig = sum(ut %chin% l$transcript), tx_gene_has_dml_exact = sum(ut %chin% e$transcript),
    genes_of_unique_tx = uniqueN(tx2gene[ut]), tx_without_gene_in_map = sum(is.na(tx2gene[ut])),
    pairs_sig = sum(pairs$sig), sig_direct_in_rows = sum(pairs$sig & pairs$in_legacy_rows),
    sig_direct_lost_by_reader = sum(pairs$sig & !pairs$in_legacy_rows & pairs$in_exact_rows),
    sig_direct_map_not_in_table = sum(pairs$sig & pairs$direct_in_map & !pairs$in_exact_rows),
    sig_indirect = sum(pairs$sig & !pairs$direct_in_map),
    rows_direct_not_in_map = sum(pairs$in_legacy_rows & !pairs$direct_in_map),
    unique_genes_panel = genes_direct_legacy, unique_genes_exact_reader = genes_direct_exact,
    genes_with_sig_tx_in_top10 = genes_any_sig, genes_in_fig9_list = uniqueN(r$methdiff$geneID),
    genes_with_sig_tx_and_row_of_ont = genes_any_sig_ontrows, genes_any_exact_reader = genes_any_exact,
    gp_rows = nrow(an), gp_na_rows = n_na, gp_matched_rows = n_match, gp_matched_pairs = nrow(mult),
    gp_rows_equal_dml_x_slim = all(mult$rows == mult$dml * mult$slim), gp_max_rows_per_pair = if (nrow(mult)) max(mult$rows) else 0L)
}
res5 <- rbindlist(res5)
options(width = 250)
print(res5[, .(set, ont, pairs, closure_ok, unique_tx, tx_sig, tx_gene_has_dml_exact, genes_of_unique_tx, tx_without_gene_in_map)])
print(res5[, .(set, ont, pairs_sig, sig_direct_in_rows, sig_direct_lost_by_reader, sig_direct_map_not_in_table, sig_indirect, rows_direct_not_in_map)])
print(res5[, .(set, ont, unique_genes_panel, unique_genes_exact_reader, genes_with_sig_tx_in_top10, genes_with_sig_tx_and_row_of_ont, genes_in_fig9_list, genes_any_exact_reader)])
print(res5[, .(set, ont, gp_rows, gp_na_rows, gp_matched_rows, gp_matched_pairs, gp_rows_equal_dml_x_slim, gp_max_rows_per_pair)])
fwrite(res5, file.path(scratch, "B_decomposition.csv"))

# ---------------------------------------------------------------------------------------------
hdr("6. Figure 9 / Table 1 lists: DML counted twice, overlap between ontologies, effect of the reader")
res6 <- list(); allm <- list()
for (s in names(sets)) for (o in names(onts)) {
  k <- paste(s, o); r <- out[[k]]; l <- L[[s]]; e <- E[[s]]
  ut <- unique(unlist(r$inTerm, use.names = FALSE))
  m <- as.data.table(r$methdiff)
  # same definition with every line of the table read (same ten terms)
  mx <- unique(e[transcript %chin% ut & GOcat == onts[[o]], .(geneID, chr, start, end, meth.diff)])
  # without the requirement that the transcript has a row of this ontology
  mn <- unique(l[transcript %chin% ut, .(geneID, chr, start, end, meth.diff)])
  res6[[k]] <- data.table(set = s, ont = o, rows = nrow(m), distinct_dml = uniqueN(m[, .(chr, start)]), genes = uniqueN(m$geneID),
                          rows_without_GOcat_filter = nrow(mn), genes_without_GOcat_filter = uniqueN(mn$geneID),
                          rows_exact_reader_same_terms = nrow(mx), distinct_dml_exact = uniqueN(mx[, .(chr, start)]), genes_exact = uniqueN(mx$geneID))
  allm[[k]] <- cbind(set = s, ont = o, m)
}
res6 <- rbindlist(res6); allm <- rbindlist(allm)
print(res6)
cat("Per set: sum of the three ontology lists, distinct (gene, DML) rows, distinct DML, distinct genes:\n")
print(allm[, .(sum_of_three = .N, distinct_gene_dml = uniqueN(paste(geneID, chr, start)), distinct_dml = uniqueN(paste(chr, start)),
               distinct_genes = uniqueN(geneID)), by = set])
cat("DML rows of the BP list also present in the CC or MF list of the same set:\n")
print(allm[, .(in_BP = sum(ont == "BP"), BP_also_in_CC_or_MF = sum(ont == "BP" & paste(geneID, chr, start) %chin% paste(geneID, chr, start)[ont != "BP"])), by = set])
fwrite(res6, file.path(scratch, "B_fig9_lists.csv"))

# ---------------------------------------------------------------------------------------------
hdr("7. Which contrast the meth.diff of the 'interaction' table is")
d <- file.path(repo, "08_archive/tables")
h <- fread(file.path(d, "DML-getMethylDiff-HEAT-Cov10-20.csv")); p <- fread(file.path(d, "DML-getMethylDiff-PLOIDY-Cov10-20.csv"))
i <- fread(file.path(d, "DML-getMethylDiff-INTERACTION-Cov10-20.csv"))
setnames(h, 1, "row"); setnames(p, 1, "row"); setnames(i, 1, "row")
cat("rows: heat", nrow(h), "ploidy", nrow(p), "interaction", nrow(i), "\n")
mi <- merge(merge(i[, .(row, chr, start, md_int = meth.diff)], h[, .(chr, start, row_h = row, md_heat = meth.diff)], by = c("chr", "start"), all.x = TRUE),
            p[, .(chr, start, row_p = row, md_ploidy = meth.diff)], by = c("chr", "start"), all.x = TRUE)
cat("interaction loci found in heat list:", sum(!is.na(mi$md_heat)), "; in ploidy list:", sum(!is.na(mi$md_ploidy)),
    "; intersect(heat, ploidy) by chr+start:", nrow(merge(h[, .(chr, start)], p[, .(chr, start)])), "\n")
cat("max |interaction meth.diff - heat meth.diff| :", max(abs(mi$md_int - mi$md_heat)), "\n")
cat("max |interaction meth.diff - ploidy meth.diff|:", max(abs(mi$md_int - mi$md_ploidy)), "; loci where it equals the ploidy value (1e-6):", sum(abs(mi$md_int - mi$md_ploidy) < 1e-6), "\n")
cat("row label equals heat row label:", sum(mi$row == mi$row_h), "and ploidy row label:", sum(mi$row == mi$row_p), "\n")
cat("sign classes: both > 0:", sum(mi$md_heat > 0 & mi$md_ploidy > 0), "; both < 0:", sum(mi$md_heat < 0 & mi$md_ploidy < 0),
    "; opposite:", sum(sign(mi$md_heat) != sign(mi$md_ploidy)), "; heat > 0:", sum(mi$md_heat > 0), "; heat < 0:", sum(mi$md_heat < 0), "\n")
# the COMBO annotation table: which contrast are its meth.diff values?
cmb <- unique(E$COMBO[, .(chr, start, end, meth.diff)])
cm <- merge(cmb, merge(h[, .(chr, start, md_heat = meth.diff)], p[, .(chr, start, md_ploidy = meth.diff)], by = c("chr", "start")), by = c("chr", "start"), all.x = TRUE)
cat("COMBO annotation table: distinct DML", nrow(cmb), "; matched to both lists", sum(!is.na(cm$md_heat)),
    "; max |table - heat|", max(abs(cm$meth.diff - cm$md_heat)), "; max |table - ploidy|", max(abs(cm$meth.diff - cm$md_ploidy)), "\n")
t1 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables/table_1_legacy.csv"))
t1m <- merge(merge(t1, h[, .(chromosome = chr, start, md_heat = meth.diff)], by = c("chromosome", "start"), all.x = TRUE),
             p[, .(chromosome = chr, start, md_ploidy = meth.diff)], by = c("chromosome", "start"), all.x = TRUE)
setorder(t1m, dml_set, gene_id, start)
print(t1m[, .(dml_set, gene_id, gene_name = substr(gene_name, 1, 28), chromosome, start, table1 = methylation_difference,
              heat_list = round(md_heat, 2), ploidy_list = round(md_ploidy, 2))])

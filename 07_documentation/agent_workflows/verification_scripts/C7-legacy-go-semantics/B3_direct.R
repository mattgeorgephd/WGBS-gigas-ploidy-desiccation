.libPaths(c("/home/user/gannet_cache/legacy_r_lib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(GO.db) })
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
P <- readRDS(file.path(scratch, "A_port.rds"))
out <- P$out; L <- lapply(P$DMLGOterms, as.data.table); E <- lapply(P$DMLGOterms_exact, as.data.table)
geneID2GO <- P$geneID2GO; first_line <- geneID2GO[!duplicated(names(geneID2GO))]
txg <- fread("/home/user/gannet_cache/legacy_annot/analyses_GOterm-annotation_cgigas_uk_roslin_v1_mRNA.transcriptID.geneID",
             header = FALSE, col.names = c("gene", "transcript"), colClasses = "character")
tx2gene <- setNames(txg$gene, txg$transcript)
onts <- c(BP = "P", CC = "C", MF = "F")
all_goids_in_tables <- unique(unlist(lapply(E, function(d) d$GOID)))
options(width = 250)
res <- list()
for (s in names(L)) for (o in names(onts)) {
  k <- paste(s, o); r <- out[[k]]; l <- L[[s]]; e <- E[[s]]
  lo <- l[GOcat == onts[[o]]]; eo <- e[GOcat == onts[[o]]]
  pairs <- data.table(GOID = rep(names(r$inTerm), lengths(r$inTerm)), transcript = unlist(r$inTerm, use.names = FALSE))
  pairs[, sig := transcript %chin% l$transcript]
  pairs[, direct_in_map := mapply(function(g, t) g %in% first_line[[t]], GOID, transcript)]
  pairs[, in_rows := paste(GOID, transcript) %chin% paste(lo$GOID, lo$transcript)]
  pairs[, in_exact_any_cat := paste(GOID, transcript) %chin% paste(e$GOID, e$transcript)]
  pairs[, gene := tx2gene[transcript]]
  # per term summary
  term <- pairs[, .(tx = .N, sig = sum(sig), sig_direct_map = sum(sig & direct_in_map), sig_in_rows = sum(sig & in_rows),
                    genes_sig = uniqueN(gene[sig]), genes_sig_direct_map = uniqueN(gene[sig & direct_in_map]), genes_in_rows = uniqueN(gene[sig & in_rows]),
                    goid_anywhere_in_tables = GOID[1] %chin% all_goids_in_tables), by = GOID]
  term[, term := Term(GOTERM[GOID])]
  term[, p := signif(r$pval[GOID], 3)]
  res[[k]] <- data.table(set = s, ont = o,
                         genes_sig_any = uniqueN(pairs[sig == TRUE, gene]),
                         genes_sig_direct_in_map = uniqueN(pairs[sig & direct_in_map, gene]),
                         genes_panel = uniqueN(pairs[sig & in_rows, gene]),
                         top10_terms_absent_from_all_tables = sum(!term$goid_anywhere_in_tables),
                         top10_terms_with_no_row = sum(term$sig_in_rows == 0))
  if (o == "BP") { cat("\n----", k, "\n"); print(term[, .(GOID, term = substr(term, 1, 55), p, tx, sig, sig_direct_map, sig_in_rows, genes_sig, genes_sig_direct_map, genes_in_rows, goid_anywhere_in_tables)]) 
    cat("genes behind significant transcripts:\n"); print(pairs[sig == TRUE, .(terms = uniqueN(GOID), transcripts = uniqueN(transcript), direct_terms = uniqueN(GOID[direct_in_map]), terms_in_rows = uniqueN(GOID[in_rows])), by = gene]) }
}
cat("\n==== genes: any significant transcript in a top-ten term / directly annotated in the transcript-to-GO map / counted by the panel\n")
print(rbindlist(res))
# is the missing-term effect a property of the legacy GO term table (2017 snapshot)?
lt <- fread("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/annotation/legacy_go_goslim_terms.tsv.gz", quote = "", sep = "\t")
cat("\nlegacy GO/GOslim term table: rows", nrow(lt), "distinct GO IDs", uniqueN(lt$GOID), "\n")
allgo_map <- unique(unlist(first_line))
cat("distinct GO IDs in the universe map:", length(allgo_map), "; of these in the legacy term table:", sum(allgo_map %chin% lt$GOID), "\n")
for (s in names(L)) for (o in "BP") { r <- out[[paste(s, o)]]; cat(s, o, "top10 in legacy term table:", paste(r$top10, r$top10 %chin% lt$GOID), "\n") }

suppressPackageStartupMessages(library(data.table))
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scratch, "A_port.rds")); out <- P$out
txg <- fread("/home/user/gannet_cache/legacy_annot/analyses_GOterm-annotation_cgigas_uk_roslin_v1_mRNA.transcriptID.geneID", header = FALSE, col.names = c("gene", "transcript"), colClasses = "character")
tx2gene <- setNames(txg$gene, txg$transcript)
lab <- c(PLOIDY = "ploidy", HEAT = "heat", COMBO = "shared")
all <- rbindlist(lapply(names(out), function(k) { r <- out[[k]]; s <- sub(" .*", "", k); o <- sub(".* ", "", k)
  g <- vapply(r$inTerm, function(v) uniqueN(tx2gene[v[v %chin% r$sig]]), integer(1))
  f9 <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_fisher_%s_%s.csv", lab[[s]], o)))[rank <= 10]
  data.table(set = s, ont = o, go = names(g), genes_port = g, genes_step09 = f9$genes_significant[match(names(g), f9$go_id)]) }))
cat("top-ten terms:", nrow(all), "; port == step 09 (genes behind significant transcripts):", sum(all$genes_port == all$genes_step09),
    "; resting on one gene:", sum(all$genes_port == 1), "; on two or fewer:", sum(all$genes_port <= 2), "\n")
print(all[, .(one_gene = sum(genes_port == 1), two_or_fewer = sum(genes_port <= 2), max_genes = max(genes_port)), by = .(set, ont)])

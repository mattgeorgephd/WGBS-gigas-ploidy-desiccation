source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
mode <- "legacy"
gene_go <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
gene2go_all <- split(gene_go$go_id, gene_go$gene_id)
gt_all <- fread(file.path(repo, sprintf("04_output/08_gene_annotation/tables/gene_table_%s.csv.gz", mode)), colClasses = c(gene_id = "character"))
cat("all genes:", nrow(gt_all), "; biotypes:\n"); print(table(gt_all$biotype))
cat(sprintf("genes with a tested CpG: %d (all biotypes), %d protein coding; with an eligible CpG: %d (all), %d protein coding\n", sum(gt_all$n_cpg_tested > 0), sum(gt_all$n_cpg_tested > 0 & gt_all$biotype == "protein_coding"),
            sum(gt_all$n_cpg_eligible > 0), sum(gt_all$n_cpg_eligible > 0 & gt_all$biotype == "protein_coding")))
cat(sprintf("protein-coding genes: %d; with GO terms: %d\n", sum(gt_all$biotype == "protein_coding"), sum(gt_all$biotype == "protein_coding" & gt_all$n_go_terms > 0)))
gt <- gt_all[biotype == "protein_coding"]
any_dml <- gt$n_dml_ploidy + gt$n_dml_heat + gt$n_dml_shared > 0
cat(sprintf("protein-coding genes with DML but zero eligible CpGs: %d\n", sum(any_dml & gt$n_cpg_eligible == 0)))
print(gt[any_dml & n_cpg_eligible == 0, .(gene_id, gene, n_cpg_tested, n_cpg_eligible, n_dml_ploidy, n_dml_heat, n_go_terms)])
n_cpg <- pmax(gt$n_cpg_eligible, as.integer(any_dml))
uni <- cbind(gt[n_cpg > 0], bias = n_cpg[n_cpg > 0])
cat(sprintf("eligible universe: %d protein-coding genes; with GO: %d; DML genes ploidy %d heat %d shared %d (with GO: %d, %d, %d)\n", nrow(uni), sum(uni$n_go_terms > 0),
            sum(uni$n_dml_ploidy > 0), sum(uni$n_dml_heat > 0), sum(uni$n_dml_shared > 0), sum(uni$n_dml_ploidy > 0 & uni$n_go_terms > 0), sum(uni$n_dml_heat > 0 & uni$n_go_terms > 0), sum(uni$n_dml_shared > 0 & uni$n_go_terms > 0)))
cat("bias (eligible CpGs per gene) quantiles:\n"); print(quantile(uni$bias, c(0, .1, .25, .5, .75, .9, .95, .99, 1)))
tm <- list()
for (ont in c("BP", "CC", "MF")) {
  t0 <- Sys.time(); x <- term_members(ont, uni[n_go_terms > 0, gene_id], gene2go_all)
  tm[[ont]] <- list(genes = x$genes, M = incidence(x$members, x$genes))
  cat(sprintf("%s: %d genes, %d terms (%d with >= 5 genes); %.0f s\n", ont, length(x$genes), nrow(tm[[ont]]$M), sum(rowSums(tm[[ont]]$M) >= 5), as.numeric(Sys.time() - t0, units = "secs")))
}
saveRDS(list(uni = uni, tm = tm, gene2go_all = gene2go_all), file.path(sc, "eligible_setup_legacy.rds"))

# validate my replication against the pipeline's table
for (s in c("ploidy", "heat", "shared")) {
  ref <- fread(file.path(repo, sprintf("04_output/10_go_enrichment_gene_level/tables/go_gene_level_%s_%s.csv", s, mode)))[universe == "eligible"]
  de_all <- as.integer(uni[[paste0("n_dml_", s)]] > 0); names(de_all) <- uni$gene_id
  pwf <- fit_pwf(de_all, uni$bias)
  for (ont in c("BP", "CC", "MF")) {
    g <- tm[[ont]]$genes; r <- wallenius_p(tm[[ont]]$M, de_all[g], pwf[g]); r[, p_fisher := fisher_p(tm[[ont]]$M, de_all[g])]
    mm <- merge(r, ref[ontology == ont], by = "go_id")
    stopifnot(nrow(mm) == nrow(r), nrow(mm) == nrow(ref[ontology == ont]))
    cat(sprintf("%-6s %s: terms %d; max |log10 p_wall mine - pipeline| = %.2e; max |log10 p_fisher diff| = %.2e; counts equal: %s; BH(min 5): n_sig %d, min q %.4f (pipeline min q %.4f)\n", s, ont, nrow(mm),
                max(abs(log10(mm$p_wallenius.x) - log10(mm$p_wallenius.y))), max(abs(log10(mm$p_fisher.x) - log10(mm$p_fisher.y))), all(mm$x == mm$genes_with_dml & mm$n_cat == mm$genes_annotated),
                bh_count(r$p_wallenius, r$n_cat)["n_sig"], bh_count(r$p_wallenius, r$n_cat)["min_q"], min(mm$q_wallenius, na.rm = TRUE)))
  }
}

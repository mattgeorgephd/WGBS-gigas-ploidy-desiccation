# Robustness of "no term at BH q < 0.05" to go_min_genes, from the pipeline's own tables
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
res <- rbindlist(lapply(c("legacy", "corrected"), function(m) rbindlist(lapply(c("ploidy", "heat", "shared"), function(s)
  fread(file.path(repo, sprintf("04_output/10_go_enrichment_gene_level/tables/go_gene_level_%s_%s.csv", s, m)))))))
el <- res[universe == "eligible"]
mins <- c(1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 75, 100, 200)
out <- rbindlist(lapply(mins, function(k) el[genes_annotated >= k, {q <- p.adjust(p_wallenius, "BH"); qf <- p.adjust(p_fisher, "BH")
  .(min_genes = k, terms = .N, n_q05 = sum(q < 0.05), n_q10 = sum(q < 0.10), min_q = min(q), top = term[which.min(q)], top_n = genes_annotated[which.min(q)], n_fisher_q05 = sum(qf < 0.05))}, by = .(mode, set, ontology)]))
cat("Wallenius, eligible universe: BH recomputed at each minimum term size. Rows with any q < 0.10:\n")
print(out[n_q10 > 0][order(mode, set, ontology, min_genes)], nrows = 200)
cat("\nCount of q<0.05 terms by min_genes (legacy), set x ontology:\n")
print(dcast(out[mode == "legacy"], set + ontology ~ min_genes, value.var = "n_q05"))
cat("\nSmallest q by min_genes (legacy):\n")
out[, min_q_r := round(min_q, 3)]; print(dcast(out[mode == "legacy"], set + ontology ~ min_genes, value.var = "min_q_r"))
cat("\nSmallest q by min_genes (corrected):\n")
print(dcast(out[mode == "corrected"], set + ontology ~ min_genes, value.var = "min_q_r"))
cat("\nCount of q<0.05 terms by min_genes (corrected):\n")
print(dcast(out[mode == "corrected"], set + ontology ~ min_genes, value.var = "n_q05"))
# top terms heat CC legacy
cat("\nTop 12 terms, heat CC legacy eligible, Wallenius:\n")
print(el[mode == "legacy" & set == "heat" & ontology == "CC"][order(p_wallenius)][1:12, .(go_id, term, genes_annotated, genes_with_dml, expected = round(expected, 1), p_fisher = signif(p_fisher, 3), p_wallenius = signif(p_wallenius, 3), q_wallenius = signif(q_wallenius, 3))])
cat("\nTop 8 terms, heat MF and BP, ploidy CC/MF/BP legacy eligible:\n")
for (s in c("heat", "ploidy", "shared")) for (o in c("BP", "CC", "MF")) { cat(s, o, "\n"); print(el[mode == "legacy" & set == s & ontology == o & tested == TRUE][order(p_wallenius)][1:5, .(term = substr(term, 1, 50), genes_annotated, genes_with_dml, expected = round(expected, 1), p_fisher = signif(p_fisher, 3), p_wallenius = signif(p_wallenius, 3), q_wallenius = signif(q_wallenius, 3))]) }
# p-value distribution of tested terms: share below nominal levels (uniform expectation in brackets)
cat("\nShare of tested terms (>=5 genes) with Wallenius p below 0.05 / 0.01 / 0.001 (legacy, eligible):\n")
print(el[mode == "legacy" & tested == TRUE, .(terms = .N, p05 = round(mean(p_wallenius < 0.05), 4), p01 = round(mean(p_wallenius < 0.01), 4), p001 = round(mean(p_wallenius < 0.001), 5), min_p = signif(min(p_wallenius), 3),
                                              share_p_eq_1 = round(mean(p_wallenius > 0.999), 3), f05 = round(mean(p_fisher < 0.05), 4), f01 = round(mean(p_fisher < 0.01), 4)), by = .(set, ontology)])
# pooled across the three ontologies of a set (one BH family per set)
cat("\nOne BH family per set and mode (ontologies pooled), min 5 genes:\n")
print(el[tested == TRUE, {q <- p.adjust(p_wallenius, "BH"); .(terms = .N, n_q05 = sum(q < 0.05), min_q = round(min(q), 3))}, by = .(mode, set)])

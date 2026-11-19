# go_tests.R
# Gene-level GO enrichment: one definition of the test, used by step 10 and by
# the control in 06_diagnostics/go_random_cpg_control.R.

# Messages of topGO and goseq are not needed in the logs; goseq::nullp() also
# warns about its spline start values on every call, which is harmless.
quiet_go <- function(expr) {
  withCallingHandlers(suppressMessages(expr), warning = function(w) {
    if (grepl("initial point very close", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}

# A topGOdata object for one ontology over the GO-annotated genes of a universe.
#   gene_ids  genes of the universe that carry GO terms
#   gene2go   named list, gene ID -> its GO IDs (direct annotations)
# Every term is kept (nodeSize = 1); small terms are set aside when the tests are
# adjusted, see go_gene_level_test().
go_gene_level_object <- function(ontology, gene_ids, gene2go) {
  start <- factor(c(1L, rep(0L, length(gene_ids) - 1L)), levels = 0:1)   # placeholder gene list
  names(start) <- gene_ids
  quiet_go(methods::new("topGOdata", ontology = ontology, allGenes = start, annot = topGO::annFUN.gene2GO,
                        gene2GO = gene2go[gene_ids], nodeSize = 1))
}

# Test one list of genes against one universe and ontology.
#   td         object from go_gene_level_object()
#   universe   data.table with every gene of the universe, annotated or not:
#              gene_id and bias (the number of CpGs the gene offers)
#   dml_genes  IDs of the genes that hold a DML
# For every term: genes annotated, genes with a DML, the classic Fisher p-value
# (topGO), and the Wallenius p-value of goseq, whose weighting function is
# fitted on the whole universe with `bias` as the bias variable. Terms with at
# least config$go_min_genes annotated genes are "tested" and get Benjamini-
# Hochberg q-values; the others keep their p-values and get no q-value.
# Returns the table, the updated object, the fitted weighting function, and the
# genes of each term.
go_gene_level_test <- function(td, universe, dml_genes, config) {
  annotated <- topGO::allGenes(td)
  gene_list <- factor(as.integer(annotated %in% dml_genes), levels = 0:1); names(gene_list) <- annotated
  td <- topGO::updateGenes(td, gene_list)
  fisher <- methods::new("classicCount", testStatistic = topGO::GOFisherTest, name = "Fisher test")
  p_f <- topGO::score(quiet_go(topGO::getSigGroups(td, fisher)))
  stats <- topGO::termStat(td, names(p_f))
  members <- topGO::genesInTerm(td, names(p_f))

  de <- as.integer(universe$gene_id %in% dml_genes); names(de) <- universe$gene_id
  pwf <- quiet_go(goseq::nullp(de, bias.data = universe$bias, plot.fit = FALSE))
  gene2cat <- split(rep(names(members), lengths(members)), unlist(members, use.names = FALSE))
  gs <- data.table::as.data.table(quiet_go(goseq::goseq(pwf, gene2cat = gene2cat, method = "Wallenius",
                                                        use_genes_without_cat = FALSE)))
  data.table::setkey(gs, category)
  gs <- gs[names(p_f)]
  stopifnot(all(gs$numInCat == stats$Annotated), all(gs$numDEInCat == stats$Significant))

  res <- data.table::data.table(
    go_id = names(p_f), term = unname(AnnotationDbi::Term(GO.db::GOTERM[names(p_f)])),
    genes_annotated = stats$Annotated, genes_with_dml = stats$Significant, expected = stats$Expected,
    p_fisher = unname(p_f), p_wallenius = gs$over_represented_pvalue)
  res[, tested := genes_annotated >= config$go_min_genes]
  res[tested == TRUE, `:=`(q_fisher = stats::p.adjust(p_fisher, "BH"), q_wallenius = stats::p.adjust(p_wallenius, "BH"))]
  res[, enrichment := genes_with_dml / expected]
  data.table::setorder(res, p_wallenius)
  list(table = res, td = td, pwf = pwf, members = members)
}

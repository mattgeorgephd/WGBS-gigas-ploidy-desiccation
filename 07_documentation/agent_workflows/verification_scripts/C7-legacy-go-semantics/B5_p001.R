# What the lists behind Figure 9 and Table 1 would hold under the rule the Methods state (p < 0.01) instead of the ten top-ranked terms.
.libPaths(c("/home/user/gannet_cache/legacy_r_lib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(GO.db) })
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
P <- readRDS(file.path(scratch, "A_port.rds")); out <- P$out; L <- lapply(P$DMLGOterms, as.data.table)
first_line <- P$geneID2GO[!duplicated(names(P$geneID2GO))]
anc <- list(BP = as.list(GOBPANCESTOR), CC = as.list(GOCCANCESTOR), MF = as.list(GOMFANCESTOR))
onts <- c(BP = "P", CC = "C", MF = "F")
res <- list()
for (o in names(onts)) {
  # closure of every transcript's terms within this ontology
  a <- anc[[o]]
  clos <- lapply(first_line, function(d) { d <- d[d %in% names(a)]; unique(c(d, unlist(a[d], use.names = FALSE))) })
  for (s in names(L)) {
    k <- paste(s, o); r <- out[[k]]; l <- L[[s]]
    feas <- unique(r$feasible)
    sig_terms <- names(r$pval)[r$pval < 0.01]
    top10 <- r$top10
    in_set <- function(terms) feas[vapply(clos[feas], function(v) any(v %chin% terms), logical(1))]
    t10 <- in_set(top10); stopifnot(setequal(t10, unique(unlist(r$inTerm))))
    tp <- in_set(sig_terms)
    m10 <- unique(l[transcript %chin% t10 & GOcat == onts[[o]], .(geneID, chr, start, end, meth.diff)])
    mp <- unique(l[transcript %chin% tp & GOcat == onts[[o]], .(geneID, chr, start, end, meth.diff)])
    allm <- unique(l[GOcat == onts[[o]], .(geneID, chr, start, end, meth.diff)])
    res[[k]] <- data.table(set = s, ont = o, terms_p_below_0.01 = length(sig_terms),
                           top10_rows = nrow(m10), top10_genes = uniqueN(m10$geneID),
                           p001_rows = nrow(mp), p001_hyper = sum(mp$meth.diff > 0), p001_hypo = sum(mp$meth.diff < 0), p001_genes = uniqueN(mp$geneID),
                           all_rows_of_this_ontology = nrow(allm), all_genes = uniqueN(allm$geneID))
  }
}
res <- rbindlist(res); setorder(res, set, ont); options(width = 220); print(res)
fwrite(res, file.path(scratch, "B5_p001.csv"))

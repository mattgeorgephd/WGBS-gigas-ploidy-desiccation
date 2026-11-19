.libPaths(c("/home/user/gannet_cache/legacy_r_lib", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(GO.db); library(readxl) })
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
P <- readRDS(file.path(scratch, "A_port.rds"))
out <- P$out; L <- lapply(P$DMLGOterms, as.data.table); E <- lapply(P$DMLGOterms_exact, as.data.table)
geneID2GO <- P$geneID2GO; first_line <- geneID2GO[!duplicated(names(geneID2GO))]
options(width = 250)

cat("\n==== (c) HEAT BP: the term at rank 10 and the term at rank 11 have the same p-value. What if rank 11 had been taken?\n")
r <- out[["HEAT BP"]]; l <- L$HEAT; lo <- l[GOcat == "P"]
suppressPackageStartupMessages(library(topGO))
# rebuild membership of rank 11 from the map: transcripts whose own terms or ancestors include it
anc <- as.list(GOBPANCESTOR)
alt <- r$allRes$GO.ID[11]; cur <- r$allRes$GO.ID[10]
members <- function(go) names(first_line)[vapply(first_line, function(d) go %in% d || any(vapply(d, function(x) go %in% anc[[x]], logical(1))), logical(1))]
m_cur <- intersect(members(cur), r$feasible); m_alt <- intersect(members(alt), r$feasible)
cat("rank 10", cur, Term(GOTERM[[cur]]), "| members rebuilt", length(m_cur), "vs topGO", length(r$inTerm[[cur]]), "identical:", setequal(m_cur, r$inTerm[[cur]]), "\n")
cat("rank 11", alt, Term(GOTERM[[alt]]), "| members", length(m_alt), "| same transcripts as rank 10:", setequal(m_cur, m_alt), "\n")
cat("rank 11 is ancestor of rank 10:", alt %in% anc[[cur]], "\n")
inT <- r$inTerm; inT[[cur]] <- NULL; inT[[alt]] <- m_alt
pairs <- data.table(GOID = rep(names(inT), lengths(inT)), transcript = unlist(inT, use.names = FALSE))
an <- unique(merge(pairs, lo, by = c("GOID", "transcript"), all.x = TRUE)[, !"GOcat"])
md <- unique(l[transcript %chin% unique(pairs$transcript) & GOcat == "P", .(geneID, chr, start, end, meth.diff)])
cat("with rank 11 instead of rank 10: unique transcripts", uniqueN(pairs$transcript), "(93); unique genes", uniqueN(na.omit(an$geneID)), "(3); gene products", nrow(an),
    "(1441); DML", nrow(md), "(7)\n")

cat("\n==== (d) ploidy BP: the gene whose DML the faulty reader removed from the list\n")
r <- out[["PLOIDY BP"]]; ut <- unique(unlist(r$inTerm, use.names = FALSE))
mx <- unique(E$PLOIDY[transcript %chin% ut & GOcat == "P", .(geneID, chr, start, end, meth.diff)])
ml <- unique(L$PLOIDY[transcript %chin% ut & GOcat == "P", .(geneID, chr, start, end, meth.diff)])
print(mx[!paste(geneID, start) %chin% paste(ml$geneID, ml$start)])
for (k in c("PLOIDY CC", "HEAT CC", "HEAT MF")) {
  s <- sub(" .*", "", k); o <- c(BP = "P", CC = "C", MF = "F")[[sub(".* ", "", k)]]
  ut <- unique(unlist(out[[k]]$inTerm, use.names = FALSE))
  mx <- unique(E[[s]][transcript %chin% ut & GOcat == o, .(geneID, chr, start, end, meth.diff)])
  ml <- unique(L[[s]][transcript %chin% ut & GOcat == o, .(geneID, chr, start, end, meth.diff)])
  cat(k, ": DML rows with every line read", nrow(mx), "(hyper", sum(mx$meth.diff > 0), "hypo", sum(mx$meth.diff < 0), ") against legacy", nrow(ml), "\n")
}

cat("\n==== (g) summary_table.xlsx\n")
f <- file.path(repo, "08_archive/GOterms/summary_table.xlsx")
for (sh in excel_sheets(f)) { cat("--- sheet", sh, "\n"); print(as.data.frame(suppressMessages(read_excel(f, sheet = sh, col_names = FALSE))), row.names = FALSE) }

cat("\n==== annotation_summary.csv\n")
print(fread(file.path(repo, "04_output/08_gene_annotation/tables/annotation_summary.csv")))

cat("\n==== (m) shared loci in both modes\n")
for (m in c("legacy", "corrected")) { d <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/dml_shared_%s.csv", m))); cat(m, nrow(d), "\n"); print(table(d$class)); print(table(heat = d$direction_heat)) }

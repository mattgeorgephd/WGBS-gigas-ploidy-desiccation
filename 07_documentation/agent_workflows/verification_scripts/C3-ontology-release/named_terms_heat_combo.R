suppressPackageStartupMessages(library(data.table))
.libPaths(c("/home/user/gannet_cache/godb/lib313", .libPaths()))
scr <- getwd(); source(file.path(scr, "parse_archive.R"))
R <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive/GOterms"
con <- DBI::dbConnect(RSQLite::SQLite(), "/home/user/gannet_cache/godb/lib313/GO.db/extdata/GO.sqlite")
terms <- as.data.table(DBI::dbGetQuery(con, "select go_id, term, ontology from go_term where ontology = 'BP'")); DBI::dbDisconnect(con)
named <- list(
  HEAT = c("MAPK cascade", "ER overload response", "signal transduction", "mRNA methylation",
           "mRNA export from nucleus in response to heat stress", "rRNA catabolic process", "response to osmotic stress", "response to oxidative stress"),
  COMBO = c("tRNA modification", "RNA 3'-end processing", "post-translational protein modification", "calcium-mediated signaling",
            "ER overload response", "response to osmotic stress", "positive regulation of gene expression", "gene silencing by RNA"))
files <- list(HEAT = c("BOTH-HEAT/HEAT-BP-FisherTestResults.csv", "BOTH-HEAT/HEAT-BP-GOterms.csv"),
              COMBO = c("BOTH-COMBO/COMBO-BP-FisherTestResults.csv", "BOTH-COMBO/COMBO-BP-GOterms.csv"))
for (s in names(named)) {
  fish <- as.data.table(parse_archive(file.path(R, files[[s]][1]))); dml <- fread(file.path(R, files[[s]][2]))
  x <- terms[term %in% named[[s]]]
  x <- merge(x, fish[, .(go_id, annotated, significant, classic, rank)], by = "go_id", all.x = TRUE)
  x[, attached_to_DML_transcript := go_id %in% dml$GOID]
  cat("\n", s, ": terms named in Results 3.4; Fisher result in the archived BP table (", nrow(fish), " terms)\n", sep = "")
  print(x[order(rank), .(go_id, term, annotated, significant, p = classic, rank, attached_to_DML_transcript)])
  cat("named terms not found by name in GO.db 3.13.0:", paste(setdiff(named[[s]], x$term), collapse = "; "), "\n")
}

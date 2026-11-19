S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
x <- readRDS(file.path(S, "tables.rds")); E <- x$E; L <- x$L
suppressPackageStartupMessages(library(GO.db))
cat("GO.db", as.character(packageVersion("GO.db")), "\n")
ul <- readLines(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))
nm <- gsub(" ", "", sub("\t.*$", "", ul)); go <- sub("^[^\t]*\t", "", ul); first <- !duplicated(nm)
g2g <- setNames(lapply(strsplit(go[first], ",", fixed = TRUE), function(v) gsub(" ", "", v)), nm[first])
allgo <- unique(unlist(g2g))
ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = allgo, columns = "ONTOLOGY", keytype = "GOID")); ont <- setNames(ont$ONTOLOGY, ont$GOID)
cat("GO IDs unknown to this GO.db:", sum(is.na(ont)), "of", length(allgo), "\n")
for (o in c("BP", "CC", "MF")) {
  feas <- names(g2g)[vapply(g2g, function(v) any(ont[v] == o, na.rm = TRUE), logical(1))]
  for (s in names(L)) cat(sprintf("%s %-6s feasible %d; transcripts of interest: legacy reader %d, exact reader %d\n", o, s, length(feas),
                                   sum(feas %in% L[[s]]$transcript), sum(feas %in% E[[s]]$transcript)))
}

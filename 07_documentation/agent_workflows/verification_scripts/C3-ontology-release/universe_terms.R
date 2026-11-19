# For one GO.db release: the GO terms and "Annotated" counts that topGO derives
# from a universe file. These depend only on the universe and the ontology, not
# on the gene list, so they date any archived topGO table whose universe is known.
# Usage: Rscript universe_terms.R <tag> <godb_lib|default> <universe_name=path> [more universes] [onts=CC,MF,BP]
args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]; godb_lib <- args[2]
onts <- c("CC", "MF", "BP")
if (any(grepl("^onts=", args))) { onts <- strsplit(sub("^onts=", "", args[grepl("^onts=", args)]), ",")[[1]]; args <- args[!grepl("^onts=", args)] }
uni <- args[-(1:2)]; uni <- setNames(sub("^[^=]+=", "", uni), sub("=.*", "", uni))
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C3-ontology-release"
out <- file.path(scr, "universe_terms", tag); dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (godb_lib != "default") .libPaths(c(godb_lib, .libPaths()))
suppressPackageStartupMessages(library(topGO))
info <- GO.db::GO_dbInfo()
cat("GO.db", as.character(packageVersion("GO.db")), info$value[info$name == "GOSOURCEDATE"], "topGO", as.character(packageVersion("topGO")), "\n")
for (u in names(uni)) {
  geneID2GO <- topGO::readMappings(file = uni[[u]]); geneNames <- names(geneID2GO)
  gl <- factor(c(1L, rep(0L, length(geneNames) - 1L)), levels = 0:1); names(gl) <- geneNames
  for (ont in onts) {
    t0 <- Sys.time()
    td <- suppressMessages(new("topGOdata", ontology = ont, allGenes = gl, annot = annFUN.gene2GO, gene2GO = geneID2GO))
    st <- termStat(td)
    res <- data.frame(go_id = rownames(st), annotated = st$Annotated, stringsAsFactors = FALSE)
    write.csv(res, file.path(out, sprintf("%s_%s.csv", u, ont)), row.names = FALSE, quote = FALSE)
    cat(sprintf("%s %s: %d terms, %d feasible entries (%.0f s)\n", u, ont, nrow(res), numGenes(td), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    rm(td); invisible(gc())
  }
}
cat(grep("VmHWM", readLines("/proc/self/status"), value = TRUE), "\n")

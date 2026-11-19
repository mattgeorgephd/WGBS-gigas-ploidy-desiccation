# Effect of repeated transcript names in the legacy universe on the Fisher margins.
.libPaths(c("/home/user/gannet_cache/godb/lib313", .libPaths()))
suppressPackageStartupMessages(library(topGO))
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C3-ontology-release"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
cols <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
onts <- strsplit(commandArgs(trailingOnly = TRUE)[1], ",")[[1]]
tmp <- tempfile(fileext = ".tab"); writeLines(readLines(gzfile(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))), tmp)
geneID2GO <- topGO::readMappings(file = tmp); geneNames <- names(geneID2GO)
cat("map lines", length(geneNames), "unique names", length(unique(geneNames)), "\n")
# are repeated lines of a transcript identical in content?
sp <- split(vapply(geneID2GO, function(x) paste(sort(x), collapse = ","), ""), geneNames)
cat("transcripts whose repeated lines differ in GO content:", sum(vapply(sp, function(x) length(unique(x)) > 1, logical(1))), "\n")
uniq <- geneID2GO[!duplicated(geneNames)]
fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
out <- list()
for (ont in onts) for (s in c("PLOIDY", "HEAT", "COMBO")) {
  ann <- suppressWarnings(read.table(file.path(repo, "03_input_files/annotation/legacy", sprintf("DML-Cov10-20-BOTH-%s.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz", s)),
                                     header = FALSE, sep = "\t", col.names = cols))
  gl <- factor(as.integer(geneNames %in% ann$transcript)); names(gl) <- geneNames
  td <- suppressMessages(new("topGOdata", ontology = ont, allGenes = gl, annot = annFUN.gene2GO, gene2GO = geneID2GO))
  p1 <- score(suppressMessages(getSigGroups(td, fisher)))
  un <- names(uniq); gu <- factor(as.integer(un %in% ann$transcript)); names(gu) <- un
  tu <- suppressMessages(new("topGOdata", ontology = ont, allGenes = gu, annot = annFUN.gene2GO, gene2GO = uniq))
  p2 <- score(suppressMessages(getSigGroups(tu, fisher)))
  stopifnot(setequal(names(p1), names(p2))); p2 <- p2[names(p1)]
  s1 <- termStat(td, names(p1)); s2 <- termStat(tu, names(p1))
  o1 <- names(sort(p1))[1:10]; o2 <- names(sort(p2))[1:10]
  out[[length(out) + 1]] <- data.frame(set = s, ont = ont, N_legacy = numGenes(td), N_unique = numGenes(tu), sig_legacy = numSigGenes(td), sig_unique = numSigGenes(tu),
    rate_legacy = round(numSigGenes(td) / numGenes(td), 4), rate_unique = round(numSigGenes(tu) / numGenes(tu), 4),
    annotated_same = all(s1$Annotated == s2$Annotated), significant_same = all(s1$Significant == s2$Significant),
    n_p01_legacy = sum(p1 < 0.01), n_p01_unique = sum(p2 < 0.01), n_p01_both = sum(p1 < 0.01 & p2 < 0.01),
    top10_in_common = length(intersect(o1, o2)), median_log10_ratio_p = round(median(log10(p2[p1 < 0.01] / p1[p1 < 0.01])), 2))
  print(out[[length(out)]])
  rm(td, tu); invisible(gc())
}
res <- do.call(rbind, out); write.csv(res, file.path(scr, sprintf("duplicate_margins_%s.csv", paste(onts, collapse = "_"))), row.names = FALSE)

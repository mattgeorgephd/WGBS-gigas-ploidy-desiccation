# 09_topgo.R (scratch; own code, nothing sourced from the repository)
# The legacy topGO procedure with (a) the legacy reader and (b) a complete read, on the
# legacy universe and GO.db 3.13.0. Questions:
#   1. Does (a) reproduce the nine archived Fisher tables?
#   2. Would the top-ten lists and the Figure 9 counts differ under (b)?  For which DML sets?
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
.libPaths(c("/home/user/gannet_cache/legacy_r_lib", .libPaths()))
suppressPackageStartupMessages(library(topGO))
stopifnot(as.character(packageVersion("GO.db")) == "3.13.0")
cat("topGO", as.character(packageVersion("topGO")), " GO.db", as.character(packageVersion("GO.db")), "\n")
options(width = 220)
x <- readRDS(file.path(S, "tables.rds")); E <- x$E; L <- x$L
for (s in names(E)) { E[[s]]$meth.diff <- as.numeric(E[[s]]$meth.diff); E[[s]]$geneID <- as.integer(E[[s]]$geneID)
                      E[[s]]$start <- as.integer(E[[s]]$start); E[[s]]$end <- as.integer(E[[s]]$end) }
arch <- file.path(repo, "08_archive/GOterms")
adir <- c(ploidy = "BOTH-PLOIDY", heat = "BOTH-HEAT", shared = "BOTH-COMBO"); apre <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")
cats <- c(BP = "P", CC = "C", MF = "F")

geneID2GO <- suppressWarnings(topGO::readMappings(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz")))
geneNames <- names(geneID2GO)
cat("map lines:", length(geneNames), " distinct names:", length(unique(geneNames)), "\n")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

read_archived <- function(f) {                      # written with quote = FALSE: term names may hold commas
  l <- readLines(f)[-1]; p <- strsplit(l, ",", fixed = TRUE)
  data.frame(GO.ID = vapply(p, `[`, "", 1), Annotated = as.integer(vapply(p, function(v) v[length(v) - 3], "")),
             Significant = as.integer(vapply(p, function(v) v[length(v) - 2], "")), Expected = vapply(p, function(v) v[length(v) - 1], ""),
             classic = vapply(p, function(v) v[length(v)], ""), stringsAsFactors = FALSE)
}
fig9 <- function(td, top_ids, tab, cat_letter) {    # legacy lines 354-407: DML "in genes with enriched GO terms"
  tx <- unique(unlist(genesInTerm(td, whichGO = top_ids)))
  d <- tab[tab$transcript %in% tx & !is.na(tab$GOcat) & tab$GOcat == cat_letter, c("geneID", "chr", "start", "end", "meth.diff")]
  d <- unique(d)
  c(total = nrow(d), hyper = sum(d$meth.diff >= 0), hypo = sum(d$meth.diff <= 0), genes = length(unique(d$geneID)))
}

res <- list(); tops <- list()
for (o in names(cats)) {
  td <- NULL
  for (s in names(L)) for (reader in c("legacy", "exact")) {
    tab <- if (reader == "legacy") L[[s]] else E[[s]]
    gl <- factor(as.integer(geneNames %in% tab$transcript)); names(gl) <- geneNames
    td <- if (is.null(td)) suppressMessages(new("topGOdata", ontology = o, allGenes = gl, annot = annFUN.gene2GO, gene2GO = geneID2GO)) else updateGenes(td, gl)
    rf <- suppressMessages(getSigGroups(td, test.stat)); p <- score(rf)
    gt <- GenTable(td, classic = rf, ranksOf = "classic", orderBy = "classic", topNodes = length(p))
    key <- paste(s, o, reader)
    tops[[key]] <- gt$GO.ID[1:10]
    f9 <- fig9(td, gt$GO.ID[1:10], tab, cats[[o]])
    row <- data.frame(set = s, ontology = o, reader = reader, terms = nrow(gt), sig_transcripts = length(sigGenes(td)),
                      p_lt_0.01 = sum(p < 0.01), fig9_total = f9[["total"]], fig9_hyper = f9[["hyper"]], fig9_hypo = f9[["hypo"]], fig9_genes = f9[["genes"]],
                      archived_same_ids = NA, archived_same_counts = NA, archived_same_p_strings = NA, archived_same_order = NA, archived_top10_same = NA)
    if (reader == "legacy") {
      old <- read_archived(file.path(arch, adir[[s]], sprintf("%s-%s-FisherTestResults.csv", apre[[s]], o)))
      row$archived_same_ids <- setequal(old$GO.ID, gt$GO.ID) && nrow(old) == nrow(gt)
      m <- merge(old, gt[, c("GO.ID", "Annotated", "Significant", "Expected", "classic")], by = "GO.ID", suffixes = c(".old", ".new"))
      row$archived_same_counts <- nrow(m) == nrow(old) && all(m$Annotated.old == m$Annotated.new) && all(m$Significant.old == m$Significant.new)
      row$archived_same_p_strings <- nrow(m) == nrow(old) && all(trimws(m$classic.old) == trimws(as.character(m$classic.new)))
      row$archived_same_order <- identical(old$GO.ID, gt$GO.ID)
      row$archived_top10_same <- identical(old$GO.ID[1:10], gt$GO.ID[1:10])
      saveRDS(gt, file.path(S, sprintf("gentable_%s_%s_legacy.rds", s, o)))
    } else saveRDS(gt, file.path(S, sprintf("gentable_%s_%s_exact.rds", s, o)))
    res[[length(res) + 1]] <- row
    cat(sprintf("%s done\n", key))
  }
}
res <- do.call(rbind, res)
print(res, row.names = FALSE)
write.csv(res, file.path(S, "topgo_reader_comparison.csv"), row.names = FALSE)

cat("\nTop-ten lists, legacy reader against complete read\n")
cmp <- list()
for (s in names(L)) for (o in names(cats)) {
  a <- tops[[paste(s, o, "legacy")]]; b <- tops[[paste(s, o, "exact")]]
  gl <- readRDS(file.path(S, sprintf("gentable_%s_%s_legacy.rds", s, o))); ge <- readRDS(file.path(S, sprintf("gentable_%s_%s_exact.rds", s, o)))
  rank_in_exact <- match(a, ge$GO.ID)
  cmp[[length(cmp) + 1]] <- data.frame(set = s, ontology = o, top10_shared = length(intersect(a, b)), top10_identical_order = identical(a, b),
                                       ranks_of_legacy_top10_under_complete_read = paste(rank_in_exact, collapse = " "))
}
cmp <- do.call(rbind, cmp); print(cmp, row.names = FALSE)
write.csv(cmp, file.path(S, "topgo_top10_comparison.csv"), row.names = FALSE)
saveRDS(tops, file.path(S, "topgo_top10_ids.rds"))

# Independent re-run of the legacy topGO enrichment, following
# 08_archive/code/5_WGBS_Functional-enrichement-analysis-BOTH.R line by line
# (lines 72-84, 186-206, 210-246, 255-317), from the ORIGINAL legacy inputs.
# Nothing is written to the repository.
#
# Usage: Rscript rerun_legacy_topgo.R <tag> <godb_lib|default> <topgo_lib|default> [onts=BP,CC,MF] [sets=PLOIDY,HEAT,COMBO] [reader=legacy|exact]
args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]; godb_lib <- args[2]; topgo_lib <- args[3]
onts <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else c("BP", "CC", "MF")
sets <- if (length(args) >= 5) strsplit(args[5], ",")[[1]] else c("PLOIDY", "HEAT", "COMBO")
reader <- if (length(args) >= 6) args[6] else "legacy"
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C3-ontology-release"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
out  <- file.path(scr, "runs", tag); dir.create(out, recursive = TRUE, showWarnings = FALSE)

libs <- character(0)
if (topgo_lib != "default") libs <- c(libs, topgo_lib)
if (godb_lib  != "default") libs <- c(libs, godb_lib)
if (length(libs)) .libPaths(c(libs, .libPaths()))
stopifnot(!"GO.db" %in% loadedNamespaces(), !"topGO" %in% loadedNamespaces())
suppressPackageStartupMessages(library(topGO))
info <- GO.db::GO_dbInfo()
ver <- c(R = R.version.string, topGO = as.character(packageVersion("topGO")), topGO_path = find.package("topGO"),
         GO.db = as.character(packageVersion("GO.db")), GO.db_path = find.package("GO.db"),
         GO_source_date = info$value[info$name == "GOSOURCEDATE"], reader = reader)
print(ver)
writeLines(paste(names(ver), ver, sep = "\t"), file.path(out, "versions.tsv"))

# ---- legacy lines 72-84: the DML annotation tables, read.table defaults ----
cols <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
ann_file <- function(s) file.path(repo, "03_input_files/annotation/legacy",
                                  sprintf("DML-Cov10-20-BOTH-%s.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz", s))
read_ann <- function(s) {
  if (reader == "legacy") suppressWarnings(read.table(ann_file(s), header = FALSE, sep = "\t", col.names = cols))
  else read.table(ann_file(s), header = FALSE, sep = "\t", col.names = cols, quote = "", comment.char = "")
}
ann <- lapply(setNames(sets, sets), read_ann)
for (s in sets) cat(sprintf("%s: rows read %d, transcripts %d\n", s, nrow(ann[[s]]), length(unique(ann[[s]]$transcript))))

# ---- legacy lines 186-189: the universe ----
tmp <- tempfile(fileext = ".tab")
writeLines(readLines(gzfile(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))), tmp)
geneID2GO <- topGO::readMappings(file = tmp)
geneNames <- names(geneID2GO)
cat(sprintf("universe: %d map lines, %d unique names\n", length(geneNames), length(unique(geneNames))))

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
summ <- list()
for (ont in onts) for (s in sets) {
  t0 <- Sys.time()
  Genes <- ann[[s]]$transcript
  GeneList <- factor(as.integer(geneNames %in% Genes)); names(GeneList) <- geneNames       # legacy 193-206
  GOdata <- suppressMessages(new("topGOdata", ontology = ont, allGenes = GeneList,
                                 annot = annFUN.gene2GO, gene2GO = geneID2GO))                # legacy 210-246
  resultFisher <- suppressMessages(getSigGroups(GOdata, test.stat))                          # legacy 257-303
  pvalFis <- score(resultFisher)
  allRes <- GenTable(GOdata, classic = resultFisher, ranksOf = "classic", orderBy = "classic", topNodes = length(pvalFis))
  f <- file.path(out, sprintf("%s-%s-FisherTestResults.csv", s, ont))
  write.csv(allRes, f, quote = FALSE, row.names = FALSE)                                     # legacy 307-317
  saveRDS(pvalFis, file.path(out, sprintf("%s-%s-pvalues.rds", s, ont)))
  summ[[length(summ) + 1]] <- data.frame(set = s, ont = ont, terms = length(pvalFis), genes_in_object = numGenes(GOdata),
                                         sig_genes = numSigGenes(GOdata), feasible = length(feasible(GOdata)),
                                         secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  print(summ[[length(summ)]])
  rm(GOdata, resultFisher, allRes); invisible(gc())
}
write.csv(do.call(rbind, summ), file.path(out, "run_summary.csv"), row.names = FALSE)
writeLines(capture.output(print(sessionInfo())), file.path(out, "session_info.txt"))
cat("peak memory, gc max used (MB):", sum(gc()[, 6]), "\n")
cat(grep("VmHWM|VmPeak", readLines("/proc/self/status"), value = TRUE), sep = "\n")

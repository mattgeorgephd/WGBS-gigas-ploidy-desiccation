# 06_archive_check.R (scratch)
# Do the outputs that the legacy script wrote in 2021 agree with what read.table() defaults
# return today (R 4.3.3) for the same files?  Two independent traces:
#   (a) the per-ontology GO ID count tables  GOterms/BOTH-*/*-GOterms.csv
#   (b) the root-term rows of the nine archived Fisher tables (Annotated, Significant)
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
x <- readRDS(file.path(S, "tables.rds")); E <- x$E; L <- x$L; universe <- x$universe
arch <- file.path(repo, "08_archive/GOterms")
options(width = 200)

# ---------------- (a) GO ID count tables ----------------
files <- list(
  ploidy = c(P = "BOTH-PLOIDY/ploidy-BP-GOterms.csv", C = "BOTH-PLOIDY/PLOIDY-CC-GOterms.csv", F = "BOTH-PLOIDY/PLOIDY-MF-GOterms.csv"),
  heat   = c(P = "BOTH-HEAT/HEAT-BP-GOterms.csv",     C = "BOTH-HEAT/HEAT-CC-GOterms.csv",     F = "BOTH-HEAT/HEAT-MF-GOterms.csv"),
  shared = c(P = "BOTH-COMBO/COMBO-BP-GOterms.csv",   C = "BOTH-COMBO/COMBO-CC-GOterms.csv",   F = "BOTH-COMBO/COMBO-MF-GOterms.csv"))
cmp <- function(tab, old) {
  new <- as.data.frame(table(tab), stringsAsFactors = FALSE); names(new) <- c("GOID", "Freq")
  m <- merge(old, new, by = "GOID", all = TRUE, suffixes = c(".archived", ".now"))
  c(ids_archived = nrow(old), ids_now = nrow(new), ids_in_both = sum(!is.na(m$Freq.archived) & !is.na(m$Freq.now)),
    counts_equal = sum(m$Freq.archived == m$Freq.now, na.rm = TRUE), sum_archived = sum(old$Freq), sum_now = sum(new$Freq))
}
out <- list()
for (s in names(files)) for (k in names(files[[s]])) {
  old <- read.csv(file.path(arch, files[[s]][[k]]), stringsAsFactors = FALSE)
  l <- L[[s]]; e <- E[[s]]
  a <- cmp(l$GOID[which(l$GOcat == k)], old); b <- cmp(e$GOID[e$GOcat == k], old)
  out[[length(out) + 1]] <- data.frame(set = s, GOcat = k, archived_file = files[[s]][[k]],
    archived_ids = a[["ids_archived"]], archived_rows = a[["sum_archived"]],
    legacy_reader_ids = a[["ids_now"]], legacy_reader_rows = a[["sum_now"]],
    legacy_identical = a[["ids_archived"]] == a[["ids_now"]] && a[["ids_in_both"]] == a[["ids_now"]] && a[["counts_equal"]] == a[["ids_now"]],
    exact_reader_ids = b[["ids_now"]], exact_reader_rows = b[["sum_now"]],
    exact_identical = b[["ids_archived"]] == b[["ids_now"]] && b[["ids_in_both"]] == b[["ids_now"]] && b[["counts_equal"]] == b[["ids_now"]])
}
out <- do.call(rbind, out)
cat("(a) archived GO ID count tables against today's reads\n"); print(out, row.names = FALSE)
write.csv(out, file.path(S, "archive_check_goid_counts.csv"), row.names = FALSE)

# ---------------- (b) root rows of the archived Fisher tables ----------------
suppressPackageStartupMessages(library(GO.db, lib.loc = "/home/user/gannet_cache/legacy_r_lib"))
cat("\nGO.db version loaded:", as.character(packageVersion("GO.db", lib.loc = "/home/user/gannet_cache/legacy_r_lib")), "\n")
# universe map: transcript -> GO IDs (first line per name, as gene2GO[intersect(names, allGenes)] does)
ul <- readLines(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))
nm <- gsub(" ", "", sub("\t.*$", "", ul)); go <- sub("^[^\t]*\t", "", ul)
first <- !duplicated(nm)
# are the repeated lines of a transcript identical? (so that taking the first is harmless)
cat("transcripts whose repeated map lines differ:", sum(tapply(go, nm, function(v) length(unique(v))) > 1), "\n")
g2g <- setNames(lapply(strsplit(go[first], ",", fixed = TRUE), function(v) gsub(" ", "", v)), nm[first])
allgo <- unique(unlist(g2g))
ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = allgo, columns = "ONTOLOGY", keytype = "GOID"))
ont <- setNames(ont$ONTOLOGY, ont$GOID)
cat(sprintf("GO IDs in the universe map: %d; known to GO.db 3.13.0: %d (BP %d, CC %d, MF %d); unknown: %d\n", length(allgo), sum(!is.na(ont)),
            sum(ont == "BP", na.rm = TRUE), sum(ont == "CC", na.rm = TRUE), sum(ont == "MF", na.rm = TRUE), sum(is.na(ont))))
feasible <- sapply(c("BP", "CC", "MF"), function(o) names(g2g)[vapply(g2g, function(v) any(ont[v] == o, na.rm = TRUE), logical(1))], simplify = FALSE)
roots <- c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
ftab <- list(
  ploidy = c(BP = "BOTH-PLOIDY/PLOIDY-BP-FisherTestResults.csv", CC = "BOTH-PLOIDY/PLOIDY-CC-FisherTestResults.csv", MF = "BOTH-PLOIDY/PLOIDY-MF-FisherTestResults.csv"),
  heat   = c(BP = "BOTH-HEAT/HEAT-BP-FisherTestResults.csv",     CC = "BOTH-HEAT/HEAT-CC-FisherTestResults.csv",     MF = "BOTH-HEAT/HEAT-MF-FisherTestResults.csv"),
  shared = c(BP = "BOTH-COMBO/COMBO-BP-FisherTestResults.csv",   CC = "BOTH-COMBO/COMBO-CC-FisherTestResults.csv",   MF = "BOTH-COMBO/COMBO-MF-FisherTestResults.csv"))
rb <- list()
for (s in names(ftab)) for (o in names(roots)) {
  ln <- grep(paste0("^", roots[[o]], ","), readLines(file.path(arch, ftab[[s]][[o]])), value = TRUE)
  p <- strsplit(ln, ",", fixed = TRUE)[[1]]
  n <- length(p)                                   # term names may hold commas; counts are the 4th and 3rd fields from the end
  rb[[length(rb) + 1]] <- data.frame(set = s, ontology = o, archived_annotated = as.integer(p[n - 3]), archived_significant = as.integer(p[n - 2]),
    now_annotated = length(feasible[[o]]),
    now_significant_legacy_reader = sum(feasible[[o]] %in% L[[s]]$transcript),
    now_significant_exact_reader = sum(feasible[[o]] %in% E[[s]]$transcript))
}
rb <- do.call(rbind, rb)
rb$legacy_matches <- rb$archived_annotated == rb$now_annotated & rb$archived_significant == rb$now_significant_legacy_reader
rb$exact_matches <- rb$archived_significant == rb$now_significant_exact_reader
cat("\n(b) root-term rows of the archived Fisher tables against counts recomputed from the tables and the universe map\n"); print(rb, row.names = FALSE)
write.csv(rb, file.path(S, "archive_check_fisher_roots.csv"), row.names = FALSE)

# unique transcripts of interest in the universe (any ontology)
for (s in names(L)) cat(sprintf("%-6s transcripts of interest in universe: legacy reader %d, exact reader %d\n", s,
                                 length(intersect(unique(L[[s]]$transcript), universe)), length(intersect(unique(E[[s]]$transcript), universe))))

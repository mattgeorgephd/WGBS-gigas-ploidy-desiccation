# Are the margins of the legacy Fisher tests counted in transcripts or in map lines (duplicated names)?
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
x <- readRDS(file.path(S, "tables.rds")); L <- x$L; E <- x$E
suppressPackageStartupMessages(library(GO.db, lib.loc = "/home/user/gannet_cache/legacy_r_lib"))
ul <- readLines(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))
nm <- gsub(" ", "", sub("\t.*$", "", ul)); go <- sub("^[^\t]*\t", "", ul); first <- !duplicated(nm)
g2g <- setNames(lapply(strsplit(go[first], ",", fixed = TRUE), function(v) gsub(" ", "", v)), nm[first])
allgo <- unique(unlist(g2g))
ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = allgo, columns = "ONTOLOGY", keytype = "GOID")); ont <- setNames(ont$ONTOLOGY, ont$GOID)
arch <- file.path(repo, "08_archive/GOterms")
chk <- list(c("ploidy", "BP", "BOTH-PLOIDY/PLOIDY-BP-FisherTestResults.csv"), c("heat", "MF", "BOTH-HEAT/HEAT-MF-FisherTestResults.csv"),
            c("shared", "CC", "BOTH-COMBO/COMBO-CC-FisherTestResults.csv"))
for (k in chk) {
  s <- k[1]; o <- k[2]
  feas <- names(g2g)[vapply(g2g, function(v) any(ont[v] == o, na.rm = TRUE), logical(1))]
  N_u <- length(feas); M_u <- sum(feas %in% L[[s]]$transcript)
  N_d <- sum(nm %in% feas); M_d <- sum(nm %in% feas & nm %in% L[[s]]$transcript)        # counted per map line
  l <- readLines(file.path(arch, k[3]))[2:6]
  cat(sprintf("\n%s %s: margins in distinct transcripts N = %d, M = %d; in map lines N = %d, M = %d\n", s, o, N_u, M_u, N_d, M_d))
  for (ln in l) {
    p <- strsplit(ln, ",", fixed = TRUE)[[1]]; n <- length(p); a <- as.integer(p[n - 3]); sg <- as.integer(p[n - 2])
    cat(sprintf("  %s annotated %d significant %d | archived expected %s p %s | distinct-transcript margins: expected %.2f p %.2g | map-line margins: expected %.2f p %.2g\n",
                p[1], a, sg, p[n - 1], p[n], a * M_u / N_u, phyper(sg - 1, a, N_u - a, M_u, lower.tail = FALSE), a * M_d / N_d, phyper(sg - 1, a, N_d - a, M_d, lower.tail = FALSE)))
  }
}
cat(sprintf("\nmap lines per transcript: mean %.2f, max %d; transcripts on more than one line: %d of %d\n", mean(table(nm)), max(table(nm)), sum(table(nm) > 1), length(unique(nm))))

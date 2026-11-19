# 07_scope.R (scratch)
# (i) how many genes of interest take part in each ontology's test (GO.db 3.13.0)
# (ii) where on the genome the tested foreground sits, against the universe
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
x <- readRDS(file.path(S, "tables.rds")); E <- x$E; L <- x$L; universe <- x$universe; tx2gene <- x$tx2gene
options(width = 220)
suppressPackageStartupMessages(library(GO.db, lib.loc = "/home/user/gannet_cache/legacy_r_lib"))
stopifnot(as.character(packageVersion("GO.db", lib.loc = "/home/user/gannet_cache/legacy_r_lib")) == "3.13.0")

ul <- readLines(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))
nm <- gsub(" ", "", sub("\t.*$", "", ul)); go <- sub("^[^\t]*\t", "", ul); first <- !duplicated(nm)
g2g <- setNames(lapply(strsplit(go[first], ",", fixed = TRUE), function(v) gsub(" ", "", v)), nm[first])
allgo <- unique(unlist(g2g))
ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = allgo, columns = "ONTOLOGY", keytype = "GOID")); ont <- setNames(ont$ONTOLOGY, ont$GOID)
feas <- sapply(c("BP", "CC", "MF"), function(o) names(g2g)[vapply(g2g, function(v) any(ont[v] == o, na.rm = TRUE), logical(1))], simplify = FALSE)
anyfeas <- unique(unlist(feas))
cat(sprintf("universe transcripts: %d; feasible in BP %d, CC %d, MF %d; in at least one ontology %d; in none %d\n",
            length(universe), length(feas$BP), length(feas$CC), length(feas$MF), length(anyfeas), length(setdiff(universe, anyfeas))))

cat("\n(i) genes of interest that take part in a test (gene = geneID of the table row)\n")
rows <- list()
for (s in names(L)) {
  l <- unique(L[[s]][, c("transcript", "geneID")]); e <- unique(E[[s]][, c("transcript", "geneID")])
  n_all <- length(unique(e$geneID))
  f <- function(d, set) length(unique(d$geneID[d$transcript %in% set]))
  rows[[s]] <- data.frame(set = s, genes_with_GO = n_all, kept_by_reader = length(unique(l$geneID)),
    in_universe = f(l, universe), BP = f(l, feas$BP), CC = f(l, feas$CC), MF = f(l, feas$MF), any_ontology = f(l, anyfeas),
    tx_in_universe = sum(unique(l$transcript) %in% universe), tx_BP = sum(unique(l$transcript) %in% feas$BP),
    tx_CC = sum(unique(l$transcript) %in% feas$CC), tx_MF = sum(unique(l$transcript) %in% feas$MF),
    exact_in_universe = f(e, universe), exact_BP = f(e, feas$BP), exact_CC = f(e, feas$CC), exact_MF = f(e, feas$MF))
}
rows <- do.call(rbind, rows); print(rows, row.names = FALSE)
write.csv(rows, file.path(S, "genes_entering_each_test.csv"), row.names = FALSE)

# ---- (ii) genome position ----
gff <- read.table(file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_mRNA.gff.gz"), sep = "\t", quote = "", comment.char = "#",
                  colClasses = c("character", "NULL", "NULL", "integer", "integer", "NULL", "NULL", "NULL", "character"), col.names = c("chr", "a", "b", "start", "end", "c", "d", "e", "attr"))
gff$transcript <- sub("^.*transcript_id=([^;]+).*$", "\\1", gff$attr)
tx_chr <- setNames(gff$chr, gff$transcript); tx_start <- setNames(gff$start, gff$transcript)
lgs <- sprintf("NC_0475%02d.1", 59:68)
seqof <- function(t) { c <- unname(tx_chr[t]); ifelse(is.na(c), "not an mRNA of the track", ifelse(c %in% lgs, c, "scaffolds")) }
lev <- c(lgs, "scaffolds", "not an mRNA of the track")
u_seq <- table(factor(seqof(universe), levels = lev))
cat("\n(ii) transcripts by sequence: legacy universe, and the transcripts of interest that were in it, under each reader\n")
tabs <- data.frame(sequence = lev, universe = as.integer(u_seq))
for (s in names(L)) {
  tl <- intersect(unique(L[[s]]$transcript), universe); te <- intersect(unique(E[[s]]$transcript), universe)
  tabs[[paste0(s, "_legacy")]] <- as.integer(table(factor(seqof(tl), levels = lev)))
  tabs[[paste0(s, "_exact")]] <- as.integer(table(factor(seqof(te), levels = lev)))
}
print(tabs, row.names = FALSE)
first5 <- lgs[1:5]
for (s in c("ploidy", "heat")) {
  a <- tabs[[paste0(s, "_legacy")]]; b <- tabs[[paste0(s, "_exact")]]
  cat(sprintf("%-6s share of tested transcripts of interest on NC_047559 to NC_047563: legacy reader %.1f%% (%d of %d); exact reader %.1f%% (%d of %d); universe %.1f%%\n", s,
              100 * sum(a[1:5]) / sum(a), sum(a[1:5]), sum(a), 100 * sum(b[1:5]) / sum(b), sum(b[1:5]), sum(b), 100 * sum(tabs$universe[1:5]) / sum(tabs$universe)))
  cat(sprintf("%-6s ... on NC_047564 to NC_047568: legacy reader %d of %d (%.1f%%); exact reader %d of %d (%.1f%%); universe %.1f%%\n", s,
              sum(a[6:10]), sum(a), 100 * sum(a[6:10]) / sum(a), sum(b[6:10]), sum(b), 100 * sum(b[6:10]) / sum(b), 100 * sum(tabs$universe[6:10]) / sum(tabs$universe)))
}
write.csv(tabs, file.path(S, "foreground_by_chromosome.csv"), row.names = FALSE)

# The kept accession block of the ploidy table and its share of the universe
blk <- c("XM_034462509.1", "XM_034479694.1")
inblk <- function(t) t >= blk[1] & t <= blk[2]
kp <- unique(L$ploidy$transcript)
cat(sprintf("\nploidy: kept transcripts %d, of which inside the accession run %s to %s: %d (%.1f%%)\n", length(kp), blk[1], blk[2], sum(inblk(kp)), 100 * mean(inblk(kp))))
cat(sprintf("universe transcripts inside that accession run: %d of %d (%.1f%%)\n", sum(inblk(universe)), length(universe), 100 * mean(inblk(universe))))
tp <- intersect(kp, universe)
cat(sprintf("tested ploidy transcripts of interest inside the run: %d of %d (%.1f%%)\n", sum(inblk(tp)), length(tp), 100 * mean(inblk(tp))))
ub <- universe[inblk(universe)]
cat("sequence of the universe transcripts inside the run:\n"); print(table(factor(seqof(ub), levels = lev)))
# how do XM_034 accessions relate to position?
x34 <- gff[grepl("^XM_034", gff$transcript), ]; x34$acc <- as.numeric(sub("^XM_([0-9]+)\\..*$", "\\1", x34$transcript))
agg <- aggregate(acc ~ chr, data = x34[x34$chr %in% lgs, ], FUN = function(v) c(n = length(v), min = min(v), median = median(v), max = max(v)))
cat("\nXM_034 accession numbers by linkage group (n, min, median, max):\n"); print(agg)

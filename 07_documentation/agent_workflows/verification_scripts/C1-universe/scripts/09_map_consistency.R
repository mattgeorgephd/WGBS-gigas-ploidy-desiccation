suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; dl <- commandArgs(trailingOnly = TRUE)[1]
full <- fread(cmd = sprintf("zcat %s/03_input_files/annotation/geneid2go.tab.gz", repo), header = FALSE, sep = "\t", quote = "", col.names = c("id", "go"))
for (f in c("data_gene_annotations_BOTH_geneid2go-union1x-BOTH.tab", "data_gene_annotations_RONIT_geneid2go-union1x-RONIT.tab", "data_gene_annotations_YAAMINI_geneid2go-union1x-YAAMINI.tab")) {
  u <- fread(file.path(dl, f), header = FALSE, sep = "\t", quote = "", col.names = c("id", "go"))
  s <- full[id %in% u$id]
  cat(sprintf("%-62s lines %6d ; equals the full-map lines of its ids, same order: %s\n", f, nrow(u), identical(s$id, u$id) && identical(s$go, u$go)))
}
# repo copy of the BOTH universe equals the gannet copy?
a <- readLines(gzfile(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz"))); b <- readLines(file.path(dl, "data_gene_annotations_BOTH_geneid2go-union1x-BOTH.tab"))
cat("repo copy identical to gannet copy (data/gene_annotations/BOTH):", identical(a, b), "\n")

# Is a tested CpG at 1-based positions start+1 .. start+2 (0-based start)?
suppressPackageStartupMessages({library(data.table); library(Biostrings)})
scr <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C5-feature-odds"
cpg <- readRDS(file.path(scr, "cpg_features.rds"))[, .(chr, start, end)]
g <- readDNAStringSet("/home/user/gannet_cache/genome/cgigas_uk_roslin_v1_genomic-mito.fa")
names(g) <- sub(" .*", "", names(g))
cat("sequences in genome:", length(g), " total bp:", sum(as.numeric(width(g))), "\n")
cat("all CpG sequences present:", all(unique(cpg$chr) %in% names(g)), "\n")
res <- cpg[, {
  s <- g[[chr[1]]]
  d0 <- as.character(Views(s, start = start + 1L, width = 2))   # 0-based start convention
  d1 <- as.character(Views(s, start = pmax(start, 1L), width = 2))  # if start were 1-based
  .(n = .N, cg_0based = sum(toupper(d0) == "CG"), cg_1based = sum(toupper(d1) == "CG"))
}, by = chr]
cat("CpGs:", sum(res$n), " 'CG' at [start+1,start+2]:", sum(res$cg_0based),
    " 'CG' at [start,start+1]:", sum(res$cg_1based), "\n")

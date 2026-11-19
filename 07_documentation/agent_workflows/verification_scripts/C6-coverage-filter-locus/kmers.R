s <- readLines("block.fa")[2]
rc <- function(x) paste(rev(strsplit(chartr("ACGT", "TGCA", x), "")[[1]]), collapse = "")
st <- seq(1, nchar(s) - 29, by = 60)
k <- substring(s, st, st + 29)
d <- data.frame(pos = 39633066 + st - 1, kmer = k, rc = vapply(k, rc, ""), stringsAsFactors = FALSE)
write.table(d, "kmers_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(unique(c(d$kmer, d$rc)), "kmers.txt")
cat(nrow(d), "30-mers, one every 60 bp of the block\n")

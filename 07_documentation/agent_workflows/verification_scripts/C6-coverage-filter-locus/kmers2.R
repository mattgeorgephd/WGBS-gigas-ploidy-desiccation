d <- read.delim("kmers_table.tsv", stringsAsFactors = FALSE)
cnt <- read.table("kmer_counts.txt", col.names = c("n", "kmer"), stringsAsFactors = FALSE)
g <- function(k) { i <- match(k, cnt$kmer); ifelse(is.na(i), 0L, cnt$n[i]) }
d$fwd <- g(d$kmer); d$rev <- g(d$rc); d$total <- d$fwd + d$rev
print(d[, c("pos", "kmer", "fwd", "rev", "total")], row.names = FALSE)
cat("\n30-mers with exactly one genome-wide copy:", sum(d$total == 1), "of", nrow(d), "\n")
cat("with 2 copies:", sum(d$total == 2), "; 3 or more:", sum(d$total >= 3), "\n")
print(d[d$total > 1, c("pos", "total")], row.names = FALSE)

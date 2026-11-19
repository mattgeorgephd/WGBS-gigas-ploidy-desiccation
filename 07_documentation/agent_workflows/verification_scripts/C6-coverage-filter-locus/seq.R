s <- readLines("locus_seq_39625001_39645000.txt", warn = FALSE)
a <- 39625001L
sub1 <- function(from, to) substr(s, from - a + 1, to - a + 1)     # 1-based genome coordinates
# check CpG positions: coverage-file start is 0-based, so the C is at start + 1
for (st in c(39633323L, 39635520L, 39633065L, 39636662L)) cat(st, ":", sub1(st + 1, st + 2), "\n")
prof <- function(from, to, label) {
  x <- sub1(from, to); ch <- strsplit(x, "")[[1]]
  up <- toupper(ch)
  cat(sprintf("%-28s %d-%d len %d | GC %.1f%% | lowercase (soft-masked) %.1f%% | N %d | CpG sites %d | CpG o/e %.2f\n", label, from, to, length(ch),
              100 * mean(up %in% c("G", "C")), 100 * mean(ch %in% c("a", "c", "g", "t")), sum(up == "N"),
              lengths(regmatches(toupper(x), gregexpr("CG", toupper(x)))),
              lengths(regmatches(toupper(x), gregexpr("CG", toupper(x)))) * length(ch) / (sum(up == "C") * sum(up == "G"))))
}
prof(39633066, 39636664, "high-coverage block")
prof(39633324, 39635522, "tested-CpG span")
prof(39628066, 39633065, "upstream 5 kb")
prof(39636665, 39641664, "downstream 5 kb")
prof(39625001, 39645000, "whole 20 kb window")
# 500 bp bins: masked fraction and N
cat("\n500 bp bins: start, %soft-masked, N count, CpG sites\n")
for (st in seq(39630323L, 39639822L, by = 500L)) {
  x <- sub1(st + 1L, st + 500L); ch <- strsplit(x, "")[[1]]
  cat(sprintf("%d  masked %5.1f%%  N %3d  CpG %2d\n", st, 100 * mean(ch %in% c("a", "c", "g", "t")), sum(toupper(ch) == "N"),
              lengths(regmatches(toupper(x), gregexpr("CG", toupper(x))))))
}
writeLines(c(">NC_047563.1:39633066-39636664", toupper(sub1(39633066, 39636664))), "block.fa")
writeLines(c(">NC_047563.1:39632324-39637323", toupper(sub1(39632324, 39637323))), "block_wide.fa")

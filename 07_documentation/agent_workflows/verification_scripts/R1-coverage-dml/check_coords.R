suppressPackageStartupMessages(library(data.table))
x <- readLines("two_seqs.fa"); h <- grep("^>", x)
seqs <- setNames(vapply(seq_along(h), function(i) toupper(paste(x[(h[i] + 1):(if (i < length(h)) h[i + 1] - 1 else length(x))], collapse = "")), ""), sub("^>(\\S+).*", "\\1", x[h]))
print(nchar(seqs))
for (f in c("zr3534_1", "zr3644_13")) {
  d <- fread(cmd = sprintf("LC_ALL=C awk -F'\\t' '$1==\"NC_001276.1\" || $1==\"NW_022994895.1\"' /home/user/gannet_cache/cov_files/%s_R1.CpG_report.merged_CpG_evidence.cov", f), header = FALSE)
  for (s in names(seqs)) {
    dd <- d[V1 == s]; cg0 <- as.integer(gregexpr("CG", seqs[[s]], fixed = TRUE)[[1]]) - 1L     # 0-based position of each C of a CG
    cat(sprintf("%s %s: %d rows in the coverage file ; CG sites in the reference %d ; rows whose start is the 0-based position of a C followed by G: %d ; rows where start is 1-based instead: %d ; end == start+2: %s\n",
                f, s, nrow(dd), length(cg0), sum(dd$V2 %in% cg0), sum((dd$V2 - 1L) %in% cg0), all(dd$V3 == dd$V2 + 2L)))
  }
}

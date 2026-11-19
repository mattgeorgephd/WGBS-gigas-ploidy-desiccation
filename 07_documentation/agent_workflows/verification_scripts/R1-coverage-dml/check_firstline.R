suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
lcv <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
fp <- unique(lcv$first_line_position)
cov <- sapply(1:22, function(j) { r <- fread(sprintf("rawcheck/%s.sel", lib$cov_file[j]), header = FALSE); k <- paste0(r$V1, ":", r$V2); v <- (r$V5 + r$V6)[match(fp, k)]; v[is.na(v)] <- 0; v })
rownames(cov) <- fp; colnames(cov) <- 1:22
isfirst <- sapply(1:22, function(j) fp == lcv$first_line_position[j])
cat("first-line CpGs:", length(fp), "\n")
only_first <- rowSums(cov >= 1 & !isfirst) == 0
cat("first-line CpGs covered (>=1 read) only in libraries where they are the first line:", sum(only_first), "->", paste(fp[only_first], collapse = ", "), "\n")
cat("so the legacy 'any library' count at 1x should be 12918252 -", sum(only_first), "=", 12918252 - sum(only_first), "\n")
cat("first-line CpGs with >= 10 reads in all 22 libraries (complete read):", sum(rowSums(cov >= 10) == 22), "\n")
cat("first-line CpGs with >= 1 read in all 22 libraries (complete read):", sum(rowSums(cov >= 1) == 22), "; still in all 22 after legacy loss:", sum(rowSums(cov >= 1 & !isfirst) == 22), "\n")
for (x in c(1, 3, 5, 10, 20)) cat(sprintf("  %dx: first-line CpGs in all 22 with complete read: %d ; lost in legacy: %d\n", x, sum(rowSums(cov >= x) == 22), sum(rowSums(cov >= x) == 22 & rowSums(cov >= x & !isfirst) < 22)))
print(cov)

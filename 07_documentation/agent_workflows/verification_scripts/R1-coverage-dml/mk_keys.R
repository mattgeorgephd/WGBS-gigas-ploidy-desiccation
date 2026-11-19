suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tdir <- file.path(repo, "04_output/05_differential_methylation/tables")
u <- fread(file.path(tdir, "united_cpg_legacy.csv.gz"), select = c("chr", "start"))
set.seed(1)
idx <- sort(unique(c(sample(nrow(u), 3000), fread(file.path(tdir, "dml_shared_legacy.csv"))$united_row,
                     fread(file.path(tdir, "dml_ploidy_legacy.csv"))$united_row[1:200], 1:20, nrow(u) - 0:19)))
lcv <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
fp <- unique(lcv$first_line_position)
keys <- rbind(u[idx], data.table(chr = sub(":.*", "", fp), start = as.integer(sub(".*:", "", fp))))
fwrite(keys, "rawcheck/keys.tsv", sep = "\t", col.names = FALSE)
fwrite(data.table(united_row = idx), "rawcheck/idx.csv")
cat(nrow(keys), "keys\n")

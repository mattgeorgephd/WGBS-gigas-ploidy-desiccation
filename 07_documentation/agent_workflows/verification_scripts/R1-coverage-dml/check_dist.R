suppressPackageStartupMessages(library(data.table))
a <- fread("dist_awk.tsv", col.names = c("library_id", "k", "n1", "n10"))
labs <- c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
d <- fread("/home/user/WGBS-gigas-ploidy-desiccation/04_output/04_cpg_coverage/tables/methylation_distribution_by_library.csv")[mode == "legacy"]
d[, k := match(bin, labs)]
m <- merge(a, dcast(d, library_id + k ~ min_coverage, value.var = "n_cpg"), by = c("library_id", "k"))
cat("bins compared:", nrow(m), "; 1x bins equal:", sum(m$n1 == m$`1`), "; 10x bins equal:", sum(m$n10 == m$`10`), "\n")
print(m[n1 != `1` | n10 != `10`])

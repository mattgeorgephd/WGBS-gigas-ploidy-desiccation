suppressPackageStartupMessages({library(readxl); library(data.table)})
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; f <- file.path(repo, "08_archive/analyses/coverage_summary.xlsx")
print(excel_sheets(f)); options(width = 200)
tc <- as.data.table(read_excel(f, sheet = "total_cpg")); print(tc)
cc <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/cpg_counts_by_coverage.csv"))[mode == "legacy"]
pp <- as.data.table(read_excel(f, sheet = "ploidy_plot")); print(head(pp, 12))
tp <- as.data.table(read_excel(f, sheet = "trt_plot")); print(head(tp, 12))
dg <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/methylation_distribution_by_group.csv"))[mode == "legacy" & min_coverage == 1]
m <- merge(pp, dg[grouping %in% c("ploidy", "all"), .(ploidy = sub("all libraries", "both", group), bin, mean_percent, se_percent)], by = c("ploidy", "bin"))
cat("ploidy_plot vs regenerated: max |percent - mean_percent| =", max(abs(m$percent - m$mean_percent)), " ; max |se diff| =", max(abs(m$se - m$se_percent)), "\n")
print(m[order(-abs(percent - mean_percent))][1:6])

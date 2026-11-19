suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
tdir <- file.path(repo, "04_output/05_differential_methylation/tables")
ul <- fread(file.path(tdir, "united_cpg_legacy.csv.gz")); uc <- fread(file.path(tdir, "united_cpg_corrected.csv.gz"))
fc <- fread(file.path(tdir, "library_filtering_corrected.csv"))
covc <- paste0("coverage", 1:22)
cl <- as.matrix(ul[, ..covc])
lim <- round(fc$scaling_factor * (ceiling(fc$upper_cutoff) - 1))   # largest scaled coverage a kept CpG can have
pass <- rowSums(sweep(cl, 2, lim, ">")) == 0
cat("legacy united rows passing every library's upper cutoff:", sum(pass), " corrected united rows:", nrow(uc), "\n")
kl <- paste(ul$chr, ul$start)[pass]; kc <- paste(uc$chr, uc$start)
cat("same CpG set:", setequal(kl, kc), " same order:", identical(kl, kc), "\n")
cat("corrected CpGs not in legacy united:", sum(!kc %in% paste(ul$chr, ul$start)), "\n")
if (identical(kl, kc)) cat("identical count columns for shared CpGs:", isTRUE(all.equal(as.data.frame(ul[pass, -(1:4)]), as.data.frame(uc[, -(1:4)]), check.attributes = FALSE)), "\n")
# how many libraries exceed the cutoff at removed CpGs
ex <- rowSums(sweep(cl, 2, lim, ">"))
cat("removed CpGs by number of libraries at or above their cutoff:\n"); print(table(ex[!pass]))
# I-09 region
reg <- ul$chr == "NC_047563.1" & ul$start >= 39633323 & ul$start <= 39635520
cat("legacy united CpGs in NC_047563.1:39,633,323-39,635,520:", sum(reg), " of which survive in corrected:", sum(reg & pass), "\n")
# first-line CpGs
lcv <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
fp <- unique(lcv$first_line_position)
cat("first-line CpGs present in legacy united:", sum(fp %in% paste0(ul$chr, ":", ul$start)), " in corrected united:", sum(fp %in% paste0(uc$chr, ":", uc$start)), " of", length(fp), "distinct\n")
# DML lost / gained between modes
for (ct in c("ploidy", "heat")) {
  dl <- fread(file.path(tdir, sprintf("dml_%s_legacy.csv", ct))); dc <- fread(file.path(tdir, sprintf("dml_%s_corrected.csv", ct)))
  a <- paste(dl$chr, dl$start); b <- paste(dc$chr, dc$start)
  cat(sprintf("%s DML: legacy %d, corrected %d, common %d, legacy-only %d (of which CpG removed by filter: %d), corrected-only %d\n", ct, length(a), length(b), length(intersect(a, b)),
              length(setdiff(a, b)), sum(!setdiff(a, b) %in% kc), length(setdiff(b, a))))
}

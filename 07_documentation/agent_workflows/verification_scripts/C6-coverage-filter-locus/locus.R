suppressPackageStartupMessages(library(data.table))
options(width = 250)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))[, .(library_id, cohort, ploidy, treatment, acute_heat_shock)]
lc  <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/library_coverage.csv"))
cut <- fread("cutoffs_mine.csv")
lf  <- fread(t5("library_filtering_legacy.csv"))
reg <- rbindlist(lapply(lib$library_id, function(l) fread(sprintf("pass/region_%s.tsv", l), header = FALSE,
          col.names = c("library_id", "chr", "start", "end", "pct", "nC", "nT"))))
reg[, cov := nC + nT]
L0 <- 39633323L; L1 <- 39635520L
# the tested CpGs of the locus, from the pipeline's legacy united table
u <- fread(t5("united_cpg_legacy.csv.gz"), select = c("chr", "start"))
tested <- u[chr == "NC_047563.1" & start >= L0 & start <= L1, start]
cat("tested CpGs in NC_047563.1:", L0, "-", L1, ":", length(tested), " (first", min(tested), ", last", max(tested), ")\n")
uc <- fread(t5("united_cpg_corrected.csv.gz"), select = c("chr", "start"))
cat("of these, present in the corrected united table:", uc[chr == "NC_047563.1" & start %in% tested, .N], "\n")
cat("tested CpGs within 5 kb either side in legacy:", u[chr == "NC_047563.1" & start >= L0 - 5000 & start < L0, .N], "upstream,",
    u[chr == "NC_047563.1" & start > L1 & start <= L1 + 5000, .N], "downstream; in corrected:",
    uc[chr == "NC_047563.1" & start >= L0 - 5000 & start < L0, .N], uc[chr == "NC_047563.1" & start > L1 & start <= L1 + 5000, .N], "\n")
x <- merge(reg[start %in% tested], cut[, .(library_id, cutoff = cutoff_mine, median_all, mean_all)], by = "library_id")
x <- merge(x, lc[, .(library_id, upper_cutoff)], by = "library_id")
x <- merge(x, lf[, .(library_id, median_ge10 = median_coverage)], by = "library_id")
x <- merge(x, lib, by = "library_id")
stopifnot(nrow(x) == length(tested) * 22)
x[, over := cov >= cutoff]; stopifnot(all(x$over == (x$cov >= x$upper_cutoff)))
# ---- per CpG: in how many libraries at or above the cutoff ----
per_cpg <- x[, .(n_lib_over = sum(over), n_heat_over = sum(over & treatment == "heat"), n_control_over = sum(over & treatment == "control"),
                 n_dip_over = sum(over & ploidy == "diploid"), n_trip_over = sum(over & ploidy == "triploid"),
                 min_cov = min(cov), max_cov = max(cov)), by = start][order(start)]
cat("\nEvery one of the", nrow(per_cpg), "CpGs at or above the cutoff in at least one library:", all(per_cpg$n_lib_over >= 1), "\n")
cat("number of libraries at or above the cutoff, per CpG: min", min(per_cpg$n_lib_over), " median", median(per_cpg$n_lib_over), " max", max(per_cpg$n_lib_over), "\n")
cat("distribution (libraries over : number of CpGs):\n"); print(table(per_cpg$n_lib_over))
cat("per CpG, libraries over the cutoff among 10 heat: range", range(per_cpg$n_heat_over), "; among 12 control: range", range(per_cpg$n_control_over),
    "; among 11 diploid:", range(per_cpg$n_dip_over), "; among 11 triploid:", range(per_cpg$n_trip_over), "\n")
print(per_cpg)
# ---- per library ----
per_lib <- x[, .(n_cpg = .N, n_over = sum(over), median_cov = as.numeric(median(cov)), min_cov = min(cov), max_cov = max(cov),
                 cutoff = cutoff[1], median_all = median_all[1], mean_all = round(mean_all[1], 2), median_ge10 = median_ge10[1],
                 pooled_meth = round(100 * sum(nC) / sum(cov), 1), median_meth = round(median(100 * nC / cov), 1)),
             by = .(library_id, cohort, ploidy, treatment, acute_heat_shock)]
per_lib[, `:=`(ratio_to_median_all = round(median_cov / median_all, 2), ratio_to_mean_all = round(median_cov / mean_all, 2),
               ratio_to_median_ge10 = round(median_cov / median_ge10, 2), ratio_to_cutoff = round(median_cov / cutoff, 2))]
per_lib <- per_lib[match(lib$library_id, library_id)]
print(per_lib[, .(library_id, ploidy, treatment, hs = acute_heat_shock, n_over, median_cov, min_cov, max_cov, cutoff, median_all, median_ge10,
                  ratio_to_median_all, ratio_to_mean_all, ratio_to_median_ge10, ratio_to_cutoff, pooled_meth)])
cat("\nlibraries with at least one of the 54 CpGs at or above the cutoff:", sum(per_lib$n_over > 0), "of 22\n")
cat("libraries with all 54 at or above:", sum(per_lib$n_over == 54), "; with none:", sum(per_lib$n_over == 0), "\n")
print(per_lib[, .(n_lib = .N, lib_with_any_over = sum(n_over > 0), cpg_lib_pairs_over = sum(n_over), pairs = sum(n_cpg),
                  pct_pairs_over = round(100 * sum(n_over) / sum(n_cpg), 1),
                  ratio_median_all_min = min(ratio_to_median_all), ratio_median_all_median = median(ratio_to_median_all), ratio_median_all_max = max(ratio_to_median_all),
                  ratio_ge10_median = median(ratio_to_median_ge10), meth_min = min(pooled_meth), meth_median = median(pooled_meth), meth_max = max(pooled_meth)),
              by = .(treatment, ploidy)])
print(per_lib[, .(n_lib = .N, lib_with_any_over = sum(n_over > 0), pct_pairs_over = round(100 * sum(n_over) / sum(n_cpg), 1),
                  ratio_median_all_min = min(ratio_to_median_all), ratio_median_all_median = median(ratio_to_median_all), ratio_median_all_max = max(ratio_to_median_all),
                  meth_median = median(pooled_meth)), by = treatment])
print(per_lib[, .(n_lib = .N, lib_with_any_over = sum(n_over > 0), pct_pairs_over = round(100 * sum(n_over) / sum(n_cpg), 1),
                  ratio_median_all_min = min(ratio_to_median_all), ratio_median_all_median = median(ratio_to_median_all), ratio_median_all_max = max(ratio_to_median_all),
                  meth_median = median(pooled_meth)), by = ploidy])
cat("overall: ratio of region median coverage to library median (all covered CpGs): range", range(per_lib$ratio_to_median_all), " median", median(per_lib$ratio_to_median_all), "\n")
cat("overall: ratio to library median among CpGs with >= 10 reads: range", range(per_lib$ratio_to_median_ge10), " median", median(per_lib$ratio_to_median_ge10), "\n")
cat("correlation across 22 libraries between coverage ratio and pooled methylation in the region: Pearson",
    round(cor(per_lib$ratio_to_median_all, per_lib$pooled_meth), 3), " Spearman", round(cor(per_lib$ratio_to_median_all, per_lib$pooled_meth, method = "spearman"), 3), "\n")
# which single library removes most
cat("CpGs of the 54 that only ONE library removes:", sum(per_cpg$n_lib_over == 1), "\n")
ov <- dcast(x, start ~ library_id, value.var = "over")
m <- as.matrix(ov[, -1]); 
cat("smallest set of libraries needed: CpGs covered by the top library alone:", max(colSums(m)), "(", names(which.max(colSums(m))), ")\n")
fwrite(per_lib, "locus_per_library.csv"); fwrite(per_cpg, "locus_per_cpg.csv")
saveRDS(list(reg = reg, x = x, tested = tested), "locus.rds")

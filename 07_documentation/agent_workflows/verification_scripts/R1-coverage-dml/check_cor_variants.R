suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"; t5 <- file.path(repo, "04_output/05_differential_methylation/tables")
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv")); setorder(lib, legacy_sample_index)
u <- fread(file.path(t5, "united_cpg_legacy.csv.gz"))
nc <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); nt <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); mm <- nc / (nc + nt)
rng <- function(r, idx = 1:22) { o <- r[idx, idx][upper.tri(r[idx, idx])]; sprintf("%.3f to %.3f", min(o), max(o)) }
rp <- cor(mm); rs <- cor(mm, method = "spearman")
h <- which(lib$treatment == "heat"); c0 <- which(lib$treatment == "control"); d <- which(lib$ploidy == "diploid"); t3 <- which(lib$ploidy == "triploid")
cat("10x, 480,579 CpGs. Pearson all:", rng(rp), "| heat only:", rng(rp, h), "| control only:", rng(rp, c0), "| diploid:", rng(rp, d), "| triploid:", rng(rp, t3), "\n")
cat("Spearman all:", rng(rs), "| heat only:", rng(rs, h), "| control only:", rng(rs, c0), "\n")
el <- rowMeans(100 * mm) >= 10 & rowMeans(100 * mm) <= 90
cat("Pearson on CpGs with mean 10-90%:", rng(cor(mm[el, ])), "\n")

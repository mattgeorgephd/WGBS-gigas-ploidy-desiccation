# How do the exploratory rows of part 6 change if "high" is cut at the second largest gap (two libraries) instead of the largest (five)?
# Exact integer arithmetic, all CpGs, both modes.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE)
bis <- read.csv(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"), stringsAsFactors = FALSE)
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG); o <- order(L); g <- diff(L[o])
cut2 <- order(g, decreasing = TRUE)[2]; high2 <- seq_len(22) %in% o[(cut2 + 1):22]
cat("second largest gap:", round(sort(g, decreasing = TRUE)[2], 4), "; libraries above it:", paste(lib$library_id[high2], collapse = " "), "\n")
for (mode in c("legacy", "corrected")) {
  u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
  Cl <- lapply(1:22, function(j) as.numeric(u[[paste0("numCs", j)]])); Nl <- lapply(1:22, function(j) Cl[[j]] + u[[paste0("numTs", j)]])
  add <- function(lst, idx) Reduce(`+`, lst[idx])
  for (coh in c("heat", "control")) {
    a0 <- which(lib$treatment == coh & lib$ploidy == "diploid" & !high2); b0 <- which(lib$treatment == coh & lib$ploidy == "triploid" & !high2)
    libs <- c(a0, b0); ctot <- add(Cl, libs); ntot <- add(Nl, libs)
    gr <- combn(libs, length(a0), simplify = FALSE)
    cnt <- vapply(gr, function(A) { ca <- add(Cl, A); na <- add(Nl, A); cb <- ctot - ca; nb <- ntot - na; sum(5 * abs(cb * na - ca * nb) > na * nb) }, numeric(1))
    obs <- cnt[vapply(gr, function(A) setequal(A, a0), logical(1))]
    cat(sprintf("%s, %s cohort without %s: %d diploid against %d triploid; observed %d; %d groupings; median %s; largest other %d; at or above %d; p = %.4f\n", mode, coh,
                paste(lib$library_id[lib$treatment == coh & high2], collapse = " "), length(a0), length(b0), obs, length(cnt), format(median(cnt)),
                max(cnt[cnt != obs | !vapply(gr, function(A) setequal(A, a0), logical(1))][!vapply(gr, function(A) setequal(A, a0), logical(1))]), sum(cnt >= obs), mean(cnt >= obs)))
  }
  rm(u, Cl, Nl); invisible(gc())
}

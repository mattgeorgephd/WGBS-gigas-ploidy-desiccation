# 00_cache_counts.R
# Read the united count table of one mode once and keep compact integer
# matrices (C = methylated reads, N = all reads) in the scratch folder, so the
# later scripts do not have to parse the 480,579-row table again.
# Independent code: does not source anything from the repository.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "legacy"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
u <- fread(file.path(repo, "04_output/05_differential_methylation/tables", sprintf("united_cpg_%s.csv.gz", mode)))
stopifnot(ncol(u) == 4 + 3 * 22)
C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); storage.mode(C) <- "integer"
Tm <- as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); storage.mode(Tm) <- "integer"
cov <- as.matrix(u[, paste0("coverage", 1:22), with = FALSE]); storage.mode(cov) <- "integer"
N <- C + Tm
cat("rows:", nrow(u), "\n")
cat("coverage == numCs + numTs everywhere:", all(cov == N), "\n")
cat("minimum coverage in the table:", min(cov), "\n")
dimnames(C) <- dimnames(N) <- NULL
saveRDS(list(chr = u$chr, start = u$start, C = C, N = N), file.path(scratch, sprintf("counts_%s.rds", mode)), compress = FALSE)
cat("saved\n")

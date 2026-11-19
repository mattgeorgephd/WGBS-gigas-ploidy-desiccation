# Independent replicate of 06_diagnostics/go_random_cpg_control.R (different seed, my own implementation),
# recording p-value calibration as well as counts.
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages(library(GenomicRanges))
set.seed(991); n_draws <- as.integer(commandArgs(TRUE)[1])
CP <- readRDS(file.path(sc, "cpg_table_legacy.rds")); st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); tmf <- readRDS(file.path(sc, "tested_setup_legacy.rds"))
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
elig_rows <- which(CP$mean_pm >= 10 & CP$mean_pm <= 90)
h <- suppressWarnings(findOverlaps(GRanges(CP$chr[elig_rows], IRanges(CP$start[elig_rows] + 1L, CP$start[elig_rows] + 2L)), GRanges(gt$chr, IRanges(gt$start, gt$end)), ignore.strand = TRUE))
pairs <- data.table(cpg = queryHits(h), gene_id = gt$gene_id[subjectHits(h)]); setkey(pairs, cpg)
uniE <- cbind(gt[n_cpg_eligible > 0], bias = gt[n_cpg_eligible > 0, n_cpg_eligible]); uniT <- gt[n_cpg_tested > 0]
sizes <- c(ploidy = 1038, heat = 2351)
one <- function(genes) rbindlist(lapply(names(tmf), function(ont) {
  gT <- tmf[[ont]]$genes; MT <- tmf[[ont]]$M; MT <- MT[rowSums(MT) >= 5, ]; deT <- as.numeric(gT %in% genes)
  A <- sum(p.adjust(fisher_p(MT, deT), "BH") < 0.05)
  de <- as.integer(uniE$gene_id %in% genes); names(de) <- uniE$gene_id; pwf <- fit_pwf(de, uniE$bias)
  g <- st$tm[[ont]]$genes; M <- st$tm[[ont]]$M; M <- M[rowSums(M) >= 5, ]
  r <- wallenius_p(M, de[g], pwf[g]); B <- sum(p.adjust(fisher_p(M, de[g]), "BH") < 0.05)
  data.table(ontology = ont, genes_hit = length(genes), A = A, B = B, C = sum(p.adjust(r$p_wallenius, "BH") < 0.05), C_min_q = min(p.adjust(r$p_wallenius, "BH")), w_p05 = mean(r$p_wallenius < 0.05), w_p01 = mean(r$p_wallenius < 0.01),
             w_p05_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.05), w_p01_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.01), w_p001_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.001))
}))
out <- rbindlist(lapply(names(sizes), function(s) rbindlist(lapply(seq_len(n_draws), function(d) { drawn <- sample.int(length(elig_rows), sizes[[s]]); cbind(list_size = s, draw = d, one(unique(pairs[.(drawn), nomatch = NULL]$gene_id))) }))))
fwrite(out, file.path(sc, "random_cpg_replicate.csv"))
cat("Replicate of the random-CpG control (", n_draws, "draws per list size):\n")
print(out[, .(draws = .N, genes_hit = median(genes_hit), A_median = median(A), A_min = min(A), A_max = max(A), B_median = median(B), B_min = min(B), B_max = max(B), C_median = median(C), C_max = max(C), C_any = mean(C > 0), C_min_q_p10 = round(quantile(C_min_q, .1), 3),
              w_p05 = round(mean(w_p05), 4), w_p01 = round(mean(w_p01), 4), w_p05_big = round(mean(w_p05_big), 4), w_p01_big = round(mean(w_p01_big), 4), w_p001_big = round(mean(w_p001_big), 5)), by = .(list_size, ontology)])

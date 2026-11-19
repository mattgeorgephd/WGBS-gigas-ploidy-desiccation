# A harder version of the repository's random-CpG control: CpGs are drawn from the eligible set with probability
# proportional to their label-permutation probability of passing the 20-point rule (so the draws follow the
# mean- and variance-dependent propensity of real DML calls), then tested like real DML (analysis C).
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages(library(GenomicRanges))
set.seed(1515); n_draws <- as.integer(commandArgs(TRUE)[1])
ps <- readRDS(file.path(sc, "perm_sets_ploidy_legacy.rds")); CP <- readRDS(file.path(sc, "cpg_table_legacy.rds")); st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); uni <- st$uni; tm <- st$tm
sets <- ps$sets[lengths(ps$sets) > 0]; null_rate <- tabulate(unlist(sets), nbins = nrow(CP)) / length(sets)
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
elig <- which(CP$eligible); w <- null_rate[elig]
cat(sprintf("eligible CpGs %d; with non-zero propensity %d\n", length(elig), sum(w > 0)))
h <- suppressWarnings(findOverlaps(GRanges(CP$chr[elig], IRanges(CP$start[elig] + 1L, CP$start[elig] + 2L)), GRanges(gt$chr, IRanges(gt$start, gt$end)), ignore.strand = TRUE))
pairs <- data.table(cpg = queryHits(h), gene_id = gt$gene_id[subjectHits(h)]); setkey(pairs, cpg)
prep <- lapply(tm, function(t) list(genes = t$genes, M = t$M[rowSums(t$M) >= 5, ]))
one <- function(genes) { de <- as.integer(uni$gene_id %in% genes); names(de) <- uni$gene_id; pwf <- fit_pwf(de, uni$bias)
  rbindlist(lapply(names(prep), function(ont) { P <- prep[[ont]]; r <- wallenius_p(P$M, de[P$genes], pwf[P$genes]); q <- p.adjust(r$p_wallenius, "BH")
    data.table(ontology = ont, de = sum(de[P$genes]), n_q05 = sum(q < 0.05), min_q = min(q), n_q05_min15 = sum(p.adjust(r$p_wallenius[r$n_cat >= 15], "BH") < 0.05), p_organelle = if (ont == "CC") r[go_id == "GO:0043226", p_wallenius] else NA_real_,
               w_p05_big = mean(r$p_wallenius[r$n_cat >= 50] < 0.05), fisher_q05 = sum(p.adjust(fisher_p(P$M, de[P$genes]), "BH") < 0.05)) })) }
sizes <- c(ploidy = 1038, heat = 2351)
out <- rbindlist(lapply(names(sizes), function(nm) rbindlist(lapply(seq_len(n_draws), function(d) { drawn <- sample.int(length(elig), sizes[[nm]], prob = w); gg <- unique(pairs[.(drawn), nomatch = NULL]$gene_id); cbind(list_size = nm, draw = d, genes_hit = length(gg), one(gg)) }))))
fwrite(out, file.path(sc, "propensity_control.csv"))
print(out[, .(draws = .N, genes_hit = median(genes_hit), any_q05 = round(mean(n_q05 > 0), 3), max_terms = max(n_q05), any_q05_min15 = round(mean(n_q05_min15 > 0), 3), median_min_q = round(median(min_q), 3), q10_min_q = round(quantile(min_q, .1), 3),
              w_p05_big = round(mean(w_p05_big), 4), fisher_median = median(fisher_q05), organelle_p_le_real = round(mean(p_organelle <= 1.174e-4), 3), organelle_p_median = signif(median(p_organelle), 2)), by = .(list_size, ontology)])
cat("draws with a hit in any ontology:\n"); print(out[, .(any = any(n_q05 > 0)), by = .(list_size, draw)][, .(share = mean(any)), by = list_size])

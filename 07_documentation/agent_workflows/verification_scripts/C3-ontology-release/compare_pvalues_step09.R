# Raw p-values: my re-run from legacy inputs (fresh topGOdata per set, as in the legacy script)
# against the committed step-09 CSV (updateGenes on one object per ontology, regenerated inputs).
suppressPackageStartupMessages(library(data.table))
scr <- getwd(); repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
sets <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")
out <- list()
for (s in names(sets)) for (o in c("BP", "CC", "MF")) {
  p <- readRDS(file.path(scr, "runs/godb313_topgo254", sprintf("%s-%s-pvalues.rds", sets[[s]], o)))
  n <- fread(file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_fisher_%s_%s.csv", s, o)))
  stopifnot(setequal(names(p), n$go_id))
  d <- abs(p[n$go_id] - n$p_fisher) / pmax(n$p_fisher, 1e-300)
  out[[length(out) + 1]] <- data.table(set = s, ont = o, terms = nrow(n), max_rel_diff = max(d), n_rel_diff_gt_1e12 = sum(d > 1e-12),
                                       min_p = min(n$p_fisher), n_p_exact_equal = sum(p[n$go_id] == n$p_fisher))
}
print(rbindlist(out))

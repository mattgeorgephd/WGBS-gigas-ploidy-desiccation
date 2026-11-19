# (d) Does the conclusion depend on how well goseq's spline fits the bias? Alternative weighting functions.
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages(library(mgcv))
st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); uni <- st$uni; tm <- st$tm
pwfs <- function(de, bias) {
  d <- data.frame(de = de, lb = log2(bias))
  list(goseq_raw   = fit_pwf(setNames(de, uni$gene_id), bias),
       goseq_log2  = fit_pwf(setNames(de, uni$gene_id), log2(bias) + 1),          # same goseq spline on a log scale (+1 keeps it positive)
       gam_log2    = setNames(as.numeric(predict(gam(de ~ s(lb, k = 8), binomial, d), type = "response")), uni$gene_id),
       isotonic    = { o <- order(bias); f <- isoreg(bias[o], de[o]); p <- numeric(length(de)); p[o] <- f$yf; setNames(pmax(p, 1e-3), uni$gene_id) },
       exact_count = { b <- pmin(bias, 30); bb <- ifelse(bias > 100, 101, ifelse(bias > 50, 51, ifelse(bias > 30, 31, b))); setNames(ave(de, bb), uni$gene_id) })
}
res <- list()
for (s in c("ploidy", "heat", "shared")) {
  de <- as.integer(uni[[paste0("n_dml_", s)]] > 0); P <- pwfs(de, uni$bias)
  d <- data.table(bias = uni$bias, de = de, as.data.table(P)); d[, bin := cut(bias, c(0, 1, 2, 3, 5, 10, 20, 50, 100, Inf))]
  cat("\n", s, ": observed share with a DML and mean fitted weight, by number of eligible CpGs\n")
  print(d[, c(list(genes = .N, observed = round(mean(de), 3)), lapply(.SD, function(x) round(mean(x), 3))), by = bin, .SDcols = names(P)][order(bin)])
  for (nm in names(P)) for (ont in names(tm)) {
    g <- tm[[ont]]$genes; M <- tm[[ont]]$M
    r <- wallenius_p(M, setNames(de, uni$gene_id)[g], P[[nm]][g])
    for (k in c(5, 10, 15, 20, 50, 100)) { b <- bh_count(r$p_wallenius, r$n_cat, k); kk <- r$n_cat >= k
      res[[length(res) + 1]] <- data.table(set = s, pwf = nm, ontology = ont, min_genes = k, terms = b[["tested"]], n_q05 = b[["n_sig"]], min_q = round(b[["min_q"]], 4), top = r[kk][which.min(p_wallenius), go_id], top_n = r[kk][which.min(p_wallenius), n_cat]) }
  }
}
res <- rbindlist(res); fwrite(res, file.path(sc, "pwf_alternatives.csv"))
cat("\nTerms with BH q < 0.05 (eligible universe, Wallenius) by weighting function and minimum term size:\n")
print(dcast(res, set + ontology + pwf ~ min_genes, value.var = "n_q05"))
cat("\nSmallest q:\n"); print(dcast(res, set + ontology + pwf ~ min_genes, value.var = "min_q"))
cat("\nTop term at min_genes = 5:\n"); tt <- res[min_genes == 5]; tt[, term := AnnotationDbi::Term(GO.db::GOTERM[top])]; print(tt[, .(set, ontology, pwf, min_q, top_n, term = substr(term, 1, 45))])

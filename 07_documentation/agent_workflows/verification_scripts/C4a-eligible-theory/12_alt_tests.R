# (d) An independent test that does not use the Wallenius two-class approximation:
#     logistic regression of "gene holds a DML" on a spline of log2(eligible CpGs), with a score test for each GO term
#     (large terms) and an exact Poisson-binomial tail (small terms).
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages(library(splines))
st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); uni <- st$uni; tm <- st$tm
pb_tail <- function(p, x) { if (x <= 0) return(1); pmf <- 1; for (pi in p) pmf <- c(pmf * (1 - pi), 0) + c(0, pmf * pi); sum(pmf[(x + 1):length(pmf)]) }
out <- list()
for (s in c("ploidy", "heat", "shared")) {
  de_all <- as.integer(uni[[paste0("n_dml_", s)]] > 0); names(de_all) <- uni$gene_id; pwf <- fit_pwf(de_all, uni$bias)
  for (ont in names(tm)) {
    g <- tm[[ont]]$genes; M <- tm[[ont]]$M; y <- de_all[g]; lb <- log2(uni$bias[match(g, uni$gene_id)])
    X <- cbind(1, ns(lb, df = 5)); f <- glm.fit(X, y, family = binomial()); ph <- f$fitted.values; w <- ph * (1 - ph)
    n_cat <- as.vector(rowSums(M)); x <- as.vector(M %*% y); U <- x - as.vector(M %*% ph)
    MWX <- as.matrix(M %*% (w * X)); V <- as.vector(M %*% w) - rowSums((MWX %*% solve(crossprod(X, w * X))) * MWX)
    p_score <- pnorm(U / sqrt(pmax(V, 1e-12)), lower.tail = FALSE)
    Mt <- as(M, "TsparseMatrix"); idx <- split(Mt@j + 1L, Mt@i + 1L)
    small <- which(n_cat < 100 & n_cat >= 1)
    p_pb <- rep(NA_real_, nrow(M)); p_pb[small] <- vapply(small, function(k) pb_tail(ph[idx[[as.character(k)]]], x[k]), numeric(1))
    p_alt <- ifelse(n_cat < 100, p_pb, p_score)
    wal <- wallenius_p(M, y, pwf[g])
    d <- data.table(set = s, ontology = ont, go_id = rownames(M), n_cat = n_cat, x = x, expected_logistic = x - U, p_alt = p_alt, p_score = p_score, p_wallenius = wal$p_wallenius)
    out[[length(out) + 1]] <- d
  }
}
out <- rbindlist(out); fwrite(out, file.path(sc, "alt_tests.csv"))
cat("BH q<0.05 counts by minimum term size: logistic/Poisson-binomial test against the Wallenius test\n")
res <- rbindlist(lapply(c(5, 10, 15, 20, 50, 100), function(k) out[n_cat >= k, .(min_genes = k, terms = .N, alt_q05 = sum(p.adjust(p_alt, "BH") < 0.05), alt_min_q = round(min(p.adjust(p_alt, "BH")), 4),
   wall_q05 = sum(p.adjust(p_wallenius, "BH") < 0.05), wall_min_q = round(min(p.adjust(p_wallenius, "BH")), 4), top_alt = go_id[which.min(p_alt)], top_alt_n = n_cat[which.min(p_alt)]), by = .(set, ontology)]))
res[, top_term := substr(AnnotationDbi::Term(GO.db::GOTERM[top_alt]), 1, 40)]
print(res[order(set, ontology, min_genes)], nrows = 200)
cat("\nAgreement of the two tests (terms with >= 5 genes): Spearman correlation of p-values, and ratio of medians\n")
print(out[n_cat >= 5, .(terms = .N, spearman = round(cor(p_alt, p_wallenius, method = "spearman"), 3), share_alt_smaller = round(mean(p_alt < p_wallenius), 3), alt_p05 = round(mean(p_alt < 0.05), 4), wall_p05 = round(mean(p_wallenius < 0.05), 4),
              alt_p001 = round(mean(p_alt < 0.001), 5), wall_p001 = round(mean(p_wallenius < 0.001), 5)), by = .(set, ontology)])
cat("\norganelle (GO:0043226), heat, CC:\n"); print(out[set == "heat" & go_id == "GO:0043226"])

# Calibration of the logistic/Poisson-binomial alternative test under (i) the label-permutation null and (ii) the model null.
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
suppressPackageStartupMessages({library(splines); library(GenomicRanges)})
set.seed(1313); nperm <- as.integer(commandArgs(TRUE)[1]); nsim <- as.integer(commandArgs(TRUE)[2])
st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); uni <- st$uni; tm <- st$tm
ps <- readRDS(file.path(sc, "perm_sets_ploidy_legacy.rds")); CP <- readRDS(file.path(sc, "cpg_table_legacy.rds"))
gt <- fread(file.path(repo, "04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz"), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
h <- suppressWarnings(findOverlaps(GRanges(CP$chr, IRanges(CP$start + 1L, CP$start + 2L)), GRanges(gt$chr, IRanges(gt$start, gt$end)), ignore.strand = TRUE))
pairs <- data.table(cpg = queryHits(h), gene_id = gt$gene_id[subjectHits(h)]); setkey(pairs, cpg)
pb_tail <- function(p, x) { if (x <= 0) return(1); pmf <- 1; for (pi in p) pmf <- c(pmf * (1 - pi), 0) + c(0, pmf * pi); sum(pmf[(x + 1):length(pmf)]) }
prep <- lapply(tm, function(t) { M <- t$M[rowSums(t$M) >= 5, ]; Mt <- as(M, "TsparseMatrix"); list(genes = t$genes, M = M, idx = split(Mt@j + 1L, Mt@i + 1L), n_cat = as.vector(rowSums(M)),
                                                       X = suppressWarnings(cbind(1, ns(log2(uni$bias[match(t$genes, uni$gene_id)]), df = 5)))) })
alt_run <- function(de_all) rbindlist(lapply(names(prep), function(ont) {
  P <- prep[[ont]]; y <- de_all[P$genes]; f <- glm.fit(P$X, y, family = binomial()); ph <- f$fitted.values; w <- ph * (1 - ph)
  x <- as.vector(P$M %*% y); U <- x - as.vector(P$M %*% ph); MWX <- as.matrix(P$M %*% (w * P$X)); V <- as.vector(P$M %*% w) - rowSums((MWX %*% solve(crossprod(P$X, w * P$X))) * MWX)
  p <- pnorm(U / sqrt(pmax(V, 1e-12)), lower.tail = FALSE); small <- which(P$n_cat < 100)
  p[small] <- vapply(small, function(k) pb_tail(ph[P$idx[[as.character(k)]]], x[k]), numeric(1))
  pwf <- fit_pwf(de_all, uni$bias); wal <- wallenius_p(P$M, y, pwf[P$genes])$p_wallenius
  cnt <- function(pp, k) sum(p.adjust(pp[P$n_cat >= k], "BH") < 0.05)
  data.table(ontology = ont, de = sum(y), alt_5 = cnt(p, 5), alt_10 = cnt(p, 10), alt_50 = cnt(p, 50), alt_100 = cnt(p, 100), wal_5 = cnt(wal, 5), wal_10 = cnt(wal, 10), wal_15 = cnt(wal, 15), wal_50 = cnt(wal, 50), wal_100 = cnt(wal, 100),
             alt_p05_big = mean(p[P$n_cat >= 50] < 0.05), alt_p001_big = mean(p[P$n_cat >= 50] < 0.001), alt_p05_small = mean(p[P$n_cat < 50] < 0.05))
}))
ok <- which(lengths(ps$sets) > 0)[seq_len(nperm)]
rp <- rbindlist(lapply(ok, function(b) { de <- as.integer(uni$gene_id %in% unique(pairs[.(ps$sets[[b]]), nomatch = NULL]$gene_id)); names(de) <- uni$gene_id; cbind(perm = b, alt_run(de)) }))
fwrite(rp, file.path(sc, "alt_calibration_perm.csv"))
cat("Label-permutation null (ploidy), share of permutations with at least one BH q<0.05 term, by test and minimum term size:\n")
print(rp[, c(list(perms = .N), lapply(.SD, function(v) round(mean(v > 0), 3))), by = ontology, .SDcols = c("alt_5", "alt_10", "alt_50", "alt_100", "wal_5", "wal_10", "wal_15", "wal_50", "wal_100")])
print(rp[, .(alt_p05_big = round(mean(alt_p05_big), 4), alt_p001_big = round(mean(alt_p001_big), 5), alt_p05_small = round(mean(alt_p05_small), 4), max_alt_50 = max(alt_50), q95_alt_50 = quantile(alt_50, .95), q95_alt_100 = quantile(alt_100, .95), q95_wal_100 = quantile(wal_100, .95)), by = ontology])
pw <- lapply(c(ploidy = "ploidy", heat = "heat"), function(s) { de <- as.integer(uni[[paste0("n_dml_", s)]] > 0); names(de) <- uni$gene_id; fit_pwf(de, uni$bias) })
rm_ <- rbindlist(lapply(names(pw), function(s) rbindlist(lapply(seq_len(nsim), function(i) { de <- setNames(rbinom(nrow(uni), 1, pw[[s]]), uni$gene_id); cbind(set = s, sim = i, alt_run(de)) }))))
fwrite(rm_, file.path(sc, "alt_calibration_model.csv"))
cat("\nModel null (independent genes, fitted PWF as probability), share of runs with at least one BH q<0.05 term:\n")
print(rm_[, c(list(sims = .N), lapply(.SD, function(v) round(mean(v > 0), 3))), by = .(set, ontology), .SDcols = c("alt_5", "alt_10", "alt_50", "alt_100", "wal_5", "wal_10", "wal_15", "wal_50", "wal_100")])
print(rm_[, .(alt_p05_big = round(mean(alt_p05_big), 4), alt_p001_big = round(mean(alt_p001_big), 5), alt_p05_small = round(mean(alt_p05_small), 4)), by = .(set, ontology)])

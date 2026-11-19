suppressPackageStartupMessages(library(methylKit))
for (md in c("legacy", "corrected")) {
  x <- readRDS(sprintf("pq_%s.rds", md))
  for (ct in c("p", "h")) {
    p <- x[[paste0("p", ct)]]; q <- x[[paste0("q", ct)]]
    q1 <- methylKit:::p.adjusted(p, method = "SLIM"); q2 <- methylKit:::p.adjusted(p, method = "SLIM")
    cat(sprintf("[%s %s] SLIM from CSV p-values: max rel diff vs committed q = %.3g ; deterministic (two calls identical): %s ; pi0 = %.4f ; calls q<0.01 differing: %d\n", md, ct,
                max(abs(q1 - q) / pmax(q, 1e-300)), identical(q1, q2), methylKit:::SLIMfunc(p)$pi0_Est, sum((q1 < 0.01) != (q < 0.01))))
  }
}

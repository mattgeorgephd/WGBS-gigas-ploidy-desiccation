# (d) Calibration of goseq's Wallenius approximation under its own model, and power (positive control).
source("/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C4a-eligible-theory/lib_go.R")
args <- commandArgs(TRUE); what <- args[1]; nsim <- as.integer(args[2]); set.seed(as.integer(args[3]))
st <- readRDS(file.path(sc, "eligible_setup_legacy.rds")); uni <- st$uni; tm <- st$tm
pwf_real <- lapply(c(ploidy = "ploidy", heat = "heat"), function(s) { de <- as.integer(uni[[paste0("n_dml_", s)]] > 0); names(de) <- uni$gene_id; fit_pwf(de, uni$bias) })
# fit quality of the PWF: fitted against observed share by bias bin
for (s in names(pwf_real)) {
  d <- data.table(bias = uni$bias, de = uni[[paste0("n_dml_", s)]] > 0, pwf = pwf_real[[s]])
  d[, bin := cut(bias, c(0, 1, 2, 5, 10, 20, 50, 100, Inf))]
  cat(s, ": observed share of genes with a DML against mean fitted PWF, by number of eligible CpGs\n"); print(d[, .(genes = .N, observed = round(mean(de), 3), fitted = round(mean(pwf), 3)), by = bin][order(bin)])
  cat(sprintf("   sum(pwf) = %.1f, observed DML genes = %d; max pwf %.3f at bias %d; pwf at bias 743: %.3f\n", sum(d$pwf), sum(d$de), max(d$pwf), d$bias[which.max(d$pwf)], d[bias == max(bias), pwf]))
}
one <- function(de_all, ont, spiked = NULL) {
  pwf <- fit_pwf(de_all, uni$bias); g <- tm[[ont]]$genes; M <- tm[[ont]]$M; keep <- rowSums(M) >= 5; M <- M[keep, ]
  r <- wallenius_p(M, de_all[g], pwf[g]); r[, q := p.adjust(p_wallenius, "BH")]; r[, pf := fisher_p(M, de_all[g])]; r[, qf := p.adjust(pf, "BH")]
  r
}
if (what == "time") { t0 <- Sys.time(); de <- setNames(rbinom(nrow(uni), 1, pwf_real$heat), uni$gene_id); r <- one(de, "BP"); cat("one BP run:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n"); t0 <- Sys.time(); r <- one(de, "CC"); cat("one CC run:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n") }
if (what == "null") {
  out <- list()
  for (s in c("heat", "ploidy")) for (ont in c("BP", "CC")) for (i in seq_len(nsim)) {
    de <- setNames(rbinom(nrow(uni), 1, pwf_real[[s]]), uni$gene_id); r <- one(de, ont)
    cls <- function(k) c(n = sum(k), p05 = mean(r$p_wallenius[k] < 0.05), p01 = mean(r$p_wallenius[k] < 0.01), p001 = mean(r$p_wallenius[k] < 0.001))
    out[[length(out) + 1]] <- data.table(set = s, ontology = ont, sim = i, n_q05 = sum(r$q < 0.05), min_q = min(r$q), n_q05_fisher = sum(r$qf < 0.05),
      class = c("all (>=5)", "5-19", "20-99", "100-499", ">=500", ">=50 & weight>1.15", ">=50 & weight<0.87"),
      rbind(cls(rep(TRUE, nrow(r))), cls(r$n_cat < 20), cls(r$n_cat >= 20 & r$n_cat < 100), cls(r$n_cat >= 100 & r$n_cat < 500), cls(r$n_cat >= 500), cls(r$n_cat >= 50 & r$weight > 1.15), cls(r$n_cat >= 50 & r$weight < 0.87)))
  }
  out <- rbindlist(out); fwrite(out, file.path(sc, sprintf("calib_null_%s.csv", args[3])))
  cat("\nType I error of the Wallenius p-value when genes are drawn independently with the fitted PWF as their probability:\n")
  print(out[, .(sims = .N, terms = round(mean(n)), p05 = round(mean(p05, na.rm = TRUE), 4), p01 = round(mean(p01, na.rm = TRUE), 4), p001 = round(mean(p001, na.rm = TRUE), 5)), by = .(set, ontology, class)])
  cat("\nRuns with at least one BH q<0.05 term (Wallenius; Fisher for comparison):\n")
  print(unique(out[, .(set, ontology, sim, n_q05, min_q, n_q05_fisher)])[, .(sims = .N, any_wallenius = mean(n_q05 > 0), max_wallenius = max(n_q05), median_min_q = round(median(min_q), 3), q10_min_q = round(quantile(min_q, .1), 3), any_fisher = mean(n_q05_fisher > 0), median_fisher = median(n_q05_fisher)), by = .(set, ontology)])
}
if (what == "power") {
  ont <- "BP"; g <- tm[[ont]]$genes; M <- tm[[ont]]$M; sz <- rowSums(M)
  targets <- c(10, 30, 100, 300, 1000)
  out <- list()
  for (s in c("heat", "ploidy")) for (target in targets) {
    cand <- names(sz)[sz >= target * 0.9 & sz <= target * 1.1]
    for (or in c(1.5, 2, 3)) for (i in seq_len(nsim)) {
      term <- sample(cand, 1); in_term <- uni$gene_id %in% g[M[term, ] > 0]
      p <- pwf_real[[s]]; p[in_term] <- or * p[in_term] / (1 - p[in_term] + or * p[in_term])
      de <- setNames(rbinom(nrow(uni), 1, p), uni$gene_id); r <- one(de, ont)
      out[[length(out) + 1]] <- data.table(set = s, target = target, or = or, sim = i, term = term, size = sz[[term]], x = r[go_id == term, x], p = r[go_id == term, p_wallenius], q = r[go_id == term, q],
                                           any_q05 = sum(r$q < 0.05), mean_pwf_in_term = mean(pwf_real[[s]][in_term]))
    }
    cat(s, target, "done\n")
  }
  out <- rbindlist(out); fwrite(out, file.path(sc, sprintf("power_%s.csv", args[3])))
  cat("\nPower: share of runs in which the spiked BP term itself has BH q < 0.05 (eligible universe, Wallenius, min 5 genes), and in which any term has:\n")
  print(dcast(out[, .(power = round(mean(q < 0.05), 2)), by = .(set, target, or)], set + target ~ or, value.var = "power"))
  print(dcast(out[, .(any = round(mean(any_q05 > 0), 2)), by = .(set, target, or)], set + target ~ or, value.var = "any"))
}

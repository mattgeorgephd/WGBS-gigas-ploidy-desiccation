# 05_alt4_gap_rule.R
# Task (b)4: do other label-free rules pick the same five libraries, and does
# "five libraries stand apart" survive the 24 libraries of the control project?
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
options(width = 200)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S2-I49-evidence-and-alternatives"
bis <- fread(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"))
L <- setNames(100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG), bis$library_id)
der <- readRDS(file.path(scratch, "derived_legacy.rds")); me <- der$me
five <- c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15")

rules <- function(x, label) {
  cat("\n----", label, "(n =", length(x), ") ----\n")
  s <- sort(x); g <- diff(s)
  cat("largest gap:", round(max(g), 3), "between", names(s)[which.max(g)], "and", names(s)[which.max(g) + 1], "; second:", round(sort(g, decreasing = TRUE)[2], 3),
      "; third:", round(sort(g, decreasing = TRUE)[3], 3), "\n")
  cat("  above the largest gap:", names(s)[(which.max(g) + 1):length(s)], "\n")
  o <- order(-g)[2]; cat("  above the second largest gap:", names(s)[(o + 1):length(s)], "\n")
  # k-means, k = 2 and 3 (1-D; deterministic optimum by trying every split of the sorted values)
  best_split <- function(k) {
    n <- length(s)
    if (k == 2) { ss <- sapply(1:(n - 1), function(i) sum((s[1:i] - mean(s[1:i]))^2) + sum((s[(i + 1):n] - mean(s[(i + 1):n]))^2)); i <- which.min(ss); list(cut = i, ss = ss[i]) }
    else { best <- list(ss = Inf); for (i in 1:(n - 2)) for (j in (i + 1):(n - 1)) { v <- sum((s[1:i] - mean(s[1:i]))^2) + sum((s[(i + 1):j] - mean(s[(i + 1):j]))^2) + sum((s[(j + 1):n] - mean(s[(j + 1):n]))^2); if (v < best$ss) best <- list(ss = v, i = i, j = j) }; best }
  }
  b2 <- best_split(2); cat("k-means, k = 2 (exact optimum): upper cluster =", names(s)[(b2$cut + 1):length(s)], "\n")
  b3 <- best_split(3); cat("k-means, k = 3 (exact optimum): middle =", names(s)[(b3$i + 1):b3$j], "; upper =", names(s)[(b3$j + 1):length(s)], "\n")
  tot <- sum((s - mean(s))^2)
  cat("  share of variance between clusters: k = 2", round(1 - b2$ss / tot, 3), "; k = 3", round(1 - b3$ss / tot, 3), "\n")
  # top quartile
  q3 <- quantile(x, 0.75); cat("above the 75th percentile (", round(q3, 3), "):", names(x)[x > q3], "\n")
  # Tukey fence and MAD rule
  iqr <- IQR(x); fence <- q3 + 1.5 * iqr; cat("above the Tukey fence Q3 + 1.5 IQR (", round(fence, 3), "):", names(x)[x > fence], "\n")
  md <- median(x) + 3 * mad(x); cat("above median + 3 MAD (", round(md, 3), "):", names(x)[x > md], "\n")
  # Gaussian mixtures by EM, 1 to 3 components, unequal variances (variance floor to avoid degeneracy)
  em <- function(k, equal = FALSE, nstart = 200) {
    best <- NULL; n <- length(x); floorv <- (0.05 * sd(x))^2
    for (st in 1:nstart) {
      set.seed(st); mu <- sort(sample(x, k)); sg <- rep(var(x) / k, k); pi_ <- rep(1 / k, k); ll_old <- -Inf
      for (it in 1:500) {
        d <- sapply(1:k, function(j) pi_[j] * dnorm(x, mu[j], sqrt(sg[j])))
        if (k == 1) d <- matrix(d, ncol = 1)
        ll <- sum(log(rowSums(d))); r <- d / rowSums(d)
        nk <- colSums(r); pi_ <- nk / n; mu <- colSums(r * x) / nk
        sg <- sapply(1:k, function(j) max(sum(r[, j] * (x - mu[j])^2) / nk[j], floorv))
        if (equal) sg <- rep(sum(nk * sg) / n, k)
        if (abs(ll - ll_old) < 1e-9) break; ll_old <- ll
      }
      if (is.finite(ll) && (is.null(best) || ll > best$ll)) best <- list(ll = ll, mu = mu, sg = sg, pi = pi_, r = r)
    }
    npar <- if (equal) 2 * k else 3 * k - 1
    best$bic <- -2 * best$ll + npar * log(n); best
  }
  m1 <- em(1, nstart = 1); m2e <- em(2, TRUE); m2 <- em(2); m3e <- em(3, TRUE)
  cat("Gaussian mixture BIC (lower is better): 1 component", round(m1$bic, 1), "; 2 equal variance", round(m2e$bic, 1), "; 2 unequal variance", round(m2$bic, 1), "; 3 equal variance", round(m3e$bic, 1), "\n")
  up <- which.max(m2e$mu); cat("  2 components, equal variance: means", round(sort(m2e$mu), 2), "; libraries with posterior > 0.5 in the upper one:", names(x)[m2e$r[, up] > 0.5], "\n")
  up <- which.max(m2$mu); cat("  2 components, unequal variance: means", round(sort(m2$mu), 2), " sd", round(sqrt(m2$sg[order(m2$mu)]), 2), "; upper:", names(x)[m2$r[, up] > 0.5], "\n")
  o3 <- order(m3e$mu); cat("  3 components, equal variance: means", round(m3e$mu[o3], 2), "; middle:", names(x)[apply(m3e$r, 1, which.max) == o3[2]], "; upper:", names(x)[apply(m3e$r, 1, which.max) == o3[3]], "\n")
  invisible(NULL)
}
rules(L, "genome-wide level, 22 libraries")
rules(me, "mean over eligible CpGs, 22 libraries")
rules(L[grepl("zr3534", names(L))], "genome-wide level, heat cohort alone")
rules(L[grepl("zr3644", names(L))], "genome-wide level, control cohort alone (12)")

cat("\n==== the 24 libraries of project zr3644 (Bismark summary report) ====\n")
sm <- fread(file.path(repo, "03_input_files/bismark_reports/061021-big_bismark_summary_report.txt"), sep = "\t")
setnames(sm, gsub(" ", "_", names(sm)))
sm[, lib := sub("_R1_bismark_bt2_pe.bam", "", File)]
z <- sm[grepl("^zr3644_", lib)]
z[, level := 100 * Methylated_CpGs / (Methylated_CpGs + Unmethylated_CpGs)]
z[, chh := 100 * Methylated_CHHs / (Methylated_CHHs + Unmethylated_CHHs)]
z[, num := as.integer(sub("zr3644_", "", lib))]
z[, in_study := num %in% c(1:6, 13:18)]
setorder(z, level)
print(z[, .(lib, in_study, level = round(level, 3), chh = round(chh, 3), read_pairs_M = round(Total_Reads / 1e6, 1))])
cat("agreement with bismark_by_library.csv for the 12 libraries of the study: max abs difference", max(abs(z[in_study == TRUE][order(num), level] - L[paste0("zr3644_", c(1:6, 13:18))])), "\n")
x24 <- setNames(z$level, z$lib)
rules(x24, "genome-wide level, all 24 libraries of project zr3644")
cat("values of the 12 libraries of the project that are not in the study:", round(sort(z[in_study == FALSE, level]), 2), "\n")
cat("libraries of the project between 10.41 and 11.02:", z[level > 10.4105 & level < 11.018, paste(lib, round(level, 2))], "\n")
all34 <- c(L[grepl("zr3534", names(L))], x24)
rules(all34, "heat cohort (10) plus all 24 libraries of project zr3644")
# other projects in the same summary report, for scale (other tissues or species are possible; shown only as context)
sm[, level := 100 * Methylated_CpGs / (Methylated_CpGs + Unmethylated_CpGs)]
cat("\nother rows of the summary report (different projects; context only):\n")
print(sm[!grepl("^zr3644_|^zr3534_", lib), .(lib, level = round(level, 2))])

cat("\n==== how unusual is a largest gap of this size under a unimodal distribution with a tail? ====\n")
# Parametric bootstrap: fit a shifted log-normal and a normal to the 22 values, draw samples of 22, and record how
# often the largest gap is at least the observed one with at least 5 values above it (and with exactly k above).
x <- L; n <- length(x); obs_gap <- max(diff(sort(x)))
sim <- function(rdraw, B = 20000) {
  set.seed(42); hit <- 0; hit5 <- 0; ratio <- numeric(B)
  for (b in 1:B) { s <- sort(rdraw(n)); g <- diff(s); k <- which.max(g); nab <- n - k
    if (max(g) >= obs_gap) hit <- hit + 1
    if (max(g) >= obs_gap && nab >= 3) hit5 <- hit5 + 1
    ratio[b] <- max(g) / sort(g, decreasing = TRUE)[2] }
  c(p_gap = hit / B, p_gap_and_3_or_more_above = hit5 / B, p_ratio = mean(ratio >= obs_gap / sort(diff(sort(x)), decreasing = TRUE)[2]))
}
shift <- min(x) - 0.5; lx <- log(x - shift)
cat("normal(mean, sd of the 22):", round(sim(function(n) rnorm(n, mean(x), sd(x))), 4), "\n")
cat("shifted log-normal fitted to the 22:", round(sim(function(n) shift + exp(rnorm(n, mean(lx), sd(lx)))), 4), "\n")
cat("observed largest gap", round(obs_gap, 3), "; ratio to the second largest", round(obs_gap / sort(diff(sort(x)), decreasing = TRUE)[2], 2), "\n")
# skewness and a normality test of the 22 values
cat("Shapiro-Wilk on the 22 levels: p =", signif(shapiro.test(x)$p.value, 3), "; on log(level): p =", signif(shapiro.test(log(x))$p.value, 3), "\n")

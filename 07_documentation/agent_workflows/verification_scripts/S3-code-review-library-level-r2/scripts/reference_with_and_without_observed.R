# Reference statistics of part 6 with and without the observed grouping (and its mirror image), and under the rule without margin.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
a <- rbindlist(lapply(c("legacy","corrected"), function(m) fread(file.path(S, "out", sprintf("relabel_all_groupings_%s.csv", m)))))
a[, `:=`(n_exact = as.numeric(n_exact), n_S = as.numeric(n_S))]
options(width=250)
md <- function(x) as.numeric(median(x)); q <- function(x, p) as.numeric(quantile(x, p, names = FALSE))
print(a[, .(groupings = .N, median_all = md(n_exact), median_others = md(n_exact[!is_observed & !is_mirror]),
            p05_all = q(n_exact, 0.05), p05_others = q(n_exact[!is_observed & !is_mirror], 0.05),
            p95_all = q(n_exact, 0.95), p95_others = q(n_exact[!is_observed & !is_mirror], 0.95),
            observed = n_exact[is_observed], ratio_all = round(n_exact[is_observed]/md(n_exact), 5), ratio_others = round(n_exact[is_observed]/md(n_exact[!is_observed & !is_mirror]), 5),
            strictly_above = sum(n_exact > n_exact[is_observed]), equal_incl_obs = sum(n_exact == n_exact[is_observed])), by = .(mode, comparison = substr(comparison, 1, 58))])
print(a[, .(observed_exact = n_exact[is_observed], observed_S = n_S[is_observed], median_exact = md(n_exact), median_S = md(n_S), min_exact = min(n_exact), min_S = min(n_S),
            p05_exact = q(n_exact, .05), p05_S = q(n_S, .05), p95_exact = q(n_exact, .95), p95_S = q(n_S, .95),
            max_other_exact = max(n_exact[!is_observed & !is_mirror]), max_other_S = max(n_S[!is_observed & !is_mirror]),
            at_or_above_exact = sum(n_exact >= n_exact[is_observed]), at_or_above_S = sum(n_S >= n_S[is_observed]), ties_observed = n_ties[is_observed], max_ties = max(n_ties)),
        by = .(mode, comparison = substr(comparison, 1, 58))])

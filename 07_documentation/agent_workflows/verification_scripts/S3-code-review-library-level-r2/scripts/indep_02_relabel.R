# indep_02_relabel.R  <mode>
# Independent recomputation of part 6 of library_methylation_level.R (relabelling references) and
# of the comparison with dml_sensitivity/permutation_reference.csv (task 4).
# Different code: groupings enumerated from bit patterns, group sums by adding columns (no matrix
# product), every threshold decision by exact integer arithmetic on ALL CpGs (no screening).
# For each grouping three counts:
#   exact   5 * |cb*na - ca*nb| > na*nb            (|difference| > 20 as rational numbers)
#   S       |100 * (cb/nb - ca/na)| > 20            (floating point, no margin: dml_sensitivity.R and methylKit)
#   L       |100 * (cb/nb - ca/na)| > 20 + 1e-9     (floating point with the margin: library_methylation_level.R)
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE); mode <- args[1]
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
out <- file.path(S, "out"); committed <- file.path(repo, "06_diagnostics/library_methylation_level")
t0 <- proc.time()
lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE)
bis <- read.csv(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"), stringsAsFactors = FALSE)
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG)
o <- order(L); high <- seq_len(22) %in% o[(which.max(diff(L[o])) + 1):22]
u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
n_cpg <- nrow(u)
Cl <- lapply(1:22, function(j) as.numeric(u[[paste0("numCs", j)]]))
Nl <- lapply(1:22, function(j) Cl[[j]] + u[[paste0("numTs", j)]])
PM <- lapply(1:22, function(j) 100 * Cl[[j]] / Nl[[j]])
w <- function(cond) which(cond)
dh <- w(lib$ploidy == "diploid" & lib$treatment == "heat"); th <- w(lib$ploidy == "triploid" & lib$treatment == "heat")
dc <- w(lib$ploidy == "diploid" & lib$treatment == "control"); tc <- w(lib$ploidy == "triploid" & lib$treatment == "control")
ac <- w(lib$acute_heat_shock == "yes"); hi <- w(high)
comparisons <- list(
  "ploidy inside the heat cohort, all 10 libraries" = list(list(a = dh, b = th)),
  "ploidy inside the control cohort, all 12 libraries" = list(list(a = dc, b = tc)),
  "ploidy inside the heat cohort, without the high libraries" = list(list(a = setdiff(dh, hi), b = setdiff(th, hi))),
  "ploidy inside the control cohort, without the high libraries" = list(list(a = setdiff(dc, hi), b = setdiff(tc, hi))),
  "ploidy, both cohorts, without the high libraries" = list(list(a = setdiff(dh, hi), b = setdiff(th, hi)), list(a = setdiff(dc, hi), b = setdiff(tc, hi))),
  "acute heat shock inside the heat cohort" = list(list(a = intersect(dh, ac), b = setdiff(dh, ac)), list(a = intersect(th, ac), b = setdiff(th, ac))))
popcount <- function(x) { n <- 0L; while (any(x > 0L)) { n <- n + bitwAnd(x, 1L); x <- bitwShiftR(x, 1L) }; n }
subsets_of_stratum <- function(members, k) {            # every subset of size k, from bit patterns
  members <- sort(members); n <- length(members); masks <- 0:(2^n - 1)
  masks <- masks[popcount(masks) == k]
  lapply(masks, function(mk) members[bitwAnd(mk, bitwShiftL(1L, 0:(n - 1))) != 0L])
}
add <- function(lst, idx) Reduce(`+`, lst[idx])
keep_detail <- mode == "legacy"
all_rows <- list(); detail <- list(); screen_rows <- list()
for (nm in names(comparisons)) {
  strata <- comparisons[[nm]]
  libs <- sort(unlist(lapply(strata, function(s) c(s$a, s$b))))
  per <- lapply(strata, function(s) subsets_of_stratum(c(s$a, s$b), length(s$a)))
  grid <- expand.grid(lapply(per, seq_along))
  groupings <- lapply(seq_len(nrow(grid)), function(i) sort(unlist(lapply(seq_along(per), function(k) per[[k]][[grid[i, k]]]))))
  a_obs <- sort(unlist(lapply(strata, `[[`, "a"))); b_obs <- sort(unlist(lapply(strata, `[[`, "b")))
  ctot <- add(Cl, libs); ntot <- add(Nl, libs)
  # the screen of the script: range of the library values of the comparison no more than 20 points
  span <- do.call(pmax, PM[libs]) - do.call(pmin, PM[libs]); dropped <- !(span > 20)
  # the same screen decided exactly: largest and smallest library value as fractions
  imax <- max.col(do.call(cbind, PM[libs]), ties.method = "first"); imin <- max.col(-do.call(cbind, PM[libs]), ties.method = "first")
  Cmat <- do.call(cbind, Cl[libs]); Nmat <- do.call(cbind, Nl[libs]); ii <- seq_len(n_cpg)
  cmx <- Cmat[cbind(ii, imax)]; nmx <- Nmat[cbind(ii, imax)]; cmn <- Cmat[cbind(ii, imin)]; nmn <- Nmat[cbind(ii, imin)]
  span_exact_above <- 5 * (cmx * nmn - cmn * nmx) > nmx * nmn
  rm(Cmat, Nmat)
  res <- matrix(0, length(groupings), 7, dimnames = list(NULL, c("n_exact", "n_ties", "n_S", "n_L", "n_dropped_exact", "n_dropped_S", "max_product")))
  for (g in seq_along(groupings)) {
    A <- groupings[[g]]
    ca <- add(Cl, A); na <- add(Nl, A); cb <- ctot - ca; nb <- ntot - na
    x <- cb * na - ca * nb; lhs <- 5 * abs(x); rhs <- na * nb
    d <- 100 * (cb / nb - ca / na)
    ex <- lhs > rhs; tie <- lhs == rhs; sS <- abs(d) > 20; sL <- abs(d) > 20 + 1e-9
    res[g, ] <- c(sum(ex), sum(tie), sum(sS), sum(sL), sum(ex & dropped), sum(sS & dropped), max(rhs))
    if (keep_detail && nm %in% names(comparisons)[1:2] && any(sS != ex)) {
      k <- which(sS != ex)
      detail[[length(detail) + 1]] <- data.table(comparison = nm, grouping = g, libraries_in_group_a = paste(lib$library_id[A], collapse = " "),
        observed_or_mirror = identical(A, a_obs) || identical(A, b_obs), chr = u$chr[k], start = u$start[k],
        c_a = ca[k], n_a = na[k], c_b = cb[k], n_b = nb[k], exact_tie = tie[k], d_float = sprintf("%.17g", d[k]), span_float = sprintf("%.17g", span[k]), dropped_by_screen = dropped[k])
    }
  }
  stopifnot(max(res[, "max_product"]) < 2^53 / 5)
  is_obs <- vapply(groupings, function(a) identical(a, a_obs), logical(1)); is_mir <- vapply(groupings, function(a) identical(a, b_obs), logical(1))
  stopifnot(sum(is_obs) == 1)
  all_rows[[nm]] <- data.table(mode = mode, comparison = nm, grouping = seq_along(groupings),
    libraries_in_group_a = vapply(groupings, function(a) paste(lib$library_id[a], collapse = " "), character(1)), is_observed = is_obs, is_mirror = is_mir, as.data.table(res))
  screen_rows[[nm]] <- data.table(mode = mode, comparison = nm, n_libraries = length(libs), groupings = length(groupings),
    cpgs_kept_by_screen = sum(!dropped), cpgs_kept_by_exact_screen = sum(span_exact_above), screen_float_differs_from_exact = sum(dropped == span_exact_above),
    cpgs_with_span_exactly_20 = sum(5 * (cmx * nmn - cmn * nmx) == nmx * nmn),
    largest_count_of_dropped_cpgs_above_exactly = max(res[, "n_dropped_exact"]), largest_count_of_dropped_cpgs_above_without_margin = max(res[, "n_dropped_S"]),
    groupings_where_L_differs_from_exact = sum(res[, "n_L"] != res[, "n_exact"]), groupings_where_S_differs_from_exact = sum(res[, "n_S"] != res[, "n_exact"]),
    largest_S_minus_exact = max(res[, "n_S"] - res[, "n_exact"]), smallest_S_minus_exact = min(res[, "n_S"] - res[, "n_exact"]),
    largest_product_of_pooled_reads = max(res[, "max_product"]), smallest_possible_distance_from_20 = 20 / max(res[, "max_product"]))
  cat(sprintf("%s: %d groupings done at %.0f s\n", nm, length(groupings), (proc.time() - t0)[3]))
}
allg <- rbindlist(all_rows); fwrite(allg, file.path(out, sprintf("relabel_all_groupings_%s.csv", mode)))
scr <- rbindlist(screen_rows); fwrite(scr, file.path(out, sprintf("relabel_screen_and_margin_%s.csv", mode)))
if (length(detail)) fwrite(rbindlist(detail), file.path(out, sprintf("relabel_cpgs_where_no_margin_differs_%s.csv", mode)))
options(width = 250); print(scr)

# ---- summary rows as the script defines them, from the exact counts -------------------------------
ref_committed <- fread(file.path(committed, "relabelling_reference.csv"))[mode == get("mode", envir = globalenv())]
summ <- function(cnt, is_obs, is_mir) {
  obs <- cnt[is_obs]
  data.table(n_above_threshold_observed = obs, groupings = length(cnt), reference_median = median(cnt), reference_p05 = unname(quantile(cnt, 0.05)),
             reference_p95 = unname(quantile(cnt, 0.95)), reference_min = min(cnt), largest_other_grouping = max(cnt[!is_obs & !is_mir]),
             observed_over_reference_median = obs / median(cnt), groupings_at_or_above_observed = sum(cnt >= obs), p_enumeration = sum(cnt >= obs) / length(cnt),
             n_mirror = sum(is_mir), reference_max = max(cnt), second_largest_distinct = sort(unique(cnt), decreasing = TRUE)[2])
}
mine <- allg[, summ(n_exact, is_observed, is_mirror), by = comparison]
mm <- merge(mine, ref_committed, by = "comparison", suffixes = c("", ".committed"), sort = FALSE)
stopifnot(nrow(mm) == 6)
cat("\nrelabelling_reference.csv (", mode, ") against my exact counts:\n", sep = "")
for (col in c("n_above_threshold_observed", "groupings", "reference_median", "reference_p05", "reference_p95", "reference_min", "largest_other_grouping",
              "observed_over_reference_median", "groupings_at_or_above_observed", "p_enumeration")) {
  a <- signif(as.numeric(mm[[col]]), 6); b <- as.numeric(mm[[paste0(col, ".committed")]])
  cat(sprintf("  %-34s %s\n", col, if (isTRUE(all.equal(a, b, tolerance = 1e-12))) "all 6 rows equal" else paste("DIFFERENT:", paste(a, collapse = " "), "against", paste(b, collapse = " "))))
}
print(mine)
# group sizes
sizes <- rbindlist(lapply(names(comparisons), function(nm) data.table(comparison = nm, n_a = sum(sapply(comparisons[[nm]], function(s) length(s$a))), n_b = sum(sapply(comparisons[[nm]], function(s) length(s$b))))))
ms <- merge(sizes, ref_committed[, .(comparison, n_a_c = n_a, n_b_c = n_b)], by = "comparison")
cat("  n_a and n_b:", if (all(ms$n_a == ms$n_a_c, ms$n_b == ms$n_b_c)) "all 6 rows equal" else "DIFFERENT", "\n")

# ---- every grouping of the small sets ---------------------------------------------------------------
small_committed <- fread(file.path(committed, "relabelling_small_sets.csv"))[mode == get("mode", envir = globalenv())]
sm <- merge(small_committed, allg[, .(comparison, libraries_in_group_a, n_exact, is_observed, is_mirror)], by = c("comparison", "libraries_in_group_a"), all.x = TRUE, sort = FALSE)
cat("\nrelabelling_small_sets.csv (", mode, "): ", nrow(small_committed), " rows; matched to my groupings: ", sum(!is.na(sm$n_exact)),
    "; counts equal: ", sum(sm$n_above_threshold == sm$n_exact), "; observed flag equal: ", sum(sm$observed_grouping == sm$is_observed),
    "; mirror flag equal: ", sum(sm$observed_grouping_with_labels_exchanged == sm$is_mirror), "\n", sep = "")
cat("  rows per comparison in the committed table:\n"); print(small_committed[, .N, by = comparison])
cat("  sorted by decreasing count within comparison in the committed table:", all(small_committed[, !is.unsorted(rev(n_above_threshold)), by = comparison]$V1), "\n")

# ---- task 4: the counts without margin against dml_sensitivity/permutation_reference.csv -----------
if (mode == "legacy") {
  pr <- fread(file.path(repo, "06_diagnostics/dml_sensitivity/permutation_reference.csv"))[mode == "legacy" & libraries %in% c("heat_only", "control_only")]
  for (nm in names(comparisons)[1:2]) {
    x <- allg[comparison == nm]; cnt <- x$n_S; obs <- cnt[x$is_observed]
    q <- quantile(cnt, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), names = FALSE)
    row <- pr[libraries == (if (grepl("heat", nm)) "heat_only" else "control_only")]
    theirs <- unlist(row[, .(n_above_cutoff_min, n_above_cutoff_q05, n_above_cutoff_q25, n_above_cutoff_median, n_above_cutoff_q75, n_above_cutoff_q95, n_above_cutoff_max)])
    cat("\n", nm, "\n  counts without margin (S): observed ", obs, " (permutation_reference.csv: ", row$n_effect_above_cutoff, "); quantiles ", paste(q, collapse = " "),
        "\n  permutation_reference.csv:                 quantiles ", paste(theirs, collapse = " "),
        "\n  identical: ", isTRUE(all.equal(q, unname(theirs))) && obs == row$n_effect_above_cutoff, "; at least observed: ", sum(cnt >= obs), " (", row$n_labellings_at_least_observed, ")\n", sep = "")
    qe <- quantile(x$n_exact, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), names = FALSE)
    cat("  exact counts:                              observed ", x$n_exact[x$is_observed], "; quantiles ", paste(qe, collapse = " "), "; at least observed: ", sum(x$n_exact >= x$n_exact[x$is_observed]), "\n", sep = "")
    cat("  groupings in which S counts more than the exact rule: ", sum(x$n_S > x$n_exact), " of ", nrow(x), "; by 1: ", sum(x$n_S - x$n_exact == 1), "; by 2: ", sum(x$n_S - x$n_exact == 2),
        "; by 3 or more: ", sum(x$n_S - x$n_exact >= 3), "; fewer: ", sum(x$n_S < x$n_exact), "\n", sep = "")
    cat("  exact ties per grouping: median ", median(x$n_ties), ", range ", min(x$n_ties), " to ", max(x$n_ties), "; share of ties that the rule without margin counts: ",
        round(sum(x$n_S - x$n_exact) / sum(x$n_ties), 4), "\n", sep = "")
  }
  dd <- rbindlist(detail)
  cat("\nCpG by grouping cases where the rule without margin differs from the exact rule:", nrow(dd), "; all exact ties:", all(dd$exact_tie),
      "; distinct CpGs:", uniqueN(dd[, .(chr, start)]), "; dropped by the screen of the script:", sum(dd$dropped_by_screen), "\n")
  print(dd[, .N, by = .(comparison, d_float)][order(comparison, -N)])
}
cat("\nelapsed seconds:", round((proc.time() - t0)[3]), "\n"); print(gc())

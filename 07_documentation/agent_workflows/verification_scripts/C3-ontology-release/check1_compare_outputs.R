# Check 1: committed step-09 outputs against the nine archived tables.
# Read-only on the repository. Usage: Rscript check1_compare_outputs.R
suppressPackageStartupMessages(library(data.table))
scr  <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C3-ontology-release"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(scr, "parse_archive.R"))
sets <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")
res <- list()
for (s in names(sets)) for (o in c("BP", "CC", "MF")) {
  fa <- file.path(repo, "08_archive/GOterms", paste0("BOTH-", sets[[s]]), sprintf("%s-%s-FisherTestResults.csv", sets[[s]], o))
  fn <- file.path(repo, "04_output/09_go_enrichment_legacy/tables", sprintf("legacy_fisher_%s_%s.csv", s, o))
  a <- as.data.table(parse_archive(fa)); n <- fread(fn)
  stopifnot(!anyDuplicated(a$go_id), !anyDuplicated(n$go_id))
  m <- merge(a, n, by = "go_id", suffixes = c("_a", "_n"))
  same_ids <- setequal(a$go_id, n$go_id)
  # order
  same_order <- identical(a$go_id, n$go_id)
  n_pos_diff <- if (nrow(a) == nrow(n)) sum(a$go_id != n$go_id) else NA
  # p-values: as printed by GenTable for the new vector, in the new row order
  n[, printed := gentable_p(p_fisher)]
  mm <- merge(a[, .(go_id, classic, rank_a = rank)], n[, .(go_id, printed, p_fisher, rank_n = rank)], by = "go_id")
  norm <- function(x) gsub(" ", "", x)
  p_string_same <- sum(norm(mm$classic) == norm(mm$printed))
  pa <- classic_numeric(mm$classic)
  parseable <- !is.na(pa)
  # relative difference where parseable: archive is printed to 2 significant digits
  rel <- abs(pa[parseable] - mm$p_fisher[parseable]) / pmax(pa[parseable], .Machine$double.xmin)
  # a 2-significant-digit print can be off by up to 5% of the leading-digit unit; test by re-rounding
  resig <- signif(mm$p_fisher[parseable], 2); resig3 <- signif(mm$p_fisher[parseable], 3)
  ok2 <- abs(resig - pa[parseable]) <= 1e-12 * pmax(pa[parseable], 1e-300) | abs(resig3 - pa[parseable]) <= 1e-12 * pmax(pa[parseable], 1e-300)
  lt <- grepl("<", mm$classic)
  res[[length(res) + 1]] <- data.table(
    set = s, ont = o, terms_archive = nrow(a), terms_new = nrow(n), same_id_set = same_ids,
    only_archive = length(setdiff(a$go_id, n$go_id)), only_new = length(setdiff(n$go_id, a$go_id)),
    annotated_same = sum(m$annotated_a == m$annotated_n), significant_same = sum(m$significant_a == m$significant_n),
    both_counts_same = sum(m$annotated_a == m$annotated_n & m$significant_a == m$significant_n),
    expected_same = sum(abs(m$expected_num - m$expected_n) < 1e-9),
    term_same = sum(m$term_a == short_term(m$term_n)),
    same_row_order = same_order, rows_at_other_position = n_pos_diff,
    p_parseable = sum(parseable), p_lt_1e30_archive = sum(lt),
    p_lt_1e30_new = sum(mm$p_fisher < 1e-30), p_lt_agree = sum(lt & mm$p_fisher < 1e-30),
    p_string_identical = p_string_same,
    p_resignif_equal = sum(ok2), p_max_rel_diff = if (any(parseable)) max(rel) else NA_real_,
    sum_annot_archive = sum(a$annotated), sum_sig_archive = sum(a$significant),
    top10_same_order = identical(a$go_id[1:10], n$go_id[1:10]),
    n_p_lt_001_archive = sum(pa < 0.01, na.rm = TRUE) + sum(lt), n_p_lt_001_new = sum(n$p_fisher < 0.01))
  if (p_string_same < nrow(mm)) {
    cat(sprintf("\n-- %s %s: first printed-p mismatches --\n", s, o))
    print(head(mm[norm(classic) != norm(printed)], 8))
  }
  if (sum(m$term_a == short_term(m$term_n)) < nrow(m)) {
    cat(sprintf("\n-- %s %s: first term-name mismatches --\n", s, o))
    print(head(m[term_a != short_term(term_n), .(go_id, term_a, term_n)], 8))
  }
}
res <- rbindlist(res)
options(width = 250)
print(res[, 1:14])
print(res[, c(1, 2, 15:ncol(res)), with = FALSE])
fwrite(res, file.path(scr, "check1_results.csv"))
cat("\nheader of archive files:", attr(parse_archive(fa), "header"), "\n")

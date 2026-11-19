# The script's own relabel_reference() (taken from the parsed file) on the real count tables: every grouping of the six
# comparisons against my exact counts (relabel_all_groupings_<mode>.csv), grouping by grouping.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
ex <- parse(file.path(repo, "06_diagnostics/library_methylation_level.R"), keep.source = FALSE); env <- new.env()
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% c("above", "relabel_reference")) eval(e, env)
assign("cutoff", 20, env); assign("tie_tolerance", 1e-9, env); assign("n_lib", 22L, env)
lib <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE)
bis <- read.csv(file.path(repo, "04_output/02_alignment_qc/tables/bismark_by_library.csv"), stringsAsFactors = FALSE)
L <- 100 * bis$meth_CpG / (bis$meth_CpG + bis$unmeth_CpG); o <- order(L); hi <- which(seq_len(22) %in% o[(which.max(diff(L[o])) + 1):22])
w <- function(cond) which(cond)
dh <- w(lib$ploidy == "diploid" & lib$treatment == "heat"); th <- w(lib$ploidy == "triploid" & lib$treatment == "heat")
dc <- w(lib$ploidy == "diploid" & lib$treatment == "control"); tc <- w(lib$ploidy == "triploid" & lib$treatment == "control"); ac <- w(lib$acute_heat_shock == "yes")
comparisons <- list(
  "ploidy inside the heat cohort, all 10 libraries" = list(list(a = dh, b = th)),
  "ploidy inside the control cohort, all 12 libraries" = list(list(a = dc, b = tc)),
  "ploidy inside the heat cohort, without the high libraries" = list(list(a = setdiff(dh, hi), b = setdiff(th, hi))),
  "ploidy inside the control cohort, without the high libraries" = list(list(a = setdiff(dc, hi), b = setdiff(tc, hi))),
  "ploidy, both cohorts, without the high libraries" = list(list(a = setdiff(dh, hi), b = setdiff(th, hi)), list(a = setdiff(dc, hi), b = setdiff(tc, hi))),
  "acute heat shock inside the heat cohort" = list(list(a = intersect(dh, ac), b = setdiff(dh, ac)), list(a = intersect(th, ac), b = setdiff(th, ac))))
for (mode in c("legacy", "corrected")) {
  u <- fread(file.path(repo, sprintf("04_output/05_differential_methylation/tables/united_cpg_%s.csv.gz", mode)))
  C <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]); N <- C + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); dimnames(C) <- dimnames(N) <- NULL
  m <- list(C = C, N = N); rm(u)
  mine <- fread(file.path(S, "out", sprintf("relabel_all_groupings_%s.csv", mode)))
  for (nm in names(comparisons)) {
    r <- env$relabel_reference(m, comparisons[[nm]])
    theirs <- data.table(libraries_in_group_a = vapply(r$group_a, function(a) paste(lib$library_id[a], collapse = " "), character(1)), n_script = r$counts, obs = r$is_observed, mir = r$is_mirror)
    x <- merge(theirs, mine[comparison == nm], by = "libraries_in_group_a")
    cat(sprintf("%s, %s: %d groupings by the function, %d matched to mine; count equal to the exact count in %d; equal to the count without margin in %d; observed and mirror flags equal: %s; CpGs kept by the screen %d\n",
                mode, nm, length(r$counts), nrow(x), sum(x$n_script == x$n_exact), sum(x$n_script == x$n_S), all(x$obs == x$is_observed, x$mir == x$is_mirror), r$n_cpg_screened))
  }
  rm(m, C, N); invisible(gc())
}

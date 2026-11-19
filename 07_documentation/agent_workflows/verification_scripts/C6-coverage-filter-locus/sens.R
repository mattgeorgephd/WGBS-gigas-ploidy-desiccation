suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
t5 <- function(f) file.path(repo, "04_output/05_differential_methylation/tables", f)
z <- readRDS("cand_matrix.rds"); cand <- z$cand; cov <- z$cov
cut <- fread("cutoffs_mine.csv"); ord <- colnames(cov)
cm <- function(col) cut[match(ord, library_id)][[col]]
over <- function(cutoffs, strict = FALSE) { m <- matrix(cutoffs, nrow(cov), 22, byrow = TRUE); rowSums(if (strict) cov > m else cov >= m) }
dml <- list(ploidy = fread(t5("dml_ploidy_legacy.csv")), heat = fread(t5("dml_heat_legacy.csv")), shared = fread(t5("dml_shared_legacy.csv")))
key <- paste(cand$chr, cand$start)
L0 <- 39633323L; L1 <- 39635520L; loc <- cand$chr == "NC_047563.1" & cand$start >= L0 & cand$start <= L1
row <- function(label, keep) {
  kk <- key[keep]
  data.table(rule = label, cpg_tested = sum(keep), removed = sum(!keep), pct_removed = round(100 * mean(!keep), 2),
             ploidy_dml = sum(paste(dml$ploidy$chr, dml$ploidy$start) %in% kk), heat_dml = sum(paste(dml$heat$chr, dml$heat$start) %in% kk),
             shared = sum(paste(dml$shared$chr, dml$shared$start) %in% kk), locus_cpg_kept = sum(keep & loc))
}
n2 <- over(cm("cutoff_mine"))
res <- rbind(
  row("no upper filter (legacy)", rep(TRUE, nrow(cov))),
  row("99.9th pct, base >= 2 reads, remove >= cutoff in any library (pipeline corrected)", n2 == 0),
  row("same, remove only > cutoff", over(cm("cutoff_mine"), strict = TRUE) == 0),
  row("99.9th pct, base >= 1 read (mincov = 1 import)", over(cm("cutoff_base1")) == 0),
  row("99.9th pct, base >= 10 reads (methRead default mincov = 10)", over(cm("cutoff_base10")) == 0),
  row("base >= 2, remove only if >= cutoff in at least 2 libraries", n2 < 2),
  row("base >= 2, remove only if >= cutoff in more than half (12+) of libraries", n2 < 12),
  row("base >= 2, remove only if >= cutoff in all 22 libraries", n2 < 22))
options(width = 250); print(res)
cat("\nNote: DML counts for the alternative rules are legacy DML at retained CpGs; valid because p-values and differences do not depend on the CpG set and the q-value criterion is not binding (see qbind.R).\n")
cat("check against the pipeline's corrected lists:", res[2, ploidy_dml], res[2, heat_dml], res[2, shared], "\n")
fwrite(res, "filter_rule_sensitivity.csv")

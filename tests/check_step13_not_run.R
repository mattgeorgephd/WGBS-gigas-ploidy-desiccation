# tests/check_step13_not_run.R
#
# T4 of the settings layout (settings_build_spec.md, section 8): step 13 reports
# the values of optional analyses as "not run" when their outputs are absent,
# and never stops on them; a missing table of the paper still stops the paper part.
# Run in a SCRATCH COPY of the repository (it deletes files), never in the
# working repository:
#   cp -a <repository> <scratch> && cd <scratch> && CT_TEST_SCRATCH=1 Rscript tests/check_step13_not_run.R
# Needs the folders 04_output/not_adopted/corrected/ and legacy/ as committed (with corrected/10_go_enrichment_gene_level/).

suppressPackageStartupMessages({ library(here); library(data.table) })
if (file.exists(file.path(here(), ".git")) && !nzchar(Sys.getenv("CT_TEST_SCRATCH"))) {
  stop("Run this test in a scratch copy of the repository (cp -a), and set CT_TEST_SCRATCH=1 to confirm.", call. = FALSE)
}
rscript <- file.path(R.home("bin"), "Rscript")
results <- list()
ok <- function(name, pass) { results[[length(results) + 1]] <<- isTRUE(pass); cat(sprintf("%-4s %s\n", if (isTRUE(pass)) "ok" else "FAIL", name)) }
run13 <- function(env = character(0)) system2(rscript, "01_pipeline/13_manuscript_values.R", stdout = FALSE, stderr = FALSE, env = env)
values <- function() fread(here("05_manuscript", "not_adopted", "manuscript_values_by_analysis.csv"))
map <- as.data.table(utils::read.csv(here("03_input_files", "manuscript", "manuscript_values_map.csv"), stringsAsFactors = FALSE,
                                     na.strings = "", colClasses = "character", encoding = "UTF-8"))
with_expr <- map[!is.na(expression_optional), id]
reads_10 <- map[!is.na(expression_optional) & grepl("10_go_enrichment_gene_level", expression_optional), id]

# 0. The state after the migration: only the folders corrected/ and legacy/ (other analysis folders are removed here).
for (d in setdiff(list.dirs(here("04_output", "not_adopted"), recursive = FALSE, full.names = FALSE),
                  c("corrected", "legacy", "comparisons", "ct_variants", "13_manuscript_values"))) {
  unlink(here("04_output", "not_adopted", d), recursive = TRUE)
}
ok("step 13, optional part, exits with status 0", run13("WGBS_ANALYSIS=corrected") == 0)
v0 <- values()
for (A in c("cov5", "mnchisq", "ctexcl")) {
  ok(sprintf("%s: %d rows with an expression, all 'not run: analysis not run'", A, length(with_expr)),
     v0[analysis == A & id %in% with_expr, all(status == "not run: analysis not run") && .N == length(with_expr)])
}

# 1. Step 10 of 'corrected' removed: exactly the rows that read step 10 become not run.
unlink(here("04_output", "not_adopted", "corrected", "10_go_enrichment_gene_level"), recursive = TRUE)
ok("without corrected/10_go_enrichment_gene_level/, step 13 exits with status 0", run13("WGBS_ANALYSIS=corrected") == 0)
v1 <- values()
a0 <- v0[analysis == "corrected" & id %in% with_expr]; a1 <- v1[analysis == "corrected" & id %in% with_expr]
changed <- a1$id[a1$status != a0$status | a1$regenerated %in% NA != a0$regenerated %in% NA]
# Of the rows whose expression names step 10, those that read it become not run; two rows ("a || b") do not
# reach their step 10 part when the first part is TRUE, and keep their value.
ok(sprintf("of the %d rows whose expression names step 10, the %d changed ones are 'not run: 10_go_enrichment_gene_level/<file> absent'",
           length(reads_10), length(changed)),
   all(changed %in% reads_10) && a1[id %in% changed, all(startsWith(status, "not run: 10_go_enrichment_gene_level/") & endsWith(status, " absent"))])
cat(sprintf("     rows naming step 10 that keep their value: %s\n", paste(setdiff(reads_10, changed), collapse = ", ")))
ok(sprintf("the other %d rows are unchanged", length(with_expr) - length(changed)),
   isTRUE(all.equal(a0[!id %in% changed], a1[!id %in% changed], check.attributes = FALSE)))

# 2. A stamp that describes another analysis: every row of that id is 'not run: folder stale'.
s <- fread(here("04_output", "not_adopted", "corrected", "analysis.csv"), colClasses = "character")
s$min_coverage <- "5"; fwrite(s, here("04_output", "not_adopted", "corrected", "analysis.csv"))
ok("with a stale stamp, step 13 exits with status 0", run13("WGBS_ANALYSIS=corrected") == 0)
ok("every row of 'corrected' is 'not run: folder stale'", values()[analysis == "corrected", all(status == "not run: folder stale")])

# 3. A table of the paper removed: the paper part stops with an error.
unlink(here("04_output", "05_differential_methylation", "tables", "dml_counts_by_mode.csv"))
ok("a missing paper table stops the paper part (exit status not 0)", run13() != 0)
cat(sprintf("\nT4: %d of %d checks pass\n", sum(unlist(results)), length(results)))
quit(status = if (all(unlist(results))) 0 else 1)

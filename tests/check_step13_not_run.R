# tests/check_step13_not_run.R
#
# T4 of the settings layout (settings_build_spec.md, section 8): step 13 reports
# the values of optional analyses as "not run" when their outputs are absent,
# and never stops on them; a missing table of the paper still stops the paper part.
# Run it in a SCRATCH COPY of the repository WITHOUT .git (it deletes files),
# never in the original checkout; tests/scratch_guard.R refuses to run it there:
#   mkdir <scratch> && tar --exclude=./.git --exclude=./04_output/not_adopted/ct_variants -cf - . | tar -xf - -C <scratch>
#     (from the repository root; or cp -a, then delete <scratch>/.git and <scratch>/04_output/not_adopted/ct_variants/)
#   cd <scratch> && CT_TEST_SCRATCH=1 Rscript tests/check_step13_not_run.R
# Needs the folders 04_output/not_adopted/corrected/ and legacy/ as committed (with corrected/10_go_enrichment_gene_level/).

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("tests", "scratch_guard.R")); scratch_guard("tests/check_step13_not_run.R")
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
n_reads_10 <- 53     # rows of 'corrected' that read a table of step 10 (2026-11-20; two more name step 10 and do not reach it)
ok(sprintf("of the %d rows whose expression names step 10, %d (expected %d) changed, all 'not run: 10_go_enrichment_gene_level/<file> absent'",
           length(reads_10), length(changed), n_reads_10),
   length(changed) == n_reads_10 && all(changed %in% reads_10) &&
     a1[id %in% changed, all(startsWith(status, "not run: 10_go_enrichment_gene_level/") & endsWith(status, " absent"))])
cat(sprintf("     rows naming step 10 that keep their value: %s\n", paste(setdiff(reads_10, changed), collapse = ", ")))
ok(sprintf("the other %d rows are unchanged", length(with_expr) - length(changed)),
   isTRUE(all.equal(a0[!id %in% changed], a1[!id %in% changed], check.attributes = FALSE)))

# 1b. A step that the stamp does not list: the rows that read its tables become 'not run: <step> not done (stamp)',
#     although the tables are on disk (a table left by an earlier run counts only while the stamp lists its step).
stamp_path <- here("04_output", "not_adopted", "corrected", "analysis.csv")
s <- fread(stamp_path, colClasses = "character")
s1 <- copy(s); s1$steps_done <- paste(setdiff(strsplit(s$steps_done, ";")[[1]], "08"), collapse = ";"); fwrite(s1, stamp_path)
reads_08 <- map[!is.na(expression_optional) & grepl('(^|[^_])tab\\("08_gene_annotation"', expression_optional), id]
ok("without step 08 in the stamp, step 13 exits with status 0", run13("WGBS_ANALYSIS=corrected") == 0)
v2 <- values()[analysis == "corrected" & id %in% with_expr]
changed_08 <- v2$id[v2$status != a1$status]
ok(sprintf("of the %d rows that read a table of step 08 with tab(), %d changed, all to 'not run: 08_gene_annotation not done (stamp)'",
           length(reads_08), length(changed_08)),
   length(changed_08) > 0 && all(changed_08 %in% reads_08) &&
     v2[id %in% changed_08, all(status == "not run: 08_gene_annotation not done (stamp)")])
cat(sprintf("     rows reading step 08 that keep their status (a table of the paper folder, or step 10 absent first): %s\n",
            paste(setdiff(reads_08, changed_08), collapse = ", ")))
fwrite(s, stamp_path)

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

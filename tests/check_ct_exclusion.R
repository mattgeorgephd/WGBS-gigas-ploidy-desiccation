# tests/check_ct_exclusion.R
#
# T3-d0 of the settings layout (settings_build_spec.md, section 8): unit test of
# the C-to-T exclusion of step 05 (option ctexcl) with synthetic variant lists,
# run in a SCRATCH COPY of the repository WITHOUT .git and without
# 04_output/not_adopted/ct_variants/, never in the original checkout: it writes a
# synthetic list (and a synthetic run record, ct_variant_run.csv) into
# 04_output/not_adopted/ct_variants/tables/, where the real list lives.
# tests/scratch_guard.R refuses to run it in the original checkout.
#   mkdir <scratch> && tar --exclude=./.git --exclude=./04_output/not_adopted/ct_variants -cf - . | tar -xf - -C <scratch>
#     (from the repository root; or cp -a, then delete <scratch>/.git and <scratch>/04_output/not_adopted/ct_variants/)
#     (the copy needs 04_output/05_differential_methylation/cache/methylkit_legacy.rds)
#   cd <scratch> && CT_TEST_SCRATCH=1 WGBS_CACHE_DIR=<download cache> Rscript tests/check_ct_exclusion.R
# Three runs of step 05 with WGBS_ANALYSIS=ctexcl (fit from the paper's cache, no refit):
#   1. an empty list: every table equals the paper's, except the value of the
#      column mode and the file names;
#   2. a list of 23 CpGs (10 ploidy DML, 5 heat DML of which 2 are shared, 5
#      other tested CpGs, 5 untested positions), 3 of the kept DML marked as not
#      assessable: the DML lists equal the paper's minus the listed loci, the
#      p-values, q-values and differences of every kept row equal the paper's,
#      the united table loses exactly the listed tested CpGs, united_row is
#      renumbered, ct_exclusion_summary_ctexcl.csv adds up, and
#      united_at_dml_ctexcl.csv holds the united rows at the DML;
#   3. the list changed after the stamp was written: step 05 stops.
# The assessability table is synthetic as well: its coordinates are rebuilt from
# the reference genome as the variant step does, every CpG with 10 evidence reads
# except the three.

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("tests", "scratch_guard.R")); scratch_guard("tests/check_ct_exclusion.R")
rscript <- file.path(R.home("bin"), "Rscript")
ct_dir <- here("04_output", "not_adopted", "ct_variants", "tables")
dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(synthetic_run_marker, file.path(ct_dir, "ct_variant_run.csv"))     # marks the folder as this test's
p05 <- function(f) here("04_output", "05_differential_methylation", "tables", f)
c05 <- function(f) here("04_output", "not_adopted", "ctexcl", "05_differential_methylation", "tables", f)
results <- list()
ok <- function(name, pass) { results[[length(results) + 1]] <<- pass; cat(sprintf("%-4s %s\n", if (isTRUE(pass)) "ok" else "FAIL", name)) }

# ---- synthetic assessability (all CpGs of the reference) and its index ------------------------
source(here("run_config.R")); invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
genome <- Biostrings::readDNAStringSet(fetch_gannet_file("cgigas_uk_roslin_v1_genomic-mito.fa", config))
names(genome) <- sub(" .*", "", names(genome))
coords <- rbindlist(lapply(seq_along(genome), function(j)
  data.table(chr = names(genome)[j], start = BiocGenerics::start(Biostrings::matchPattern("CG", genome[[j]])) - 1L)))
rm(genome); invisible(gc())
md5_text <- function(starts) { tmp <- tempfile(); on.exit(unlink(tmp)); writeLines(format(starts, scientific = FALSE, trim = TRUE), tmp); unname(tools::md5sum(tmp)) }
index <- coords[, .(n_cpg = .N, first_start = start[1], last_start = start[.N], md5_starts = md5_text(start)), by = chr]
index[, first_row := cumsum(c(1L, head(n_cpg, -1L)))]
write_table(index[, .(chr, n_cpg, first_row, first_start, last_start, md5_starts)], file.path(ct_dir, "ct_variant_cpg_index.csv"))

paper <- list(united = fread(p05("united_cpg_legacy.csv.gz")), ploidy = fread(p05("dml_ploidy_legacy.csv")),
              heat = fread(p05("dml_heat_legacy.csv")), shared = fread(p05("dml_shared_legacy.csv")),
              dp = fread(p05("diffmeth_ploidy_legacy.csv.gz")), dh = fread(p05("diffmeth_heat_legacy.csv.gz")),
              counts = fread(p05("dml_counts_by_mode.csv")), sweep = fread(p05("dml_counts_by_threshold_legacy.csv")))
write_assess <- function(not_assessable) {
  a <- data.table(crick_pooled = rep(10L, nrow(coords)), n_assessable_zr3534 = 10L, n_assessable_zr3644 = 12L)
  if (nrow(not_assessable)) a[coords[not_assessable, on = c("chr", "start"), which = TRUE], crick_pooled := 0L]
  write_table(a, file.path(ct_dir, "ct_variant_assessability.csv.gz"))
}
write_list <- function(v) write_table(cbind(v, pos1 = v$start + 1L, listed = rep(TRUE, nrow(v))), file.path(ct_dir, "ct_variant_cpgs.csv.gz"))
stamp <- function() {        # the stamp that step 04 would write for the analysis ctexcl, with the md5 of the list
  dir.create(here("04_output", "not_adopted", "ctexcl"), recursive = TRUE, showWarnings = FALSE)
  write_table(data.table(id = "ctexcl", filter = "legacy", min_coverage = "10", overdispersion = "none", test = "not passed",
                         ct_excluded = "TRUE", fit_id = "legacy", fit_fingerprint_md5 = "", ct_variant_list_md5 = ct_variant_list_md5(),
                         git_commit = "", written_utc = "", steps_done = "04"), here("04_output", "not_adopted", "ctexcl", "analysis.csv"))
}
run05 <- function() {
  log <- tempfile()
  status <- system2(rscript, "01_pipeline/05_differential_methylation.R", stdout = log, stderr = log, env = "WGBS_ANALYSIS=ctexcl")
  l <- readLines(log); cat(tail(l, 4), sep = "\n")
  list(status = status, log = l)
}

# ---- 1. empty list -------------------------------------------------------------------------------
write_list(data.table(chr = character(0), start = integer(0))); write_assess(data.table(chr = character(0), start = integer(0))); stamp()
r <- run05()
ok("empty list: step 05 finishes, from the cached paper fit (no refit)", r$status == 0 && any(grepl("using the cached fit methylkit_legacy.rds", r$log)))
same_as_paper <- function(f_paper, f_ct, drop_mode = TRUE) {
  a <- fread(p05(f_paper)); b <- fread(c05(f_ct))
  if (drop_mode && "mode" %in% names(a)) { a[, mode := NULL]; b[, mode := NULL] }
  isTRUE(all.equal(a, b, check.attributes = FALSE))
}
for (f in c("united_cpg_%s.csv.gz", "diffmeth_ploidy_%s.csv.gz", "diffmeth_heat_%s.csv.gz", "dml_ploidy_%s.csv", "dml_heat_%s.csv",
            "dml_shared_%s.csv", "dml_counts_by_threshold_%s.csv", "library_filtering_%s.csv")) {
  ok(sprintf("empty list: %s equals the paper's", sprintf(f, "ctexcl")), same_as_paper(sprintf(f, "legacy"), sprintf(f, "ctexcl")))
}
ok("empty list: dml_counts_by_mode.csv equals the paper's except mode",
   same_as_paper("dml_counts_by_mode.csv", "dml_counts_by_mode.csv"))
s <- fread(c05("ct_exclusion_summary_ctexcl.csv"))
ok("empty list: summary removes nothing", all(s$cpg_removed == 0) && all(s$dml_removed == 0) && all(s$dml_before == s$dml_after))

# ---- 2. synthetic list ---------------------------------------------------------------------------
set.seed(1)
sh <- paper$shared[sample(.N, 2), .(chr, start)]
ploidy_only <- paper$ploidy[!paste(chr, start) %in% paste(paper$shared$chr, paper$shared$start)][sample(.N, 8), .(chr, start)]
heat_only <- paper$heat[!paste(chr, start) %in% paste(paper$shared$chr, paper$shared$start)][sample(.N, 3), .(chr, start)]
dml_keys <- paste(c(paper$ploidy$chr, paper$heat$chr), c(paper$ploidy$start, paper$heat$start))
other <- paper$united[!paste(chr, start) %in% dml_keys][sample(.N, 5), .(chr, start)]
untested <- coords[!paste(chr, start) %in% paste(paper$united$chr, paper$united$start)][sample(.N, 5)]
listed <- unique(rbind(sh, ploidy_only, heat_only, other, untested))
kept_dml <- paper$ploidy[!paste(chr, start) %in% paste(listed$chr, listed$start)][1:3, .(chr, start)]
write_list(listed); write_assess(kept_dml); stamp()
r <- run05()
ok("list of 23: step 05 finishes from the cached fit", r$status == 0 && any(grepl("using the cached fit", r$log)))
key <- function(d) paste(d$chr, d$start)
lk <- key(listed)
u <- fread(c05("united_cpg_ctexcl.csv.gz"))
ok("united table: the paper's minus the 18 listed tested CpGs", nrow(u) == nrow(paper$united) - 18 &&
     identical(key(u), key(paper$united)[!key(paper$united) %in% lk]))
for (set in c("ploidy", "heat", "shared")) {
  d <- fread(c05(sprintf("dml_%s_ctexcl.csv", set))); p <- paper[[set]][!key(paper[[set]]) %in% lk]
  cols <- setdiff(names(p), c("united_row", "united_row_ploidy", "united_row_heat"))
  ok(sprintf("dml_%s: the paper's list minus the listed loci, same p, q, difference (%d -> %d)", set, nrow(paper[[set]]), nrow(d)),
     identical(key(d), key(p)) && isTRUE(all.equal(d[, ..cols], p[, ..cols], check.attributes = FALSE)))
}
dp <- fread(c05("diffmeth_ploidy_ctexcl.csv.gz"))
ok("diffmeth_ploidy: united_row renumbered 1..n, p and q of every kept row equal the paper's",
   identical(dp$united_row, seq_len(nrow(dp))) &&
     isTRUE(all.equal(dp[, .(chr, start, pvalue, qvalue, meth.diff)], paper$dp[!key(paper$dp) %in% lk, .(chr, start, pvalue, qvalue, meth.diff)],
                      check.attributes = FALSE)))
s <- fread(c05("ct_exclusion_summary_ctexcl.csv"))
ok("summary adds up: 18 tested CpGs removed; ploidy DML 10 removed, heat DML 5 removed; 3 kept DML not assessable",
   all(s$cpg_removed == 18) && all(s$cpg_tested_before - s$cpg_removed == s$cpg_tested_after) &&
     s[contrast == "ploidy", dml_removed] == 10 && s[contrast == "heat", dml_removed] == 5 &&
     all(s$dml_before - s$dml_removed == s$dml_after) && s[contrast == "ploidy", dml_kept_not_assessable] == 3 &&
     s[contrast == "ploidy", cpg_kept_not_assessable] == 3)
ex <- fread(c05("ct_excluded_cpgs_ctexcl.csv.gz"))
ok("ct_excluded_cpgs: the 18 tested CpGs with their DML flags", nrow(ex) == 18 && sum(ex$dml_ploidy) == 10 && sum(ex$dml_heat) == 5)
counts <- fread(c05("dml_counts_by_mode.csv"))
ok("dml_counts_by_mode: cpg_analysed and DML counts after the removal",
   counts$cpg_analysed == paper$counts$cpg_analysed - 18 && counts$ploidy_dml == paper$counts$ploidy_dml - 10 &&
     counts$heat_dml == paper$counts$heat_dml - 5 && counts$shared == paper$counts$shared - 2)
sw <- fread(c05("dml_counts_by_threshold_ctexcl.csv"))
ok("counts by threshold at 20 points equal the DML lists after the removal",
   sw[difference == 20 & set == "ploidy", n] == counts$ploidy_dml && sw[difference == 20 & set == "heat", n] == counts$heat_dml)
uad <- fread(c05("united_at_dml_ctexcl.csv"))
at <- sort(unique(c(fread(c05("dml_ploidy_ctexcl.csv"))$united_row, fread(c05("dml_heat_ctexcl.csv"))$united_row)))
ok("united_at_dml_ctexcl.csv: the rows of the united table at the DML of either list, with united_row",
   identical(uad$united_row, at) && isTRUE(all.equal(uad[, !"united_row"], u[at], check.attributes = FALSE)))

# ---- 3. a list changed after step 04 ----------------------------------------------------------------
write_list(listed[-1])                     # the stamp keeps the md5 of the list of 23
r <- run05()
ok("a list changed since step 04 (md5 in the stamp differs): step 05 stops", r$status != 0 && any(grepl("is not the list that step 04 used", r$log)))
cat(sprintf("\nT3-d0: %d of %d checks pass\n", sum(unlist(results)), length(results)))
quit(status = if (all(unlist(results))) 0 else 1)

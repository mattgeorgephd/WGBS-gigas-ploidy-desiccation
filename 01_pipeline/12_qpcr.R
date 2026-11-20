# 12_qpcr.R
#
# Purpose   The qPCR data behind the HSP90 and "mitochondrial abundance"
#           statements of the Discussion (issue I-29). The Methods and Results of
#           the draft do not describe qPCR, so this step also records what was
#           measured and how it was analysed. It does five things:
#             1. Builds the per-animal table from the replicate-level files and
#                compares it with the consolidated file, which was built by hand.
#             2. Repeats the legacy analysis: actin as the only reference gene,
#                ln(2^-(Cq target - Cq actin)) per animal, aov() on Ploidy,
#                Desiccation (the legacy name of the heat treatment) and their
#                interaction, then TukeyHSD(). Both legacy variants are run:
#                non-detects dropped, and non-detects set to
#                config$qpcr$legacy_nondetect. The result is compared with the
#                knitted legacy output and with the audit note.
#             3. Fits delta Cq ~ ploidy * treatment per gene in two ways:
#                (a) complete cases, (b) a censored regression. Effects are
#                reported as fold changes with confidence intervals, with a
#                Benjamini-Hochberg correction across genes for each term.
#             4. Examines the reference gene: actin Cq by group, samples whose
#                actin Cq marks a failed sample, and how far each result moves
#                with the way the data are normalized.
#             5. Reports COX1 as what it is, a raw Cq without a reference, and
#                sets every qPCR statement of the draft beside the numbers.
#
# Run       Rscript 01_pipeline/12_qpcr.R
#           (needs no other step; under a minute, most of it the simulation
#           that checks the test of the censored fit. The step does not depend
#           on the optional analyses: it is part of the paper run only.)
#
# Inputs    03_input_files/qpcr/raw/qpcr_data_<gene>.csv   two wells per animal, 40 animals
#           03_input_files/qpcr/raw/qpcr_data_consolidated.csv   one row per animal, 32 animals
#           03_input_files/qpcr/legacy/dct_analysis_NA-to-45_COX1convert.md   printed legacy ANOVA and TukeyHSD
#             tables (a copy of 08_archive/ELISA/dct_analysis_NA-to-45_COX1convert.md)
# Outputs   04_output/12_qpcr/tables/
#             qpcr_per_animal.csv               one row per animal: Cq, detected wells, replicate difference
#             qpcr_delta_cq.csv                 one row per animal and target: delta Cq as analysed
#             qpcr_consolidated_check.csv       the consolidated file against the table built here
#             qpcr_nondetects.csv               non-detects per gene and group; which genes are analysed
#             qpcr_anova_legacy.csv             legacy ANOVA, both variants, checked against the legacy output
#             qpcr_tukey_legacy.csv             legacy TukeyHSD comparisons, checked against the legacy output;
#                                               for two groups also the Student test of those groups alone
#             qpcr_group_summary.csv            relative expression by group
#             qpcr_models.csv                   complete-case and censored fits: fold changes, intervals, p, q,
#                                               leave-one-out ranges; fits that repeat an earlier fit are marked
#             qpcr_censored_test_size.csv       simulation: how often the test of the censored fit rejects a true null
#             qpcr_sensitivity.csv              each gene and term under every analysis, side by side
#             qpcr_reference_gene.csv           actin by group, group differences, failed samples
#             qpcr_normalization_check.csv      each effect with and without normalization to actin
#             qpcr_cox1.csv                     COX1 raw Cq by group and the tests that can be made
#             qpcr_manuscript_statements.csv    each qPCR statement of the draft, the numbers, a verdict
#           04_output/12_qpcr/figures/figure_qpcr_expression.(png|pdf), figure_qpcr_cox1_actin.(png|pdf)
#           04_output/12_qpcr/session_info.txt
#
# Manuscript items: Discussion 4.2, 4.5, 4.7 and Conclusions 4.9 (HSP90 expression;
#           "relative mitochondrial abundance"). No Methods, Results, figure or table.
# Notes     - Sample IDs. Replicate files: D_N_N01 is ploidy (D diploid, T
#             triploid), treatment (N control, D heat), then the acute heat shock
#             (N none, H one hour at 45 C) joined to the animal number. The
#             consolidated file writes the same animal as D_C_N_01 (C control,
#             D heat).
#           - Acute heat shock. The 8 animals given it (numbers 09, 10, 19, 20
#             of each ploidy) are in the replicate files and were left out of the
#             legacy analysis. They are kept in qpcr_per_animal.csv, with their
#             flag, and in the figures, and are left out of every model and
#             test. Two animals per group are too few to form groups of their
#             own, and counted as control or heat animals they would change
#             what those groups are.
#           - Scale. Delta Cq is Cq target minus Cq actin; a lower delta Cq is a
#             higher expression. A difference is the first-named group minus the
#             second, in cycles, and its fold change is 2^-difference. This
#             assumes that target and reference double every cycle (Livak and
#             Schmittgen 2001, Methods 25:402, doi:10.1006/meth.2001.1262). No
#             primer efficiencies are in the repository.
#           - Technical replicates. An animal's Cq is the mean of its detected
#             wells, as in the consolidated file. An animal with one detected
#             well is represented by that well; the other well is treated as a
#             failed reaction. This overstates the expression of such an animal
#             when the cause was a low template, so the number of these animals
#             is reported per gene (qpcr_nondetects.csv), and a second censored
#             fit takes the other reading: the undetected well is above the
#             cycle limit, so the mean of the two wells is above the mean of the
#             detected Cq and the limit. The data cannot say which reading is
#             right; where the two fits differ, the gene is too close to the
#             detection limit for a firm estimate.
#           - Why a non-detect is not replaced by a number. The legacy scripts
#             set every non-detect to Cq 45 (config$qpcr$legacy_nondetect). The
#             runs had config$qpcr$max_cycles cycles, so no well can have that
#             Cq, and the value enters means and variances as if it had been
#             measured. After normalization the "expression" of a substituted
#             value is set by the reference gene alone: a sample with a late
#             actin Cq becomes a high expresser of a gene that was not detected
#             in it. McCall et al. (2014, Bioinformatics 30:2310,
#             doi:10.1093/bioinformatics/btu239) show that the common treatments
#             of non-detects, substituting the cycle limit and dropping them,
#             give biased estimates of expression and of differences in
#             expression, because non-detects are not missing at random.
#           - Censored regression. What is known about an animal with no detected
#             well is that its Cq is above config$qpcr$max_cycles, so its delta
#             Cq is above max_cycles minus its actin Cq. The model uses exactly
#             that: a Gaussian regression with right-censored observations
#             (Tobin 1958, Econometrica 26:24, doi:10.2307/1907382), fitted with
#             survival::survreg. It assumes that a non-detect is caused by a low
#             template, that delta Cq is normally distributed within a group,
#             and that an animal's bound says nothing about its delta Cq. The
#             last assumption is doubtful here: the bound falls as the actin Cq
#             rises, and so does delta Cq wherever a target Cq rises by less
#             than one cycle per cycle of actin (see "One reference gene"). The
#             fit with the actin Cq as covariate in
#             qpcr_normalization_check.csv has one bound, max_cycles, for every
#             animal and does not need that assumption.
#           - Small-sample correction of the censored regression. Maximum
#             likelihood standard errors are too small in a sample of this size,
#             so they are multiplied by sqrt(n / (n - k)), k being the number of
#             mean parameters, and tests use the t distribution with n - k
#             degrees of freedom. Without censoring this is the linear model of
#             (a): a gene without non-detects gives the same estimates,
#             intervals and p values, and the script stops if it does not. With
#             censoring the correction has no derivation, so the step checks it
#             by simulation: values are drawn under the null hypothesis for the
#             animals of the primary analysis, with their censoring bounds, at
#             shares of censored animals up to the largest share among the fits
#             of the step. qpcr_censored_test_size.csv and the log give how
#             often the corrected and the uncorrected test reject. The step
#             stops if the corrected test rejects more than 1.5 times as often
#             as its nominal level (the liberal criterion of Bradley 1978,
#             British Journal of Mathematical and Statistical Psychology 31:144,
#             doi:10.1111/j.2044-8317.1978.tb00581.x).
#           - Which genes are analysed. The limit of detection of a qPCR assay is
#             the concentration at which 95% of wells are detected (Bustin et al.
#             2009, Clinical Chemistry 55:611, doi:10.1373/clinchem.2008.112797,
#             section 7.4.3). A gene with fewer detected wells is labelled "near
#             the detection limit": it is analysed, and its result depends on the
#             treatment of non-detects. A gene detected in fewer than half of the
#             wells of any group is not analysed: the centre of that group is not
#             observed, and an estimate would rest on the assumed distribution
#             alone. Both thresholds are settings (lod_detected and min_detected
#             of config$qpcr in run_config.R); the class of each gene
#             is written to the log and to qpcr_nondetects.csv.
#           - Terms. The main effect of a factor is its effect averaged over the
#             two levels of the other factor. With unequal group sizes this
#             differs slightly from the sequential tests of the legacy tables.
#             Simple effects (one factor within one level of the other) use the
#             residual variance of all four groups; for the complete-case fits
#             the table adds Student and Welch tests on the two groups alone
#             (Welch 1947, Biometrika 34:28; Ruxton 2006, Behavioral Ecology
#             17:688, doi:10.1093/beheco/ark016). Simple effects are not adjusted
#             for the four simple effects of a gene. qpcr_tukey_legacy.csv
#             carries the same Student test, on the legacy response, beside
#             each Tukey comparison of two groups; the audit note quotes it.
#           - Multiple testing. Benjamini and Hochberg (1995, Journal of the
#             Royal Statistical Society B 57:289) across the analysed genes,
#             separately for each term, analysis and set of animals.
#           - Leave one out. With 8 animals per group one animal can decide a
#             test, so every fit is repeated with each animal left out in turn;
#             the table gives the lowest and the highest p value, the animal
#             whose removal produced it, and the number of these fits with p
#             below alpha. The tests of two groups alone are repeated in the
#             same way. The legacy tests are reproduced as run and have no
#             leave-one-out.
#           - Repeated fits. Two analyses of a gene can use the same animals with
#             the same values and the same censoring. Excluding a failed sample
#             changes nothing in a complete-case fit of a gene in which that
#             sample has no value; a gene without replicate-level data has no
#             known one-well animals; and without a censored animal the
#             censored regression is the linear model. qpcr_models.csv names
#             the earlier fit in the column same_fit_as. Counts of analyses in
#             qpcr_sensitivity.csv and in the statement table count results
#             that differ: p values that agree to six significant digits are
#             one result.
#           - One reference gene. Normalization against a single reference gene
#             needs evidence that the gene is stable under the conditions
#             compared (Bustin et al. 2009, section 8.1; Vandesompele et al.
#             2002, Genome Biology 3:research0034). No such evidence exists
#             here, so qpcr_normalization_check.csv fits each gene three ways:
#             delta Cq, raw Cq, and raw Cq with the actin Cq as a covariate. If
#             actin tracked the amount of template and nothing else, the slope
#             of a target Cq on the actin Cq would be 1; a smaller slope means
#             that part of the variation in actin is not shared by the target.
#           - Failed samples. A sample is called failed when its actin Cq lies
#             above the upper outer fence of the actin Cq of the animals in the
#             models: third quartile plus three interquartile ranges (Tukey
#             1977, Exploratory Data Analysis, Addison-Wesley). The analysis
#             this step treats as primary is the censored regression without
#             failed samples; every fit is also made with them. A sample
#             without any detected actin well is failed as well; it has no
#             delta Cq and is in no model of delta Cq, with or without the
#             failed samples.
#           - COX1 has no replicate-level file and no nuclear reference measured
#             on the same template. Mitochondrial abundance is the ratio of a
#             mitochondrial to a nuclear target in the same DNA (Quiros et al.
#             2017, Current Protocols in Mouse Biology 7:47,
#             doi:10.1002/cpmo.21), so it cannot be estimated from these data.
#           - Verdicts in qpcr_manuscript_statements.csv are derived from the
#             tables by rules that are written out in section 9.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "survival", "ggplot2", "patchwork", "ggrepel", "scales"))

step <- "12_qpcr"
t0 <- start_step(step, "qPCR: legacy analysis, censored regression, reference gene, COX1")
require_run_kind("paper", "this step does not depend on the analysis and writes the paper's outputs")

# ---- Settings ------------------------------------------------------------------------
# All settings of this step are in run_config.R, list config$qpcr: the reference
# gene, the cycle limit, the thresholds of the detection classes and of the
# verdicts, the seed of the simulation, and, as config$qpcr$design, the facts
# about the files that the step checks the data against.
opt <- config$qpcr
design <- opt$design

group_levels <- c("diploid control", "diploid heat", "triploid control", "triploid heat")
animal_sets <- c(all = "all animals", clean = "failed samples excluded")
ci_label <- sprintf("%g%% CI", 100 * opt$conf_level)
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
# Tables are written with 10 significant digits, which removes binary rounding noise such as 0.350000000000001.
write_tab <- function(x, name) {
  x <- as.data.table(x)
  write_table(x[, lapply(.SD, function(v) if (is.double(v)) signif(v, 10) else v)], output_path(step, "tables", name))
}

# =============================================================================
# 1. Read and check the inputs
# =============================================================================
raw_dir <- input_path("qpcr", "raw")
consolidated_file <- "qpcr_data_consolidated.csv"
replicate_files <- setdiff(list.files(raw_dir, pattern = "^qpcr_data_.+\\.csv$"), consolidated_file)
if (length(replicate_files) == 0) stop("No replicate-level file qpcr_data_<gene>.csv in ", raw_dir, call. = FALSE)

# Split sample IDs with a pattern; stop on any ID that does not fit.
id_parts <- function(id, pattern, where) {
  if (length(id) == 0) stop(where, ": no rows.", call. = FALSE)
  m <- regmatches(id, regexec(pattern, id))
  bad <- unique(id[lengths(m) == 0])
  if (length(bad) > 0) stop(where, ": sample ID(s) not in the expected form: ", paste(bad, collapse = ", "), call. = FALSE)
  do.call(rbind, m)
}

# ---- the consolidated file (built by hand; one row per animal) ----
cons <- as.data.table(read.csv(file.path(raw_dir, consolidated_file), check.names = FALSE, strip.white = TRUE,
                               stringsAsFactors = FALSE))
if (names(cons)[1] != "Sample" || anyDuplicated(cons$Sample) > 0) {
  stop("Consolidated file: the first column must be 'Sample', with one row per animal.")
}
cons_genes <- setdiff(names(cons), "Sample")
if (!all(vapply(cons[, ..cons_genes], is.numeric, logical(1)))) stop("Consolidated file: a gene column holds text other than NA.")
ref <- opt$reference_gene
if (!ref %in% cons_genes) stop("Reference gene '", ref, "' is not a column of the consolidated file.")
pc <- id_parts(cons$Sample, "^([DT])_([CD])_(N)_([0-9]{2})$", "consolidated file")
cons[, animal_id := paste0(pc[, 2], "_", c(C = "N", D = "D")[pc[, 3]], "_", pc[, 4], pc[, 5])]
legacy_code <- data.frame(Ploidy = pc[, 2], Desiccation = pc[, 3])      # legacy factor codes, in the row order of the file

# ---- the replicate-level files (two wells per animal) ----
gene_key  <- sub("^qpcr_data_(.+)\\.csv$", "\\1", replicate_files)
rep_genes <- cons_genes[match(gene_key, tolower(cons_genes))]            # gene names as spelled in the consolidated file
rep_genes[is.na(rep_genes)] <- toupper(gene_key[is.na(rep_genes)])
if (!ref %in% rep_genes) stop("No replicate-level file for the reference gene '", ref, "'.")
in_cons_order <- order(match(rep_genes, cons_genes), rep_genes)          # genes in the column order of the consolidated file
replicate_files <- replicate_files[in_cons_order]; rep_genes <- rep_genes[in_cons_order]
read_wells <- function(file, gene) {
  x <- read.csv(file.path(raw_dir, file), colClasses = "character", check.names = FALSE, strip.white = TRUE,
                na.strings = character(0))
  names(x) <- trimws(names(x))                                           # the actin header is "Sample " with a trailing space
  if (!identical(names(x), c("Sample", "Cq Value"))) stop(file, ": expected the columns 'Sample' and 'Cq Value'.", call. = FALSE)
  data.table(gene = gene, animal_id = x[["Sample"]], cq_text = x[["Cq Value"]])
}
wells <- rbindlist(Map(read_wells, replicate_files, rep_genes))
wells[, cq := suppressWarnings(as.numeric(cq_text))]
wells[, detected := !is.na(cq)]
unreadable <- wells[detected == FALSE & cq_text != design$nondetect_text]
if (nrow(unreadable) > 0) stop("Cq value(s) that are neither a number nor '", design$nondetect_text, "': ",
                               paste(unique(paste(unreadable$gene, unreadable$animal_id, unreadable$cq_text)), collapse = "; "))
if (wells[detected == TRUE, any(cq <= 0 | cq > opt$max_cycles)]) {
  stop("A detected Cq is outside 0 to config$qpcr$max_cycles (", opt$max_cycles, "). A non-detect is defined by that limit; check it.")
}

# ---- the animals, from the sample IDs ----
animals <- unique(wells[, .(animal_id)])
pr <- id_parts(animals$animal_id, "^([DT])_([ND])_([NH])([0-9]{2})$", "replicate files")
animals[, `:=`(ploidy = factor(c(D = "diploid", T = "triploid")[pr[, 2]], levels = c("diploid", "triploid")),
               treatment = factor(c(N = "control", D = "heat")[pr[, 3]], levels = c("control", "heat")),
               acute_heat_shock = pr[, 4] == "H", animal_number = as.integer(pr[, 5]))]
animals[, group := factor(paste(ploidy, treatment), levels = group_levels)]
# Numbering of the project (root README): 01 to 10 control, 11 to 20 heat; 09, 10, 19 and 20 had the acute heat shock.
if (animals[, any((animal_number > 10) != (treatment == "heat")) || any((animal_number %% 10 %in% c(9, 0)) != acute_heat_shock)]) {
  stop("Replicate files: the treatment or acute heat shock code of a sample ID contradicts its animal number.")
}
sizes <- animals[, .N, by = .(group, acute_heat_shock)]
if (nrow(sizes) != 2 * length(group_levels) || sizes[acute_heat_shock == FALSE, any(N != design$animals_per_group)] ||
    sizes[acute_heat_shock == TRUE, any(N != design$acute_per_group)]) {
  stop("Replicate files: expected ", design$animals_per_group, " animals per group plus ", design$acute_per_group,
       " with the acute heat shock. Found:\n", paste(capture.output(print(sizes[order(group, acute_heat_shock)])), collapse = "\n"))
}
per_file <- wells[, .N, by = .(gene, animal_id)]
if (nrow(per_file) != length(rep_genes) * nrow(animals) || any(per_file$N != design$wells_per_animal)) {
  stop("Replicate files: every file must hold ", design$wells_per_animal, " wells for each of the ", nrow(animals), " animals.")
}
if (!setequal(cons$animal_id, animals[acute_heat_shock == FALSE, animal_id])) {
  stop("The consolidated file does not hold exactly the animals without acute heat shock of the replicate files.")
}
n_animals_read <- nrow(animals)
animals <- merge(animals, cons[, .(animal_id, consolidated_id = Sample)], by = "animal_id", all.x = TRUE)
stopifnot(nrow(animals) == n_animals_read)
animals[, in_models := !acute_heat_shock]                                # see the header note on the acute heat shock
cons_only <- setdiff(cons_genes, rep_genes)                              # genes without a replicate-level file
log_msg("%d replicate-level files (%s); %d animals, %d wells; %d non-detect wells",
        length(rep_genes), paste(rep_genes, collapse = ", "), nrow(animals), nrow(wells), sum(!wells$detected))
log_msg("consolidated file: %d animals, %d genes; without a replicate-level file: %s",
        nrow(cons), length(cons_genes), if (length(cons_only) > 0) paste(cons_only, collapse = ", ") else "none")
mito <- opt$mito_gene
if (!mito %in% c(rep_genes, cons_genes)) {
  stop("The gene named in mito_gene ('", mito, "') is in neither the replicate files nor the consolidated file.")
}

# =============================================================================
# 2. Per-animal values, and the consolidated file checked against them
# =============================================================================
per_gene <- wells[, .(n_detected = sum(detected),
                      cq = if (any(detected)) mean(cq[detected]) else NA_real_,
                      replicate_difference = if (sum(detected) == 2) abs(diff(cq[detected])) else NA_real_),
                  by = .(gene, animal_id)]
cons_long <- melt(cons, id.vars = c("Sample", "animal_id"), measure.vars = cons_genes, variable.name = "gene",
                  value.name = "cq_consolidated", variable.factor = FALSE)

chk <- merge(cons_long, per_gene[, .(gene, animal_id, cq_built = cq)], by = c("gene", "animal_id"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(chk) == nrow(cons_long))
chk[, has_file := gene %chin% rep_genes]
chk[, differs := has_file & (is.na(cq_built) != is.na(cq_consolidated) |
                               (!is.na(cq_built) & !is.na(cq_consolidated) & abs(cq_built - cq_consolidated) > opt$cq_tolerance))]
consolidated_check <- chk[, .(
  replicate_file = has_file[1], cells = .N, consolidated_missing = sum(is.na(cq_consolidated)),
  cells_compared = if (has_file[1]) sum(!is.na(cq_built) & !is.na(cq_consolidated)) else NA_integer_,
  missing_in_both = if (has_file[1]) sum(is.na(cq_built) & is.na(cq_consolidated)) else NA_integer_,
  largest_absolute_difference = if (has_file[1]) round(max(c(0, abs(cq_built - cq_consolidated)), na.rm = TRUE), 6) else NA_real_,
  cells_that_differ = if (has_file[1]) sum(differs) else NA_integer_,
  which_cells = paste(Sample[differs], collapse = "; ")), by = gene]
consolidated_check[, result := fifelse(!replicate_file, "cannot be checked: no replicate-level file exists",
                               fifelse(cells_that_differ == 0, "consolidated values are the means of the detected wells",
                                       "differs from the replicate-level file"))]
consolidated_check <- consolidated_check[match(cons_genes, gene)]
write_tab(consolidated_check, "qpcr_consolidated_check.csv")
log_msg(paste("consolidated file against the replicate files: %d genes compared, %d cells differ; %d missing cells in the file",
              "(%d without a detected well in the replicate files; %d in genes without a replicate-level file)"),
        sum(consolidated_check$replicate_file), sum(consolidated_check$cells_that_differ, na.rm = TRUE),
        sum(consolidated_check$consolidated_missing), sum(consolidated_check$missing_in_both, na.rm = TRUE),
        consolidated_check[replicate_file == FALSE, sum(consolidated_missing)])

# =============================================================================
# 3. Reference gene and failed samples; the per-animal table
# =============================================================================
animals <- merge(animals, per_gene[gene == ref, .(animal_id, reference_cq = cq)], by = "animal_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(animals) == n_animals_read)
# An animal without a detected reference well has no delta Cq: it is a failed sample and is in no model of delta Cq.
no_reference <- animals[is.na(reference_cq), animal_id]
if (length(no_reference) > 0) {
  log_msg("no %s Cq for: %s; counted as failed sample(s) and left out of every model of delta Cq", ref, paste(no_reference, collapse = ", "))
}
if (animals[in_models == TRUE, sum(!is.na(reference_cq))] < 2) stop("Fewer than two animals in the models have a Cq of the reference gene ", ref, ".")
ref_quartiles <- quantile(animals[in_models == TRUE & !is.na(reference_cq), reference_cq], c(0.25, 0.75), names = FALSE)
ref_fence <- ref_quartiles[2] + opt$failed_fence_iqr * diff(ref_quartiles)
animals[, failed_sample := is.na(reference_cq) | reference_cq > ref_fence]
setorder(animals, ploidy, animal_number)
failed_ids <- animals[in_models == TRUE & failed_sample == TRUE, animal_id]            # failed samples among the animals in the models
failed_text <- if (length(failed_ids) == 0) "none" else
  paste(sprintf("%s (%s)", failed_ids, animals[match(failed_ids, animal_id), fifelse(is.na(reference_cq), "no Cq", sprintf("%.2f", reference_cq))]),
        collapse = ", ")
log_msg("%s Cq of the %d animals in the models (%d with a Cq): quartiles %.2f and %.2f, outer fence %.2f; failed sample(s): %s",
        ref, sum(animals$in_models), animals[in_models == TRUE, sum(!is.na(reference_cq))], ref_quartiles[1], ref_quartiles[2], ref_fence,
        failed_text)
# The two sets of animals of the delta Cq models. "All animals" are the animals in the models that have a reference Cq.
in_set <- list(all = animals[in_models == TRUE & !is.na(reference_cq), animal_id],
               clean = animals[in_models == TRUE & failed_sample == FALSE, animal_id])

wide <- dcast(per_gene, animal_id ~ gene, value.var = c("cq", "n_detected", "replicate_difference"))
for (gene_i in rep_genes) setnames(wide, paste0(c("cq_", "n_detected_", "replicate_difference_"), gene_i),
                                   paste0(gene_i, c("_cq", "_n_detected", "_replicate_difference")))
setcolorder(wide, c("animal_id", unlist(lapply(rep_genes, function(gene_i) paste0(gene_i, c("_cq", "_n_detected", "_replicate_difference"))))))
per_animal <- merge(animals[, .(animal_id, consolidated_id, ploidy, treatment, acute_heat_shock, animal_number, in_models, failed_sample)],
                    wide, by = "animal_id", sort = FALSE)
if (length(cons_only) > 0) {
  from_cons <- cons[, c("animal_id", cons_only), with = FALSE]
  setnames(from_cons, cons_only, paste0(cons_only, "_cq"))
  per_animal <- merge(per_animal, from_cons, by = "animal_id", all.x = TRUE, sort = FALSE)
}
per_animal[, genes_without_replicate_data := paste(cons_only, collapse = "; ")]   # their Cq are copied from the consolidated file
stopifnot(nrow(per_animal) == n_animals_read)
write_tab(per_animal, "qpcr_per_animal.csv")

# =============================================================================
# 4. Delta Cq per animal and target; non-detects; which genes are analysed
# =============================================================================
targets_rep  <- setdiff(rep_genes, c(ref, mito))
targets_cons <- setdiff(cons_only, mito)        # targets known from the consolidated file only
targets <- c(targets_rep, targets_cons)
dq <- rbind(per_gene[gene %chin% targets_rep, .(gene, animal_id, n_detected, cq)],
            cons_long[gene %chin% targets_cons, .(gene, animal_id, n_detected = NA_integer_, cq = cq_consolidated)])
n_dq <- nrow(dq)
dq <- merge(dq, animals[, .(animal_id, ploidy, treatment, group, acute_heat_shock, animal_number, in_models, failed_sample, reference_cq)],
            by = "animal_id", sort = FALSE)
stopifnot(nrow(dq) == n_dq)
dq[, replicate_data := gene %chin% targets_rep]
dq[, detection := fifelse(!replicate_data, fifelse(is.na(cq), "no value in the consolidated file", "value in the consolidated file"),
                  fifelse(n_detected == 0, "no well detected", fifelse(n_detected == 1, "one well detected", "both wells detected")))]
dq[, censored := is.na(cq)]
dq[, delta_cq := cq - reference_cq]
dq[, delta_cq_above := fifelse(censored, opt$max_cycles - reference_cq, NA_real_)]   # all that is known when no well was detected
dq[, relative_expression := 2^-delta_cq]
dq[, relative_expression_below := 2^-delta_cq_above]
dq[, gene := factor(gene, levels = targets)]
setorder(dq, gene, ploidy, animal_number)
write_tab(dq[, .(gene, replicate_data, animal_id, ploidy, treatment, acute_heat_shock, in_models, failed_sample, detection, n_detected,
                 cq, reference_cq, delta_cq, delta_cq_above, relative_expression, relative_expression_below)], "qpcr_delta_cq.csv")

# ---- non-detects per gene and group ----
wells_by_group <- merge(wells, animals[, .(animal_id, group)], by = "animal_id", sort = FALSE)
stopifnot(nrow(wells_by_group) == nrow(wells))
count_wells <- function(ids, scope, with_consolidated = TRUE) {
  w <- wells_by_group[animal_id %chin% ids]
  by_group <- function(by_cols) {
    a <- w[, .(n_det = sum(detected)), by = c("gene", "animal_id", by_cols)]
    of_wells <- w[, .(n_wells = .N, n_wells_not_detected = sum(!detected), median_cq_detected_wells = median(cq[detected]),
                      latest_cq_detected = if (any(detected)) max(cq[detected]) else NA_real_), by = c("gene", by_cols)]
    of_animals <- a[, .(n_animals = .N, n_animals_both_wells = sum(n_det == design$wells_per_animal),
                        n_animals_one_well = sum(n_det > 0 & n_det < design$wells_per_animal), n_animals_no_well = sum(n_det == 0)),
                    by = c("gene", by_cols)]
    both <- merge(of_wells, of_animals, by = c("gene", by_cols))
    stopifnot(nrow(both) == nrow(of_wells), nrow(both) == nrow(of_animals))
    both
  }
  out <- rbind(by_group("group"), cbind(by_group(character(0)), group = "all groups"), use.names = TRUE)
  cons_part <- cons_long[gene %chin% targets_cons & animal_id %chin% ids]
  if (with_consolidated && nrow(cons_part) > 0) {     # no wells are known; an animal without a value counts as "no well"
    n_cons_part <- nrow(cons_part)
    cons_part <- merge(cons_part, animals[, .(animal_id, group)], by = "animal_id", sort = FALSE)
    stopifnot(nrow(cons_part) == n_cons_part)
    cons_counts <- function(d) d[, .(n_wells = NA_integer_, n_wells_not_detected = NA_integer_, median_cq_detected_wells = NA_real_,
                                     latest_cq_detected = NA_real_, n_animals = .N, n_animals_both_wells = NA_integer_,
                                     n_animals_one_well = NA_integer_, n_animals_no_well = sum(is.na(cq_consolidated)))]
    out <- rbind(out, cons_part[, cons_counts(.SD), by = .(gene, group)],
                 cbind(cons_part[, cons_counts(.SD), by = gene], group = "all groups"), use.names = TRUE)
  }
  out[, animals_included := scope]
  out[]
}
every_animal <- "every animal in the replicate files"
nondetects <- rbind(count_wells(animals$animal_id, every_animal, with_consolidated = FALSE),
                    count_wells(in_set$all, animal_sets[["all"]]), count_wells(in_set$clean, animal_sets[["clean"]]))
# Share detected: of wells where wells are known, of animals for a gene known from the consolidated file only.
nondetects[, share_detected := fifelse(is.na(n_wells), 1 - n_animals_no_well / n_animals, 1 - n_wells_not_detected / n_wells)]
classes <- nondetects[group != "all groups" & animals_included != every_animal & gene %chin% targets,
                      .(lowest_share_in_a_group = min(share_detected)), by = .(gene, animals_included)]
classes <- merge(classes, nondetects[group == "all groups", .(gene, animals_included, share_all = share_detected)],
                 by = c("gene", "animals_included"))
classes[, detection_class := fifelse(lowest_share_in_a_group < opt$min_detected, "not analysed",
                             fifelse(share_all >= opt$lod_detected, "quantifiable", "near the detection limit"))]
stopifnot(nrow(classes) == length(targets) * length(animal_sets))
n_nd <- nrow(nondetects)
nondetects <- merge(nondetects, classes[, .(gene, animals_included, lowest_share_in_a_group, detection_class)],
                    by = c("gene", "animals_included"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(nondetects) == n_nd)
nondetects[, replicate_data := gene %chin% rep_genes]
nondetects[, `:=`(gene = factor(gene, levels = unique(c(rep_genes, targets_cons))), group = factor(group, levels = c(group_levels, "all groups")),
                  animals_included = factor(animals_included, levels = c(every_animal, animal_sets)))]
setorder(nondetects, animals_included, gene, group)
setcolorder(nondetects, c("animals_included", "gene", "replicate_data", "group", "n_wells", "n_wells_not_detected", "share_detected", "n_animals",
                          "n_animals_both_wells", "n_animals_one_well", "n_animals_no_well", "median_cq_detected_wells", "latest_cq_detected",
                          "lowest_share_in_a_group", "detection_class"))
write_tab(nondetects, "qpcr_nondetects.csv")
class_of <- function(gene_i, set_i) classes[gene == gene_i & animals_included == animal_sets[[set_i]], detection_class]
for (set_i in names(animal_sets)) {
  cl <- classes[animals_included == animal_sets[[set_i]]][match(targets, gene)]
  genes_of <- function(class_i) if (any(cl$detection_class == class_i)) paste(cl[detection_class == class_i, gene], collapse = " ") else "none"
  log_msg("%s (%d): quantifiable %s; near the detection limit %s; not analysed %s", animal_sets[[set_i]], length(in_set[[set_i]]),
          genes_of("quantifiable"), genes_of("near the detection limit"), genes_of("not analysed"))
}
if (all(classes$detection_class == "not analysed")) stop("No gene was detected in enough wells to be analysed.")

# =============================================================================
# 5. The legacy analysis, as run
# =============================================================================
# Input: the consolidated file. Response: ln(2^-(Cq target - Cq actin)); for COX1 ln(2^-Cq) in the R Markdown file and
# ln(Cq) in the R script. Model and term names are the legacy ones; "Desiccation" is the heat treatment.
# The legacy script without substitution (08_archive/gene_expression/scripts/dct_analysis.R) held ANOVA calls for seven
# targets; it has none for MBD2, SOD, ATP synthetase or COX1. Both variants are run here for every gene of the file.
legacy_terms <- c(ploidy = "Ploidy", treatment = "Desiccation", interaction = "Ploidy:Desiccation")
legacy_variants <- c(dropped = "dropped", substituted = sprintf("set to Cq %g", opt$legacy_nondetect))
legacy_response <- sprintf("ln(2^-(Cq - %s Cq))", ref)
legacy_targets <- setdiff(cons_genes, c(ref, mito))          # the legacy scripts read the consolidated file only
legacy_words <- function(label) {       # "T:C-D:C", as printed by TukeyHSD(), in plain words
  side <- function(s) {
    p <- strsplit(s, ":", fixed = TRUE)[[1]]
    if (length(p) == 2) paste(c(D = "diploid", T = "triploid")[p[1]], c(C = "control", D = "heat")[p[2]]) else NA_character_
  }
  two <- strsplit(label, "-", fixed = TRUE)[[1]]
  paste(side(two[1]), "vs", side(two[2]))
}
legacy_input <- function(variant) {
  d <- as.data.frame(cons[, ..cons_genes])
  if (variant == "substituted") d[is.na(d)] <- opt$legacy_nondetect       # the legacy line: dCt[is.na(dCt)] <- 45
  cbind(d, legacy_code)
}
# Student test of two groups alone on a legacy response, not adjusted. It is not part of the legacy analysis: the audit
# note quotes it for the two control groups. pair is a comparison as TukeyHSD() prints it, "T:C-D:C".
legacy_pair_p <- function(pair, response_values, cell_of_row) {
  cells <- strsplit(pair, "-", fixed = TRUE)[[1]]
  a <- response_values[cell_of_row == cells[1] & !is.na(response_values)]
  b <- response_values[cell_of_row == cells[2] & !is.na(response_values)]
  if (length(a) < 2 || length(b) < 2) NA_real_ else t.test(a, b, var.equal = TRUE)$p.value
}
legacy_anova <- list(); legacy_tukey <- list()
for (variant_i in names(legacy_variants)) {
  ld <- legacy_input(variant_i)
  legacy_cell <- paste(ld$Ploidy, ld$Desiccation, sep = ":")
  responses <- c(
    lapply(legacy_targets, function(gene_i) list(gene = gene_i, response = legacy_response, y = log(2^-(ld[[gene_i]] - ld[[ref]])))),
    if (mito %in% cons_genes) list(
      list(gene = mito, response = "ln(2^-Cq)", y = log(2^-ld[[mito]])),
      list(gene = mito, response = "ln(Cq)", y = log(ld[[mito]]))))
  for (r_i in responses) {
    ld$y <- r_i$y
    fit <- aov(y ~ Ploidy + Desiccation + Ploidy:Desiccation, data = ld)
    s <- summary(fit)[[1]]; rn <- trimws(rownames(s)); at <- match(legacy_terms, rn); at_res <- match("Residuals", rn)
    id_cols <- data.table(gene = r_i$gene, response = r_i$response, nondetects = legacy_variants[[variant_i]])
    legacy_anova[[length(legacy_anova) + 1]] <- cbind(id_cols, data.table(
      n_animals = nobs(fit), term = names(legacy_terms), legacy_term = unname(legacy_terms), df = s[at, "Df"], df_residual = s[at_res, "Df"],
      sum_sq = s[at, "Sum Sq"], sum_sq_residual = s[at_res, "Sum Sq"], F_value = s[at, "F value"], p_value = s[at, "Pr(>F)"]))
    tk <- TukeyHSD(fit)
    legacy_tukey[[length(legacy_tukey) + 1]] <- cbind(id_cols, rbindlist(lapply(names(tk), function(k_i) {
      m <- tk[[k_i]]
      data.table(legacy_term = k_i, legacy_comparison = rownames(m),
                 comparison = if (k_i == "Ploidy") "triploid vs diploid" else if (k_i == "Desiccation") "heat vs control" else
                   vapply(rownames(m), legacy_words, ""),
                 difference = m[, "diff"], lower = m[, "lwr"], upper = m[, "upr"], p_adjusted = m[, "p adj"],
                 p_two_groups_student = if (k_i != legacy_terms[["interaction"]]) NA_real_ else
                   vapply(rownames(m), legacy_pair_p, 0, response_values = ld$y, cell_of_row = legacy_cell))
    })))
  }
}
legacy_anova <- rbindlist(legacy_anova); legacy_tukey <- rbindlist(legacy_tukey)
# On the ln(2^-delta Cq) scale a difference of means is the log of a ratio of geometric means.
legacy_tukey[, ratio := fifelse(response == "ln(Cq)", NA_real_, exp(difference))]

# ---- the printed legacy output: ANOVA and TukeyHSD tables of the knitted R Markdown file ----
knitted_file <- input_path("qpcr", "legacy", "dct_analysis_NA-to-45_COX1convert.md")
if (!file.exists(knitted_file)) stop("Legacy output not found: ", knitted_file)
read_knitted <- function(path) {
  l <- readLines(path, warn = FALSE)
  at_tukey <- grep("^TukeyHSD\\([A-Za-z0-9]+\\)\\s*$", l); at_summary <- grep("^summary\\([A-Za-z0-9]+\\)\\s*$", l)
  chunk_start <- grep("^``` r", l)
  bad_form <- function(what) stop("Legacy output: ", what, " not in the expected form: ", path, call. = FALSE)
  if (length(at_summary) == 0 || length(at_tukey) != length(at_summary) || any(at_tukey >= at_summary)) {
    bad_form("the TukeyHSD() and summary() calls are")
  }
  fields <- function(lines, pattern) strsplit(trimws(sub("^\\s*##\\s+", "", grep(pattern, lines, value = TRUE))), "\\s+")
  parts <- lapply(seq_along(at_summary), function(k) {
    i <- at_summary[k]; g <- sub("^summary\\((.+)\\)\\s*$", "\\1", l[i])
    stop_at <- min(c(chunk_start[chunk_start > i], length(l) + 1)) - 1
    f <- fields(l[i:stop_at], "^\\s*##\\s+(Ploidy|Desiccation|Ploidy:Desiccation)\\s")       # Df, Sum Sq, Mean Sq, F value, Pr(>F)
    at_aov <- grep(sprintf("^%s<-aov\\(", g), l)
    ss <- if (length(at_aov) > 0) fields(l[at_aov[1]:i], "^\\s*##\\s+Sum of Squares") else list()   # "Sum" "of" "Squares", then four values
    tk <- fields(l[at_tukey[k]:i], "^\\s*##\\s+[A-Z](:[A-Z])?-[A-Z](:[A-Z])?\\s")              # comparison, diff, lwr, upr, p adj
    if (length(f) != length(legacy_terms) || any(lengths(f) < 6) || length(ss) != 1 || length(ss[[1]]) != length(legacy_terms) + 4) {
      bad_form(paste("the ANOVA table of", g, "is"))
    }
    if (length(tk) == 0 || any(lengths(tk) != 5)) bad_form(paste("the TukeyHSD table of", g, "is"))
    pick_field <- function(x, j) vapply(x, `[`, "", j)
    list(anova = data.table(gene = g, legacy_term = pick_field(f, 1), printed_sum_sq = ss[[1]][3 + seq_along(legacy_terms)],
                            printed_F = pick_field(f, 5), printed_p = pick_field(f, 6)),
         tukey = data.table(gene = g, legacy_comparison = pick_field(tk, 1), printed_difference = pick_field(tk, 2),
                            printed_lower = pick_field(tk, 3), printed_upper = pick_field(tk, 4), printed_p_adjusted = pick_field(tk, 5)))
  })
  out <- list(anova = rbindlist(lapply(parts, `[[`, "anova")), tukey = rbindlist(lapply(parts, `[[`, "tukey")))
  printed <- unlist(c(out$anova[, .(printed_sum_sq, printed_F, printed_p)], out$tukey[, !c("gene", "legacy_comparison")]))
  if (anyNA(suppressWarnings(as.numeric(printed)))) bad_form("a printed value is not a number; the tables are")
  # The knitted file ran the substitution; its COX1 response was ln(2^-Cq).
  for (x in out) x[, `:=`(nondetects = legacy_variants[["substituted"]], response = fifelse(gene == mito, "ln(2^-Cq)", legacy_response))]
  out
}
# A regenerated value agrees with a printed one when it rounds to it at the printed number of decimals.
agrees <- function(value, printed) {
  decimals <- nchar(sub("^[^.]*\\.?", "", printed))
  abs(value - as.numeric(printed)) <= 0.5 * 10^-decimals * (1 + 1e-6)
}
knitted <- read_knitted(knitted_file)
n_anova <- nrow(legacy_anova); n_tukey <- nrow(legacy_tukey)
legacy_anova <- merge(legacy_anova, knitted$anova, by = c("gene", "response", "nondetects", "legacy_term"), all.x = TRUE, sort = FALSE)
legacy_tukey <- merge(legacy_tukey, knitted$tukey, by = c("gene", "response", "nondetects", "legacy_comparison"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(legacy_anova) == n_anova, sum(!is.na(legacy_anova$printed_F)) == nrow(knitted$anova),
          nrow(legacy_tukey) == n_tukey, sum(!is.na(legacy_tukey$printed_difference)) == nrow(knitted$tukey))
legacy_anova[, agrees_with_legacy_output := fifelse(is.na(printed_F), NA,
  agrees(sum_sq, printed_sum_sq) & agrees(F_value, printed_F) & agrees(p_value, printed_p))]
legacy_tukey[, agrees_with_legacy_output := fifelse(is.na(printed_difference), NA,
  agrees(difference, printed_difference) & agrees(lower, printed_lower) & agrees(upper, printed_upper) & agrees(p_adjusted, printed_p_adjusted))]

# ---- the audit note ----
# Values quoted from 07_documentation/audit/elisa_qpcr_audit.md, section 3.3. Used for comparison only.
audit_note <- data.table(
  gene = rep(c("HSP90", "HSP90", "COX1"), each = 3), response = rep(c(legacy_response, legacy_response, "ln(2^-Cq)"), each = 3),
  nondetects = rep(legacy_variants[c("substituted", "dropped", "substituted")], each = 3), term = rep(names(legacy_terms), 3),
  audit_note_F = c("4.90", "6.39", NA, "3.61", "5.22", NA, NA, NA, "3.56"),
  audit_note_p = c("0.035", "0.017", "0.68", "0.068", "0.030", "0.84", "0.69", "0.69", "0.070"))
legacy_anova <- merge(legacy_anova, audit_note, by = c("gene", "response", "nondetects", "term"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(legacy_anova) == n_anova, sum(!is.na(legacy_anova$audit_note_p)) == nrow(audit_note))
legacy_anova[, agrees_with_audit_note := fifelse(is.na(audit_note_p), NA,
  (is.na(audit_note_F) | agrees(F_value, audit_note_F)) & agrees(p_value, audit_note_p))]
legacy_order <- c(legacy_targets, mito)
legacy_anova <- legacy_anova[order(nondetects, match(gene, legacy_order), response, match(term, names(legacy_terms)))]
legacy_tukey <- legacy_tukey[order(nondetects, match(gene, legacy_order), response)]
setcolorder(legacy_anova, c("gene", "response", "nondetects", "n_animals", "term", "legacy_term"))
setcolorder(legacy_tukey, c("gene", "response", "nondetects", "legacy_term", "legacy_comparison", "comparison"))
write_tab(legacy_anova, "qpcr_anova_legacy.csv")
write_tab(legacy_tukey, "qpcr_tukey_legacy.csv")
by_gene <- legacy_anova[!is.na(agrees_with_legacy_output), .(agree = all(agrees_with_legacy_output)), by = gene]
log_msg("legacy ANOVA (non-detects %s) against the knitted output: %d of %d genes agree at the printed precision%s",
        legacy_variants[["substituted"]], sum(by_gene$agree), nrow(by_gene),
        if (all(by_gene$agree)) "" else paste0("; not: ", paste(by_gene[agree == FALSE, gene], collapse = ", ")))
log_msg("legacy TukeyHSD against the knitted output: %d of %d printed comparisons agree",
        sum(legacy_tukey$agrees_with_legacy_output, na.rm = TRUE), sum(!is.na(legacy_tukey$agrees_with_legacy_output)))
log_msg("legacy ANOVA against the audit note: %d of %d quoted terms agree",
        sum(legacy_anova$agrees_with_audit_note, na.rm = TRUE), sum(!is.na(legacy_anova$agrees_with_audit_note)))

# =============================================================================
# 6. Models of delta Cq: complete cases and censored regression
# =============================================================================
term_info <- data.table(
  term = c("ploidy", "treatment", "interaction", "ploidy_in_control", "ploidy_in_heat", "treatment_in_diploid", "treatment_in_triploid"),
  contrast = c("triploid vs diploid, averaged over the two treatments", "heat vs control, averaged over the two ploidies",
               "heat vs control in triploids, relative to heat vs control in diploids",
               "triploid vs diploid, control animals", "triploid vs diploid, heat animals",
               "heat vs control, diploids", "heat vs control, triploids"))
# Weights on the four group means, in the order of group_levels.
K <- rbind(c(-1, -1, 1, 1) / 2, c(-1, 1, -1, 1) / 2, c(1, -1, -1, 1), c(-1, 0, 1, 0), c(0, -1, 0, 1), c(-1, 1, 0, 0), c(0, 0, -1, 1))
dimnames(K) <- list(term_info$term, group_levels)
model_terms <- term_info$term[1:3]; simple_terms <- term_info$term[4:7]
contrast_of <- function(term_i) term_info$contrast[match(term_i, term_info$term)]
methods <- c("complete cases", "censored regression", "censored regression, one-well animals as bounds")
n_groups <- length(group_levels)
t_half <- function(df_i, se_i) qt(1 - (1 - opt$conf_level) / 2, df_i) * se_i      # half width of a confidence interval

# One row per animal for a fit: y is the value, or the bound when no well was detected (exact = FALSE).
# With one_well_as_bound an animal with some but not all wells detected is a bound as well: each of its other wells is
# above the cycle limit, so the mean of its wells is above the mean of the detected Cq and the limit.
model_data <- function(gene_i, set_i, normalization = "delta Cq", one_well_as_bound = FALSE) {
  d <- dq[gene == gene_i & animal_id %chin% in_set[[set_i]] & !is.na(reference_cq)]
  shift <- if (normalization == "delta Cq") d$reference_cq else 0
  partly <- one_well_as_bound & !is.na(d$n_detected) & d$n_detected > 0 & d$n_detected < design$wells_per_animal
  mean_with_limit <- (d$cq * d$n_detected + (design$wells_per_animal - d$n_detected) * opt$max_cycles) / design$wells_per_animal
  data.frame(animal_id = d$animal_id, group = d$group, exact = !d$censored & !partly,
             y = ifelse(d$censored, opt$max_cycles, ifelse(partly, mean_with_limit, d$cq)) - shift, x = d$reference_cq)
}
# Group means and their covariance matrix. With a covariate the means are adjusted to the mean of x.
fit_groups <- function(d, method, with_covariate = FALSE) {
  X <- model.matrix(~ 0 + group, d); colnames(X) <- group_levels
  if (with_covariate) X <- cbind(X, slope = d$x - mean(d$x))
  k <- ncol(X)
  if (method == "complete cases") {
    use <- d$exact; n <- sum(use)
    f <- lm.fit(X[use, , drop = FALSE], d$y[use])
    if (f$rank < k || n <= k) stop("Model not estimable: too few animals with a value in a group.")
    V <- chol2inv(f$qr$qr[seq_len(k), seq_len(k), drop = FALSE]) * sum(f$residuals^2) / (n - k)
    list(b = unname(f$coefficients), V = V, df = n - k, n = n, n_censored = 0L, residual_sd = sqrt(sum(f$residuals^2) / (n - k)),
         residuals = f$residuals, group = d$group[use])
  } else {
    n <- nrow(d)
    f <- tryCatch(survreg(Surv(d$y, d$exact) ~ 0 + X, dist = "gaussian"),
                  warning = function(w) stop("Censored regression did not fit cleanly: ", conditionMessage(w), call. = FALSE))
    list(b = unname(coef(f)), V = f$var[seq_len(k), seq_len(k), drop = FALSE] * n / (n - k), df = n - k, n = n, n_censored = sum(!d$exact),
         residual_sd = f$scale * sqrt(n / (n - k)))
  }
}
# What goes into a fit: the animals used, their values, and which values are bounds. Two fits of a gene with the same
# content are one fit; without a bound the censored regression is the linear model.
fit_content <- function(d, method) {
  used <- if (method == "complete cases") d[d$exact, ] else d
  used <- used[order(used$animal_id), ]
  paste(used$animal_id, used$y, used$exact, collapse = ";")
}
contrast_table <- function(f) {
  est <- drop(K %*% f$b[seq_len(n_groups)])
  se <- sqrt(diag(K %*% f$V[seq_len(n_groups), seq_len(n_groups)] %*% t(K)))
  data.table(term = term_info$term, contrast = term_info$contrast, n_animals = f$n, n_censored = f$n_censored,
             difference_cycles = est, se_cycles = se, ci_low_cycles = est - t_half(f$df, se), ci_high_cycles = est + t_half(f$df, se),
             fold_change = 2^-est, fold_low = 2^-(est + t_half(f$df, se)), fold_high = 2^-(est - t_half(f$df, se)),
             t_value = est / se, df_residual = f$df, p_value = 2 * pt(-abs(est / se), f$df))
}
# Brown and Forsythe (1974, Journal of the American Statistical Association 69:364): one-way ANOVA of the absolute
# deviations from the group medians.
brown_forsythe_p <- function(y, g) anova(lm(abs(y - ave(y, g, FUN = median)) ~ g))[1, "Pr(>F)"]
# Student and Welch tests of a simple effect on its two groups alone, and both tests again with each animal of the two
# groups left out in turn. The tests need two animals with a value in each group, the leave-one-out runs three.
two_group_tests <- function(d, term_i) {
  w <- K[term_i, ]
  used <- d[d$exact & d$group %in% names(w)[w != 0], ]
  in_first <- used$group == names(w)[w > 0]
  both_tests <- function(keep) c(t.test(used$y[in_first & keep], used$y[!in_first & keep], var.equal = TRUE)$p.value,
                                 t.test(used$y[in_first & keep], used$y[!in_first & keep])$p.value)
  smallest_group <- min(sum(in_first), sum(!in_first))
  p_full <- if (smallest_group >= 2) both_tests(TRUE) else c(NA_real_, NA_real_)
  p_left_out <- if (smallest_group >= 3) vapply(seq_len(nrow(used)), function(j) both_tests(seq_len(nrow(used)) != j), numeric(2)) else NULL
  range_of <- function(i) if (is.null(p_left_out)) list(NA_real_, NA_real_, NA_integer_) else
    list(min(p_left_out[i, ]), max(p_left_out[i, ]), sum(p_left_out[i, ] < opt$alpha))
  student <- range_of(1); welch <- range_of(2)
  data.table(n_two_groups = nrow(used), p_two_groups_student = p_full[1], p_two_groups_welch = p_full[2],
             p_student_leave_one_out_low = student[[1]], p_student_leave_one_out_high = student[[2]],
             n_student_leave_one_out_p_below_alpha = student[[3]],
             p_welch_leave_one_out_low = welch[[1]], p_welch_leave_one_out_high = welch[[2]], n_welch_leave_one_out_p_below_alpha = welch[[3]])
}
# The fit repeated with one animal left out at a time: lowest and highest p of each contrast, the animal left out, and
# the number of these fits with p below alpha.
leave_one_out <- function(d, method) {
  rows_i <- if (method == "complete cases") which(d$exact) else seq_len(nrow(d))
  p <- vapply(rows_i, function(j) contrast_table(fit_groups(d[-j, ], method))$p_value, numeric(nrow(term_info)))
  data.table(p_leave_one_out_low = apply(p, 1, min), left_out_for_low = d$animal_id[rows_i[apply(p, 1, which.min)]],
             p_leave_one_out_high = apply(p, 1, max), left_out_for_high = d$animal_id[rows_i[apply(p, 1, which.max)]],
             n_leave_one_out = length(rows_i), n_leave_one_out_p_below_alpha = as.integer(rowSums(p < opt$alpha)))
}

models <- list(); group_fit <- list(); primary_sd <- list()
for (set_i in names(animal_sets)) for (gene_i in targets) {
  d_i <- model_data(gene_i, set_i)
  data_for <- setNames(list(d_i, d_i, model_data(gene_i, set_i, one_well_as_bound = TRUE)), methods)
  id_cols <- data.table(gene = gene_i, replicate_data = gene_i %in% targets_rep, detection_class = class_of(gene_i, set_i),
                        animals_included = animal_sets[[set_i]], n_animals_without_value = sum(!d_i$exact))
  for (m_i in methods) {
    d_m <- data_for[[m_i]]
    # A fit needs two animals with a value in every group; a gene that is not analysed gets empty rows as well.
    if (id_cols$detection_class == "not analysed" || min(table(d_m$group[d_m$exact])) < 2) {
      models[[length(models) + 1]] <- cbind(id_cols, analysis = m_i, term_info)
      next
    }
    fit_m <- fit_groups(d_m, m_i)
    rows <- cbind(id_cols, analysis = m_i, contrast_table(fit_m), leave_one_out(d_m, m_i), fit_content = fit_content(d_m, m_i))
    if (m_i == "complete cases") {
      rows[, `:=`(residual_shapiro_p = shapiro.test(fit_m$residuals)$p.value,
                  brown_forsythe_p = brown_forsythe_p(fit_m$residuals, fit_m$group))]
      two_groups <- rbindlist(lapply(simple_terms, function(term_i) two_group_tests(d_m, term_i)))
      rows[match(simple_terms, term), (names(two_groups)) := two_groups]
    }
    if (m_i == methods[2] && set_i == "clean") primary_sd[[gene_i]] <- fit_m$residual_sd
    models[[length(models) + 1]] <- rows
    group_fit[[length(group_fit) + 1]] <- data.table(
      gene = gene_i, animals_included = animal_sets[[set_i]], analysis = m_i, group = group_levels,
      mean_delta_cq = fit_m$b[seq_len(n_groups)], se_mean = sqrt(diag(fit_m$V)[seq_len(n_groups)]))
  }
}
models <- rbindlist(models, fill = TRUE); group_fit <- rbindlist(group_fit)
if (!"p_value" %in% names(models) || length(primary_sd) == 0) {
  stop("No gene could be fitted in the primary analysis: every analysed gene has a group with fewer than two animals with a value.")
}
bh <- function(p) { q <- rep(NA_real_, length(p)); q[!is.na(p)] <- p.adjust(p[!is.na(p)], method = "BH"); q }
models[, `:=`(q_value = bh(p_value), q_two_groups_student = bh(p_two_groups_student), q_two_groups_welch = bh(p_two_groups_welch),
              genes_in_correction = sum(!is.na(p_value))), by = .(analysis, animals_included, term)]
models[, primary := analysis == methods[2] & animals_included == animal_sets[["clean"]]]
models <- models[order(match(animals_included, animal_sets), match(analysis, methods), match(gene, targets), match(term, term_info$term))]
# A fit with the content of an earlier fit of the same gene is named after that fit. Such fits must give the same numbers;
# this includes the linear model against the censored regression where nothing is censored.
models[, fit_label := paste(analysis, animals_included, sep = ", ")]
models[!is.na(fit_content), same_fit_as := fit_label[match(fit_content, fit_content)], by = gene]
models[same_fit_as == fit_label, same_fit_as := NA_character_]
repeats <- models[!is.na(fit_content), .(p_gap = diff(range(p_value)), p_largest = max(p_value), estimate_gap = diff(range(difference_cycles)),
                                         se_gap = diff(range(se_cycles))), by = .(gene, term, fit_content)]
if (repeats[, any(p_gap > opt$same_p_tolerance * p_largest | estimate_gap > 1e-6 | se_gap > 1e-6)]) {
  stop("Fits of the same animals and values give different numbers for: ",
       paste(repeats[p_gap > opt$same_p_tolerance * p_largest | estimate_gap > 1e-6 | se_gap > 1e-6, unique(gene)], collapse = ", "))
}
models[, c("fit_content", "fit_label") := NULL]
setcolorder(models, c("gene", "replicate_data", "detection_class", "analysis", "animals_included", "primary", "same_fit_as", "term", "contrast",
                      "n_animals", "n_censored", "n_animals_without_value", "difference_cycles", "se_cycles", "ci_low_cycles", "ci_high_cycles",
                      "fold_change", "fold_low", "fold_high", "t_value", "df_residual", "p_value", "q_value", "genes_in_correction",
                      "p_leave_one_out_low", "left_out_for_low", "p_leave_one_out_high", "left_out_for_high", "n_leave_one_out",
                      "n_leave_one_out_p_below_alpha", "n_two_groups", "p_two_groups_student", "q_two_groups_student",
                      "p_student_leave_one_out_low", "p_student_leave_one_out_high", "n_student_leave_one_out_p_below_alpha",
                      "p_two_groups_welch", "q_two_groups_welch", "p_welch_leave_one_out_low", "p_welch_leave_one_out_high",
                      "n_welch_leave_one_out_p_below_alpha", "residual_shapiro_p", "brown_forsythe_p"))
write_tab(models, "qpcr_models.csv")
# Analyses of the same animals and values give the same p value. Counts of analyses count such repeats once:
# first_result() is TRUE for a p value that equals none of the values before it, within the tolerance.
first_result <- function(p) vapply(seq_along(p), function(i) {
  earlier <- p[seq_len(i - 1)]
  !any(abs(earlier - p[i]) <= opt$same_p_tolerance * pmax(earlier, p[i]))
}, logical(1))
distinct_p <- function(p) { p <- p[!is.na(p)]; p[first_result(p)] }
pick <- function(gene_i, term_i, analysis_i = methods[2], set_i = "clean") {
  r <- models[gene == gene_i & term == term_i & analysis == analysis_i & animals_included == animal_sets[[set_i]]]
  stopifnot(nrow(r) == 1); r
}
smallest_q <- vapply(model_terms, function(term_i) models[primary == TRUE & term == term_i, min(q_value, na.rm = TRUE)], numeric(1))
log_msg("delta Cq models: %d genes analysed; primary analysis (censored regression, %s): smallest q per term: %s",
        models[primary == TRUE & term == "ploidy", sum(!is.na(p_value))], animal_sets[["clean"]],
        paste(sprintf("%s %.2f", model_terms, smallest_q), collapse = ", "))

# ---- how often the test of the censored regression rejects a true null hypothesis (simulation) ----
# The small-sample correction in fit_groups() is exact when nothing is censored. Here it is checked with censoring.
# One delta Cq per animal of the primary analysis is drawn from one normal distribution, so that no group differs. A
# value above the animal's bound (max_cycles minus its reference Cq) becomes a bound, as for a non-detect, and the data
# set is fitted and tested as in the analysis. The standard deviation is the median residual SD of the primary fits.
# The mean is set so that the expected share of censored animals rises in steps of 0.1 to the largest share among the
# fits above. A test is taken to hold its level when it rejects in 0.5 to 1.5 times the nominal share of data sets
# (Bradley 1978, British Journal of Mathematical and Statistical Psychology 31:144, the liberal criterion). The step
# stops if the corrected test rejects more often than that. A test that rejects less often is conservative: its p
# values are too large, not too small, and the table says so.
censored_share <- models[!is.na(p_value) & n_censored > 0, n_censored / n_animals]
test_size <- data.table(share_censored_expected = numeric(0), mean_delta_cq = numeric(0), sd_delta_cq = numeric(0), n_animals = integer(0),
                        data_sets = integer(0), data_sets_fitted = integer(0), share_censored_simulated = numeric(0), term = character(0),
                        contrast = character(0), nominal_level = numeric(0), rejection_rate = numeric(0),
                        rejection_rate_uncorrected = numeric(0), monte_carlo_se = numeric(0), level_check = character(0))
if (length(censored_share) > 0) {
  sim <- animals[match(in_set$clean, animal_id), .(group, bound = opt$max_cycles - reference_cq)]
  sd_null <- median(unlist(primary_sd))
  steps <- seq(0.1, 0.9, by = 0.1)
  shares <- c(steps[steps < max(censored_share) - 1e-9], max(censored_share))
  n_terms <- length(model_terms)
  set.seed(opt$size_check_seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  one_data_set <- function(mean_null) {       # p of the three model terms, corrected then uncorrected, and the share censored
    drawn <- rnorm(nrow(sim), mean_null, sd_null)
    exact <- drawn < sim$bound
    not_fitted <- c(rep(NA_real_, 2 * n_terms), mean(!exact))
    if (min(table(sim$group[exact])) < 2) return(not_fitted)                 # the step does not fit such a gene
    f <- tryCatch(fit_groups(data.frame(group = sim$group, y = pmin(drawn, sim$bound), exact = exact), methods[2]), error = function(e) NULL)
    if (is.null(f)) return(not_fitted)
    ct <- contrast_table(f)[match(model_terms, term)]
    c(ct$p_value, 2 * pnorm(-abs(ct$t_value) * sqrt(f$n / f$df)), mean(!exact))    # without the correction: z test, maximum likelihood SE
  }
  test_size <- rbindlist(lapply(shares, function(share_i) {
    mean_null <- uniroot(function(m) mean(pnorm(sim$bound, m, sd_null, lower.tail = FALSE)) - share_i,
                         range(sim$bound) + c(-10, 10) * sd_null, tol = 1e-9)$root
    runs <- t(replicate(opt$size_check_runs, one_data_set(mean_null)))
    fitted <- !is.na(runs[, 1])
    if (!any(fitted)) stop("Check of the censored test: no simulated data set could be fitted at a censored share of ", share_i, ".")
    rate <- function(columns) colMeans(runs[fitted, columns, drop = FALSE] < opt$alpha)
    data.table(share_censored_expected = share_i, mean_delta_cq = mean_null, sd_delta_cq = sd_null, n_animals = nrow(sim),
               data_sets = nrow(runs), data_sets_fitted = sum(fitted), share_censored_simulated = mean(runs[, 2 * n_terms + 1]),
               term = model_terms, contrast = contrast_of(model_terms), nominal_level = opt$alpha,
               rejection_rate = rate(seq_len(n_terms)), rejection_rate_uncorrected = rate(n_terms + seq_len(n_terms)))
  }))
  test_size[, monte_carlo_se := sqrt(rejection_rate * (1 - rejection_rate) / data_sets_fitted)]
  test_size[, level_check := fifelse(rejection_rate > 1.5 * nominal_level, "rejects too often",
                             fifelse(rejection_rate < 0.5 * nominal_level, "rejects too seldom (conservative)",
                                     "within 0.5 to 1.5 times the nominal level"))]
  log_msg(paste("censored regression under the null, %d data sets per share of censored animals (%s; %d to %d of them fitted):",
                "the corrected test rejects at %g in %.1f%% to %.1f%%, the uncorrected test in %.1f%% to %.1f%%"),
          opt$size_check_runs, paste(sprintf("%.2f", shares), collapse = ", "), min(test_size$data_sets_fitted), max(test_size$data_sets_fitted),
          opt$alpha, 100 * min(test_size$rejection_rate), 100 * max(test_size$rejection_rate),
          100 * min(test_size$rejection_rate_uncorrected), 100 * max(test_size$rejection_rate_uncorrected))
  if (any(test_size$level_check == "rejects too seldom (conservative)")) {
    log_msg("the corrected test is conservative at a censored share of %s",
            paste(sprintf("%.2f", test_size[level_check == "rejects too seldom (conservative)", unique(share_censored_expected)]), collapse = ", "))
  }
} else log_msg("no fit has a censored animal: the censored regression is the linear model throughout, and no simulation is needed")
write_tab(test_size, "qpcr_censored_test_size.csv")
if (any(test_size$level_check == "rejects too often")) {
  stop("The corrected test of the censored regression rejects a true null hypothesis too often (qpcr_censored_test_size.csv):\n",
       paste(capture.output(print(test_size[level_check == "rejects too often", .(share_censored_expected, term, rejection_rate)])), collapse = "\n"))
}

# ---- relative expression by group ----
group_summary <- list()
for (set_i in names(animal_sets)) for (variant_i in names(legacy_variants)) {
  d_v <- dq[animal_id %chin% in_set[[set_i]] & !is.na(reference_cq)]
  d_v[, value := if (variant_i == "substituted") fifelse(censored, opt$legacy_nondetect - reference_cq, delta_cq) else delta_cq]
  group_summary[[length(group_summary) + 1]] <- d_v[, .(
    animals_included = animal_sets[[set_i]], nondetects = legacy_variants[[variant_i]], n_animals = .N, n_animals_with_value = sum(!censored),
    mean_delta_cq = mean(value, na.rm = TRUE), sd_delta_cq = sd(value, na.rm = TRUE), geometric_mean_expression = 2^-mean(value, na.rm = TRUE),
    median_expression = median(2^-value, na.rm = TRUE), arithmetic_mean_expression = mean(2^-value, na.rm = TRUE)), by = .(gene, group)]
}
group_n <- rbindlist(lapply(names(animal_sets), function(set_i) {
  dq[animal_id %chin% in_set[[set_i]] & !is.na(reference_cq),
     .(animals_included = animal_sets[[set_i]], n_animals = .N, n_animals_with_value = sum(!censored)),
     by = .(gene = as.character(gene), group = as.character(group))]
}))
censored_part <- merge(group_fit[analysis == "censored regression"], group_n, by = c("gene", "group", "animals_included"), sort = FALSE)
stopifnot(nrow(censored_part) == sum(group_fit$analysis == "censored regression"))
group_summary <- rbind(rbindlist(group_summary),
                       censored_part[, .(gene, group, animals_included, nondetects = "censored regression", n_animals, n_animals_with_value,
                                         mean_delta_cq, se_mean_delta_cq = se_mean, geometric_mean_expression = 2^-mean_delta_cq)], fill = TRUE)
num_cols <- names(group_summary)[vapply(group_summary, is.double, logical(1))]
group_summary[, (num_cols) := lapply(.SD, function(v) fifelse(is.nan(v), NA_real_, v)), .SDcols = num_cols]
group_summary <- group_summary[order(match(animals_included, animal_sets), match(nondetects, c(legacy_variants, "censored regression")),
                                     match(gene, targets), match(group, group_levels))]
setcolorder(group_summary, c("gene", "group", "animals_included", "nondetects", "n_animals", "n_animals_with_value", "mean_delta_cq", "sd_delta_cq",
                             "se_mean_delta_cq"))
write_tab(group_summary, "qpcr_group_summary.csv")

# ---- every analysis side by side, for the three model terms ----
wide_p <- dcast(models[term %chin% model_terms], gene + term ~ analysis + animals_included, value.var = "p_value")
setnames(wide_p, c(outer(methods, animal_sets, paste, sep = "_")),
         c("p_complete_all", "p_censored_all", "p_censored_bounds_all",
           "p_complete_without_failed", "p_censored_without_failed", "p_censored_bounds_without_failed"))
legacy_p <- dcast(legacy_anova[gene %chin% legacy_targets], gene + term ~ nondetects, value.var = "p_value")
setnames(legacy_p, unname(legacy_variants), c("p_legacy_nondetects_dropped", "p_legacy_nondetects_substituted"))
sensitivity <- merge(legacy_p, wide_p, by = c("gene", "term"), all.y = TRUE)       # every target; legacy p where the legacy input holds the gene
sensitivity <- merge(sensitivity, models[primary == TRUE & term %chin% model_terms,
                                         .(gene, term, detection_class, fold_change_primary = fold_change, fold_low_primary = fold_low,
                                           fold_high_primary = fold_high, q_value_primary = q_value,
                                           p_primary_leave_one_out_low = p_leave_one_out_low, p_primary_leave_one_out_high = p_leave_one_out_high)],
                     by = c("gene", "term"))
stopifnot(nrow(sensitivity) == length(targets) * length(model_terms))
p_cols <- c("p_legacy_nondetects_substituted", "p_legacy_nondetects_dropped", "p_complete_all", "p_complete_without_failed", "p_censored_all",
            "p_censored_without_failed", "p_censored_bounds_all", "p_censored_bounds_without_failed")
# Several of the eight analyses can be one fit of the same animals (same_fit_as in qpcr_models.csv; the two legacy variants
# of a gene without a missing value). The counts below are counts of results that differ. Every analysis of this table
# normalizes to the reference gene; qpcr_normalization_check.csv holds the other normalizations.
p_by_row <- as.matrix(sensitivity[, ..p_cols])
sensitivity[, analyses_run := rowSums(!is.na(p_by_row))]
sensitivity[, distinct_results := apply(p_by_row, 1, function(p) length(distinct_p(p)))]
sensitivity[, distinct_results_p_below_alpha := apply(p_by_row, 1, function(p) sum(distinct_p(p) < opt$alpha))]
of_this_table <- sprintf("distinct results of this table (all relative to %s)", tolower(ref))
sensitivity[, reading := fifelse(detection_class == "not analysed", "legacy analyses only: the gene is not analysed here",
  fifelse(distinct_results_p_below_alpha == 0, sprintf("p >= %g in each of the %d %s", opt$alpha, distinct_results, of_this_table),
  fifelse(distinct_results_p_below_alpha == distinct_results, sprintf("p < %g in each of the %d %s", opt$alpha, distinct_results, of_this_table),
          sprintf("p < %g in %d of the %d %s: the conclusion depends on the analysis", opt$alpha, distinct_results_p_below_alpha,
                  distinct_results, of_this_table))))]
sensitivity <- sensitivity[order(match(gene, targets), match(term, model_terms))]
setcolorder(sensitivity, c("gene", "term", "detection_class", p_cols))
write_tab(sensitivity, "qpcr_sensitivity.csv")
cat("Gene and term with p <", opt$alpha, "in at least one analysis:\n")
print(sensitivity[distinct_results_p_below_alpha > 0, c("gene", "term", p_cols, "q_value_primary"), with = FALSE][
  , lapply(.SD, function(v) if (is.double(v)) round(v, 3) else v)])

# =============================================================================
# 7. Reference gene: actin by group, and what normalization does to each result
# =============================================================================
ref_rows <- list()
describe_ref <- function(d) d[, .(n = .N, mean_cq = mean(reference_cq), sd_cq = sd(reference_cq), median_cq = median(reference_cq),
                                  min_cq = min(reference_cq), max_cq = max(reference_cq))]
for (set_i in names(animal_sets)) {
  a_i <- animals[animal_id %chin% in_set[[set_i]] & !is.na(reference_cq)]
  ref_rows[[length(ref_rows) + 1]] <- cbind(section = sprintf("1 %s Cq by group", ref), animals_included = animal_sets[[set_i]],
    rbind(a_i[, describe_ref(.SD), by = .(item = as.character(group))][match(group_levels, item)], cbind(item = "all groups", describe_ref(a_i))))
  ct <- contrast_table(fit_groups(data.frame(group = a_i$group, y = a_i$reference_cq, exact = TRUE), "complete cases"))[term %chin% model_terms]
  ref_rows[[length(ref_rows) + 1]] <- ct[, .(section = sprintf("2 difference in %s Cq between groups", ref), animals_included = animal_sets[[set_i]],
                                             item = contrast, n = n_animals, difference_cycles, ci_low = ci_low_cycles, ci_high = ci_high_cycles,
                                             template_ratio = fold_change, t_value, df_residual, p_value)]
}
ref_rows[[length(ref_rows) + 1]] <- if (length(failed_ids) > 0) {
  failed_rows <- merge(
    animals[match(failed_ids, animal_id), .(section = "3 failed samples", item = animal_id, n = 1L, mean_cq = reference_cq,
                                            animals_included = fifelse(is.na(reference_cq), "in no model of delta Cq", animal_sets[["all"]]))],
    dq[animal_id %chin% failed_ids & replicate_data == TRUE, .(targets_with_no_detected_well = sum(detection == "no well detected"),
                                                               targets_measured = .N), by = .(item = animal_id)],
    by = "item", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(failed_rows) == length(failed_ids))
  failed_rows
} else data.table(section = "3 failed samples", animals_included = animal_sets[["all"]], item = "none", n = 0L)
reference_gene <- rbindlist(ref_rows, fill = TRUE)
reference_gene[section == "3 failed samples",
               `:=`(outer_fence_cq = ref_fence, first_quartile_cq = ref_quartiles[1], third_quartile_cq = ref_quartiles[2])]
reference_gene[, what_it_is := fifelse(
  grepl("^1", section), sprintf("%s Cq per animal (mean of the detected wells); a later Cq means less %s template in the well", ref, ref),
  fifelse(grepl("^2", section), paste("contrast of group means of the", ref, "Cq, in cycles; template_ratio is 2^-difference, the ratio of", ref,
                                      "template that the difference implies; every delta Cq inherits this difference with the opposite sign"),
  fifelse(item %chin% no_reference, sprintf("animal without a detected %s well: it has no delta Cq", ref),
          sprintf("animal whose %s Cq (in mean_cq) is above the outer fence: third quartile + %g interquartile ranges of the animals in the models",
                  ref, opt$failed_fence_iqr))))]
setcolorder(reference_gene, c("section", "animals_included", "item"))
write_tab(reference_gene, "qpcr_reference_gene.csv")
cat(sprintf("Difference in %s Cq between groups (cycles):\n", ref))
print(reference_gene[grepl("^2", section), .(animals_included, item, difference_cycles = round(difference_cycles, 2), ci_low = round(ci_low, 2),
                                             ci_high = round(ci_high, 2), p_value = round(p_value, 3))])

# ---- each effect under three normalizations (censored regression) ----
# All seven contrasts, the simple effects included, so that a claim about one group pair can be checked on that pair.
normalizations <- c("delta Cq", "raw Cq", sprintf("raw Cq with the %s Cq as covariate", ref))
slope_term <- sprintf("slope of the target Cq on the %s Cq", ref)
norm_rows <- list()
for (set_i in names(animal_sets)) for (gene_i in targets) {
  if (class_of(gene_i, set_i) == "not analysed") next
  for (norm_i in normalizations) {
    covariate_i <- norm_i == normalizations[3]
    f_i <- fit_groups(model_data(gene_i, set_i, norm_i), "censored regression", with_covariate = covariate_i)
    rows <- contrast_table(f_i)[, .(term, contrast, null_value = 0, estimate = difference_cycles, ci_low = ci_low_cycles,
                                    ci_high = ci_high_cycles, fold_change, p_value)]
    if (covariate_i) {      # the slope is 1 when the target Cq follows the reference Cq cycle for cycle
      slope <- f_i$b[n_groups + 1]; slope_se <- sqrt(f_i$V[n_groups + 1, n_groups + 1])
      rows <- rbind(rows, data.table(term = slope_term, contrast = sprintf("cycles of target Cq per cycle of %s Cq, within groups", ref),
                                     null_value = 1, estimate = slope, ci_low = slope - t_half(f_i$df, slope_se),
                                     ci_high = slope + t_half(f_i$df, slope_se), fold_change = NA_real_,
                                     p_value = 2 * pt(-abs((slope - 1) / slope_se), f_i$df)))
    }
    norm_rows[[length(norm_rows) + 1]] <- cbind(gene = gene_i, animals_included = animal_sets[[set_i]], detection_class = class_of(gene_i, set_i),
                                                normalization = norm_i, rows, n_animals = f_i$n, n_censored = f_i$n_censored)
  }
}
normalization_check <- rbindlist(norm_rows)
write_tab(normalization_check, "qpcr_normalization_check.csv")
slopes <- normalization_check[term == slope_term & animals_included == animal_sets[["clean"]]]
log_msg("slope of the target Cq on the %s Cq (%s; 1 if %s tracks the template): %s", ref, animal_sets[["clean"]], ref,
        paste(sprintf("%s %.2f", slopes$gene, slopes$estimate), collapse = ", "))

# =============================================================================
# 8. COX1: a raw Cq without a reference
# =============================================================================
mito_cq <- if (mito %in% rep_genes) per_gene[gene == mito, .(animal_id, cq)] else cons_long[gene == mito, .(animal_id, cq = cq_consolidated)]
mito_d <- merge(animals[in_models == TRUE], mito_cq, by = "animal_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(mito_d) == sum(animals$in_models))
if (anyNA(mito_d$cq)) {
  log_msg("%s: no value for %d of the %d animals in the models (%s); they are in no %s test", mito, sum(is.na(mito_d$cq)), nrow(mito_d),
          paste(mito_d[is.na(cq), animal_id], collapse = ", "), mito)
}
mito_d <- mito_d[!is.na(cq)]
mito_frame <- data.frame(group = mito_d$group, y = mito_d$cq, exact = TRUE)
describe_mito <- function(d) d[, .(n = .N, mean_cq = mean(cq), sd_cq = sd(cq), median_cq = median(cq), min_cq = min(cq), max_cq = max(cq))]
mito_ct <- contrast_table(fit_groups(mito_frame, "complete cases"))
mito_two <- rbindlist(lapply(simple_terms, function(term_i) two_group_tests(mito_frame, term_i)))
mito_tukey <- legacy_tukey[gene == mito & response == "ln(2^-Cq)" & nondetects == legacy_variants[["dropped"]] &
                             legacy_term == legacy_terms[["interaction"]]]
mito_r <- cor.test(mito_d$cq, mito_d$reference_cq)
r_with_ref <- function(v_cq, v_ref) { ok <- !is.na(v_cq) & !is.na(v_ref); if (sum(ok) > 2) cor(v_cq[ok], v_ref[ok]) else NA_real_ }
analysed_all <- targets[vapply(targets, class_of, "", set_i = "all") != "not analysed"]
other_r <- dq[animal_id %chin% in_set$all & gene %in% analysed_all, .(r = r_with_ref(cq, reference_cq)), by = gene]
what <- list(       # what each kind of row is
  by_group = sprintf("%s Cq per animal, not normalized; a lower Cq means more %s template in the well", mito, mito),
  contrast = paste("contrast of group means of the raw Cq, in cycles; the ratio is 2^-estimate. The legacy ANOVA of ln(2^-Cq) has the same F and p,",
                   "because that response is a linear function of the Cq. Simple effects also carry Student and Welch tests on the two groups alone"),
  tukey = "legacy TukeyHSD() on ln(2^-Cq)",
  r_mito = sprintf(paste("across animals; the Cq of the genes measured on the same template as %s follows the %s Cq, the %s Cq does not;",
                         "this points to another template preparation for %s (an inference, no record exists)"), ref, ref, mito, mito),
  r_targets = "estimate is the median over genes; ci_low and ci_high hold the lowest and the highest correlation; n is the number of genes",
  abundance = sprintf("no nuclear gene was measured on the template used for %s, so no ratio of mitochondrial to nuclear copies can be formed", mito))
cox1 <- rbind(
  cbind(quantity = "raw Cq by group", what_it_is = what$by_group,
        rbind(mito_d[, describe_mito(.SD), by = .(level = as.character(group))][match(group_levels, level)],
              cbind(level = "all groups", describe_mito(mito_d)))),
  mito_ct[, .(quantity = "difference in raw Cq, two-way model", what_it_is = what$contrast, level = contrast, n = n_animals,
              estimate = difference_cycles, ci_low = ci_low_cycles, ci_high = ci_high_cycles, ratio_2_to_minus_estimate = fold_change,
              F_value = t_value^2, df_residual, p_value,
              p_two_groups_student = c(rep(NA_real_, length(model_terms)), mito_two$p_two_groups_student),
              p_two_groups_welch = c(rep(NA_real_, length(model_terms)), mito_two$p_two_groups_welch))],
  if (nrow(mito_tukey) > 0) data.table(
    quantity = "smallest Tukey-adjusted p among the six pairs of groups", what_it_is = what$tukey,
    level = mito_tukey[which.min(p_adjusted), comparison], n = nrow(mito_d), p_value = mito_tukey[, min(p_adjusted)]),
  data.table(quantity = sprintf("Pearson correlation of the %s Cq with the %s Cq", mito, ref), what_it_is = what$r_mito, level = "all groups",
             n = sum(!is.na(mito_d$reference_cq)), estimate = unname(mito_r$estimate), ci_low = mito_r$conf.int[1], ci_high = mito_r$conf.int[2],
             p_value = mito_r$p.value),
  data.table(quantity = sprintf("Pearson correlation of the analysed genes with the %s Cq: median, lowest, highest", ref),
             what_it_is = what$r_targets, level = "all groups", n = nrow(other_r), estimate = median(other_r$r), ci_low = min(other_r$r),
             ci_high = max(other_r$r)),
  data.table(quantity = "mitochondrial abundance relative to the nuclear genome", what_it_is = what$abundance, level = "not estimable", n = 0L),
  fill = TRUE)
setcolorder(cox1, c("quantity", setdiff(names(cox1), c("quantity", "what_it_is")), "what_it_is"))
write_tab(cox1, "qpcr_cox1.csv")
log_msg("%s raw Cq: %s; interaction %s; correlation with the %s Cq r = %.3f (analysed targets: median r = %.2f)", mito,
        paste(sprintf("%s %.2f", group_levels, cox1[quantity == "raw Cq by group"][match(group_levels, level), mean_cq]), collapse = ", "),
        fmt_p(cox1[level == contrast_of("interaction"), p_value]), ref, unname(mito_r$estimate), median(other_r$r))

# =============================================================================
# 9. The qPCR statements of the draft, beside the numbers
# =============================================================================
# Each sentence is paired with the claim it makes about these data. Numbers and verdicts are built from the tables
# above; the rule behind each verdict is written out where the claim is evaluated.
hsp <- "HSP90"                                   # the one expression gene the draft names
# Counts of tests in the text are counts of results that differ (distinct_p): a repeat of one fit is not counted again.
of_n <- function(p) { p <- distinct_p(p); sprintf("%d of %d", sum(p < opt$alpha), length(p)) }
p_span <- function(p) if (diff(range(p)) < 5e-4) fmt_p(p[1]) else sprintf("p = %.3f to %.3f", min(p), max(p))
fold_ci <- function(r) sprintf("%.2f-fold (%s %.2f to %.2f)", r$fold_change, ci_label, r$fold_low, r$fold_high)
with_conditions <- function(verdict, conditions) {
  if (verdict == "partly supported") paste0(verdict, ": ", paste(conditions, collapse = "; ")) else verdict
}
claims <- list()

# ---- HSP90 ----
if (!hsp %in% targets || any(vapply(names(animal_sets), function(set_i) class_of(hsp, set_i), "") == "not analysed")) {
  for (k_i in c("hsp_control", "hsp_heat")) claims[[k_i]] <- list(
    numbers = sprintf("%s is not analysed; see qpcr_nondetects.csv.", hsp), verdict = "cannot be evaluated",
    reasoning = sprintf("%s was detected in too few wells for a model to be fitted.", hsp))
} else {
  h <- function(term_i) models[gene == hsp & term == term_i & !is.na(p_value)]
  h1 <- function(term_i) pick(hsp, term_i)                                    # the primary analysis
  legacy_of <- function(term_i) legacy_anova[gene == hsp & term == term_i]
  legacy_text <- function(d) sprintf("p = %.3f (non-detects %s) and %.3f (non-detects %s)",
                                     d$p_value[1], d$nondetects[1], d$p_value[2], d$nondetects[2])
  norm_of <- function(term_i) {
    r <- normalization_check[gene == hsp & animals_included == animal_sets[["clean"]] & term == term_i]
    stopifnot(nrow(r) == length(normalizations)); r
  }
  by_norm <- function(d) paste(sprintf("%s %.2f-fold (%s)", d$normalization, d$fold_change, fmt_p(d$p_value)), collapse = "; ")
  # The tests of a simple effect on its two groups alone: one row per test and set of animals.
  two <- function(term_i) {
    r <- h(term_i)[analysis == "complete cases"]
    rbind(r[, .(test = "Student", animals_included, n = n_two_groups, p = p_two_groups_student, q = q_two_groups_student,
                loo_low = p_student_leave_one_out_low, loo_high = p_student_leave_one_out_high,
                loo_below = n_student_leave_one_out_p_below_alpha)],
          r[, .(test = "Welch", animals_included, n = n_two_groups, p = p_two_groups_welch, q = q_two_groups_welch,
                loo_low = p_welch_leave_one_out_low, loo_high = p_welch_leave_one_out_high,
                loo_below = n_welch_leave_one_out_p_below_alpha)])
  }
  loo <- function(r) sprintf("p = %.3f to %.3f with one animal left out, below %g in %d of %d fits", r$p_leave_one_out_low,
                             r$p_leave_one_out_high, opt$alpha, r$n_leave_one_out_p_below_alpha, r$n_leave_one_out)
  two_text <- function(t) sprintf("%s %s (q = %.2f; %s)", t$test, fmt_p(t$p), t$q, ifelse(
    is.na(t$loo_high), "no leave-one-out run: a group has two animals",
    sprintf("p = %.3f to %.3f with one animal left out, below %g in %d of %d tests", t$loo_low, t$loo_high, opt$alpha, t$loo_below, t$n)))

  # Claim: HSP90 is higher in triploid than in diploid control animals.
  #   supported         the direction holds in every fit and under each normalization; every test of the contrast has
  #                     p < alpha and q < fdr; and no test of this step loses p < alpha when one animal is left out
  #   partly supported  the direction holds and at least one test has p < alpha
  #   not supported     otherwise
  # "Partly supported" lists each requirement of "supported" that is not met. The leave-one-out condition looks at the
  # tests with p < alpha, which carry the verdict. The legacy tests are reproduced as run and have no leave-one-out.
  cc <- h("ploidy_in_control"); p1 <- h1("ploidy_in_control"); av <- h1("ploidy"); tw <- two("ploidy_in_control")
  tw1 <- tw[animals_included == animal_sets[["clean"]]]
  tk <- legacy_tukey[gene == hsp & legacy_comparison == "T:C-D:C"]            # one row per legacy treatment of non-detects
  nrm <- norm_of("ploidy_in_control"); nrm_av <- norm_of("ploidy")
  p_substituted <- tk[nondetects == legacy_variants[["substituted"]], p_two_groups_student]
  stopifnot(nrow(tk) == length(legacy_variants), nrow(tw1) == 2, length(p_substituted) == 1)
  tests <- rbind(
    cc[, .(family = "fits of the two-way model", label = paste0("fit by ", analysis, ", ", animals_included), p = p_value,
           loo_high = p_leave_one_out_high, loo_below = n_leave_one_out_p_below_alpha, loo_n = n_leave_one_out, legacy = FALSE)],
    tw[, .(family = "tests of the two control groups alone", label = paste(test, "test"), p, loo_high, loo_below, loo_n = n, legacy = FALSE)],
    tk[, .(family = "legacy Tukey tests", label = paste("non-detects", nondetects), p = p_adjusted, loo_high = NA_real_,
           loo_below = NA_integer_, loo_n = NA_integer_, legacy = TRUE)])[!is.na(p)]
  tests[, first_of_its_result := first_result(p), by = family]
  tests <- tests[first_of_its_result == TRUE]                                 # a repeat of one fit or test is one result
  q_all <- c(cc$q_value, tw$q)
  direction_ok <- all(cc$fold_change > 1) && all(nrm$fold_change > 1)
  carrying <- tests[p < opt$alpha & legacy == FALSE]
  stable <- all(!is.na(carrying$loo_high) & carrying$loo_high < opt$alpha)
  verdict <- if (!direction_ok || !any(tests$p < opt$alpha)) "not supported" else
    if (all(tests$p < opt$alpha) && all(q_all < opt$fdr) && stable) "supported" else "partly supported"
  significant_in <- tests[, .(some = any(p < opt$alpha)), by = family][some == TRUE, family]
  loo_condition <- if (any(carrying$loo_high >= opt$alpha, na.rm = TRUE)) "not in every test that leaves one animal out" else
    "leave-one-out could not be run for every test"
  conditions <- c(if (!all(tests$p < opt$alpha)) sprintf("p < %g only in %s", opt$alpha, paste(significant_in, collapse = " and ")),
                  if (!all(q_all < opt$fdr)) sprintf("not after correction across %d genes", p1$genes_in_correction),
                  if (!stable) loo_condition)
  counts <- tests[, .(text = sprintf("%s %s", of_n(p), family[1])), by = family]$text
  claims$hsp_control <- list(
    numbers = paste(
      sprintf("%s, triploid vs diploid, control animals: %s, %s in the two-way model (primary analysis; %s over the %d distinct fits; %s),",
              hsp, fold_ci(p1), fmt_p(p1$p_value), p_span(cc$p_value), length(distinct_p(cc$p_value)), loo(p1)),
      sprintf("q = %.2f across %d genes; the two control groups alone, %d animals: %s;", p1$q_value, p1$genes_in_correction, tw1$n[1],
              paste(two_text(tw1), collapse = ", ")),
      sprintf("the same Student test with the legacy substitution of Cq %g: %s;", opt$legacy_nondetect, fmt_p(p_substituted)),
      sprintf("legacy Tukey p = %.2f (non-detects %s) and %.2f (non-detects %s).",
              tk$p_adjusted[1], tk$nondetects[1], tk$p_adjusted[2], tk$nondetects[2]),
      sprintf("The same contrast by normalization: %s.", by_norm(nrm)),
      sprintf("Ploidy averaged over treatments: %s, %s (%s), q = %.2f; legacy ANOVA %s; by normalization: %s.", fold_ci(av),
              fmt_p(av$p_value), loo(av), av$q_value, legacy_text(legacy_of("ploidy")), by_norm(nrm_av))),
    verdict = with_conditions(verdict, conditions),
    reasoning = paste0(
      sprintf("The primary estimate is %.1f-fold and the direction %s; ", p1$fold_change,
              if (direction_ok) "is the stated one in every fit and under each normalization" else
                "is not the stated one in every fit and under each normalization"),
      sprintf("counting results that differ, p is below %g in %s; ", opt$alpha, paste(counts, collapse = ", ")),
      if (nrow(carrying) > 0) sprintf("with one animal left out, p stays below %g in %s; ", opt$alpha, paste(sprintf(
        "%d of %d runs of the %s", carrying$loo_below, carrying$loo_n, carrying$label), collapse = " and ")) else "",
      sprintf("the smallest q across the %d genes analysed is %.2f.", p1$genes_in_correction, min(q_all))))

  # Claim: heat raises HSP90 in diploids and in triploids.
  #   supported         the direction holds in both ploidies; p < alpha for the heat effect within each ploidy in every
  #                     test; q < fdr; the effect has p < alpha under each normalization; and no fit of the effect
  #                     averaged over the ploidies loses p < alpha when one animal is left out
  #   partly supported  the direction holds in both ploidies and p < alpha for the heat effect averaged over the
  #                     ploidies in at least one analysis
  #   not supported     otherwise
  avg <- h("treatment"); p1 <- h1("treatment"); dip <- h("treatment_in_diploid"); tri <- h("treatment_in_triploid")
  d1 <- h1("treatment_in_diploid"); t1 <- h1("treatment_in_triploid"); nrm <- norm_of("treatment"); lg <- legacy_of("treatment")
  tests_avg <- c(avg$p_value, lg$p_value)
  tests_each <- list(diploids = c(dip$p_value, two("treatment_in_diploid")$p), triploids = c(tri$p_value, two("treatment_in_triploid")$p))
  ref_shift <- reference_gene[grepl("^2", section) & animals_included == animal_sets[["clean"]] & item == contrast_of("treatment")]
  stopifnot(nrow(ref_shift) == 1, !anyNA(tests_avg), !anyNA(unlist(tests_each)))
  robust <- all(nrm$fold_change > 1 & nrm$p_value < opt$alpha)
  stable <- all(avg[p_value < opt$alpha, p_leave_one_out_high] < opt$alpha)
  direction_ok <- all(dip$fold_change > 1) && all(tri$fold_change > 1)
  verdict <- if (!direction_ok || !any(tests_avg < opt$alpha)) "not supported" else
    if (all(unlist(tests_each) < opt$alpha) && all(avg$q_value < opt$fdr) && robust && stable) "supported" else "partly supported"
  conditions <- c(if (!all(unlist(tests_each) < opt$alpha)) "only when the two ploidies are averaged",
                  if (!all(avg$q_value < opt$fdr)) sprintf("not after correction across %d genes", p1$genes_in_correction),
                  if (!robust) sprintf("only if %s is a valid reference for the heat comparison", tolower(ref)),
                  if (!stable) "not in every fit that leaves one animal out")
  claims$hsp_heat <- list(
    numbers = paste(
      sprintf("%s, heat vs control averaged over ploidies: %s, %s (primary analysis; %s over the %d distinct fits; %s), q = %.2f across %d genes;",
              hsp, fold_ci(p1), fmt_p(p1$p_value), p_span(avg$p_value), length(distinct_p(avg$p_value)), loo(p1), p1$q_value,
              p1$genes_in_correction),
      sprintf("legacy ANOVA %s.", legacy_text(lg)),
      sprintf("Within diploids: %s, %s; within triploids: %s, %s; interaction %s.", fold_ci(d1), fmt_p(d1$p_value), fold_ci(t1), fmt_p(t1$p_value),
              fmt_p(h1("interaction")$p_value)),
      sprintf("Heat effect averaged over ploidies, by normalization: %s.", by_norm(nrm)),
      sprintf("%s Cq, heat minus control animals: %.2f cycles (%s %.2f to %.2f, %s).", ref, ref_shift$difference_cycles, ci_label,
              ref_shift$ci_low, ref_shift$ci_high, fmt_p(ref_shift$p_value))),
    verdict = with_conditions(verdict, conditions),
    reasoning = paste(
      sprintf("Relative to %s, and counting results that differ, the heat effect has p < %g in %s analyses of the same data when the",
              tolower(ref), opt$alpha, of_n(tests_avg)),
      sprintf("ploidies are averaged, in %s tests within diploids and in %s within triploids;", of_n(tests_each$diploids),
              of_n(tests_each$triploids)),
      sprintf("with one animal left out the primary analysis has p < %g in %d of %d fits;", opt$alpha, p1$n_leave_one_out_p_below_alpha,
              p1$n_leave_one_out),
      sprintf("heat minus control is %.2f cycles for the raw %s Cq and %.2f cycles for the %s Cq, and the effect %s.",
              nrm[normalization == "raw Cq", estimate], hsp, ref_shift$difference_cycles, tolower(ref),
              if (robust) "holds under each normalization" else sprintf("does not hold under each normalization, so it rests on %s", tolower(ref)))))
}

# ---- COX1 ----
# Claim: heat changes the "relative mitochondrial abundance" of a ploidy in the stated direction. Without a nuclear
# reference on the same template no relative abundance exists, so these data cannot fully support the claim.
#   partly supported  the raw Cq moves in the stated direction with p < alpha in at least one test
#   not supported     otherwise
mito_claim <- function(term_i, more_template) {
  groups_i <- names(K[term_i, ])[order(K[term_i, ])]; groups_i <- groups_i[K[term_i, groups_i] != 0]     # control group, then heat group
  r <- cox1[level == contrast_of(term_i)]; g <- cox1[quantity == "raw Cq by group"][match(groups_i, level)]
  stopifnot(nrow(r) == 1, nrow(g) == 2)
  tests <- c("two-way model" = r$p_value, "Student t test of the two groups" = r$p_two_groups_student,
             "Welch t test of the two groups" = r$p_two_groups_welch)
  direction_ok <- if (more_template) r$estimate < 0 else r$estimate > 0
  list(numbers = paste(
         sprintf("%s raw Cq, %s: mean %.2f (SD %.2f, n = %d); %s: mean %.2f (SD %.2f, n = %d).",
                 mito, g$level[1], g$mean_cq[1], g$sd_cq[1], g$n[1], g$level[2], g$mean_cq[2], g$sd_cq[2], g$n[2]),
         sprintf("Heat minus control: %.2f cycles (%s %.2f to %.2f), which is %.2f-fold in template if every well held the same amount of template;",
                 r$estimate, ci_label, r$ci_low, r$ci_high, r$ratio_2_to_minus_estimate),
         sprintf("%s. Ploidy by treatment interaction: %s.", paste(sprintf("%s %s", names(tests), fmt_p(tests)), collapse = ", "),
                 fmt_p(cox1[level == contrast_of("interaction"), p_value])),
         sprintf("Correlation of the %s Cq with the %s Cq: r = %.2f. Nuclear reference on the same template: none.", mito, tolower(ref),
                 unname(mito_r$estimate))),
       verdict = if (direction_ok && any(tests < opt$alpha)) "partly supported: only as a difference in raw Cq" else "not supported",
       reasoning = paste(
         "Relative mitochondrial abundance cannot be estimated without a nuclear reference on the same template;",
         sprintf("the raw %s Cq of heat animals is %.2f cycles %s that of control animals, which is %s the stated direction,",
                 mito, abs(r$estimate), if (r$estimate < 0) "below" else "above", if (direction_ok) "in" else "against"),
         sprintf("and p is below %g in %s tests.", opt$alpha, of_n(tests))))
}
claims$mito_diploid  <- mito_claim("treatment_in_diploid", more_template = TRUE)
claims$mito_triploid <- mito_claim("treatment_in_triploid", more_template = FALSE)
claim_text <- c(hsp_control = sprintf("%s expression is higher in triploid than in diploid control animals", hsp),
                hsp_heat = sprintf("heat raises %s expression in diploids and in triploids", hsp),
                mito_diploid = "relative mitochondrial abundance rises with heat in diploids",
                mito_triploid = "relative mitochondrial abundance falls with heat in triploids")

# Sentences quoted from the manuscript draft (Google Doc of 2026-11-18), word for word.
quoted <- c(
  d42_s2 = paste0("Our qPCR results complement this finding: relative mitochondrial abundance increased in diploids following stress exposure, ",
                  "suggesting that diploids induce a mitochondrial biogenesis response to meet increasing cellular energy demands ",
                  "(Sokolova et al., 2012)."),
  d42_s3 = paste0("Conversely, relative mitochondrial abundance decreased in triploids following stress exposure, suggesting that triploids are ",
                  "unable to launch an effective mitochondrial biogenesis response and thereby cannot optimally meet the energy demands of the ",
                  "cell under stress conditions."),
  d45_p2_s1 = paste0("Second, our physiological measurements indicate that diploids induce a mitochondrial biogenesis response (increased ",
                     "relative mitochondrial abundance) following stress exposure to meet the increasing energy needs of the cell."),
  d45_p2_s2 = paste0("Conversely, triploids showed decreased relative mitochondrial abundance following stress exposure, suggesting an ",
                     "inability to launch an effective mitochondrial biogenesis response."),
  d45_p2_s3 = paste0("This finding is consistent with prior transcriptomic and metabolic evidence that triploids experience metabolic ",
                     "depression following heatwave exposure (George et al., 2023)."),
  d45_p3_s1 = paste0("Third, qPCR results from companion experiments indicated that triploids expressed HSP90 at high levels even in control ",
                     "conditions, suggesting non-optimal allocation of energy towards maintaining a heightened immune/stress response when ",
                     "not needed."),
  d45_p3_s2 = paste0("While upregulation of HSP90 in both diploids and triploids following heat exposure indicates a functional heat shock ",
                     "response in both ploidy types, the constitutively elevated HSP90 in triploids may represent an energetic cost that ",
                     "limits their capacity to mount additional stress responses."),
  d45_p4_s1 = paste0("Overall, these findings indicate that diploids are able to launch a quicker and potentially more effective stress ",
                     "response (i.e., 5-mC DNA methylation changes, increased mitochondrial abundance) than triploids."),
  d47_item3 = "The constitutively elevated HSP90 expression in triploids suggests potential for pre-emptive stress response pathway activation.",
  c49_p2_s1 = paste0("The lack of global 5-mC methylation change in triploids following stress exposure, combined with their decreased ",
                     "mitochondrial abundance, suggests that triploids may have an impaired capacity for adaptive epigenetic responses ",
                     "compared to diploids."))
# One row per sentence and claim, in the order of the draft. A sentence that only refers back to another sentence is
# listed with the claim of that sentence.
statements <- data.table(
  sentence = c("d42_s2", "d42_s3", "d45_p2_s1", "d45_p2_s2", "d45_p2_s3", "d45_p3_s1", "d45_p3_s2", "d45_p3_s2", "d45_p4_s1", "d47_item3",
               "c49_p2_s1"),
  location = c("Discussion 4.2, paragraph 3, sentence 2", "Discussion 4.2, paragraph 3, sentence 3", "Discussion 4.5, paragraph 2, sentence 1",
               "Discussion 4.5, paragraph 2, sentence 2", "Discussion 4.5, paragraph 2, sentence 3", "Discussion 4.5, paragraph 3, sentence 1",
               "Discussion 4.5, paragraph 3, sentence 2", "Discussion 4.5, paragraph 3, sentence 2", "Discussion 4.5, paragraph 4, sentence 1",
               "Discussion 4.7, third item of the list on summer mortality", "Conclusions 4.9, paragraph 2, sentence 1"),
  claim_key = c("mito_diploid", "mito_triploid", "mito_diploid", "mito_triploid", "mito_triploid", "hsp_control", "hsp_heat", "hsp_control",
                "mito_diploid", "hsp_control", "mito_triploid"))
stopifnot(setequal(statements$sentence, names(quoted)), all(statements$claim_key %in% names(claims)))
from_claim <- function(keys, field) vapply(keys, function(k_i) claims[[k_i]][[field]], "")
statements[, `:=`(statement_id = sprintf("Q%02d", .I), statement = quoted[sentence], claim_tested = claim_text[claim_key],
                  numbers = from_claim(claim_key, "numbers"), verdict = from_claim(claim_key, "verdict"),
                  reasoning = from_claim(claim_key, "reasoning"))]
# Notes: what else a reader of the sentence needs to know.
id_of <- function(sentence_i, claim_i) statements[sentence == sentence_i & claim_key == claim_i, statement_id]
notes <- c(
  d45_p2_s1 = sprintf("The draft calls the data 'physiological measurements'. The qPCR data hold one mitochondrial target, %s.", mito),
  d45_p2_s3 = sprintf("'This finding' is the statement before it (%s). The sentence depends on that statement and reports no qPCR result of its own.",
                      id_of("d45_p2_s2", "mito_triploid")),
  d45_p3_s1 = "The draft says 'companion experiments'. Methods and Results do not describe the qPCR or its animals.",
  d45_p4_s1 = "The sentence also names 5-mC changes, which step 03 covers. The sentence after it builds on this one.",
  c49_p2_s1 = "The sentence also names the global 5-mC result, which step 03 covers. The sentence after it builds on this one.")
statements[, note := fifelse(sentence %chin% names(notes), unname(notes[sentence]), "")]
other_claim <- c(hsp_heat = id_of("d45_p3_s2", "hsp_control"), hsp_control = id_of("d45_p3_s2", "hsp_heat"))
statements[sentence == "d45_p3_s2", note := sprintf("The sentence makes two claims; the other one is %s.", other_claim[claim_key])]
stopifnot(!anyNA(statements$statement), !anyNA(statements$claim_tested), !anyNA(statements$note))
write_tab(statements[, .(statement_id, location, statement, claim_tested, numbers, verdict, reasoning, note)], "qpcr_manuscript_statements.csv")
for (k_i in names(claims)) log_msg("%-13s %s", k_i, claims[[k_i]]$verdict)

# =============================================================================
# 10. Figures
# =============================================================================
# Horizontal position: treatment on the axis, the two ploidies side by side, the animals of a group in a fixed order.
# An animal therefore sits at the same place in every panel, and nothing is random.
x_labels <- c("control", "control\n+ acute", "heat", "heat\n+ acute")
animals[, x_axis := 2 * (treatment == "heat") + acute_heat_shock + 1]
animals[, x_pos := x_axis + fifelse(ploidy == "diploid", -0.21, 0.21) + 0.04 * (rank(animal_number) - (.N + 1) / 2), by = .(group, acute_heat_shock)]
centres <- animals[acute_heat_shock == FALSE, .(x_mid = mean(x_pos)), by = group]
marks <- c("both wells detected" = 16, "one of two wells detected" = 1, "no well detected: expression is below the mark" = 6)
ring <- function(d) geom_point(data = d, shape = 21, size = 4.4, stroke = 0.45, colour = "black", fill = NA)
bar <- function(d) geom_segment(data = d, aes(x = x_mid - 0.19, xend = x_mid + 0.19, y = y, yend = y), linewidth = 0.75, colour = "black",
                                inherit.aes = FALSE)

# ---- relative expression, every gene with replicate-level data ----
# An animal without a reference Cq has no relative expression and no place in panel A of the second figure.
not_shown <- if (length(no_reference) == 0) "" else
  sprintf("Not shown: %s (no %s Cq). ", paste(no_reference, collapse = ", "), tolower(ref))
fig_rows <- dq[replicate_data == TRUE & !is.na(reference_cq)]
fig <- merge(fig_rows, animals[, .(animal_id, x_pos)], by = "animal_id", sort = FALSE)
stopifnot(nrow(fig) == nrow(fig_rows))
fig[, y := fifelse(censored, relative_expression_below, relative_expression)]
fig[, mark := factor(names(marks)[fifelse(n_detected == design$wells_per_animal, 1L, fifelse(n_detected > 0, 2L, 3L))], levels = names(marks))]
class_levels <- c("quantifiable", "near the detection limit", "not analysed")
panel_class <- vapply(targets_rep, class_of, "", set_i = "clean")
panel_levels <- sprintf("%s\n%s", targets_rep, panel_class)[order(match(panel_class, class_levels))]   # quantifiable genes first
panel_of <- function(g) factor(sprintf("%s\n%s", g, panel_class[as.character(g)]), levels = panel_levels)
fig[, panel := panel_of(gene)]
bar_rows <- group_fit[analysis == "censored regression" & animals_included == animal_sets[["clean"]] & gene %chin% targets_rep]
bars <- merge(bar_rows, centres, by = "group")
bars[, `:=`(y = 2^-mean_delta_cq, panel = panel_of(gene))]
stopifnot(nrow(bars) == nrow(bar_rows), !anyNA(bars$y), !anyNA(fig$y))
p_expr <- ggplot(fig, aes(x_pos, y)) +
  bar(bars) + ring(fig[failed_sample == TRUE]) +
  geom_point(aes(colour = ploidy, shape = mark), size = 1.9, stroke = 0.7) +
  facet_wrap(~panel, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = config$colours$ploidy, name = NULL) +
  scale_shape_manual(values = marks, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = seq_along(x_labels), labels = x_labels, limits = c(0.45, length(x_labels) + 0.55)) +
  scale_y_continuous(trans = "log2", breaks = scales::trans_breaks("log2", function(v) 2^v, n = 5),
                     labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  labs(x = NULL, y = sprintf("Expression relative to %s, 2^-(delta Cq)", tolower(ref)), caption = paste0(
    "One mark per animal. Bars: group means of the censored regression without failed samples; none where a gene is not analysed. ",
    "Ring: failed sample (", tolower(ref), " Cq above the outer fence).\nAnimals given the acute heat shock (+ acute) are shown and ",
    "are in no model. ", not_shown, "Under each gene: its detection class.")) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_wgbs(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0), axis.text.x = element_text(size = 7.5),
        plot.caption = element_text(size = 8, hjust = 0), strip.text = element_text(lineheight = 0.95))
save_figure(p_expr, output_path(step, "figures", "figure_qpcr_expression"), width = 12.5, height = 7.6)

# ---- raw Cq of the reference gene and of COX1 ----
ref_fig_rows <- animals[!is.na(reference_cq)]
ref_fig <- merge(ref_fig_rows, per_gene[gene == ref, .(animal_id, n_detected)], by = "animal_id", sort = FALSE)
ref_fig[, mark := factor(names(marks)[fifelse(n_detected == design$wells_per_animal, 1L, 2L)], levels = names(marks)[1:2])]
ref_bar_rows <- animals[animal_id %chin% in_set$clean, .(y = mean(reference_cq)), by = group]
ref_bars <- merge(ref_bar_rows, centres, by = "group")
stopifnot(nrow(ref_fig) == nrow(ref_fig_rows), nrow(ref_bars) == nrow(ref_bar_rows), !anyNA(ref_bars$y))
p_ref <- ggplot(ref_fig, aes(x_pos, reference_cq)) +
  geom_hline(yintercept = ref_fence, linewidth = 0.3, colour = "grey45") +
  annotate("text", x = length(x_labels) + 0.5, y = ref_fence, label = "outer fence", hjust = 1, vjust = -0.5, size = 2.9, colour = "grey30") +
  bar(ref_bars) + ring(ref_fig[failed_sample == TRUE]) +
  geom_point(aes(colour = ploidy, shape = mark), size = 2.1, stroke = 0.7) +
  ggrepel::geom_text_repel(data = ref_fig[failed_sample == TRUE], aes(label = animal_id), size = 2.9, nudge_x = 0.45, seed = 1,
                           min.segment.length = 0) +
  scale_colour_manual(values = config$colours$ploidy, name = NULL) + scale_shape_manual(values = marks[1:2], name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = seq_along(x_labels), labels = x_labels, limits = c(0.45, length(x_labels) + 0.55)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  coord_cartesian(ylim = range(ref_fig$reference_cq) + c(-0.6, 0.9)) +      # a fence far above the data must not stretch the axis
  labs(x = NULL, y = sprintf("%s Cq (a later cycle means less template)", ref)) + theme_wgbs(base_size = 11)
mito_fig <- merge(mito_d, animals[, .(animal_id, x_pos, x_axis)], by = "animal_id", sort = FALSE)
stopifnot(nrow(mito_fig) == nrow(mito_d))
x_used <- sort(unique(mito_fig$x_axis))                                    # only the columns that hold data
mito_fig[, x_two := match(x_axis, x_used) + x_pos - x_axis]
mito_bars <- mito_fig[, .(y = mean(cq), x_mid = mean(x_two)), by = group]   # the animals in the models, as in qpcr_cox1.csv
stopifnot(!anyNA(mito_bars$y))
p_mito <- ggplot(mito_fig, aes(x_two, cq)) +
  bar(mito_bars) + ring(mito_fig[failed_sample == TRUE]) +
  geom_point(aes(colour = ploidy), shape = 16, size = 2.1, show.legend = FALSE) +   # panel A carries the colour legend
  scale_colour_manual(values = config$colours$ploidy, name = NULL) +
  scale_x_continuous(breaks = seq_along(x_used), labels = x_labels[x_used], limits = c(0.45, length(x_used) + 0.55)) +
  labs(x = NULL, y = sprintf("%s Cq, not normalized", mito)) + theme_wgbs(base_size = 11)
p_raw <- p_ref + p_mito + plot_layout(guides = "collect", widths = c(1.7, 1)) +
  plot_annotation(tag_levels = "A", caption = paste0(
    "One mark per animal. Bars: group means (A: control and heat animals without failed samples; B: all animals shown).\n",
    "Ring: failed sample (", tolower(ref), " Cq above the outer fence", if (length(no_reference) > 0) sprintf(" or no %s Cq", tolower(ref)) else "",
    ").", if (nzchar(not_shown)) paste0(" ", trimws(sub("Not shown", "Not in A", not_shown))) else "", "\n",
    mito, if (mito %in% rep_genes) " values are means of the detected wells; no nuclear reference exists for them." else
      " values come from the consolidated file; no replicate-level data and no nuclear reference exist for them.")) &
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0), plot.caption = element_text(size = 8, hjust = 0))
save_figure(p_raw, output_path(step, "figures", "figure_qpcr_cox1_actin"), width = 10, height = 5.4)

finish_step(step, t0)

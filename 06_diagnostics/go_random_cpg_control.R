# go_random_cpg_control.R
#
# Purpose   A negative control for the GO enrichment. CpGs are drawn at random
#           from the eligible CpGs (mean methylation 10% to 90%), as many as
#           there are ploidy DML or heat DML, and treated as if they were DML.
#           They have nothing to do with ploidy or heat, so a sound enrichment
#           analysis should call next to no term enriched for them.
#             A  gene level, tested universe, Fisher, Benjamini-Hochberg
#             B  gene level, eligible universe, Fisher, Benjamini-Hochberg
#             C  gene level, eligible universe, Wallenius, Benjamini-Hochberg
#                (the primary analysis of step 10)
#           If A calls many terms enriched for random CpGs, enrichment against
#           the tested universe is a property of where methylation is, and says
#           nothing about a treatment. C should stay near zero.
#
# Run       Rscript 06_diagnostics/go_random_cpg_control.R     (about 10 minutes)
#
# Reads     04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz, dml_{ploidy,heat}_legacy.csv
#           04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz, gene_go_map.csv.gz
#           03_input_files/genome_features/ gene track
# Writes    06_diagnostics/go_random_cpg_control.csv          one row per draw and ontology
#           06_diagnostics/go_random_cpg_control_summary.csv  medians and ranges, next to the real DML

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges", "topGO", "goseq"))

n_draws <- 10        # random draws per list size
set.seed(20261119)
mode <- "legacy"
onts <- c("BP", "CC", "MF")
dml_in <- function(f) diagnostic_input("05_differential_methylation", "tables", f)
ann_in <- function(f) diagnostic_input("08_gene_annotation", "tables", f)

# ---- Eligible CpGs and the protein-coding genes they lie in --------------------------
united <- fread(dml_in(sprintf("united_cpg_%s.csv.gz", mode)))
eligible <- which(dml_eligible(mean_percent_methylation(united), config))
cpg <- united[eligible, .(chr, start)]; rm(united); gc(verbose = FALSE)
gene_tab <- fread(ann_in(sprintf("gene_table_%s.csv.gz", mode)), colClasses = c(gene_id = "character"))[biotype == "protein_coding"]
pairs <- overlap_pairs(cpg_ranges(cpg$chr, cpg$start), GRanges(gene_tab$chr, IRanges(gene_tab$start, gene_tab$end)))
pairs[, gene_id := gene_tab$gene_id[feature]]
sizes <- c(ploidy = nrow(fread(dml_in(sprintf("dml_ploidy_%s.csv", mode)))), heat = nrow(fread(dml_in(sprintf("dml_heat_%s.csv", mode)))))
cat(sprintf("%d eligible CpGs; list sizes: %s\n", nrow(cpg), paste(names(sizes), sizes, collapse = ", ")))

# ---- Universes and topGO objects, as in step 10 ----------------------------------------
gene_go <- fread(ann_in("gene_go_map.csv.gz"), colClasses = c(gene_id = "character"))
gene2go_all <- split(gene_go$go_id, gene_go$gene_id)
universe <- list(tested = cbind(gene_tab[n_cpg_tested > 0], bias = gene_tab[n_cpg_tested > 0, n_cpg_tested]),
                 eligible = cbind(gene_tab[n_cpg_eligible > 0], bias = gene_tab[n_cpg_eligible > 0, n_cpg_eligible]))
objects <- lapply(universe, function(u) lapply(setNames(onts, onts), function(o) go_gene_level_object(o, u[n_go_terms > 0, gene_id], gene2go_all)))

count_enriched <- function(tab, column) tab[tested == TRUE, sum(get(column) < config$go_fdr)]
one_list <- function(dml_genes, label, draw) {
  rbindlist(lapply(onts, function(o) {
    a <- go_gene_level_test(objects$tested[[o]], universe$tested, dml_genes, config)$table
    b <- go_gene_level_test(objects$eligible[[o]], universe$eligible, dml_genes, config)$table
    data.table(list_size = label, draw = draw, ontology = o, genes_hit = length(dml_genes),
               A_tested_fisher = count_enriched(a, "q_fisher"), A_strongest_term = a[tested == TRUE][order(p_fisher)][1, term],
               B_eligible_fisher = count_enriched(b, "q_fisher"), C_eligible_wallenius = count_enriched(b, "q_wallenius"))
  }))
}

rows <- list()
for (label in names(sizes)) for (d in seq_len(n_draws)) {
  drawn <- sample.int(nrow(cpg), sizes[[label]])
  rows[[length(rows) + 1]] <- one_list(unique(pairs[cpg %in% drawn, gene_id]), label, d)
  cat(sprintf("%s-sized list, draw %2d: terms called enriched A %s | B %s | C %s\n", label, d,
              paste(tail(rows, 1)[[1]]$A_tested_fisher, collapse = "/"), paste(tail(rows, 1)[[1]]$B_eligible_fisher, collapse = "/"),
              paste(tail(rows, 1)[[1]]$C_eligible_wallenius, collapse = "/")))
}
random <- rbindlist(rows)
fwrite(random, here("06_diagnostics", "go_random_cpg_control.csv"))

# ---- The real DML lists through the same code, for comparison ---------------------------
real <- rbindlist(lapply(names(sizes), function(label) {
  one_list(gene_tab[get(paste0("n_dml_", label)) > 0 & (n_cpg_tested > 0), gene_id], label, 0L)
}))
summary_tab <- merge(
  random[, .(random_draws = .N, random_genes_hit_median = median(genes_hit),
             A_random_median = median(A_tested_fisher), A_random_min = min(A_tested_fisher), A_random_max = max(A_tested_fisher),
             B_random_median = median(B_eligible_fisher), B_random_max = max(B_eligible_fisher),
             C_random_median = median(C_eligible_wallenius), C_random_max = max(C_eligible_wallenius)), by = .(list_size, ontology)],
  real[, .(list_size, ontology, real_genes = genes_hit, A_real = A_tested_fisher, B_real = B_eligible_fisher, C_real = C_eligible_wallenius)],
  by = c("list_size", "ontology"))
fwrite(summary_tab, here("06_diagnostics", "go_random_cpg_control_summary.csv"))
print(summary_tab)

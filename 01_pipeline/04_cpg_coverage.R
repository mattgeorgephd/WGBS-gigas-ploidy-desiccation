# 04_cpg_coverage.R
#
# Purpose   Describe the methylation data before any group comparison:
#             - how many CpGs are covered at 1x, 3x, 5x, 10x, 20x in every
#               library of a sample set (all, diploid, triploid, control, heat),
#               and how many are covered in at least one library
#             - the genome-wide number of CpG sites, counted in the reference
#             - the distribution of percent methylation per CpG in each library
#             - pairwise correlation between libraries and a principal component
#               analysis of all libraries
#           and draw Figure S1.
#
# Inputs    03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files and the reference genome listed in
#           03_input_files/gannet_manifest.csv (cache)
# Outputs   04_output/04_cpg_coverage/tables/
#             library_coverage.csv               per-library CpG counts, coverage, upper cutoff
#             genome_cpg_total.csv               CpG sites in the reference, by sequence class
#             cpg_counts_by_coverage.csv         counts per mode, sample set, and threshold
#             methylation_distribution_by_library.csv, methylation_distribution_by_group.csv
#             library_correlation_<level>_<mode>.csv, library_correlation_summary.csv
#             pca_scores_<level>_<mode>.csv, pca_variance.csv
#           04_output/04_cpg_coverage/figures/figure_S1_<mode>.(png|pdf)
#
# Manuscript items: Results 3.1 (CpG counts, correlations, PCA); Figure S1.
# Notes     - "Covered in all libraries of a set" is the definition behind the
#             5x, 10x, 20x counts in the manuscript; the manuscript's 1x count
#             uses "covered in at least one library" instead (issue I-16). Both
#             are written here, labelled.
#           - Modes differ in the upper coverage filter (issue I-05) and in
#             whether the first line of each coverage file is read (issue I-38).
#             Each file is read once, in full; legacy mode then ignores the CpG
#             on the first line of a library's file for that library.
#           - Percent methylation for the correlation and the PCA is computed
#             from counts scaled between libraries the way methylKit scales them
#             (median normalization, with rounding), as in the legacy analysis.
#           - Needs about 8 GB of memory.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "ggplot2", "patchwork", "Biostrings"))

step <- "04_cpg_coverage"
t0 <- start_step(step, "CpG coverage, methylation distributions, library similarity")
lib <- read_libraries(); sets <- sample_sets(lib); n_lib <- nrow(lib)
levels_x <- config$coverage_levels
tab <- function(name) output_path(step, "tables", name)

# =============================================================================
# 1. Read each library once: positions, counts, per-library summaries
# =============================================================================
chr_levels <- character(0); per_lib <- vector("list", n_lib); universe <- numeric(0)
first_key <- numeric(n_lib)                               # the CpG on the first line of each file
lib_rows <- list(); dist_rows <- list()
bins <- seq(0, 100, by = 10)
bin_labels <- c("0-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")

for (j in seq_len(n_lib)) {
  d <- read_cov(fetch_gannet_file(lib$cov_file[j], config))
  d[, coverage := nC + nT]
  stopifnot(all(d$coverage >= 1))                         # Bismark writes covered CpGs only
  chr_levels <- union(chr_levels, unique(d$chr))          # indices of earlier names never change
  key <- cpg_key(d$chr, d$start, chr_levels)
  first_key[j] <- key[1]
  cutoff <- upper_coverage_cutoff(d$coverage, config$import_min_coverage, config$upper_coverage_percentile)
  lib_rows[[j]] <- data.table(
    library_id = lib$library_id[j], n_cpg_1x = nrow(d),
    first_line_position = sprintf("%s:%d", d$chr[1], d$start[1]), first_line_coverage = d$coverage[1],
    n_cpg_10x = sum(d$coverage >= config$min_coverage),
    mean_coverage = mean(d$coverage), median_coverage = median(d$coverage),
    upper_cutoff = cutoff, n_cpg_at_or_above_cutoff = sum(d$coverage >= cutoff))
  pm <- 100 * d$nC / d$coverage
  for (lv in c(1, config$min_coverage)) for (mode in config$modes) {
    keep <- d$coverage >= lv & (mode == "legacy" | d$coverage < cutoff)
    if (first_line_lost(mode)) keep[1] <- FALSE
    h <- tabulate(cut(pm[keep], breaks = bins, include.lowest = TRUE, labels = FALSE), nbins = 10)
    dist_rows[[length(dist_rows) + 1]] <- data.table(
      library_id = lib$library_id[j], mode = mode, min_coverage = lv, bin = bin_labels,
      n_cpg = h, percent_of_cpg = 100 * h / sum(h))
  }
  per_lib[[j]] <- list(key = key, coverage = d$coverage, nC = d$nC)
  universe <- unique(c(universe, key))
  log_msg("%-10s %9d CpGs with reads; median coverage %g; 99.9th percentile cutoff %g",
          lib$library_id[j], nrow(d), median(d$coverage), cutoff)
  rm(d, key, pm); gc(verbose = FALSE)
}
library_coverage <- merge(as.data.table(lib[, c("library_id", "cohort", "ploidy", "treatment")]),
                          rbindlist(lib_rows), by = "library_id", sort = FALSE)
write_table(library_coverage, tab("library_coverage.csv"))
cutoffs <- library_coverage$upper_cutoff

# Coverage and methylated-count matrices: one row per CpG seen in any library,
# one column per library (legacy order); zero where a library has no reads.
universe <- sort(universe); n_cpg <- length(universe)
cov_mat <- matrix(0L, n_cpg, n_lib); nC_mat <- matrix(0L, n_cpg, n_lib)
for (j in seq_len(n_lib)) {
  i <- match(per_lib[[j]]$key, universe)
  cov_mat[i, j] <- per_lib[[j]]$coverage; nC_mat[i, j] <- per_lib[[j]]$nC
  per_lib[j] <- list(NULL)
}
rm(per_lib); gc(verbose = FALSE)
first_row <- match(first_key, universe)                   # row of each library's first-line CpG
log_msg("%d CpGs have at least one read in at least one library", n_cpg)

# =============================================================================
# 2. Genome-wide number of CpG sites, from the reference sequence
# =============================================================================
genome <- readDNAStringSet(fetch_gannet_file("cgigas_uk_roslin_v1_genomic-mito.fa", config))
seq_name <- sub(" .*", "", names(genome))
cg <- vcountPattern("CG", genome)
seq_class <- fifelse(seq_name %in% config$linkage_groups, "linkage groups",
             fifelse(seq_name == config$mito_sequence, "mitochondrial genome", "unplaced scaffolds"))
genome_cpg <- data.table(sequence_class = seq_class, bp = width(genome), cpg = cg)[
  , .(n_sequences = .N, bp = sum(as.numeric(bp)), cpg_sites = sum(cpg)), by = sequence_class]
genome_cpg <- rbind(genome_cpg, data.table(sequence_class = "total", n_sequences = sum(genome_cpg$n_sequences),
                                           bp = sum(genome_cpg$bp), cpg_sites = sum(genome_cpg$cpg_sites)))
write_table(genome_cpg, tab("genome_cpg_total.csv"))
total_cpg <- genome_cpg[sequence_class == "total", cpg_sites]
log_msg("reference: %d sequences, %.1f Mb, %d CpG sites", length(genome), sum(as.numeric(width(genome))) / 1e6, total_cpg)
rm(genome); gc(verbose = FALSE)

# =============================================================================
# 3. CpG counts by coverage threshold, sample set, and mode
# =============================================================================
# A CpG counts for a library when its coverage is at least x and, in corrected
# mode, below that library's upper cutoff. In legacy mode the CpG on the first
# line of the library's file does not count for that library (issue I-38).
covered <- function(j, x, mode) {
  cv <- cov_mat[, j]
  ok <- if (mode == "legacy") cv >= x else cv >= x & cv < cutoffs[j]
  if (first_line_lost(mode)) ok[first_row[j]] <- FALSE
  ok
}
count_rows <- list()
for (mode in config$modes) for (s in names(sets)) for (x in levels_x) {
  in_all <- rep(TRUE, n_cpg); in_any <- rep(FALSE, n_cpg)
  for (j in sets[[s]]) { ok <- covered(j, x, mode); in_all <- in_all & ok; in_any <- in_any | ok }
  count_rows[[length(count_rows) + 1]] <- data.table(
    mode = mode, sample_set = s, n_libraries = length(sets[[s]]), min_coverage = x,
    cpg_in_all_libraries = sum(in_all), cpg_in_any_library = sum(in_any))
}
cpg_counts <- rbindlist(count_rows)
cpg_counts[, `:=`(pct_of_genome_cpg_all = 100 * cpg_in_all_libraries / total_cpg,
                  pct_of_genome_cpg_any = 100 * cpg_in_any_library / total_cpg)]
write_table(cpg_counts, tab("cpg_counts_by_coverage.csv"))
print(dcast(cpg_counts[sample_set == "all"], mode ~ min_coverage, value.var = "cpg_in_all_libraries"))

# =============================================================================
# 4. Methylation distribution per library and per group
# =============================================================================
dist_lib <- merge(rbindlist(dist_rows), as.data.table(lib[, c("library_id", "ploidy", "treatment")]), by = "library_id")
dist_lib[, bin := factor(bin, levels = bin_labels)]
write_table(dist_lib, tab("methylation_distribution_by_library.csv"))
se <- function(v) sd(v) / sqrt(length(v))
dist_group <- rbind(
  dist_lib[, .(grouping = "ploidy", n = .N, mean_percent = mean(percent_of_cpg), se_percent = se(percent_of_cpg)),
           by = .(mode, min_coverage, group = as.character(ploidy), bin)],
  dist_lib[, .(grouping = "treatment", n = .N, mean_percent = mean(percent_of_cpg), se_percent = se(percent_of_cpg)),
           by = .(mode, min_coverage, group = as.character(treatment), bin)],
  dist_lib[, .(grouping = "all", n = .N, mean_percent = mean(percent_of_cpg), se_percent = se(percent_of_cpg)),
           by = .(mode, min_coverage, group = rep("all libraries", nrow(dist_lib)), bin)])
write_table(dist_group, tab("methylation_distribution_by_group.csv"))

# =============================================================================
# 5. Library correlation and PCA on CpGs covered in all 22 libraries
# =============================================================================
# Percent methylation after methylKit-style median normalization of the counts
# of the CpGs each library retains at this coverage level.
percent_methylation <- function(rows, x, mode) {
  pm <- matrix(NA_real_, length(rows), n_lib)
  med <- vapply(seq_len(n_lib), function(j) median(cov_mat[covered(j, x, mode), j]), numeric(1))
  sc <- max(med) / med
  for (j in seq_len(n_lib)) {
    cv <- cov_mat[rows, j]; nc <- nC_mat[rows, j]
    new_cov <- round(sc[j] * cv)
    new_c <- round(new_cov * nc / cv); new_t <- round(new_cov * (cv - nc) / cv)
    pm[, j] <- 100 * new_c / (new_c + new_t)
  }
  colnames(pm) <- lib$library_id
  pm
}
cor_summary <- list(); pca_var <- list(); pca_scores <- list()
for (mode in config$modes) for (x in c(1, config$min_coverage)) {
  in_all <- rep(TRUE, n_cpg)
  for (j in seq_len(n_lib)) in_all <- in_all & covered(j, x, mode)
  rows <- which(in_all)
  pm <- percent_methylation(rows, x, mode)
  pm <- pm[complete.cases(pm), , drop = FALSE]
  tag <- sprintf("%dx_%s", x, mode)

  r <- cor(pm)                                            # Pearson, as methylKit::getCorrelation
  write_table(data.table(library_id = rownames(r), r), tab(sprintf("library_correlation_%s.csv", tag)))
  off <- r[upper.tri(r)]
  same_cohort <- outer(lib$cohort, lib$cohort, "==")[upper.tri(r)]
  cor_summary[[tag]] <- data.table(mode = mode, min_coverage = x, n_cpg = nrow(pm),
                                   r_min = min(off), r_max = max(off), r_median = median(off),
                                   r_median_within_cohort = median(off[same_cohort]),
                                   r_median_between_cohorts = median(off[!same_cohort]))

  # PCA as methylKit::PCASamples: keep CpGs whose SD across libraries is above
  # the median SD, then centre and scale.
  sds <- sqrt(rowSums((pm - rowMeans(pm))^2) / (n_lib - 1))
  pc <- prcomp(t(pm[sds > quantile(sds, 0.5), , drop = FALSE]), center = TRUE, scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  pca_var[[tag]] <- data.table(mode = mode, min_coverage = x, pc = seq_along(ve), percent_variance = ve)
  pca_scores[[tag]] <- data.table(mode = mode, min_coverage = x, library_id = rownames(pc$x), pc$x[, 1:4],
                                  lib[, c("cohort", "ploidy", "treatment", "acute_heat_shock")])
  write_table(pca_scores[[tag]], tab(sprintf("pca_scores_%s.csv", tag)))
  log_msg("[%s] %d CpGs; Pearson r %.3f to %.3f; PC1 %.1f%%, PC2 %.1f%%", tag, nrow(pm), min(off), max(off), ve[1], ve[2])
  rm(pm, pc, sds); gc(verbose = FALSE)
}
write_table(rbindlist(cor_summary), tab("library_correlation_summary.csv"))
write_table(rbindlist(pca_var), tab("pca_variance.csv"))
pca_all <- rbindlist(pca_scores); pca_var_all <- rbindlist(pca_var)

# =============================================================================
# 6. Figure S1 (one per mode)
# =============================================================================
set_cols <- c(all = "chartreuse4", diploid = config$colours$ploidy[["diploid"]],
              triploid = config$colours$ploidy[["triploid"]], control = config$colours$ploidy[["diploid"]],
              heat = config$colours$ploidy[["triploid"]])
set_shapes <- c(all = 16, diploid = 16, triploid = 16, control = 15, heat = 15)
for (mode_i in config$modes) {
  cc <- cpg_counts[mode == mode_i]
  cc[, sample_set := factor(sample_set, levels = names(set_cols))]
  pA <- ggplot(cc, aes(min_coverage, cpg_in_all_libraries / 1e6, colour = sample_set, shape = sample_set)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.6) +
    scale_colour_manual(values = set_cols, name = NULL) + scale_shape_manual(values = set_shapes, name = NULL) +
    scale_x_continuous(breaks = levels_x) +
    labs(x = "Minimum coverage in every library (x)", y = "CpG loci (millions)") + theme_wgbs()
  dist_panel <- function(grouping_i, cols) {
    dd <- dist_group[mode == mode_i & min_coverage == 1 & grouping %in% c(grouping_i, "all")]
    dd[, group := factor(group, levels = names(cols))]
    ggplot(dd, aes(bin, mean_percent, fill = group)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.75) +
      geom_errorbar(aes(ymin = mean_percent - se_percent, ymax = mean_percent + se_percent),
                    position = position_dodge(width = 0.8), width = 0.25) +
      scale_fill_manual(values = cols, name = NULL) +
      labs(x = "Methylation per CpG (%)", y = "CpG loci (% of library, mean and SE)") +
      theme_wgbs() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  pB <- dist_panel("ploidy", c(diploid = set_cols[["diploid"]], triploid = set_cols[["triploid"]], `all libraries` = "chartreuse4"))
  pC <- dist_panel("treatment", c(control = "grey35", heat = "firebrick2", `all libraries` = "chartreuse4"))
  ps <- pca_all[mode == mode_i & min_coverage == 1]; pv <- pca_var_all[mode == mode_i & min_coverage == 1]
  pD <- ggplot(ps, aes(PC1, PC2, colour = ploidy, shape = treatment)) + geom_point(size = 3) +
    scale_colour_manual(values = config$colours$ploidy, name = NULL) +
    scale_shape_manual(values = c(control = 16, heat = 17), name = NULL) +
    labs(x = sprintf("PC 1 (%.1f%%)", pv$percent_variance[1]), y = sprintf("PC 2 (%.1f%%)", pv$percent_variance[2])) +
    theme_wgbs()
  fig <- (pA + pB) / (pC + pD) + plot_annotation(tag_levels = "A")
  save_figure(fig, output_path(step, "figures", sprintf("figure_S1_%s", mode_i)), width = 11, height = 9)
}

finish_step(step, t0)

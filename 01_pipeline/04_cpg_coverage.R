# 04_cpg_coverage.R
#
# Purpose   Describe the methylation data before any group comparison:
#             - how many CpGs are covered at 1x, 3x, 5x, 10x, 20x in every
#               library of a sample set (all, diploid, triploid, control, heat),
#               and how many are covered in at least one library
#             - which genome feature tracks the CpGs of each of those sets overlap
#             - the genome-wide number of CpG sites, counted in the reference
#             - the distribution of percent methylation per CpG in each library
#             - pairwise correlation between libraries and a principal component
#               analysis of all libraries
#           and draw Figure S1.
#
# Run       Rscript 01_pipeline/04_cpg_coverage.R     (after step 01; about 12 minutes, about 9.5 GB of memory)
#
# Inputs    03_input_files/metadata/wgbs_libraries.csv
#           22 coverage files and the reference genome listed in
#           03_input_files/gannet_manifest.csv (cache)
#           03_input_files/genome_features/*.gz   the 13 tracks; what each holds is in the README there
# Outputs   04_output/04_cpg_coverage/tables/
#             library_coverage.csv               per-library CpG counts, coverage, upper cutoff
#             genome_cpg_total.csv               CpG sites in the reference, by sequence class
#             cpg_counts_by_coverage.csv         counts per mode, sample set, and threshold
#             cpg_counts_by_feature_track.csv    the CpGs of each of those sets that overlap each
#                                                feature track, and the two percentages the draft
#                                                made of such counts
#             methylation_distribution_by_library.csv, methylation_distribution_by_group.csv
#             library_correlation_<level>_<mode>.csv, library_correlation_summary.csv
#             pca_scores_<level>_<mode>.csv, pca_variance.csv
#           04_output/04_cpg_coverage/figures/figure_S1_<mode>.(png|pdf)
#           04_output/04_cpg_coverage/session_info.txt
#
# Manuscript items: Results 3.1 (CpG counts, genomic distribution, correlations, PCA); Figure S1.
# Notes     - cpg_counts_by_coverage.csv holds three counts per mode, sample set
#             and threshold x:
#               cpg_in_all_libraries            every library of the set retains
#                                               the CpG with at least x reads
#               cpg_in_any_library              at least one library of the set
#                                               retains it with at least x reads
#               cpg_in_any_library_unfiltered   at least one library of the set
#                                               has at least x reads at the CpG,
#                                               every line of every file counted
#                                               and no upper filter applied; the
#                                               same number in both modes
#             What a library retains depends on the mode (next note). The first
#             count is the definition behind the 5x, 10x and 20x counts of the
#             manuscript. The manuscript's 1x count (12,794,074) is on an unknown
#             definition and could not be traced (issue I-16): it is none of the
#             three counts at 1x (5,852,939; 12,918,250; 12,918,252 in legacy mode).
#           - Modes differ in the upper coverage filter (issue I-05) and in
#             whether the first line of each coverage file is read (issue I-38).
#             Each file is read once, in full. In legacy mode a library does not
#             retain the CpG on the first line of its file, and there is no upper
#             filter. In corrected mode a library does not retain a CpG whose
#             coverage is at or above the library's upper cutoff: the 99.9th
#             percentile of coverage taken over the library's CpGs with at least
#             config$import_min_coverage (2) reads, the CpGs that methylKit's
#             import keeps in step 05. CpGs with a single read are left out of
#             the percentile and are still counted at 1x.
#           - cpg_counts_by_feature_track.csv gives, for the CpG sets behind the
#             first two counts (cpg_set "in all libraries" and "in any library"),
#             the CpGs that overlap each of the 13 feature tracks: n_overlap,
#             where a CpG is counted once per track, as by `bedtools intersect
#             -u`, and in every track it overlaps. Two columns turn the counts
#             into percentages the way the draft did. Both divide by the sum of
#             the counts of eight tracks (exon UTR, CDS, intron, upstream flank,
#             downstream flank, lncRNA, TE, intergenic) and not by the number of
#             CpGs, so they describe overlaps, not loci (issue I-18):
#               legacy_share_pct              the eight counts as they are: the
#                                             rule of Results 3.3 and of step 07
#               legacy_share_exon_twice_pct   the count of the exon track is first
#                                             added to that of the UTR track, as
#                                             in the legacy workbook ("exon UTR" =
#                                             exon + UTR), so the CpGs of coding
#                                             exons are counted under exon UTR
#                                             and again under CDS: the rule of
#                                             Results 3.1 (issue I-15)
#             Results 3.1 gives eight percentages "across all samples at 1x"
#             (exon UTR 27.3, introns 22.9, CDS 16, intergenic 14.3, transposable
#             elements 13.1, upstream flanks 2.7, downstream flanks 2, lncRNA
#             1.8); its "transposable elements" are the repeat track TE. What they match: one of the 100 CpG sets of the table, under
#             the second rule only. It is the set of legacy mode with at least 10
#             reads in every one of the ten heat libraries (1,921,856 CpGs;
#             27.30, 22.89, 16.02, 14.27, 13.05, 2.66, 2.01, 1.79), which
#             reproduces all eight values and their order at the one decimal the
#             draft prints. What they do not match: any 1x set, any set of all 22
#             libraries, and any set under the first rule. For the CpGs with a
#             read in all 22 libraries the second rule gives 20.55, 24.00, 11.69,
#             17.32, 19.03, 2.94, 2.65, 1.82, with introns first. The counts of
#             the legacy workbook are not reproduced either (832,055 for exon UTR
#             there, 830,319 here): they are the overlaps of the 1,926,243 CpGs
#             of an earlier analysis of the ten heat libraries
#             (08_archive/old_WGBS/WGBS/DML/ALL-CpG-Cov10.csv), of which 1,912,802
#             are among the 1,921,856 of this step. That comparison was made
#             once, outside the pipeline, which does not read the archive.
#             The track named TE holds repeats of every class, read as the legacy
#             analysis read the file. By the author's decision of 2026-11-20
#             (issues.md, I-50) it is labelled "repetitive elements
#             (RepeatMasker, all classes)"; the track code TE stands for the
#             repetitive elements track. Column track_label says what each track
#             holds (02_R_functions/genomic_overlaps.R).
#           - Percent methylation for the correlation and the PCA is computed
#             from counts scaled between libraries as the legacy analysis scaled
#             them, by methylKit's median normalization, which is written out
#             here without a call to methylKit. One detail is not methylKit's: a
#             scaled count is round(new_cov * n / cov) here and
#             round(new_cov * (n / cov)) in methylKit and in step 05, and the two
#             can differ where the product lies exactly halfway between two whole
#             numbers. In legacy mode that changes 224 of the 128.8 million
#             CpG-by-library values at 1x and 66 of the 10.6 million at 10x
#             (corrected mode: 167 and 41). No correlation moves by more than
#             7e-7 and no percent variance of a principal component by more than
#             3e-5 points. The code is left as it is, because
#             06_diagnostics/check_methylkit_equivalence.R records this order of
#             operations and its effect (issues.md, I-16).

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "ggplot2", "patchwork", "Biostrings", "GenomicRanges"))

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
# 3. CpG counts by coverage threshold, sample set, and mode; the same CpG sets
#    by genome feature
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

# Which feature tracks each CpG overlaps: one logical vector per track, the
# `bedtools intersect -u` of the legacy notebooks.
position <- key_to_position(universe, chr_levels)
cpg_gr <- cpg_ranges(position$chr, position$start)
tracks <- names(feature_track_files())
in_track <- lapply(setNames(tracks, tracks), function(t) overlaps_track(cpg_gr, read_feature_track(t)))
rm(position, cpg_gr); invisible(gc(verbose = FALSE))

# Overlap counts of one CpG set, and the two ways the draft turned such counts
# into percentages. Both divide the count of a track by the sum of the counts of
# eight tracks, so a CpG that lies in two tracks is counted twice (issue I-18).
# Results 3.3 takes the eight counts as they are (legacy_share_pct, as in step
# 07). Results 3.1 first adds the count of the exon track to that of the UTR
# track and calls the sum "exon UTR", which counts the CpGs of every exon under
# that name and those of coding exons again under CDS (issue I-15).
legacy_eight <- c("exon_UTR", "CDS", "intron", "upstream_flank", "downstream_flank", "lncRNA", "TE", "intergenic")
pct <- function(n, of) if (of > 0) 100 * n / of else rep(NA_real_, length(n))    # no percentage of nothing
feature_counts <- function(in_set) {
  n <- vapply(in_track, function(hit) sum(hit & in_set), integer(1))
  exon_twice <- n[legacy_eight]
  exon_twice["exon_UTR"] <- n[["exon_UTR"]] + n[["exon"]]
  data.table(n_loci = sum(in_set), track = tracks, n_overlap = unname(n), pct_of_loci = unname(pct(n, sum(in_set))),
             legacy_share_pct = ifelse(tracks %in% legacy_eight, unname(pct(n, sum(n[legacy_eight]))), NA_real_),
             legacy_share_exon_twice_pct = unname(pct(exon_twice, sum(exon_twice))[tracks]),
             track_label = unname(feature_track_labels()[tracks]))
}

count_rows <- list(); feature_rows <- list()
for (mode in config$modes) for (s in names(sets)) for (x in levels_x) {
  in_all <- rep(TRUE, n_cpg); in_any <- rep(FALSE, n_cpg)
  for (j in sets[[s]]) { ok <- covered(j, x, mode); in_all <- in_all & ok; in_any <- in_any | ok }
  set_id <- data.table(mode = mode, sample_set = s, n_libraries = length(sets[[s]]), min_coverage = x)
  count_rows[[length(count_rows) + 1]] <- cbind(set_id, cpg_in_all_libraries = sum(in_all), cpg_in_any_library = sum(in_any))
  feature_rows[[length(feature_rows) + 1]] <- rbind(cbind(set_id, cpg_set = "in all libraries", feature_counts(in_all)),
                                                    cbind(set_id, cpg_set = "in any library", feature_counts(in_any)))
}
cpg_counts <- rbindlist(count_rows)
cpg_counts[, `:=`(pct_of_genome_cpg_all = 100 * cpg_in_all_libraries / total_cpg,
                  pct_of_genome_cpg_any = 100 * cpg_in_any_library / total_cpg)]

# CpGs with at least x reads in at least one library of the set, counted on every
# line of every file and with no upper filter, so the count is the same in both
# modes. cpg_in_any_library applies the rules of the mode first.
unfiltered <- unique(cpg_counts[, .(sample_set, min_coverage)])
unfiltered[, cpg_in_any_library_unfiltered := mapply(function(s, x) {
  in_any <- rep(FALSE, n_cpg)
  for (j in sets[[s]]) in_any <- in_any | cov_mat[, j] >= x
  sum(in_any)
}, sample_set, min_coverage, USE.NAMES = FALSE)]
cpg_counts[unfiltered, on = c("sample_set", "min_coverage"), cpg_in_any_library_unfiltered := i.cpg_in_any_library_unfiltered]
cpg_counts[, pct_of_genome_cpg_any_unfiltered := 100 * cpg_in_any_library_unfiltered / total_cpg]
stopifnot(cpg_counts[sample_set == "all" & min_coverage == 1, all(cpg_in_any_library_unfiltered == n_cpg)])
write_table(cpg_counts, tab("cpg_counts_by_coverage.csv"))
print(dcast(cpg_counts[sample_set == "all"], mode ~ min_coverage, value.var = "cpg_in_all_libraries"))

cpg_features <- rbindlist(feature_rows)
write_table(cpg_features, tab("cpg_counts_by_feature_track.csv"))
# The shares as Results 3.1 computed them, for the CpGs with config$min_coverage reads in every library of a set.
# Printed only when legacy mode is configured: without it the selection is empty and dcast() stops the step.
if ("legacy" %in% config$modes) {
  print(dcast(cpg_features[mode == "legacy" & cpg_set == "in all libraries" & min_coverage == config$min_coverage & track %in% legacy_eight],
              track ~ sample_set, value.var = "legacy_share_exon_twice_pct"), digits = 3)
}
rm(in_track, cpg_features); invisible(gc(verbose = FALSE))

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

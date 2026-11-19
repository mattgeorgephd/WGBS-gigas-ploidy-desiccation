# 07_genome_features.R
#
# Purpose   Where in the genome the analysed CpGs and the DML lie.
#             - overlap of every analysed CpG with each of the 13 feature tracks
#               (the bedtools intersect -u step of the legacy notebooks, done
#               with GenomicRanges), and the same for each DML set
#             - one mutually exclusive class per CpG, so that percentages add to
#               100 and describe loci rather than overlap events (issue I-18)
#             - tests of association between DML status and genomic feature
#             - DML counts per linkage group, against gene and CpG content
#             - DML grouped into regions (neighbours within 1 kb), because
#               adjacent CpGs are not independent observations (issue I-09)
#           and draw Figure 3.
#
# Inputs    04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           04_output/05_differential_methylation/tables/dml_{ploidy,heat,shared}_<mode>.csv
#           03_input_files/genome_features/*.gz
# Outputs   04_output/07_genome_features/tables/
#             feature_overlap_by_track_<mode>.csv   n and % of loci overlapping each track, per set
#             feature_class_exclusive_<mode>.csv    one class per locus, per set
#             feature_association_tests_<mode>.csv  Fisher tests per track; chi-squared on classes
#             dml_by_linkage_group_<mode>.csv
#             dml_regions_<mode>.csv, dml_region_summary_<mode>.csv
#           04_output/07_genome_features/figures/figure_3_<mode>.(png|pdf)
#
# Manuscript items: Results 3.3; Figure 3.
# Notes     - The background is the set of CpGs that were tested for differential
#             methylation, in the same mode. The legacy background was assembled
#             from two other CpG sets (issue I-14).
#           - Two backgrounds are used for the association tests: all analysed
#             CpGs, and "eligible" CpGs, those whose mean methylation across
#             libraries is between 10% and 90%. A 20-point difference between
#             two groups almost always needs a mean in that range (all but 2 of
#             the legacy DML have one), so the other CpGs, about 86% of those
#             analysed, are almost never DML. Methylation in this genome sits
#             mostly in gene bodies, so enrichment of DML in genes relative to
#             all CpGs is expected whatever the treatment. The eligible
#             background removes most of that and not all of it: inside the
#             range the chance of passing the threshold still depends on the
#             mean and on the variation among animals, and CpG lists made with
#             shuffled ploidy labels show the same odds ratios against eligible
#             CpGs as the real DML (issues.md, I-41). Neither background is a
#             test of a treatment effect on location.
#           - "legacy_share_pct" reproduces the manuscript's percentages: the
#             overlap count of a track divided by the sum of the counts of eight
#             tracks. Those shares double count loci that lie in two tracks.
#           - Class priority for the exclusive assignment: CDS, lncRNA, UTR,
#             intron, upstream flank, downstream flank, intergenic; anything left
#             is "other". lncRNA comes before UTR and intron because the UTR track
#             (exon minus CDS) and the intron track include the exons and introns
#             of lncRNA genes. Transposable element overlap is reported alongside,
#             not as a class, because elements sit inside all the other classes.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges", "ggplot2", "patchwork"))

step <- "07_genome_features"
t0 <- start_step(step, "genomic location of analysed CpGs and DML")
dml_dir <- function(f) output_path("05_differential_methylation", "tables", f)
tab <- function(name) output_path(step, "tables", name)

tracks <- names(feature_track_files())
track_gr <- lapply(setNames(tracks, tracks), function(t) { log_msg("reading track %s", t); read_feature_track(t) })
legacy_eight <- c("exon_UTR", "CDS", "intron", "upstream_flank", "downstream_flank", "lncRNA", "TE", "intergenic")
class_priority <- c("CDS", "lncRNA", "exon_UTR", "intron", "upstream_flank", "downstream_flank", "intergenic")
class_labels <- c(CDS = "CDS", lncRNA = "lncRNA", exon_UTR = "UTR or other non-coding exon", intron = "Intron",
                  upstream_flank = "Upstream flank", downstream_flank = "Downstream flank",
                  intergenic = "Intergenic", other = "Other")

# Genes per linkage group, for the chromosome panel.
gene_gr <- read_feature_track("gene", keep_attributes = TRUE)
genes_per_seq <- as.data.table(table(sequence = as.character(GenomicRanges::seqnames(gene_gr))))
setnames(genes_per_seq, "N", "n_genes")

for (mode in config$modes) {
  log_msg("[%s] reading analysed CpGs and DML", mode)
  united <- fread(dml_dir(sprintf("united_cpg_%s.csv.gz", mode)))
  # CpGs that can reach the difference threshold (see dml_eligible()).
  eligible <- dml_eligible(mean_percent_methylation(united), config)
  cpg <- united[, .(chr, start)]
  cpg[, united_row := .I]
  rm(united); gc(verbose = FALSE)
  dml <- list(ploidy = fread(dml_dir(sprintf("dml_ploidy_%s.csv", mode))),
              heat   = fread(dml_dir(sprintf("dml_heat_%s.csv", mode))))
  shared <- fread(dml_dir(sprintf("dml_shared_%s.csv", mode)))

  # Set membership of every analysed CpG.
  sets <- list(
    `all analysed CpGs`  = rep(TRUE, nrow(cpg)),
    `eligible CpGs`      = eligible,
    `ploidy DML`         = cpg$united_row %in% dml$ploidy$united_row,
    `heat DML`           = cpg$united_row %in% dml$heat$united_row,
    `shared`             = cpg$united_row %in% shared$united_row,
    `shared, both hyper` = cpg$united_row %in% shared[class == "both_hyper", united_row],
    `shared, both hypo`  = cpg$united_row %in% shared[class == "both_hypo", united_row],
    `shared, opposite`   = cpg$united_row %in% shared[class == "opposite", united_row])

  # ---- Overlap with each track ------------------------------------------------
  gr <- cpg_ranges(cpg$chr, cpg$start)
  ov <- vapply(tracks, function(t) overlaps_track(gr, track_gr[[t]]), logical(nrow(cpg)))
  by_track <- rbindlist(lapply(names(sets), function(s) {
    n <- colSums(ov[sets[[s]], , drop = FALSE])
    data.table(mode = mode, set = s, n_loci = sum(sets[[s]]), track = tracks, n_overlap = n,
               pct_of_loci = 100 * n / sum(sets[[s]]),
               legacy_share_pct = ifelse(tracks %in% legacy_eight, 100 * n / sum(n[legacy_eight]), NA_real_))
  }))
  write_table(by_track, tab(sprintf("feature_overlap_by_track_%s.csv", mode)))

  # ---- One class per CpG ---------------------------------------------------------
  cls <- rep("other", nrow(cpg))
  for (t in rev(class_priority)) cls[ov[, t]] <- t          # highest priority assigned last
  cls <- factor(cls, levels = c(class_priority, "other"))
  exclusive <- rbindlist(lapply(names(sets), function(s) {
    n <- as.vector(table(cls[sets[[s]]]))
    data.table(mode = mode, set = s, n_loci = sum(sets[[s]]), class = levels(cls), n = n,
               pct = 100 * n / sum(sets[[s]]),
               pct_in_TE = 100 * as.vector(tapply(ov[sets[[s]], "TE"], cls[sets[[s]]], mean)))
  }))
  write_table(exclusive, tab(sprintf("feature_class_exclusive_%s.csv", mode)))

  # ---- Association tests -----------------------------------------------------------
  # Per track: is a DML more or less likely than other analysed CpGs to lie in it?
  # Loci are treated as independent, which they are not where DML cluster, so
  # these p-values are optimistic; see the region summary below.
  # Two backgrounds: every analysed CpG, and the eligible CpGs. Most CpGs in this
  # genome are unmethylated and cannot differ between groups, and methylation is
  # concentrated in gene bodies, so the second background asks the sharper
  # question.
  backgrounds <- list(`all analysed CpGs` = rep(TRUE, nrow(cpg)), `eligible CpGs` = eligible)
  tests <- rbindlist(lapply(names(backgrounds), function(b) {
    rbindlist(lapply(c("ploidy DML", "heat DML", "shared"), function(s) {
      rbindlist(lapply(c(legacy_eight, "gene"), function(t) {
        use <- backgrounds[[b]] | sets[[s]]
        m <- table(factor(sets[[s]][use], c(TRUE, FALSE)), factor(ov[use, t], c(TRUE, FALSE)))
        ft <- fisher.test(m)
        data.table(mode = mode, background = b, set = s, track = t, dml_in_track = m[1, 1], dml_total = sum(m[1, ]),
                   pct_dml = 100 * m[1, 1] / sum(m[1, ]), pct_other_cpg = 100 * m[2, 1] / sum(m[2, ]),
                   odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
                   p_fisher = ft$p.value)
      }))
    }))
  }))
  tests[, q_BH := p.adjust(p_fisher, "BH")]
  chisq <- rbindlist(lapply(c("ploidy DML", "heat DML", "shared"), function(s) {
    m <- table(sets[[s]], droplevels(cls))
    ct <- suppressWarnings(chisq.test(m))
    data.table(mode = mode, set = s, test = "chi-squared, DML status by exclusive class",
               statistic = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value,
               min_expected = min(ct$expected))
  }))
  write_table(tests, tab(sprintf("feature_association_tests_%s.csv", mode)))
  write_table(chisq, tab(sprintf("feature_class_chisq_%s.csv", mode)))

  # ---- Linkage groups ----------------------------------------------------------------
  seq_label <- function(chr) {
    i <- match(chr, config$linkage_groups)
    fifelse(!is.na(i), sprintf("LG%d", i), fifelse(chr == config$mito_sequence, "mitochondrial", "unplaced"))
  }
  lg_levels <- c(sprintf("LG%d", 1:10), "unplaced", "mitochondrial")
  genes_lg <- genes_per_seq[, .(n_genes = sum(n_genes)), by = .(sequence = seq_label(sequence))]
  per_lg <- cpg[, .(n_cpg_analysed = .N), by = .(sequence = seq_label(chr))]
  count_dir <- function(d, col, name) {
    d[, .(n = .N, n_hyper = sum(get(col) > 0), n_hypo = sum(get(col) < 0)), by = .(sequence = seq_label(chr))][
      , setnames(.SD, c("n", "n_hyper", "n_hypo"), paste0(name, c("", "_hyper", "_hypo")))]
  }
  for (x in list(count_dir(dml$ploidy, "meth.diff", "ploidy_dml"), count_dir(dml$heat, "meth.diff", "heat_dml"),
                 count_dir(shared, "meth.diff_heat", "shared"))) {
    per_lg <- merge(per_lg, x, by = "sequence", all.x = TRUE)
  }
  per_lg <- merge(per_lg, genes_lg, by = "sequence", all.x = TRUE)
  for (col in setdiff(names(per_lg), "sequence")) set(per_lg, which(is.na(per_lg[[col]])), col, 0L)
  per_lg[, sequence := factor(sequence, levels = lg_levels)]; setorder(per_lg, sequence)
  per_lg[, `:=`(mode = mode,
                pct_of_ploidy_dml = 100 * ploidy_dml / sum(ploidy_dml), pct_of_heat_dml = 100 * heat_dml / sum(heat_dml),
                ploidy_dml_per_1000_cpg = 1000 * ploidy_dml / n_cpg_analysed,
                heat_dml_per_1000_cpg = 1000 * heat_dml / n_cpg_analysed)]
  write_table(per_lg, tab(sprintf("dml_by_linkage_group_%s.csv", mode)))

  # ---- Regions: DML of one set within 1 kb of each other ------------------------------
  make_regions <- function(d, col, name) {
    d <- copy(d)[order(chr, start)]
    d[, region := cumsum(c(TRUE, diff(start) > 1000 | chr[-1] != chr[-.N]))]
    d[, .(set = name, chr = chr[1], start = min(start), end = max(end), span_bp = max(end) - min(start),
          n_dml = .N, n_hyper = sum(get(col) > 0), n_hypo = sum(get(col) < 0)), by = region][, region := NULL]
  }
  regions <- rbind(make_regions(dml$ploidy, "meth.diff", "ploidy DML"), make_regions(dml$heat, "meth.diff", "heat DML"),
                   make_regions(shared, "meth.diff_heat", "shared"))
  regions[, mode := mode]; setorder(regions, set, -n_dml)
  write_table(regions, tab(sprintf("dml_regions_%s.csv", mode)))
  region_summary <- regions[, .(n_dml = sum(n_dml), n_regions = .N, n_single_locus_regions = sum(n_dml == 1),
                                largest_region_n_dml = max(n_dml), dml_in_regions_of_5_or_more = sum(n_dml[n_dml >= 5])),
                            by = .(mode, set)]
  write_table(region_summary, tab(sprintf("dml_region_summary_%s.csv", mode)))
  print(region_summary)

  # ---- Figure 3 -------------------------------------------------------------------------
  ex <- exclusive[set != "shared"]                      # the three shared classes are shown instead
  ex[, set_label := sprintf("%s\n(n = %s)", sub(", ", ",\n", sub(" CpGs", "\nCpGs", set)),
                            formatC(n_loci, big.mark = ",", format = "d"))]
  ex[, set_label := factor(set_label, levels = unique(set_label))]
  ex[, class := factor(class, levels = rev(c(class_priority, "other")))]
  pA <- ggplot(ex, aes(set_label, pct, fill = class)) + geom_col(colour = "black", linewidth = 0.3, width = 0.8) +
    scale_fill_brewer(palette = "YlGnBu", direction = -1, labels = class_labels, name = "Genomic class") +
    labs(x = NULL, y = "CpG loci (%)") + theme_wgbs() + theme(axis.text.x = element_text(size = 9))
  lg_long <- melt(per_lg[sequence %in% sprintf("LG%d", 1:10)], id.vars = "sequence",
                  measure.vars = c("ploidy_dml_hypo", "ploidy_dml_hyper", "heat_dml_hypo", "heat_dml_hyper"))
  lg_long[, variable := factor(variable, levels = c("heat_dml_hyper", "heat_dml_hypo", "ploidy_dml_hyper", "ploidy_dml_hypo"),
                               labels = c("Heat, hyper", "Heat, hypo", "Ploidy, hyper", "Ploidy, hypo"))]
  gene_scale <- max(per_lg[sequence %in% sprintf("LG%d", 1:10), n_genes]) /
    max(lg_long[, sum(value), by = sequence]$V1)
  pB <- ggplot(lg_long, aes(sequence, value)) +
    geom_col(aes(fill = variable), colour = "black", linewidth = 0.3) +
    geom_line(data = per_lg[sequence %in% sprintf("LG%d", 1:10)], aes(sequence, n_genes / gene_scale, group = 1), linewidth = 0.7) +
    geom_point(data = per_lg[sequence %in% sprintf("LG%d", 1:10)], aes(sequence, n_genes / gene_scale), size = 2.4) +
    scale_fill_manual(values = c("orangered3", "orangered1", "royalblue3", "royalblue1"), name = NULL) +
    scale_y_continuous(sec.axis = sec_axis(~ . * gene_scale, name = "Genes")) +
    labs(x = "Linkage group", y = "DML") + theme_wgbs()
  fig <- pA / pB + plot_annotation(tag_levels = "A")
  save_figure(fig, output_path(step, "figures", sprintf("figure_3_%s", mode)), width = 10, height = 10)
  rm(ov, gr, cls); gc(verbose = FALSE)
}

finish_step(step, t0)

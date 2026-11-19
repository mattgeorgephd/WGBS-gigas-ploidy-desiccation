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
# Run       Rscript 01_pipeline/07_genome_features.R     (after step 05; under a minute, about 1 GB of memory)
#
# Inputs    04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz
#           04_output/05_differential_methylation/tables/dml_{ploidy,heat,shared}_<mode>.csv
#           03_input_files/genome_features/*.gz   the 13 tracks; what each holds is in the README there
# Outputs   04_output/07_genome_features/tables/
#             feature_overlap_by_track_<mode>.csv   loci of each set that overlap each track: n_overlap,
#                                                   pct_of_loci, legacy_share_pct; track_label says
#                                                   what the track holds
#             feature_class_exclusive_<mode>.csv    one class per locus, per set: n, pct, and pct_in_TE,
#                                                   the percent of the loci of the class that also
#                                                   overlap the repeat track (TE)
#             feature_association_tests_<mode>.csv  Fisher test per track, DML set and background;
#                                                   with track_label
#             feature_class_chisq_<mode>.csv        chi-squared test of DML status by exclusive class,
#                                                   per DML set, against all analysed CpGs only
#             dml_by_linkage_group_<mode>.csv       analysed CpGs, DML of each list (all, hyper, hypo),
#                                                   loci in both lists and genes, per linkage group
#             dml_regions_<mode>.csv                one row per region of a set: position, n_dml,
#                                                   n_hyper, n_hypo
#             dml_region_summary_<mode>.csv         per set: regions, single-locus regions, the
#                                                   largest region, DML in regions of 5 or more
#           04_output/07_genome_features/figures/figure_3_<mode>.(png|pdf)
#           04_output/07_genome_features/session_info.txt
#
# Manuscript items: Results 3.3; Figure 3.
# Notes     - The step applies no threshold of its own. A DML is a row of
#             dml_<set>_<mode>.csv as step 05 wrote it (q-value below
#             config$dml_qvalue and an absolute difference above
#             config$dml_difference, by methylKit's rule on the unrounded
#             difference).
#           - The background is the set of CpGs that were tested for differential
#             methylation, in the same mode. The legacy background was assembled
#             from two other CpG sets (issue I-14).
#           - Two backgrounds are used for the Fisher tests: all analysed CpGs,
#             and "eligible" CpGs, those whose mean methylation across libraries
#             lies between half of config$dml_difference and 100 minus that (10%
#             and 90% at the default of 20 points, both ends included). A
#             20-point difference between two groups almost always needs a mean
#             in that range (all but 2 of the legacy DML have one), so the other
#             CpGs, about 86% of those analysed, are almost never DML.
#             Methylation in this genome sits mostly in gene bodies, so
#             enrichment of DML in genes relative to all CpGs is expected
#             whatever the treatment. The eligible background removes most of
#             that and not all of it: inside the range the chance of passing the
#             threshold still depends on the mean and on the variation among
#             animals, and CpG lists made with shuffled ploidy labels show the
#             same odds ratios against eligible CpGs as the real DML (issues.md,
#             I-41). Neither background is a test of a treatment effect on
#             location. Every DML of a set is counted in the DML row of its
#             2 x 2 tables, also the few that lie outside the eligible range.
#           - The chi-squared test of DML status by exclusive class compares
#             each DML set with all analysed CpGs only (column background of
#             its table); no version against eligible CpGs is written.
#           - "legacy_share_pct" reproduces the manuscript's percentages: the
#             overlap count of a track divided by the sum of the counts of eight
#             tracks. Those shares double count loci that lie in two tracks.
#           - Class priority for the exclusive assignment: CDS, lncRNA, UTR,
#             intron, upstream flank, downstream flank, intergenic; anything left
#             is "other". These tracks cover every base of the 236 nuclear
#             sequences and no base of the mitochondrial genome, which has no
#             record in any track, so "other" is exactly the mitochondrial CpGs:
#             363 in legacy mode and none in corrected mode, where the coverage
#             filter removes them. That is why the chi-squared test has 7 degrees
#             of freedom in legacy mode and 6 in corrected mode. lncRNA comes
#             before UTR and intron because the UTR track (exon minus CDS) and
#             the intron track include the exons and introns of lncRNA genes.
#             Overlap with the repeat track is reported alongside (pct_in_TE),
#             not as a class, because repeats sit inside all the other classes.
#           - The track called TE is not a track of transposable elements and
#             not a true BED file. It is the RepeatMasker output of the assembly:
#             repeats of every class, simple repeats, low complexity, satellites
#             and rRNA included (24% of its 781,312 records, 4% of their bases),
#             with 1-based coordinates. The legacy bedtools run read it as BED,
#             which drops the first base of every interval, and this step reads
#             it the same way in both modes, because the numbers of the draft
#             come from that reading; nothing is corrected (scope decision of
#             2026-11-19). The lost base matters for few CpGs. In legacy mode
#             101 of the 480,579 analysed CpGs touch a repeat with that base
#             alone and are counted outside the track (2 of the 2,351 heat DML,
#             no ploidy DML, no locus in both lists); in corrected mode 92 of
#             454,555 (1 heat DML). The pipeline does not write these counts;
#             they come from the code review and its check
#             (07_documentation/agent_workflows/verification_scripts/
#             R2-figures-features/). The column track_label of the two tables
#             written by track says what TE holds. The track keeps the name TE,
#             and the column pct_in_TE with it, because the scripts of
#             06_diagnostics and the values map of step 13 select the track by
#             that name. See 03_input_files/genome_features/README.md, "The
#             repeat track".
#           - "hyper" and "hypo" are the sign of the difference of the list:
#             triploid minus diploid for ploidy, heat minus control for heat. A
#             locus in both lists has two differences. Its columns here,
#             shared_hyper and shared_hypo of dml_by_linkage_group and n_hyper
#             and n_hypo of the set "shared" in dml_regions, are the sign of the
#             HEAT difference alone, as in the threshold table of step 05. They
#             are not the classes both_hyper, both_hypo and opposite of
#             dml_shared_<mode>.csv: in legacy mode 106 of the 200 loci have a
#             positive and 94 a negative heat difference, while 28 are positive
#             in both contrasts, 49 negative in both and 123 of opposite sign.
#             The three classes are the sets "shared, both hyper", "shared, both
#             hypo" and "shared, opposite" of the first two tables and of
#             Figure 3A.
#           - Figure 3B draws the ten linkage groups only and stacks the ploidy
#             and the heat list, so a locus in both lists is counted twice. Its
#             caption gives the number of such loci in the panel and what lies on
#             the other sequences, which are in dml_by_linkage_group_<mode>.csv.
#           - A set can be empty (a high config$dml_difference, or a stricter
#             model). For a set without loci the percentages and the Fisher and
#             chi-squared rows hold NA, the region summary has a row of zeros,
#             and Figure 3A keeps its label without a bar; the log names every
#             such set. Nothing is tested for it. A set of one or a few loci is
#             tested as it stands; dml_total of the Fisher table and
#             min_expected of the chi-squared table show how thin such a test
#             is. With no DML on any linkage group, Figure 3B is replaced by a
#             panel that says so.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "GenomicRanges", "ggplot2", "patchwork", "RColorBrewer"))

step <- "07_genome_features"
t0 <- start_step(step, "genomic location of analysed CpGs and DML")
dml_dir <- function(f) output_path("05_differential_methylation", "tables", f)
tab <- function(name) output_path(step, "tables", name)
pct <- function(n, of) if (of > 0) 100 * n / of else rep(NA_real_, length(n))    # no percentage of nothing

tracks <- names(feature_track_files())
track_labels <- feature_track_labels()
stopifnot(setequal(names(track_labels), tracks))
track_gr <- lapply(setNames(tracks, tracks), function(t) { log_msg("reading track %s", t); read_feature_track(t) })
legacy_eight <- c("exon_UTR", "CDS", "intron", "upstream_flank", "downstream_flank", "lncRNA", "TE", "intergenic")
class_priority <- c("CDS", "lncRNA", "exon_UTR", "intron", "upstream_flank", "downstream_flank", "intergenic")
class_labels <- c(CDS = "CDS", lncRNA = "lncRNA", exon_UTR = "UTR or other non-coding exon", intron = "Intron",
                  upstream_flank = "Upstream flank", downstream_flank = "Downstream flank",
                  intergenic = "Intergenic", other = "Mitochondrial genome\n(in no track)")
# One colour per class, the same in both modes whether or not the mode has a CpG of the class.
class_levels <- rev(c(class_priority, "other"))
class_cols <- setNames(rev(brewer.pal(length(class_levels), "YlGnBu")), class_levels)
dml_sets <- c("ploidy DML", "heat DML", "shared")
linkage_group_names <- sprintf("LG%d", seq_along(config$linkage_groups))

# Genes per linkage group, for the chromosome panel.
gene_gr <- read_feature_track("gene", keep_attributes = TRUE)
genes_per_seq <- as.data.table(table(sequence = as.character(GenomicRanges::seqnames(gene_gr))))
setnames(genes_per_seq, "N", "n_genes")

for (mode in config$modes) {
  log_msg("[%s] reading analysed CpGs and DML", mode)
  united <- fread(dml_dir(sprintf("united_cpg_%s.csv.gz", mode)))
  # Eligible CpGs: mean methylation in the range that holds nearly all DML (see dml_eligible()).
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
  # A set without loci has no percentage and cannot be tested. Its rows are written with
  # NA so that every table keeps its shape, and the log says which sets are concerned.
  for (s in names(sets)[!vapply(sets, any, logical(1))]) {
    log_msg("[%s] the set '%s' has no locus: its percentages and tests are written as NA, and it has no bar in Figure 3A", mode, s)
  }

  # ---- Overlap with each track ------------------------------------------------
  gr <- cpg_ranges(cpg$chr, cpg$start)
  ov <- vapply(tracks, function(t) overlaps_track(gr, track_gr[[t]]), logical(nrow(cpg)))
  by_track <- rbindlist(lapply(names(sets), function(s) {
    n <- colSums(ov[sets[[s]], , drop = FALSE])
    data.table(mode = mode, set = s, n_loci = sum(sets[[s]]), track = tracks, n_overlap = n,
               pct_of_loci = pct(n, sum(sets[[s]])),
               legacy_share_pct = ifelse(tracks %in% legacy_eight, pct(n, sum(n[legacy_eight])), NA_real_),
               track_label = unname(track_labels[tracks]))
  }))
  write_table(by_track, tab(sprintf("feature_overlap_by_track_%s.csv", mode)))

  # ---- One class per CpG ---------------------------------------------------------
  cls <- rep("other", nrow(cpg))
  for (t in rev(class_priority)) cls[ov[, t]] <- t          # highest priority assigned last
  cls <- factor(cls, levels = c(class_priority, "other"))
  exclusive <- rbindlist(lapply(names(sets), function(s) {
    n <- as.vector(table(cls[sets[[s]]]))
    data.table(mode = mode, set = s, n_loci = sum(sets[[s]]), class = levels(cls), n = n,
               pct = pct(n, sum(sets[[s]])),
               pct_in_TE = 100 * as.vector(tapply(ov[sets[[s]], "TE"], cls[sets[[s]]], mean)))
  }))
  write_table(exclusive, tab(sprintf("feature_class_exclusive_%s.csv", mode)))

  # ---- Association tests -----------------------------------------------------------
  # Per track: is a DML more or less likely than other analysed CpGs to lie in it?
  # Loci are treated as independent, which they are not where DML cluster, so
  # these p-values are optimistic; see the region summary below.
  # Two backgrounds: every analysed CpG, and the eligible CpGs. Most CpGs in this
  # genome are unmethylated and almost never reach the difference threshold, and
  # methylation is concentrated in gene bodies, so the second background is the
  # fairer of the two. It is not a fair one either: see the header and issues.md, I-41.
  # Rows of each 2 x 2 table: every DML of the set, also the few whose mean lies outside
  # the eligible range, and the CpGs of the background that are not in the set.
  backgrounds <- list(`all analysed CpGs` = rep(TRUE, nrow(cpg)), `eligible CpGs` = eligible)
  tests <- rbindlist(lapply(names(backgrounds), function(b) {
    rbindlist(lapply(dml_sets, function(s) {
      rbindlist(lapply(c(legacy_eight, "gene"), function(t) {
        use <- backgrounds[[b]] | sets[[s]]
        m <- table(factor(sets[[s]][use], c(TRUE, FALSE)), factor(ov[use, t], c(TRUE, FALSE)))
        # fisher.test() on a set without loci would return an odds ratio of 0 with p = 1.
        ft <- if (sum(m[1, ]) > 0) fisher.test(m) else
          list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_)
        data.table(mode = mode, background = b, set = s, track = t, dml_in_track = m[1, 1], dml_total = sum(m[1, ]),
                   pct_dml = pct(m[1, 1], sum(m[1, ])), pct_other_cpg = pct(m[2, 1], sum(m[2, ])),
                   odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
                   p_fisher = ft$p.value)
      }))
    }))
  }))
  tests[, q_BH := p.adjust(p_fisher, "BH")]
  tests[, track_label := unname(track_labels[track])]
  # The chi-squared test is run against all analysed CpGs only. Without a 2 x k table (a
  # set without loci leaves one row) chisq.test() would silently run a goodness-of-fit
  # test of the class counts against equal shares; the row holds NA instead.
  chisq <- rbindlist(lapply(dml_sets, function(s) {
    m <- table(sets[[s]], droplevels(cls))
    ct <- if (all(dim(m) >= 2)) suppressWarnings(chisq.test(m)) else
      list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_, expected = NA_real_)
    data.table(mode = mode, background = "all analysed CpGs", set = s, test = "chi-squared, DML status by exclusive class",
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
  lg_levels <- c(linkage_group_names, "unplaced", "mitochondrial")
  genes_lg <- genes_per_seq[, .(n_genes = sum(n_genes)), by = .(sequence = seq_label(sequence))]
  per_lg <- cpg[, .(n_cpg_analysed = .N), by = .(sequence = seq_label(chr))]
  # n_hyper and n_hypo: sign of the difference in `col`. For the loci in both lists that
  # is the heat difference alone (see the header), as in the threshold table of step 05.
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
                pct_of_ploidy_dml = pct(ploidy_dml, sum(ploidy_dml)), pct_of_heat_dml = pct(heat_dml, sum(heat_dml)),
                ploidy_dml_per_1000_cpg = 1000 * ploidy_dml / n_cpg_analysed,
                heat_dml_per_1000_cpg = 1000 * heat_dml / n_cpg_analysed)]
  write_table(per_lg, tab(sprintf("dml_by_linkage_group_%s.csv", mode)))

  # ---- Regions: DML of one set within 1 kb of each other ------------------------------
  make_regions <- function(d, col, name) {
    if (nrow(d) == 0) {                           # a set without loci has no region
      return(data.table(set = character(), chr = character(), start = integer(), end = integer(), span_bp = integer(),
                        n_dml = integer(), n_hyper = integer(), n_hypo = integer()))
    }
    d <- copy(d)[order(chr, start)]
    d[, region := cumsum(c(TRUE, diff(start) > 1000 | chr[-1] != chr[-.N]))]
    d[, .(set = name, chr = chr[1], start = min(start), end = max(end), span_bp = max(end) - min(start),
          n_dml = .N, n_hyper = sum(get(col) > 0), n_hypo = sum(get(col) < 0)), by = region][, region := NULL]
  }
  regions <- rbind(make_regions(dml$ploidy, "meth.diff", "ploidy DML"), make_regions(dml$heat, "meth.diff", "heat DML"),
                   make_regions(shared, "meth.diff_heat", "shared"))
  regions[, mode := mode]; setorder(regions, set, -n_dml)
  write_table(regions, tab(sprintf("dml_regions_%s.csv", mode)))
  # One row per set in the order of the regions table, also for a set without loci,
  # which gets zeros.
  region_summary <- rbindlist(lapply(c("heat DML", "ploidy DML", "shared"), function(s) {
    r <- regions[set == s]
    data.table(mode = mode, set = s, n_dml = sum(r$n_dml), n_regions = nrow(r), n_single_locus_regions = sum(r$n_dml == 1),
               largest_region_n_dml = if (nrow(r) > 0) max(r$n_dml) else 0L,
               dml_in_regions_of_5_or_more = sum(r$n_dml[r$n_dml >= 5]))
  }))
  write_table(region_summary, tab(sprintf("dml_region_summary_%s.csv", mode)))
  print(region_summary)

  # ---- Figure 3 -------------------------------------------------------------------------
  ex <- exclusive[set != "shared"]                      # the three shared classes are shown instead
  ex[, set_label := sprintf("%s\nn = %s", sub(", ", ",\n", sub(" CpGs", "\nCpGs", set)),
                            formatC(n_loci, big.mark = ",", format = "d"))]
  ex[, set_label := factor(set_label, levels = unique(set_label))]
  ex[, class := factor(class, levels = class_levels)]
  # A set without loci keeps its place and its label and has no bar. A class without a CpG
  # in this mode (the mitochondrial class in corrected mode) gets no key in the legend.
  class_present <- class_levels[class_levels %in% ex[n > 0, as.character(class)]]
  pA <- ggplot(ex[n_loci > 0], aes(set_label, pct, fill = class)) + geom_col(colour = "black", linewidth = 0.3, width = 0.8) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = class_cols, breaks = class_present, labels = unname(class_labels[class_present]),
                      name = "Genomic class") +
    labs(x = NULL, y = "CpG loci (%)") + theme_wgbs() + theme(axis.text.x = element_text(size = 8))
  on_lg <- per_lg[sequence %in% linkage_group_names]
  off_lg <- per_lg[!sequence %in% linkage_group_names]      # unplaced scaffolds and the mitochondrial genome: not drawn
  lg_long <- melt(on_lg, id.vars = "sequence",
                  measure.vars = c("ploidy_dml_hypo", "ploidy_dml_hyper", "heat_dml_hypo", "heat_dml_hyper"))
  lg_long[, variable := factor(variable, levels = c("heat_dml_hyper", "heat_dml_hypo", "ploidy_dml_hyper", "ploidy_dml_hypo"),
                               labels = c("Heat, hyper", "Heat, hypo", "Ploidy, hyper", "Ploidy, hypo"))]
  if (sum(lg_long$value) == 0) {
    # Nothing to stack and no scale for the gene line: the panel says so.
    log_msg("[%s] Figure 3B: no DML on a linkage group; the panel is left empty", mode)
    pB <- ggplot() + annotate("text", x = 0, y = 0, label = "no DML on a linkage group", size = 3.6) + theme_void() +
      theme(plot.tag = element_text(face = "bold"))
  } else {
    gene_scale <- max(on_lg$n_genes) / max(lg_long[, sum(value), by = sequence]$V1)
    pB <- ggplot(lg_long, aes(sequence, value)) +
      geom_col(aes(fill = variable), colour = "black", linewidth = 0.3) +
      geom_line(data = on_lg, aes(sequence, n_genes / gene_scale, group = 1), linewidth = 0.7) +
      geom_point(data = on_lg, aes(sequence, n_genes / gene_scale), size = 2.4) +
      scale_fill_manual(values = c("orangered3", "orangered1", "royalblue3", "royalblue1"), name = NULL) +
      scale_y_continuous(sec.axis = sec_axis(~ . * gene_scale, name = "Genes")) +
      labs(x = "Linkage group", y = "DML (ploidy and heat lists stacked)",
           caption = sprintf(paste0("The two lists are stacked, so a locus in both lists is counted twice (such loci drawn here: %s).\n",
                                    "Not drawn: %s ploidy DML, %s heat DML and %s genes on sequences outside the ten linkage groups."),
                             format(sum(on_lg$shared), big.mark = ","), format(sum(off_lg$ploidy_dml), big.mark = ","),
                             format(sum(off_lg$heat_dml), big.mark = ","), format(sum(off_lg$n_genes), big.mark = ","))) +
      theme_wgbs() + theme(plot.caption = element_text(size = 9, hjust = 0))
  }
  fig <- pA / pB + plot_annotation(tag_levels = "A")
  save_figure(fig, output_path(step, "figures", sprintf("figure_3_%s", mode)), width = 10, height = 10)
  rm(ov, gr, cls); gc(verbose = FALSE)
}

finish_step(step, t0)

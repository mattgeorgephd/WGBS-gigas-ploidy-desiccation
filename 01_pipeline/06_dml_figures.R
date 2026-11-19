# 06_dml_figures.R
#
# Purpose   Figures that describe the DML sets.
#           Figure 2  A  number of DML against the difference threshold
#                     B  percent of DML that are hypermethylated, against the
#                        same thresholds; a set is drawn only where it has at
#                        least 10 DML
#                     C  PCA of the libraries on the ploidy DML
#                     D  PCA of the libraries on the heat DML
#                     E, F, G  heat maps of percent methylation at the ploidy
#                        DML, the heat DML, and the loci in both lists
#           Figure 4  for each locus in both lists, its ploidy difference and
#                     its heat difference
#
# Run       Rscript 01_pipeline/06_dml_figures.R     (after step 05; under a minute, about 0.5 GB)
# Inputs    04_output/05_differential_methylation/tables/
#             united_cpg_<mode>.csv.gz, dml_{ploidy,heat,shared}_<mode>.csv,
#             dml_counts_by_threshold_<mode>.csv
#           03_input_files/metadata/wgbs_libraries.csv
# Outputs   04_output/06_dml_figures/tables/
#             pca_dml_<mode>.csv            PC 1 and PC 2 of each library, on the ploidy DML
#                                           and on the heat DML
#             pca_dml_variance_<mode>.csv   loci in the set, loci used, percent of variance
#                                           of PC 1 to PC 4
#           04_output/06_dml_figures/figures/figure_2_<mode>.(png|pdf), figure_4_<mode>.(png|pdf)
#           04_output/06_dml_figures/session_info.txt
#
# Manuscript items: Figure 2; Figure 4.
# Notes     - The step applies no threshold of its own. A DML is a row of
#             dml_<set>_<mode>.csv as step 05 wrote it (q-value below
#             config$dml_qvalue and an absolute difference above
#             config$dml_difference, by methylKit's rule on the unrounded
#             difference), and panels A and B draw the counts of
#             dml_counts_by_threshold_<mode>.csv as they are.
#           - Panels C and D use loci that were selected because they differ
#             between the groups shown, so the separation of those groups is
#             expected and is not evidence about the loci (issue I-23). Both
#             factors are drawn in each panel so that the other grouping can be
#             read as well.
#           - Loci are matched to the methylation table by row, that is by
#             sequence and position. The legacy code matched on position alone.
#           - PCA follows methylKit::PCASamples: loci with SD above the median
#             SD, centred and scaled.
#           - "Hypermethylated" is the sign of the difference of the list: higher
#             in triploids for ploidy, higher in the heat cohort for heat. A
#             locus in both lists has two differences, and step 05 counts it by
#             the sign of its heat difference alone. The series "both lists" of
#             panel B is therefore not the share of loci that are higher in both
#             contrasts (legacy mode at 20 points: 106 of the 200 loci have a
#             positive heat difference, 28 are positive in both). The caption of
#             the panel says which sign is counted; Figure 4 shows both
#             differences of every such locus.
#           - Panel B leaves a set out at the thresholds where it has fewer than
#             10 DML, because a percentage of a handful of loci is noise; the
#             caption says so, and the axis keeps the range of panel A.
#           - The heat maps carry no names beside their two annotation strips
#             (beside the strips the names were cut by the next panel). The upper
#             strip is treatment and the lower one ploidy, as the legend of the
#             last heat map lists them.
#           - A DML set can be too small to draw (a high config$dml_difference,
#             or a stricter model): the PCA needs two loci above the median SD,
#             a clustered heat map two loci, Figure 4 one locus. Such a panel is
#             replaced by one that says so, with a message in the log; the PCA
#             tables then hold NA for that set, and the step finishes.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "ggplot2", "patchwork", "pheatmap", "RColorBrewer"))

step <- "06_dml_figures"
t0 <- start_step(step, "Figure 2 and Figure 4")
lib <- read_libraries()
dml_dir <- function(f) output_path("05_differential_methylation", "tables", f)
set_labels <- c(ploidy = "ploidy", heat = "heat", shared = "both lists")
min_dml_for_percent <- 10            # panel B: a percentage of fewer loci is not drawn

for (mode in config$modes) {
  united <- fread(dml_dir(sprintf("united_cpg_%s.csv.gz", mode)))
  n_lib <- nrow(lib)
  # Percent methylation per library from the normalized counts, as methylKit::percMethylation.
  pm <- vapply(seq_len(n_lib), function(j) 100 * united[[paste0("numCs", j)]] /
                 (united[[paste0("numCs", j)]] + united[[paste0("numTs", j)]]), numeric(nrow(united)))
  colnames(pm) <- lib$library_id
  dml <- list(ploidy = fread(dml_dir(sprintf("dml_ploidy_%s.csv", mode))),
              heat   = fread(dml_dir(sprintf("dml_heat_%s.csv", mode))),
              shared = fread(dml_dir(sprintf("dml_shared_%s.csv", mode))))
  sweep <- fread(dml_dir(sprintf("dml_counts_by_threshold_%s.csv", mode)))
  rm(united); gc(verbose = FALSE)

  # A panel that cannot be drawn because its set is too small: the reason goes to the
  # log and into the panel, so that the figure is complete and nothing is hidden.
  note_panel <- function(title, why) {
    log_msg("[%s] %s: %s; the panel is left empty", mode, title, why)
    ggplot() + annotate("text", x = 0, y = 0, label = paste(strwrap(why, width = 30), collapse = "\n"), size = 3.6) +
      labs(title = title) + theme_void() +
      theme(plot.title = element_text(size = 11, face = "italic", hjust = 0.5), plot.tag = element_text(face = "bold"))
  }

  # ---- A, B: threshold curves --------------------------------------------------
  sweep[, set := factor(set, levels = names(set_labels), labels = set_labels)]
  cols <- setNames(unname(config$colours$dml_set), set_labels)
  pA <- ggplot(sweep, aes(difference, n, colour = set, shape = set)) + geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    scale_colour_manual(values = cols, name = NULL) + scale_shape_discrete(name = NULL) +
    scale_y_continuous(labels = function(v) format(v, big.mark = ",")) +
    labs(x = "Difference threshold (%)", y = "DML") + theme_wgbs()
  # The x axis keeps the range of panel A although the sets with few DML are left out,
  # and every set keeps its colour and symbol when another set has no point at all.
  pB <- ggplot(sweep[n >= min_dml_for_percent], aes(difference, pct_hyper, colour = set, shape = set)) +
    geom_hline(yintercept = 50, linetype = 3) + geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    scale_colour_manual(values = cols, name = NULL, drop = FALSE) + scale_shape_discrete(name = NULL, drop = FALSE) +
    scale_x_continuous(limits = range(sweep$difference)) + ylim(0, 100) +
    labs(x = "Difference threshold (%)", y = "Hypermethylated DML (%)",
         caption = sprintf("Both lists: loci counted by the sign of their\nheat difference. A set is drawn where it has\nat least %d DML.",
                           min_dml_for_percent)) +
    theme_wgbs() + theme(plot.caption = element_text(size = 9, hjust = 0))

  # ---- C, D: PCA on the DML of one contrast ------------------------------------
  pca_on <- function(rows, name) {
    m <- pm[rows, , drop = FALSE]
    sds <- apply(m, 1, sd)
    used <- if (nrow(m) > 0) which(sds > quantile(sds, 0.5)) else integer(0)
    if (length(used) < 2) {                       # two components need two loci
      return(list(scores = data.table(mode = mode, loci = name, library_id = lib$library_id, PC1 = NA_real_, PC2 = NA_real_,
                                      lib[, c("ploidy", "treatment", "cohort", "acute_heat_shock")]),
                  variance = data.table(mode = mode, loci = name, n_loci = nrow(m), n_loci_used = length(used),
                                        pc = 1:4, percent_variance = NA_real_)))
    }
    pc <- prcomp(t(m[used, , drop = FALSE]), center = TRUE, scale. = TRUE)
    ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    list(scores = data.table(mode = mode, loci = name, library_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                             lib[, c("ploidy", "treatment", "cohort", "acute_heat_shock")]),
         variance = data.table(mode = mode, loci = name, n_loci = nrow(m), n_loci_used = length(used),
                               pc = 1:4, percent_variance = ve[1:4]))
  }
  pcs <- list(ploidy = pca_on(dml$ploidy$united_row, "ploidy DML"), heat = pca_on(dml$heat$united_row, "heat DML"))
  write_table(rbindlist(lapply(pcs, `[[`, "scores")), output_path(step, "tables", sprintf("pca_dml_%s.csv", mode)))
  write_table(rbindlist(lapply(pcs, `[[`, "variance")), output_path(step, "tables", sprintf("pca_dml_variance_%s.csv", mode)))
  pca_panel <- function(p, title) {
    if (anyNA(p$scores$PC1)) {
      return(note_panel(title, sprintf("too few loci for a PCA (%d in the set, %d above the median SD)",
                                       p$variance$n_loci[1], p$variance$n_loci_used[1])))
    }
    ggplot(p$scores, aes(PC1, PC2, colour = ploidy, shape = treatment)) + geom_point(size = 3) +
      scale_colour_manual(values = config$colours$ploidy, name = NULL) +
      scale_shape_manual(values = c(control = 16, heat = 17), name = NULL) +
      labs(title = title, x = sprintf("PC 1 (%.1f%%)", p$variance$percent_variance[1]),
           y = sprintf("PC 2 (%.1f%%)", p$variance$percent_variance[2])) +
      theme_wgbs() + theme(plot.title = element_text(size = 11, face = "italic"))
  }
  pC <- pca_panel(pcs$ploidy, "Ploidy DML"); pD <- pca_panel(pcs$heat, "Heat DML")
  for (p in pcs) {
    if (!anyNA(p$scores$PC1)) {
      log_msg("[%s] PCA on %s: PC1 %.1f%%, PC2 %.1f%%", mode, p$variance$loci[1],
              p$variance$percent_variance[1], p$variance$percent_variance[2])
    }
  }

  # ---- E, F, G: heat maps --------------------------------------------------------
  ann <- data.frame(ploidy = lib$ploidy, treatment = lib$treatment, row.names = lib$library_id)
  ann_cols <- list(ploidy = config$colours$ploidy, treatment = config$colours$treatment)
  # The names of the annotation strips are not drawn beside them: there the next panel cut
  # "treatment" short and the colour key ran into "ploidy". The legend names both strips.
  heat_panel <- function(rows, title, legend) {
    if (length(rows) < 2) {                       # rows cannot be clustered
      return(wrap_elements(full = note_panel(title, sprintf("too few loci for a clustered heat map (%d)", length(rows)))))
    }
    h <- pheatmap(pm[rows, , drop = FALSE], color = rev(brewer.pal(8, "RdYlBu")),
                  clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan",
                  clustering_method = "complete", show_rownames = FALSE, show_colnames = FALSE,
                  annotation_col = ann, annotation_colors = ann_cols, annotation_legend = legend, legend = legend,
                  annotation_names_col = FALSE,
                  treeheight_row = 25, treeheight_col = 20, main = title, fontsize = 9, silent = TRUE)
    wrap_elements(full = h$gtable)
  }
  # The colour key and the annotation legend go on the last heat map that can be drawn
  # (panel G unless its set is too small), and that panel gets the room for them.
  n_dml <- vapply(dml, nrow, integer(1))
  key_on <- tail(names(n_dml)[n_dml >= 2], 1)
  pE <- heat_panel(dml$ploidy$united_row, sprintf("Ploidy DML (n = %s)", format(nrow(dml$ploidy), big.mark = ",")),
                   identical(key_on, "ploidy"))
  pF <- heat_panel(dml$heat$united_row, sprintf("Heat DML (n = %s)", format(nrow(dml$heat), big.mark = ",")),
                   identical(key_on, "heat"))
  pG <- heat_panel(dml$shared$united_row, sprintf("In both lists (n = %d)", nrow(dml$shared)), identical(key_on, "shared"))

  fig2 <- ((pA + pB + pC + pD) + plot_layout(nrow = 1)) /
    ((pE + pF + pG) + plot_layout(nrow = 1, widths = ifelse(names(dml) %in% key_on, 1.25, 1))) +
    plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = "A")
  save_figure(fig2, output_path(step, "figures", sprintf("figure_2_%s", mode)), width = 16, height = 10, dpi = 300)

  # ---- Figure 4: the two differences at each locus in both lists ------------------
  s <- copy(dml$shared)
  if (nrow(s) == 0) {
    fig4 <- note_panel("Loci in both DML lists", "no locus is in both lists")
  } else {
    s[, class := factor(class, levels = c("both_hyper", "opposite", "both_hypo"),
                        labels = c("Both hypermethylated", "Opposite sign", "Both hypomethylated"))]
    setorder(s, class, meth.diff_heat)
    s[, locus := .I]
    long <- melt(s, id.vars = c("locus", "class"), measure.vars = c("meth.diff_ploidy", "meth.diff_heat"))
    long[, contrast := factor(variable, levels = c("meth.diff_ploidy", "meth.diff_heat"),
                              labels = c("triploid minus diploid", "heat minus control"))]
    fig4 <- ggplot(long, aes(value, locus, colour = contrast)) +
      geom_segment(aes(x = 0, xend = value, y = locus, yend = locus), linewidth = 0.35) + geom_point(size = 0.9) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      facet_grid(class ~ ., scales = "free_y", space = "free_y") +
      scale_colour_manual(values = c(`triploid minus diploid` = "royalblue2", `heat minus control` = "darkorange"), name = NULL) +
      labs(x = "Methylation difference (percentage points)", y = "Loci in both DML lists") +
      theme_wgbs() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top",
                           strip.text.y = element_text(angle = 0, hjust = 0))
  }
  save_figure(fig4, output_path(step, "figures", sprintf("figure_4_%s", mode)), width = 7.5, height = 8)
  rm(pm); gc(verbose = FALSE)
}

finish_step(step, t0)

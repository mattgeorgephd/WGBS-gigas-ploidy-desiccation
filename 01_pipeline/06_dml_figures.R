# 06_dml_figures.R
#
# Purpose   Figures that describe the DML sets.
#           Figure 2  A  number of DML against the difference threshold
#                     B  percent of DML that are hypermethylated, same axis
#                     C  PCA of the libraries on the ploidy DML
#                     D  PCA of the libraries on the heat DML
#                     E, F, G  heat maps of percent methylation at the ploidy
#                        DML, the heat DML, and the loci in both lists
#           Figure 4  for each locus in both lists, its ploidy difference and
#                     its heat difference
#
# Inputs    04_output/05_differential_methylation/tables/
#             united_cpg_<mode>.csv.gz, dml_{ploidy,heat,shared}_<mode>.csv,
#             dml_counts_by_threshold_<mode>.csv
#           03_input_files/metadata/wgbs_libraries.csv
# Outputs   04_output/06_dml_figures/tables/pca_dml_<mode>.csv, pca_dml_variance_<mode>.csv
#           04_output/06_dml_figures/figures/figure_2_<mode>.(png|pdf), figure_4_<mode>.(png|pdf)
#
# Manuscript items: Figure 2; Figure 4.
# Notes     - Panels C and D use loci that were selected because they differ
#             between the groups shown, so the separation of those groups is
#             expected and is not evidence about the loci (issue I-23). Both
#             factors are drawn in each panel so that the other grouping can be
#             read as well.
#           - Loci are matched to the methylation table by row, that is by
#             sequence and position. The legacy code matched on position alone.
#           - PCA follows methylKit::PCASamples: loci with SD above the median
#             SD, centred and scaled.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "ggplot2", "patchwork", "pheatmap", "RColorBrewer"))

step <- "06_dml_figures"
t0 <- start_step(step, "Figure 2 and Figure 4")
lib <- read_libraries()
dml_dir <- function(f) output_path("05_differential_methylation", "tables", f)
set_labels <- c(ploidy = "ploidy", heat = "heat", shared = "both lists")

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

  # ---- A, B: threshold curves --------------------------------------------------
  sweep[, set := factor(set, levels = names(set_labels), labels = set_labels)]
  cols <- setNames(unname(config$colours$dml_set), set_labels)
  pA <- ggplot(sweep, aes(difference, n, colour = set, shape = set)) + geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    scale_colour_manual(values = cols, name = NULL) + scale_shape_discrete(name = NULL) +
    scale_y_continuous(labels = function(v) format(v, big.mark = ",")) +
    labs(x = "Difference threshold (%)", y = "DML") + theme_wgbs()
  pB <- ggplot(sweep[n >= 10], aes(difference, pct_hyper, colour = set, shape = set)) +
    geom_hline(yintercept = 50, linetype = 3) + geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    scale_colour_manual(values = cols, name = NULL) + scale_shape_discrete(name = NULL) + ylim(0, 100) +
    labs(x = "Difference threshold (%)", y = "Hypermethylated DML (%)") + theme_wgbs()

  # ---- C, D: PCA on the DML of one contrast ------------------------------------
  pca_on <- function(rows, name) {
    m <- pm[rows, , drop = FALSE]
    sds <- apply(m, 1, sd)
    pc <- prcomp(t(m[sds > quantile(sds, 0.5), , drop = FALSE]), center = TRUE, scale. = TRUE)
    ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    list(scores = data.table(mode = mode, loci = name, library_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                             lib[, c("ploidy", "treatment", "cohort", "acute_heat_shock")]),
         variance = data.table(mode = mode, loci = name, n_loci = nrow(m), n_loci_used = sum(sds > quantile(sds, 0.5)),
                               pc = 1:4, percent_variance = ve[1:4]))
  }
  pcs <- list(ploidy = pca_on(dml$ploidy$united_row, "ploidy DML"), heat = pca_on(dml$heat$united_row, "heat DML"))
  write_table(rbindlist(lapply(pcs, `[[`, "scores")), output_path(step, "tables", sprintf("pca_dml_%s.csv", mode)))
  write_table(rbindlist(lapply(pcs, `[[`, "variance")), output_path(step, "tables", sprintf("pca_dml_variance_%s.csv", mode)))
  pca_panel <- function(p, title) {
    ggplot(p$scores, aes(PC1, PC2, colour = ploidy, shape = treatment)) + geom_point(size = 3) +
      scale_colour_manual(values = config$colours$ploidy, name = NULL) +
      scale_shape_manual(values = c(control = 16, heat = 17), name = NULL) +
      labs(title = title, x = sprintf("PC 1 (%.1f%%)", p$variance$percent_variance[1]),
           y = sprintf("PC 2 (%.1f%%)", p$variance$percent_variance[2])) +
      theme_wgbs() + theme(plot.title = element_text(size = 11, face = "italic"))
  }
  pC <- pca_panel(pcs$ploidy, "Ploidy DML"); pD <- pca_panel(pcs$heat, "Heat DML")
  log_msg("[%s] PCA on ploidy DML: PC1 %.1f%%, PC2 %.1f%%; on heat DML: PC1 %.1f%%, PC2 %.1f%%", mode,
          pcs$ploidy$variance$percent_variance[1], pcs$ploidy$variance$percent_variance[2],
          pcs$heat$variance$percent_variance[1], pcs$heat$variance$percent_variance[2])

  # ---- E, F, G: heat maps --------------------------------------------------------
  ann <- data.frame(ploidy = lib$ploidy, treatment = lib$treatment, row.names = lib$library_id)
  ann_cols <- list(ploidy = config$colours$ploidy, treatment = config$colours$treatment)
  heat_panel <- function(rows, title, legend) {
    h <- pheatmap(pm[rows, , drop = FALSE], color = rev(brewer.pal(8, "RdYlBu")),
                  clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan",
                  clustering_method = "complete", show_rownames = FALSE, show_colnames = FALSE,
                  annotation_col = ann, annotation_colors = ann_cols, annotation_legend = legend, legend = legend,
                  treeheight_row = 25, treeheight_col = 20, main = title, fontsize = 9, silent = TRUE)
    wrap_elements(full = h$gtable)
  }
  pE <- heat_panel(dml$ploidy$united_row, sprintf("Ploidy DML (n = %s)", format(nrow(dml$ploidy), big.mark = ",")), FALSE)
  pF <- heat_panel(dml$heat$united_row, sprintf("Heat DML (n = %s)", format(nrow(dml$heat), big.mark = ",")), FALSE)
  pG <- heat_panel(dml$shared$united_row, sprintf("In both lists (n = %d)", nrow(dml$shared)), TRUE)

  fig2 <- ((pA + pB + pC + pD) + plot_layout(nrow = 1)) / ((pE + pF + pG) + plot_layout(nrow = 1, widths = c(1, 1, 1.25))) +
    plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = "A")
  save_figure(fig2, output_path(step, "figures", sprintf("figure_2_%s", mode)), width = 16, height = 10, dpi = 300)

  # ---- Figure 4: the two differences at each locus in both lists ------------------
  s <- copy(dml$shared)
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
  save_figure(fig4, output_path(step, "figures", sprintf("figure_4_%s", mode)), width = 7.5, height = 8)
  rm(pm); gc(verbose = FALSE)
}

finish_step(step, t0)

# specification_curve.R
#
# A specification curve (Simonsohn, Simmons and Nelson 2020, Nature Human
# Behaviour 4:1208, doi:10.1038/s41562-020-01004-2) of every analysis in the
# comparison table: the 192 analyses (the 64 combinations of the six settings
# corrected, cov5, mnchisq, ctexcl, d50 and q05, each with the paper's contrasts,
# with ctrl and with inter), sorted by the length of the ploidy list. The upper
# panels give the four list lengths of each analysis on a log scale, each panel
# with its own axis (ploidy DML; heat DML; loci in both lists; interaction DML,
# inter only); the lower panel marks the settings of each analysis. With `inter`
# the ploidy and heat lists are those of the paper's contrasts, so those points
# repeat their neighbours. Not adopted: a description of the optional analyses.
#
# Run:     Rscript 01_pipeline/optional/specification_curve.R   (from the repository root,
#          after compare_analyses.R)
# Reads:   04_output/not_adopted/comparisons/analysis_comparison.csv
# Writes:  04_output/not_adopted/comparisons/specification_curve.png
#          04_output/not_adopted/comparisons/specification_curve.pdf

suppressPackageStartupMessages({ library(data.table); library(here); library(ggplot2); library(patchwork) })

cmp_dir <- here("04_output", "not_adopted", "comparisons")
x <- fread(file.path(cmp_dir, "analysis_comparison.csv"))
if (nrow(x) != 192) stop("Expected the 192 analyses of the full factorial; found ", nrow(x), ". Run compare_analyses.R after every run.", call. = FALSE)

tokens <- strsplit(x$analysis, "-", fixed = TRUE)
has <- function(tok) vapply(tokens, function(t) tok %in% t, logical(1))
x[, contrast := fifelse(has("ctrl"), "Control-only ploidy (ctrl)", fifelse(has("inter"), "Paper's + tested interaction (inter)", "Paper's contrasts"))]
x[, contrast := factor(contrast, levels = c("Paper's contrasts", "Control-only ploidy (ctrl)", "Paper's + tested interaction (inter)"))]

# Order: ploidy DML, then heat DML, then shared loci, then the contrast.
setorder(x, ploidy_dml, heat_dml, shared, contrast)
x[, rank := seq_len(.N)]
paper_rank <- x[analysis == "legacy", rank]
v2022_rank <- x[analysis == "corrected-cov5-mnchisq-ctexcl-d50", rank]

# Colours: the reference categorical palette, first three slots (validated all-pairs, light
# surface); shapes repeat the contrast so identity is not colour alone.
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"; off_col <- "#d9d8d4"
pal <- c("Paper's contrasts" = "#2a78d6", "Control-only ploidy (ctrl)" = "#eb6834", "Paper's + tested interaction (inter)" = "#1baf7a")
shp <- c("Paper's contrasts" = 16, "Control-only ploidy (ctrl)" = 17, "Paper's + tested interaction (inter)" = 15)

theme_spec <- function() {
  theme_minimal(base_size = 10) +
    theme(plot.background = element_rect(fill = surface, colour = NA), panel.background = element_rect(fill = surface, colour = NA),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
          axis.title.y = element_text(colour = ink, size = 9, angle = 0, vjust = 0.5, hjust = 1), legend.text = element_text(colour = ink),
          axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank())
}
marks <- list(
  annotate("rect", xmin = paper_rank - 0.5, xmax = paper_rank + 0.5, ymin = -Inf, ymax = Inf, fill = "#2a78d6", alpha = 0.10),
  annotate("rect", xmin = v2022_rank - 0.5, xmax = v2022_rank + 0.5, ymin = -Inf, ymax = Inf, fill = "#52514e", alpha = 0.10))
log_y <- scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.25),
                            breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
                            labels = c("0", "1", "10", "100", "1,000", "10,000", "100,000"))
x_scale <- scale_x_continuous(limits = c(0.5, nrow(x) + 0.5), expand = c(0, 0))

curve <- function(col, title, data = x) {
  ggplot(data[!is.na(get(col))], aes(rank, .data[[col]], colour = contrast, shape = contrast)) +
    marks + geom_point(size = 1.6) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) + scale_shape_manual(values = shp, name = NULL, drop = FALSE) +
    log_y + x_scale + labs(y = title) + theme_spec()
}
p_ploidy <- curve("ploidy_dml", "Ploidy DML\n(ctrl: among the\n12 controls)") +
  annotate("text", x = paper_rank, y = x[rank == paper_rank, ploidy_dml] * 4, label = "paper", colour = ink, size = 3, vjust = 0) +
  annotate("text", x = v2022_rank, y = 40, label = "Venkataraman\net al. 2022", colour = ink, size = 3, vjust = 0, lineheight = 0.9)
p_heat  <- curve("heat_dml", "Heat DML")
p_share <- curve("shared", "Loci in both lists\n(the draft's\n\"interaction-DML\")")
p_inter <- curve("interaction_dml", "Tested\ninteraction DML\n(inter only)")

settings <- data.table(
  label = c("5x coverage (else 10x)", "Upper coverage filter", "Overdispersion correction (MN, chi-squared)",
            "C-to-T variant loci excluded", "Threshold 50 points (else 20)", "q < 0.05 (else 0.01)"),
  token = c("cov5", "corrected", "mnchisq", "ctexcl", "d50", "q05"))
# has() reads `tokens`, which follows the row order before sorting; recompute it on the sorted table.
tokens <- strsplit(x$analysis, "-", fixed = TRUE)
m <- rbindlist(lapply(seq_len(nrow(settings)), function(i) data.table(rank = x$rank, setting = settings$label[i], on = has(settings$token[i]), contrast = x$contrast)))
m[, setting := factor(setting, levels = rev(settings$label))]
p_spec <- ggplot(m, aes(rank, setting)) + marks +
  geom_point(data = m[on == FALSE], colour = off_col, size = 0.5, shape = 16) +
  geom_point(data = m[on == TRUE], aes(colour = contrast, shape = contrast), size = 1.3) +
  scale_colour_manual(values = pal, guide = "none") + scale_shape_manual(values = shp, guide = "none") +
  x_scale + labs(y = NULL, x = "The 192 analyses, sorted by the length of the ploidy list (then heat list, shared loci)") + theme_spec() +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = ink, size = 8.5),
        axis.title.x = element_text(colour = ink2, size = 9, margin = margin(t = 6)))

fig <- (p_ploidy / p_heat / p_share / p_inter / p_spec) +
  plot_layout(heights = c(1, 1, 1, 1, 1.15), guides = "collect") +
  plot_annotation(title = "Specification curve of the 192 DML analyses (none adopted except the paper's)",
                  subtitle = paste("Each column is one analysis; counts on a log scale (0 shown at the bottom). Shaded columns: the paper (legacy)",
                                   "and the settings of Venkataraman et al. 2022 (corrected-cov5-mnchisq-ctexcl-d50).",
                                   "\nWith inter, the ploidy and heat lists equal those of the paper's contrasts. Counts: analysis_comparison.csv."),
                  theme = theme(plot.background = element_rect(fill = surface, colour = NA),
                                plot.title = element_text(colour = ink, face = "bold", size = 12),
                                plot.subtitle = element_text(colour = ink2, size = 9))) &
  theme(legend.position = "top", legend.background = element_rect(fill = surface, colour = NA))

ggsave(file.path(cmp_dir, "specification_curve.png"), fig, width = 13, height = 11, dpi = 200, bg = surface, device = ragg::agg_png)
ggsave(file.path(cmp_dir, "specification_curve.pdf"), fig, width = 13, height = 11, bg = surface)
cat("written: specification_curve.png and .pdf in", cmp_dir, "\n")

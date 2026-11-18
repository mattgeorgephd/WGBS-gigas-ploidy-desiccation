# plot_theme.R
# One ggplot2 theme for every figure.

theme_wgbs <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.text        = ggplot2::element_text(colour = "black"),
      axis.ticks       = ggplot2::element_line(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.key       = ggplot2::element_blank(),
      plot.tag         = ggplot2::element_text(face = "bold")
    )
}

# Save a plot as PNG and PDF with the same base name.
save_figure <- function(plot, path_no_ext, width, height, dpi = 400) {
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot, width = width, height = height, bg = "white")
  invisible(path_no_ext)
}

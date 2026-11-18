# 03_global_methylation.R
#
# Purpose   Global DNA methylation by two assays, and Figure 1.
#           (a) ELISA: percent 5-mC per animal computed from the raw plate
#               absorbance, the plate map, and the standard curve (N = 32).
#           (b) WGBS: percent CpG methylation per library from the Bismark
#               summary of step 02 (N = 22).
#           Both are modelled as ploidy x treatment on three scales (raw,
#           arcsine square root, log) with model diagnostics, so the dependence
#           of the conclusions on the scale is on record (issue I-10).
#
# Inputs    03_input_files/elisa/raw/20190206_{First,2nd,3rd}ReadingSTRIGG.xls
#           03_input_files/elisa/raw/20190131_RonitsSamplesforGlobalDNAMeth.xlsx
#           04_output/02_alignment_qc/tables/bismark_by_library.csv
# Outputs   04_output/03_global_methylation/tables/
#             elisa_wells.csv, elisa_standard_curve.csv, elisa_percent_5mC.csv,
#             pct_mCpG_by_library.csv, group_summary.csv,
#             anova_by_scale.csv, tukey_by_scale.csv
#           04_output/03_global_methylation/figures/figure_1.(png|pdf)
#
# Manuscript items: Results 3.1 (5-mC and percent CpG statistics); Figure 1.
# Notes     Percent 5-mC uses the kit formula with a fixed DNA input
#           (config$elisa_input_ng); see issue I-12 for why that needs
#           confirming. The WGBS comparison is between two cohorts sequenced
#           separately (issue I-01).

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "readxl", "ggplot2", "patchwork", "agricolae", "car"))

step <- "03_global_methylation"
t0 <- start_step(step, "ELISA 5-mC, percent CpG methylation, Figure 1")

# =============================================================================
# (a) ELISA: absorbance to percent 5-mC
# =============================================================================
read_plate <- function(file) {
  x <- as.data.table(read_excel(input_path("elisa", "raw", file), sheet = "plate list"))
  data.table(well = x[["Well"]], absorbance = x[[grep("^Absorbance", names(x))]])
}
reads <- lapply(c(first = "20190206_FirstReadingSTRIGG.xls", second = "20190206_2ndReadingSTRIGG.xls",
                  third = "20190206_3rdReadingSTRIGG.xls"), read_plate)
plate_map <- as.data.table(read_excel(input_path("elisa", "raw", "20190131_RonitsSamplesforGlobalDNAMeth.xlsx"),
                                      sheet = "PlateMap_list"))
setnames(plate_map, c("well", "type", "label"))

wells <- merge(plate_map, reads$first, by = "well")
wells <- merge(wells, merge(reads$second, reads$third, by = "well", suffixes = c("_second", "_third")), by = "well")
wells[, reading_sd := apply(.SD, 1, sd), .SDcols = c("absorbance", "absorbance_second", "absorbance_third")]
excluded_labels <- c("no_sample", "contam", "no_ab")       # flags typed into the plate map
wells[, used := !(label %in% excluded_labels)]
stopifnot(nrow(wells) == 96)
log_msg("plate: %d wells, %d used; largest SD among the three readings of a well: %.4f",
        nrow(wells), sum(wells$used), max(wells$reading_sd))
write_table(wells[order(well)], output_path(step, "tables", "elisa_wells.csv"))

# Mean absorbance per label, from the first reading (the three readings agree).
per_label <- wells[used == TRUE, .(mean_abs = mean(absorbance), sd_abs = sd(absorbance), n_wells = .N),
                   by = .(type, label)]

# Standards: nominal percent 5-mC is read from the label. "NC" is the negative
# control (0%). The undiluted positive control ("5 ul PC") carries no percentage
# in its label and is not used in the fit.
std <- per_label[type == "STD"]
std[, percent_5mC := fifelse(grepl("NC", label), 0,
                             suppressWarnings(as.numeric(sub(".*?([0-9.]+)\\s*%.*", "\\1", label))))]
std[!grepl("%|NC", label), percent_5mC := NA_real_]
std[, in_fit := !is.na(percent_5mC) & percent_5mC <= config$elisa_curve_max_percent]
setorder(std, percent_5mC, na.last = TRUE)
curve <- lm(mean_abs ~ percent_5mC, data = std[in_fit == TRUE])
slope <- unname(coef(curve)[2]); nc_abs <- std[percent_5mC == 0, mean_abs]
std[, `:=`(slope = slope, intercept = unname(coef(curve)[1]), r_squared = summary(curve)$r.squared)]
write_table(std, output_path(step, "tables", "elisa_standard_curve.csv"))
log_msg("standard curve on %d standards (0 to %g%%): slope %.7f, R2 %.7f, negative control OD %.7f",
        sum(std$in_fit), config$elisa_curve_max_percent, slope, summary(curve)$r.squared, nc_abs)

# Kit formula: 5-mC% = (sample OD - negative control OD) / (slope x ng DNA) x 100
elisa <- per_label[type == "Sample"]
elisa[, percent_5mC := (mean_abs - nc_abs) / (slope * config$elisa_input_ng) * 100]
elisa[, animal_id := sub("-C$", "", label)]
elisa[, ploidy := factor(fifelse(substr(animal_id, 1, 1) == "D", "diploid", "triploid"), levels = c("diploid", "triploid"))]
elisa[, animal_number := as.integer(substr(animal_id, 2, 3))]
elisa[, treatment := factor(fifelse(animal_number >= 11, "heat", "control"), levels = c("control", "heat"))]
elisa[, above_top_standard := mean_abs > max(std[in_fit == TRUE, mean_abs])]   # value is extrapolated
stopifnot(nrow(elisa) == 32, all(elisa[, .N, by = .(ploidy, treatment)]$N == 8),
          all(elisa$animal_number %in% c(1:8, 11:18)))
setorder(elisa, ploidy, treatment, animal_number)
write_table(elisa[, .(animal_id, ploidy, treatment, n_wells, mean_abs, sd_abs, percent_5mC, above_top_standard)],
            output_path(step, "tables", "elisa_percent_5mC.csv"))
log_msg("%d animals; %d with a single well; %d above the top fitted standard (extrapolated)",
        nrow(elisa), sum(elisa$n_wells == 1), sum(elisa$above_top_standard))

# =============================================================================
# (b) WGBS: percent CpG methylation per library
# =============================================================================
wgbs <- fread(output_path("02_alignment_qc", "tables", "bismark_by_library.csv"))
wgbs[, ploidy := factor(ploidy, levels = c("diploid", "triploid"))]
wgbs[, treatment := factor(treatment, levels = c("control", "heat"))]
write_table(wgbs[, .(library_id, cohort, ploidy, treatment, acute_heat_shock, pct_mCpG, pct_mCpG_predup)],
            output_path(step, "tables", "pct_mCpG_by_library.csv"))

# =============================================================================
# Models: ploidy x treatment on three scales, with diagnostics
# =============================================================================
scales <- list(raw = identity, arcsine = function(p) asin(sqrt(p / 100)), log = log)

fit_all_scales <- function(d, yvar, assay) {
  d <- as.data.frame(d); d$group <- interaction(d$treatment, d$ploidy)
  anova_rows <- list(); tukey_rows <- list(); letters_rows <- list()
  for (s in names(scales)) {
    d$response <- scales[[s]](d[[yvar]])
    fit  <- aov(response ~ ploidy * treatment, data = d)
    seq  <- as.data.frame(summary(fit)[[1]]); seq$term <- trimws(rownames(seq))
    t2   <- as.data.frame(car::Anova(fit, type = 2))
    anova_rows[[s]] <- data.frame(
      assay = assay, scale = s, term = seq$term, df = seq$Df, sum_sq = seq$`Sum Sq`, mean_sq = seq$`Mean Sq`,
      F_sequential = seq$`F value`, p_sequential = seq$`Pr(>F)`,
      F_type2 = t2[match(seq$term, rownames(t2)), "F value"], p_type2 = t2[match(seq$term, rownames(t2)), "Pr(>F)"],
      residual_shapiro_p = shapiro.test(residuals(fit))$p.value,
      brown_forsythe_p   = car::leveneTest(response ~ group, data = d, center = median)[1, "Pr(>F)"])
    one_way <- aov(response ~ group, data = d)
    tk <- as.data.frame(TukeyHSD(one_way)$group)
    tukey_rows[[s]] <- data.frame(assay = assay, scale = s, comparison = rownames(tk), difference = tk$diff,
                                  lower = tk$lwr, upper = tk$upr, p_adjusted = tk$`p adj`)
    hsd <- agricolae::HSD.test(one_way, "group", group = TRUE)$groups
    letters_rows[[s]] <- data.frame(assay = assay, scale = s, group = rownames(hsd),
                                    mean_on_scale = hsd$response, letter = hsd$groups)
  }
  list(anova = rbindlist(anova_rows), tukey = rbindlist(tukey_rows), letters = rbindlist(letters_rows))
}
res_elisa <- fit_all_scales(elisa, "percent_5mC", "ELISA percent 5-mC")
res_wgbs  <- fit_all_scales(wgbs, "pct_mCpG", "WGBS percent mCpG (deduplicated)")
res_wgbs_predup <- fit_all_scales(wgbs, "pct_mCpG_predup", "WGBS percent mCpG (before deduplication)")

anova_all <- rbind(res_elisa$anova, res_wgbs$anova, res_wgbs_predup$anova)
tukey_all <- rbind(res_elisa$tukey, res_wgbs$tukey, res_wgbs_predup$tukey)
letters_all <- rbind(res_elisa$letters, res_wgbs$letters, res_wgbs_predup$letters)
anova_all[, primary_scale := scale == config$global_methylation_scale]
write_table(anova_all, output_path(step, "tables", "anova_by_scale.csv"))
write_table(tukey_all, output_path(step, "tables", "tukey_by_scale.csv"))

summarise_groups <- function(d, yvar, assay) {
  as.data.table(d)[, .(assay = assay, n = .N, mean = mean(get(yvar)), sd = sd(get(yvar)), median = median(get(yvar)),
                       min = min(get(yvar)), max = max(get(yvar))), by = .(ploidy, treatment)]
}
group_summary <- rbind(summarise_groups(elisa, "percent_5mC", "ELISA percent 5-mC"),
                       summarise_groups(wgbs, "pct_mCpG", "WGBS percent mCpG (deduplicated)"))
group_summary <- merge(group_summary,
                       letters_all[scale == config$global_methylation_scale,
                                   .(assay, ploidy = sub(".*\\.", "", group), treatment = sub("\\..*", "", group),
                                     tukey_letter = letter)],
                       by = c("assay", "ploidy", "treatment"))
write_table(group_summary, output_path(step, "tables", "group_summary.csv"))

show <- anova_all[term != "Residuals" & !grepl("before", assay),
                  .(assay = sub(" percent.*", "", assay), scale, term, df, F = round(F_sequential, 2),
                    p = signif(p_sequential, 3), resid_normal_p = round(residual_shapiro_p, 2),
                    equal_var_p = round(brown_forsythe_p, 2))]
print(show)

# =============================================================================
# Figure 1
# =============================================================================
panel <- function(d, yvar, ylab, assay_name) {
  d <- as.data.table(d)
  lab <- group_summary[assay == assay_name]     # Tukey letters on the primary scale
  stopifnot(nrow(lab) == 4)
  lab[, y := max + 0.06 * diff(range(d[[yvar]]))]
  ggplot(d, aes(x = treatment, y = .data[[yvar]], fill = ploidy)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, position = position_dodge(width = 0.8), linewidth = 0.5) +
    geom_point(aes(group = ploidy), shape = 21, size = 1.8, alpha = 0.8,
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 1)) +
    geom_text(data = lab, aes(x = treatment, y = y, label = tukey_letter, group = ploidy),
              position = position_dodge(width = 0.8), size = 4.5, inherit.aes = FALSE) +
    scale_fill_manual(values = config$colours$ploidy, name = NULL) +
    labs(x = NULL, y = ylab) + theme_wgbs()
}
fig1 <- panel(elisa, "percent_5mC", "5-mC (% of cytosines, ELISA)", "ELISA percent 5-mC") +
  panel(wgbs, "pct_mCpG", "Methylated CpG (%, WGBS)", "WGBS percent mCpG (deduplicated)") +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
save_figure(fig1, output_path(step, "figures", "figure_1"), width = 8.5, height = 4.6)

finish_step(step, t0)

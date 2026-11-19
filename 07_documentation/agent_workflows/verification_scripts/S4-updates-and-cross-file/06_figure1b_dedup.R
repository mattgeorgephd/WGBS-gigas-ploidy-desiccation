# S4 check: the interaction contrast of Figure 1B on deduplicated counts (library_level.csv, used by I-49)
# and on counts before deduplication (what the draft's Figure 1B plots, I-13; used by I-01).
d <- read.csv("/home/user/WGBS-gigas-ploidy-desiccation/04_output/03_global_methylation/tables/pct_mCpG_by_library.csv")
high <- c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15")
f <- function(x, sub) {
  dd <- d[sub, ]; y <- dd[[x]]
  m <- tapply(y, list(dd$ploidy, dd$treatment), mean)
  con <- (m["triploid", "heat"] - m["diploid", "heat"]) - (m["triploid", "control"] - m["diploid", "control"])
  a <- anova(lm(y ~ ploidy * treatment, dd)); as <- anova(lm(asin(sqrt(y / 100)) ~ ploidy * treatment, dd))
  cat(sprintf("%-16s n = %2d : heat D - T %+.3f ; control D - T %+.3f ; interaction contrast %+.3f ; F(1,%d) = %.2f, p = %.4f (raw scale) ; F = %.2f, p = %.4f (arcsine)\n",
              x, nrow(dd), m["diploid", "heat"] - m["triploid", "heat"], m["diploid", "control"] - m["triploid", "control"], con,
              a["Residuals", "Df"], a["ploidy:treatment", "F value"], a["ploidy:treatment", "Pr(>F)"], as["ploidy:treatment", "F value"], as["ploidy:treatment", "Pr(>F)"]))
}
for (x in c("pct_mCpG", "pct_mCpG_predup")) { f(x, rep(TRUE, nrow(d))); f(x, !d$library_id %in% high) }

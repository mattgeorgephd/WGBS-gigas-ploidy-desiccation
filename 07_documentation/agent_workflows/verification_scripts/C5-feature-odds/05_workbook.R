suppressPackageStartupMessages(library(readxl))
f <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive/analyses/DML_locations-all-datasets.xlsx"
sh <- excel_sheets(f); print(sh)
options(width = 250)
for (s in sh) {
  cat("\n================ sheet:", s, "\n")
  d <- suppressMessages(read_excel(f, sheet = s, col_names = FALSE, .name_repair = "minimal"))
  cat("dim:", dim(d), "\n")
  print(as.data.frame(d[seq_len(min(nrow(d), 45)), seq_len(min(ncol(d), 22))]), digits = 6)
}

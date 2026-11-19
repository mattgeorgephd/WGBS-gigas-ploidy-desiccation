# Task (b): the help page of calculateDiffMeth and the vignette, as installed (methylKit 1.28.0).
suppressPackageStartupMessages(library(methylKit))
cat("methylKit", as.character(packageVersion("methylKit")), "\n")
rd <- tools::Rd_db("methylKit")
nm <- grep("calculateDiffMeth", names(rd), value = TRUE)
cat("Rd files:", nm, "\n")
for (n in nm) {
  cat("\n================ ", n, " ================\n")
  tools::Rd2txt(rd[[n]], options = list(underline_titles = FALSE, width = 2000))
}
cat("\n================ vignette files ================\n")
print(list.files(system.file("doc", package = "methylKit"), full.names = TRUE))
cat("\n================ p.adjusted ================\n")
print(methylKit:::p.adjusted)

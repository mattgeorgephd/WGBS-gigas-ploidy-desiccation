# Print the help page of calculateDiffMeth of the installed methylKit and the
# vignette passages on overdispersion. Read-only.
suppressMessages(library(tools))
cat("methylKit version:", as.character(packageVersion("methylKit")), "\n")
db <- Rd_db("methylKit")
rd <- db[[grep("calculateDiffMeth", names(db))[1]]]
txt <- capture.output(Rd2txt(rd, options = list(underline_titles = FALSE)))
cat(txt, sep = "\n")
vdir <- system.file("doc", package = "methylKit")
cat("\n==== vignette files:", list.files(vdir), "\n")

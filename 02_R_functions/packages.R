# packages.R
# Package loading for pipeline steps.

# Load the packages a step needs. Stops, naming every missing package, before any
# work starts; it never installs anything.
# A step that names data.table also needs R.utils: data.table::fread() reads the
# .gz tables through it, and installing data.table does not install it (it is
# only suggested). R.utils is looked up, not loaded, so a step that reads no .gz
# table runs as before.
load_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if ("data.table" %in% pkgs && !nzchar(system.file(package = "R.utils"))) {
    missing <- c(missing, "R.utils (data.table::fread() reads .gz tables with it)")
  }
  if (length(missing) > 0) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "),
         ". Install them (see README, 'Running the pipeline') and rerun.", call. = FALSE)
  }
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(pkgs)
}

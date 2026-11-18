# packages.R
# Package loading for pipeline steps.

# Load the packages a step needs. Stops, naming every missing package, before any
# work starts; it never installs anything.
load_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "),
         ". Install them (see README, 'Running the pipeline') and rerun.", call. = FALSE)
  }
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(pkgs)
}

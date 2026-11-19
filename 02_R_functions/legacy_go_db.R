# legacy_go_db.R
# The Gene Ontology release the legacy enrichment ran on.
#
# topGO takes the ontology (which terms exist, and which term is a child of
# which) from the installed GO.db package, so the same input gives different
# tables under different GO.db releases. The archived legacy tables are
# reproduced exactly, term for term, with GO.db 3.13.0 (Bioconductor 3.13, GO
# release of 2021-02-01) and with no later release. Step 09 therefore loads that
# release from a side library, without touching the GO.db every other step uses.
# The release, its address, and its md5 are set in run_config.R.

# Folder of the side library (inside the download cache).
legacy_go_db_lib <- function(config) gannet_cache_dir(config, "legacy_r_lib")

legacy_go_db_installed <- function(config) {
  d <- file.path(legacy_go_db_lib(config), "GO.db", "DESCRIPTION")
  file.exists(d) && identical(unname(read.dcf(d, fields = "Version")[1, 1]), config$legacy_go_db$version)
}

# Download the source package (verifying its md5) and install it into the side
# library. Does nothing when it is already there.
install_legacy_go_db <- function(config) {
  if (legacy_go_db_installed(config)) return(invisible(legacy_go_db_lib(config)))
  g <- config$legacy_go_db
  tarball <- file.path(gannet_cache_dir(config, "legacy_r_packages"), basename(g$url))
  if (!file.exists(tarball) || !identical(unname(tools::md5sum(tarball)), g$md5)) {
    log_msg("downloading %s", g$url)
    old <- options(timeout = config$download_timeout_s); on.exit(options(old), add = TRUE)
    status <- utils::download.file(g$url, tarball, mode = "wb", quiet = TRUE, method = "libcurl")
    if (status != 0 || !identical(unname(tools::md5sum(tarball)), g$md5)) {
      unlink(tarball)
      stop("Download failed or md5 mismatch: ", g$url, call. = FALSE)
    }
  }
  log_msg("installing GO.db %s into %s", g$version, legacy_go_db_lib(config))
  utils::install.packages(tarball, lib = legacy_go_db_lib(config), repos = NULL, type = "source",
                          INSTALL_opts = "--no-test-load", quiet = TRUE)
  if (!legacy_go_db_installed(config)) stop("GO.db ", g$version, " could not be installed", call. = FALSE)
  invisible(legacy_go_db_lib(config))
}

# Put the side library first on the search path, so that GO.db is loaded from it.
# Call before topGO or any other package that imports GO.db is loaded.
use_legacy_go_db <- function(config) {
  if ("GO.db" %in% loadedNamespaces()) {
    have <- as.character(utils::packageVersion("GO.db"))
    if (have != config$legacy_go_db$version) {
      stop("GO.db ", have, " is already loaded in this R session; the legacy enrichment needs ",
           config$legacy_go_db$version, ". Run the step in a fresh R process.", call. = FALSE)
    }
    return(invisible(have))
  }
  install_legacy_go_db(config)
  .libPaths(c(legacy_go_db_lib(config), .libPaths()))
  loadNamespace("GO.db")
  have <- as.character(utils::packageVersion("GO.db"))
  stopifnot(have == config$legacy_go_db$version)
  invisible(have)
}

# Release date of the ontology in the loaded GO.db.
go_source_date <- function() {
  info <- GO.db::GO_dbInfo()
  info$value[info$name == "GOSOURCEDATE"]
}

# legacy_go_db.R
# The Gene Ontology release the legacy enrichment ran on.
#
# topGO takes the ontology (which terms exist, and which term is a child of
# which) from the installed GO.db package, so the same input gives different
# tables under different GO.db releases. The nine archived tables behind the
# draft (dated 2022-03-07) are reproduced exactly, term for term, with GO.db
# 3.13.0 (GO release of 2021-02-01; the GO.db of Bioconductor 3.13) and with
# none of the eight other releases tested, 3.7.0 to 3.18.0 (issues.md, I-37).
# Other tables of the archive match GO.db 3.12.1 and 3.16.0. Step 09 therefore
# loads release 3.13.0 from a side library, without touching the GO.db every
# other step uses. The release, its address, and its md5 are set in run_config.R.

# Folder of the side library (inside the download cache).
legacy_go_db_lib <- function(config) gannet_cache_dir(config, "legacy_r_lib")

legacy_go_db_installed <- function(config) {
  d <- file.path(legacy_go_db_lib(config), "GO.db", "DESCRIPTION")
  file.exists(d) && identical(unname(read.dcf(d, fields = "Version")[1, 1]), config$legacy_go_db$version)
}

# Download the source package (verifying its md5) and install it into the side
# library. Does nothing when it is already there.
# An install that was killed (a lost session, a restart of the machine) leaves
# the lock folder 00LOCK-GO.db in the side library, and R then refuses every
# later install without saying why, because the installer runs quietly. The
# install takes under a minute, so a lock older than ten minutes is taken for a
# leftover and removed; a younger one may belong to an install that is running
# in another process, and the function stops and names it.
install_legacy_go_db <- function(config) {
  lib <- legacy_go_db_lib(config)
  if (legacy_go_db_installed(config)) return(invisible(lib))
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
  lock <- file.path(lib, "00LOCK-GO.db")
  if (dir.exists(lock) && difftime(Sys.time(), file.mtime(lock), units = "mins") > 10) {
    log_msg("removing the lock left by an interrupted install: %s", lock)
    unlink(lock, recursive = TRUE)
  }
  log_msg("installing GO.db %s into %s", g$version, lib)
  utils::install.packages(tarball, lib = lib, repos = NULL, type = "source",
                          INSTALL_opts = "--no-test-load", quiet = TRUE)
  if (!legacy_go_db_installed(config)) {
    stop("GO.db ", g$version, " could not be installed into ", lib, ". ",
         if (dir.exists(lock)) {
           paste0("The lock folder ", lock, " exists: another install is running, or one was interrupted less than ",
                  "ten minutes ago. If no other step is running, delete that folder and rerun.")
         } else {
           paste0("To see the messages of the installer, run: R CMD INSTALL --no-test-load -l ", shQuote(lib), " ", shQuote(tarball))
         }, call. = FALSE)
  }
  invisible(lib)
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

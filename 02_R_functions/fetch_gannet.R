# fetch_gannet.R
# Download and verify the inputs that are too large for GitHub.

read_gannet_manifest <- function() {
  utils::read.csv(input_path("gannet_manifest.csv"), stringsAsFactors = FALSE)
}

# Return the local path of one manifest file, downloading it into the cache if it
# is absent or has the wrong size. With verify_md5 = TRUE the md5 is checked
# against the manifest whether or not a download happened.
fetch_gannet_file <- function(file, config, manifest = read_gannet_manifest(), verify_md5 = FALSE) {
  row <- manifest[manifest$file == file, ]
  if (nrow(row) != 1) stop("Not in gannet_manifest.csv (or listed twice): ", file, call. = FALSE)
  dest <- file.path(gannet_cache_dir(config, row$cache_subdir), row$file)

  have <- file.exists(dest) && file.size(dest) == row$bytes
  if (!have) {
    log_msg("downloading %s (%.0f MB)", row$file, row$bytes / 1024^2)
    part <- paste0(dest, ".part")
    old <- options(timeout = config$download_timeout_s); on.exit(options(old), add = TRUE)
    status <- utils::download.file(row$url, part, mode = "wb", quiet = TRUE, method = "libcurl")
    if (status != 0 || file.size(part) != row$bytes) {
      unlink(part)
      stop("Download failed or incomplete: ", row$url, call. = FALSE)
    }
    file.rename(part, dest)
    verify_md5 <- TRUE
  }
  if (verify_md5) {
    md5 <- unname(tools::md5sum(dest))
    if (!identical(md5, row$md5)) {
      stop("md5 mismatch for ", dest, ": expected ", row$md5, ", got ", md5,
           ". Delete the file and rerun.", call. = FALSE)
    }
  }
  dest
}

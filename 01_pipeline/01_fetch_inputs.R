# 01_fetch_inputs.R
#
# Purpose   Make sure every input the pipeline needs is present and intact.
#           Inputs small enough for GitHub are already in 03_input_files/. The
#           larger ones are listed in 03_input_files/gannet_manifest.csv; this
#           step downloads those into the cache and verifies size and md5. It
#           also installs, into a side library in the cache, the GO.db release
#           that step 09 needs to reproduce the legacy enrichment.
#
# Run       Rscript 01_pipeline/01_fetch_inputs.R
#           (first; under a minute with a full cache, otherwise the time of a
#           10.5 GB download)
#
# Inputs    03_input_files/gannet_manifest.csv
#           03_input_files/metadata/wgbs_libraries.csv
#           every other file under 03_input_files/ (size and md5 only)
#           https://gannet.fish.washington.edu/panopea/ (coverage files, reference genome)
#           Bioconductor archive (GO.db; address and md5 in run_config.R)
# Outputs   04_output/01_fetch_inputs/tables/input_inventory.csv
#           cache: 03_input_files/gannet_cache/ (or $WGBS_CACHE_DIR)
#           04_output/01_fetch_inputs/session_info.txt
#
# Manuscript items: none.
# Notes     - The first run downloads 10.5 GB from gannet and the 31 MB GO.db
#             archive. Sizes are in decimal units (1 GB = 1e9 bytes), as in the
#             README files.
#           - input_inventory.csv lists every file under 03_input_files/ with its
#             size and md5, README files included, so it changes whenever one of
#             them is edited: rerun this step and commit the table with the
#             edit. Its rows are in byte order of the path, whatever the
#             collation of the locale.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table"))

step <- "01_fetch_inputs"
t0 <- start_step(step, "fetch and verify inputs")

# ---- Inputs hosted on gannet -------------------------------------------------
manifest <- read_gannet_manifest()
lib <- read_libraries()
missing_cov <- setdiff(lib$cov_file, manifest$file)
if (length(missing_cov) > 0) {
  stop("Coverage files in the sample sheet but not in the manifest: ", paste(missing_cov, collapse = ", "))
}
log_msg("%d files in the gannet manifest, %.2f GB; cache: %s", nrow(manifest),
        sum(manifest$bytes) / 1e9, gannet_cache_dir(config))

local <- vapply(manifest$file, fetch_gannet_file, character(1), config = config,
                manifest = manifest, verify_md5 = TRUE)
gannet <- data.frame(location = "gannet (cached)", file = manifest$file, bytes = manifest$bytes,
                     md5 = manifest$md5, status = "present, md5 verified")
log_msg("all %d gannet files present and verified", nrow(manifest))

# ---- Inputs kept in the repository -------------------------------------------
repo_files <- list.files(input_path(), recursive = TRUE, full.names = FALSE)
repo_files <- repo_files[!grepl("^gannet_cache/", repo_files)]
repo_files <- sort(repo_files, method = "radix")     # byte order; list.files() sorts by the collation of the locale
repo <- data.frame(location = "repository", file = repo_files,
                   bytes = file.size(input_path(repo_files)),
                   md5 = unname(tools::md5sum(input_path(repo_files))),
                   status = "present")

# The 13 genome feature tracks must be complete: later steps read all of them.
tracks <- c("gene.gff", "mRNA.gff", "exon.gff", "CDS.gff", "exonUTR.gff", "intron.bed", "nonCDS.bed",
            "lncRNA.gff", "upstream.gff", "downstream.gff", "flanks.gff", "intergenic.bed", "rm.te.bed")
track_files <- file.path("genome_features", paste0("cgigas_uk_roslin_v1_", tracks, ".gz"))
stopifnot(all(track_files %in% repo_files))

# ---- Ontology release of the legacy enrichment ---------------------------------
install_legacy_go_db(config)
go_tarball <- file.path(gannet_cache_dir(config, "legacy_r_packages"), basename(config$legacy_go_db$url))
bioc <- data.frame(location = "Bioconductor archive (cached)", file = basename(config$legacy_go_db$url),
                   bytes = if (file.exists(go_tarball)) file.size(go_tarball) else NA_real_,
                   md5 = config$legacy_go_db$md5,
                   status = sprintf("GO.db %s installed in the side library", config$legacy_go_db$version))
log_msg("GO.db %s available in %s", config$legacy_go_db$version, legacy_go_db_lib(config))

inventory <- rbind(repo, gannet, bioc)
write_table(inventory, output_path(step, "tables", "input_inventory.csv"))
log_msg("%d repository inputs (%.0f MB) and %d gannet inputs (%.2f GB) inventoried",
        nrow(repo), sum(repo$bytes) / 1e6, nrow(gannet), sum(gannet$bytes) / 1e9)

finish_step(step, t0)

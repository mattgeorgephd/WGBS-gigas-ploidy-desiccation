# build_wgbs_libraries.R
#
# Builds 03_input_files/metadata/wgbs_libraries.csv: one row per WGBS library, the
# sample sheet every pipeline step reads. Run from anywhere inside the repository:
#   Rscript 03_input_files/metadata/build_wgbs_libraries.R
#
# Inputs (03_input_files/metadata/raw/):
#   legacy_metadata.xlsx   sheet "list": the 22 libraries in the order the legacy
#                          methylKit analysis used (its sampleID 1 to 22), with
#                          ploidy (0 diploid, 1 triploid) and trt (0 control, 1 heat)
#   zr3534_wgbs_info.csv   sequencing sample sheet for the heat-treated libraries:
#                          animal ID, acute heat shock flag, SRA accessions
# Output:
#   wgbs_libraries.csv
#
# Do not edit wgbs_libraries.csv by hand; change the raw files or this script.

suppressPackageStartupMessages({
  library(readxl)
  library(here)
})

raw_dir <- here("03_input_files", "metadata", "raw")
legacy  <- as.data.frame(read_excel(file.path(raw_dir, "legacy_metadata.xlsx"), sheet = "list"))
zr3534  <- read.csv(file.path(raw_dir, "zr3534_wgbs_info.csv"), stringsAsFactors = FALSE)

stopifnot(nrow(legacy) == 22, all(legacy$sampleID == 1:22))

lib <- data.frame(
  library_id          = sub("_R1\\..*$", "", basename(legacy$path)),
  legacy_sample_index = as.integer(legacy$sampleID),
  stringsAsFactors    = FALSE
)
lib$zymo_project <- sub("_.*$", "", lib$library_id)
# Cohort: which experiment the animal came from. Legacy file names call these "ronit"
# (heat-treated animals) and "yaamini" (control animals).
lib$cohort       <- ifelse(legacy$dataset == "ronit", "heat_experiment", "control_experiment")
lib$legacy_label <- legacy$dataset
lib$tissue       <- legacy$tissue
lib$ploidy       <- ifelse(legacy$ploidy == 0, "diploid", "triploid")
lib$treatment    <- ifelse(legacy$trt == 1, "heat", "control")

# Animal identity, acute heat shock, and accessions exist only for zr3534.
m <- match(lib$library_id, zr3534$SeqID)
lib$animal_id        <- sub("-C$", "", zr3534$Library_Name[m])
lib$acute_heat_shock <- ifelse(is.na(m), "no", zr3534$Heat_Stress[m])
lib$library_kit      <- zr3534$Library_kit[m]
lib$sra_bioproject   <- zr3534$SRA_BioProject[m]
lib$sra_experiment   <- zr3534$SRA_accession[m]
lib$cov_file         <- paste0(lib$library_id, "_R1.CpG_report.merged_CpG_evidence.cov")

# Consistency checks between the two raw sources.
zr <- !is.na(m)
stopifnot(
  sum(zr) == 10,
  all(lib$treatment[zr] == "heat"),
  all(zr3534$Dessication[m[zr]] == "yes"),
  all(lib$ploidy[zr] == zr3534$Ploidy[m[zr]]),
  all(lib$treatment[!zr] == "control"),
  identical(as.vector(table(lib$treatment, lib$ploidy)), c(6L, 5L, 6L, 5L))
)

out <- here("03_input_files", "metadata", "wgbs_libraries.csv")
write.csv(lib, out, row.names = FALSE, na = "")
cat("Wrote", out, "\n")
print(table(treatment = lib$treatment, ploidy = lib$ploidy))
print(table(treatment = lib$treatment, acute_heat_shock = lib$acute_heat_shock))

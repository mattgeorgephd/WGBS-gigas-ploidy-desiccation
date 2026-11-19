suppressPackageStartupMessages(library(readxl))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
raw_dir <- file.path(repo, "03_input_files/metadata/raw")
print(excel_sheets(file.path(raw_dir, "legacy_metadata.xlsx")))
legacy <- as.data.frame(read_excel(file.path(raw_dir, "legacy_metadata.xlsx"), sheet = "list"))
options(width = 250)
print(legacy)
str(legacy)
zr3534 <- read.csv(file.path(raw_dir, "zr3534_wgbs_info.csv"), stringsAsFactors = FALSE)
lib <- data.frame(library_id = sub("_R1\\..*$", "", basename(legacy$path)), legacy_sample_index = as.integer(legacy$sampleID), stringsAsFactors = FALSE)
lib$zymo_project <- sub("_.*$", "", lib$library_id)
lib$cohort       <- ifelse(legacy$dataset == "ronit", "heat_experiment", "control_experiment")
lib$legacy_label <- legacy$dataset
lib$tissue       <- legacy$tissue
lib$ploidy       <- ifelse(legacy$ploidy == 0, "diploid", "triploid")
lib$treatment    <- ifelse(legacy$trt == 1, "heat", "control")
m <- match(lib$library_id, zr3534$SeqID)
lib$animal_id        <- sub("-C$", "", zr3534$Library_Name[m])
lib$acute_heat_shock <- ifelse(is.na(m), "no", zr3534$Heat_Stress[m])
lib$library_kit      <- zr3534$Library_kit[m]
lib$sra_bioproject   <- zr3534$SRA_BioProject[m]
lib$sra_experiment   <- zr3534$SRA_accession[m]
lib$cov_file         <- paste0(lib$library_id, "_R1.CpG_report.merged_CpG_evidence.cov")
tmp <- tempfile(fileext = ".csv"); write.csv(lib, tmp, row.names = FALSE, na = "")
cat("rebuilt sheet identical to committed wgbs_libraries.csv: ", identical(readLines(tmp), readLines(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))), "\n")
# other sheets
for (s in excel_sheets(file.path(raw_dir, "legacy_metadata.xlsx"))) { cat("\n--- sheet", s, "---\n"); print(as.data.frame(read_excel(file.path(raw_dir, "legacy_metadata.xlsx"), sheet = s)), row.names = FALSE) }

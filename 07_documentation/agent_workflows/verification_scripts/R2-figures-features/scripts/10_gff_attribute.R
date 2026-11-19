suppressPackageStartupMessages(library(data.table))
source("/home/user/WGBS-gigas-ploidy-desiccation/02_R_functions/genomic_overlaps.R")
dir <- "/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/genome_features"
strict <- function(a, key) {            # split on ';', take the field whose name is exactly `key`
  vapply(strsplit(a, ";", fixed = TRUE), function(f) { k <- sub("=.*", "", f); v <- sub("^[^=]*=", "", f); i <- match(key, k); if (is.na(i)) NA_character_ else v[i] }, "")
}
for (tr in c("gene", "mRNA")) {
  a <- fread(file.path(dir, sprintf("cgigas_uk_roslin_v1_%s.gff.gz", tr)), header = FALSE, sep = "\t", quote = "", select = 9)$V9
  for (key in c("Name", "gene_biotype", "gene", "product", "transcript_id")) {
    s <- strict(a, key); g <- gff_attribute(a, paste0(key, "="))
    cat(sprintf("%-5s %-14s present %6d; helper NA where present %d; differ %d; helper non-NA where absent %d\n", tr, key, sum(!is.na(s)), sum(is.na(g) & !is.na(s)), sum(!is.na(s) & !is.na(g) & s != g), sum(!is.na(g) & is.na(s))))
    if (sum(!is.na(s) & !is.na(g) & s != g) > 0) { w <- which(!is.na(s) & !is.na(g) & s != g)[1:3]; print(data.frame(strict = s[w], helper = g[w])) }
  }
  db <- strict(a, "Dbxref"); gid_strict <- sub(".*GeneID:([^,]+).*", "\\1", db); g <- gff_attribute(a, "GeneID:")
  cat(sprintf("%-5s GeneID: differ %d; NA %d; distinct %d\n", tr, sum(gid_strict != g, na.rm = TRUE), sum(is.na(g)), uniqueN(g)))
}

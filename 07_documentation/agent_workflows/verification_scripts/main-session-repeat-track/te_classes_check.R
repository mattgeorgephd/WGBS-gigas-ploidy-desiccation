# Independent check (main session) of the repeat-track numbers: which classes hold the analysed CpGs and DML.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); d <- a[1]; rm_out <- a[2]
x <- fread(cmd = sprintf("zcat %s | awk 'NR>3 && NF>=11 {print $5\"\\t\"$6\"\\t\"$7\"\\t\"$11}'", shQuote(rm_out)), header = FALSE, col.names = c("chr", "b", "e", "cf"))
x[, class := sub("/.*", "", cf)]
stopifnot(nrow(x) == 781312L)
u <- fread(file.path(d, "united_cpg_legacy.csv.gz"), select = c("chr", "start", "end"))
u[, row := .I]; stopifnot(nrow(u) == 480579L, all(u$end - u$start == 2L))
dp <- fread(file.path(d, "dml_ploidy_legacy.csv")); dh <- fread(file.path(d, "dml_heat_legacy.csv")); ds <- fread(file.path(d, "dml_shared_legacy.csv"))
cp <- u[, .(chr, s = start + 1L, e = start + 2L, row)]; setkey(cp, chr, s, e)
five <- c("DNA", "RC", "LTR", "LINE", "SINE"); nonte <- c("Simple_repeat", "Low_complexity", "Satellite", "rRNA")
hits <- function(shift) {            # shift = 1: read as BED (first base dropped); 0: as written
  r <- x[, .(chr, s = b + shift, e, class)][e >= s]; setkey(r, chr, s, e)
  o <- foverlaps(cp, r, type = "any", nomatch = NULL)
  o[, .(any = TRUE, five = any(class %in% five), unknown = any(class == "Unknown"), nonte = any(class %in% nonte),
        classes = paste(sort(unique(class)), collapse = "+")), by = row]
}
for (shift in c(1L, 0L)) {
  h <- hits(shift)
  cat(sprintf("\n===== reading: %s\n", if (shift) "as BED (legacy)" else "coordinates as written"))
  sets <- list(`all analysed` = u$row, `ploidy DML` = dp$united_row, `heat DML` = dh$united_row, shared = ds$united_row)
  res <- rbindlist(lapply(names(sets), function(n) { k <- h[row %in% sets[[n]]]
    data.table(set = n, n_loci = length(sets[[n]]), all_classes = nrow(k), five = sum(k$five), five_or_unknown = sum(k$five | k$unknown),
               unknown_not_five = sum(k$unknown & !k$five), only_simple_low_sat_rrna = sum(!k$five & !k$unknown)) }))
  print(res)
  cat("DML outside the five classes, by class combination:\n")
  print(rbind(h[row %in% dp$united_row & five == FALSE, .(set = "ploidy", n = .N), by = classes], h[row %in% dh$united_row & five == FALSE, .(set = "heat", n = .N), by = classes]))
}
cat("\nmean record length by class (bp):\n"); print(x[, .(records = .N, mean_bp = round(mean(e - b + 1L), 1), median_bp = as.numeric(median(e - b + 1L))), by = class][order(-records)])
cat("\nDNA records without a family:", x[cf == "DNA", .N], " Unknown:", x[cf == "Unknown", .N], "\n")

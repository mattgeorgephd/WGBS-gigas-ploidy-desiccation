suppressPackageStartupMessages(library(data.table))
options(width = 250)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
lib <- fread(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"))[, .(library_id, ploidy, treatment)]
first <- fread("first_lines.tsv", header = FALSE, col.names = c("library_id", "chr", "start", "end", "pct", "nC", "nT"))
first[, cov := nC + nT]
fp <- rbindlist(lapply(lib$library_id, function(l) { f <- sprintf("pass/fpos_%s.tsv", l)
  if (file.size(f) == 0) return(NULL)
  fread(f, header = FALSE, col.names = c("library_id", "chr", "start", "end", "pct", "nC", "nT")) }))
fp[, cov := nC + nT]
pos <- unique(first[, .(chr, start)])
cat("distinct first-line positions:", nrow(pos), "of 22 files\n")
m <- dcast(fp, chr + start ~ library_id, value.var = "cov", fill = 0L)
m <- merge(pos, m, by = c("chr", "start"), all.x = TRUE)
for (l in lib$library_id) if (!l %in% names(m)) m[, (l) := 0L]
setcolorder(m, c("chr", "start", lib$library_id))
mat <- as.matrix(m[, lib$library_id, with = FALSE]); mat[is.na(mat)] <- 0L
m[, first_line_of := vapply(seq_len(.N), function(i) paste(first[chr == m$chr[i] & start == m$start[i], library_id], collapse = ","), "")]
m[, `:=`(n_lib_1x = rowSums(mat >= 1), n_lib_3x = rowSums(mat >= 3), n_lib_5x = rowSums(mat >= 5), n_lib_10x = rowSums(mat >= 10), min_cov = apply(mat, 1, min), max_cov = apply(mat, 1, max))]
print(m[, .(chr, start, first_line_of, n_lib_1x, n_lib_3x, n_lib_5x, n_lib_10x, min_cov, max_cov)])
cat("\nfirst-line positions with >= 10 reads in all 22 libraries:", sum(m$n_lib_10x == 22), "\n")
# expected deficit of the legacy count against a complete read: positions that pass the threshold in every library
# of the set AND are the first line of at least one library of the set
sets <- list(all = lib$library_id, diploid = lib[ploidy == "diploid", library_id], triploid = lib[ploidy == "triploid", library_id],
             control = lib[treatment == "control", library_id], heat = lib[treatment == "heat", library_id])
exp <- rbindlist(lapply(names(sets), function(s) rbindlist(lapply(c(1, 3, 5, 10, 20), function(x) {
  inset <- vapply(seq_len(nrow(m)), function(i) any(first[chr == m$chr[i] & start == m$start[i], library_id] %in% sets[[s]]), TRUE)
  ok <- rowSums(mat[, sets[[s]], drop = FALSE] >= x) == length(sets[[s]])
  data.table(sample_set = s, min_coverage = x, expected_deficit = sum(ok & inset))
}))))
cc <- fread(file.path(repo, "04_output/04_cpg_coverage/tables/cpg_counts_by_coverage.csv"))
w <- dcast(cc, sample_set + min_coverage ~ mode, value.var = "cpg_in_all_libraries")
w[, deficit := corrected - legacy]
cat("\nNOTE: corrected mode also applies the upper filter, so corrected - legacy is not the first-line deficit; shown for reference only.\n")
w2 <- merge(w, exp, by = c("sample_set", "min_coverage"))
print(w2)
print(dcast(cc, sample_set + min_coverage ~ mode, value.var = "cpg_in_any_library"))

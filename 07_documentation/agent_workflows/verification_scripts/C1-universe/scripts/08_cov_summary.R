suppressPackageStartupMessages(library(data.table))
out <- commandArgs(trailingOnly = TRUE)[1]
rd <- function(f, nm) { d <- fread(file.path(out, f), header = FALSE, select = c(1, 5), col.names = c("gene_id", nm)); d[, gene_id := as.character(gene_id)]; d }
a <- rd("gene_cov_zr3644_14.tsv", "c_ctrl"); b <- rd("gene_cov_zr3534_8.tsv", "c_heat")
ann <- fread(file.path(out, "annotated_mrnas_cov.csv"), colClasses = c(gene_id = "character"))
ann <- merge(merge(ann, a, by = "gene_id"), b, by = "gene_id")
print(ann[, .(transcripts = .N, covered_ctrl_zr3644_14 = sum(c_ctrl > 0), covered_heat_zr3534_8 = sum(c_heat > 0), covered_either = sum(c_ctrl > 0 | c_heat > 0),
              pct_either = round(100 * mean(c_ctrl > 0 | c_heat > 0), 2), covered_union22_legacy = sum(c_union > 0)), by = status])
s <- fread(file.path(out, "sample40_block_genes.csv"), colClasses = c(gene_id = "character"))
s <- merge(merge(s[, .(gene_id)], a, by = "gene_id"), b, by = "gene_id")
cat("40 sampled block genes: covered in zr3644_14:", sum(s$c_ctrl > 0), "; in zr3534_8:", sum(s$c_heat > 0), "; in either:", sum(s$c_ctrl > 0 | s$c_heat > 0),
    "; range of covered CpGs zr3644_14:", paste(range(s$c_ctrl), collapse = "-"), "; zr3534_8:", paste(range(s$c_heat), collapse = "-"), "\n")
# weights: lines of the big table per transcript ~ union CpGs of the gene x GOslim rows of the transcript.
lt <- fread(cmd = "zcat /home/user/WGBS-gigas-ploidy-desiccation/03_input_files/annotation/legacy_go_goslim_terms.tsv.gz", header = TRUE, sep = "\t", quote = "")
srows <- lt[, .(s = .N), by = GOID]
long <- ann[, .(GOID = unlist(strsplit(go, ","))), by = transcript_id]
long <- merge(long, srows, by = "GOID", all.x = TRUE)
r <- long[, .(rows_known = sum(s, na.rm = TRUE), go_n = .N, go_known = sum(!is.na(s))), by = transcript_id]
ann <- merge(ann, r, by = "transcript_id")
ann[, w := as.numeric(c_union) * rows_known]
tot <- ann[hyphen == FALSE, sum(w)]
cat(sprintf("\nlower-bound line weights (union CpGs x GOslim rows of GO ids in the recovered table): all %.1f M; present %.1f M; block %.1f M (%.1f%% of all)\n",
            tot / 1e6, ann[in_uni == TRUE, sum(w)] / 1e6, ann[in_block == TRUE, sum(w)] / 1e6, 100 * ann[in_block == TRUE, sum(w)] / tot))
cat(sprintf("notebook: lines written for BOTH 157469551; ratio written / weight of present transcripts = %.3f\n", 157469551 / ann[in_uni == TRUE, sum(w)]))
cat(sprintf("expected complete BOTH table by the heat-only ratio (276062519 lines from 25779364 input lines): %.1f M lines; written %.1f%% of that\n",
            26098781 * 276062519 / 25779364 / 1e6, 100 * 157469551 / (26098781 * 276062519 / 25779364)))
cat("share of GO ids (transcript x GO pairs) found in the recovered GOslim table:", round(100 * sum(ann$go_known) / sum(ann$go_n), 1), "%\n")
# giant gene at the resume point
g <- ann[gene_id == "105323906"]
cat(sprintf("gene 105323906: %d annotated mRNAs, union CpGs %d, weight %.1f M lines = %.1f%% of all; mRNAs in universe %d, in block %d\n",
            nrow(g), g$c_union[1], sum(g$w) / 1e6, 100 * sum(g$w) / tot, sum(g$in_uni), sum(g$in_block)))
setorder(ann, acc_num); ann[, cum := cumsum(fifelse(hyphen, 0, w)) / tot]
cat(sprintf("cumulative weight at hole start (XM_011455141.3): %.3f; at hole end (XM_034461018.1): %.3f\n",
            ann[transcript_id == "XM_011455140.3", cum], ann[transcript_id == "XM_034461018.1", cum]))
cat(sprintf("observed file: junction at byte 5283561472 of 20149928521 = %.3f of the written file\n", 5283561472 / 20149928521))

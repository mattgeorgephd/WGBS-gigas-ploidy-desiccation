suppressPackageStartupMessages(library(data.table))
out <- commandArgs(trailingOnly = TRUE)[1]
ann <- fread(file.path(out, "annotated_mrnas_cov.csv"), colClasses = c(gene_id = "character"))
mem <- fread(file.path(out, "universe_membership_three_sets.csv"))
ann <- merge(ann, mem[, .(transcript_id, in_ron, in_yaa)], by = "transcript_id")
ann[, why := fifelse(in_uni, "present",
             fifelse(in_ron | in_yaa, "in the heat-only or control-only universe: lost to the failed write",
             fifelse(c_union == 0, "gene has no CpG with a read in any of 22 libraries (intended filter)",
             fifelse(hyphen, "hyphenated gene name: gene id lost by the legacy parser",
             fifelse(n_go_rec == 0, "no GO id in the recovered 2017 GOslim table (inferred cause)", "unexplained")))))]
print(ann[in_uni == FALSE, .N, by = .(in_block, why)][order(-in_block, -N)])
cat("missing in total:", ann[in_uni == FALSE, .N], "; lost to the write:", ann[in_uni == FALSE & (in_ron | in_yaa), .N],
    "; absent regardless:", ann[in_uni == FALSE & !(in_ron | in_yaa), .N], "\n")
cat("best estimate of the intended universe (union of the three legacy universes):", ann[in_uni | in_ron | in_yaa, .N], "\n")
cat("unexplained block transcripts:\n"); print(ann[why == "unexplained", .(transcript_id, gene_id, parent, c_union, n_go_rec, go = substr(go, 1, 60))])
# genes: how many genes lose every annotated transcript, and how many DML genes among them
g <- ann[, .(any_in_uni = any(in_uni), would_be_in = any(in_uni | in_ron | in_yaa)), by = gene_id]
cat("genes with annotated mRNAs:", nrow(g), "; with a transcript in the legacy universe:", sum(g$any_in_uni),
    sprintf("(%.1f%%)", 100 * mean(g$any_in_uni)), "; in the intended universe:", sum(g$would_be_in), "\n")

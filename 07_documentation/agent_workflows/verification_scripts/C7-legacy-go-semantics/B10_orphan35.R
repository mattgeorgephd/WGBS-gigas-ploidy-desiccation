suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
d <- file.path(repo, "08_archive/tables")
h <- fread(file.path(d, "DML-getMethylDiff-HEAT-Cov10-20.csv")); p <- fread(file.path(d, "DML-getMethylDiff-PLOIDY-Cov10-20.csv")); i <- fread(file.path(d, "DML-getMethylDiff-INTERACTION-Cov10-20.csv"))
mr <- fread(cmd = paste("zcat", file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_mRNA.gff.gz")), header = FALSE, sep = "\t", quote = "")
mr[, `:=`(gene_id = sub(".*GeneID:([0-9]+).*", "\\1", V9), product = sub(", transcript variant.*", "", sapply(sub(".*product=([^;]+).*", "\\1", V9), utils::URLdecode)))]
gff <- fread(cmd = paste("zcat", file.path(repo, "03_input_files/genome_features/cgigas_uk_roslin_v1_gene.gff.gz")), header = FALSE, sep = "\t", quote = "")
gff[, gene_id := sub(".*GeneID:([0-9]+).*", "\\1", V9)]
prod <- unique(mr[, .(gene_id, product)])[, .(product = paste(unique(product), collapse = "; ")), by = gene_id]
gff <- merge(gff, prod, by = "gene_id", all.x = TRUE)
genes <- GRanges(gff$V1, IRanges(gff$V4, gff$V5), gene_id = gff$gene_id, product = gff$product)
cat("genes whose product name contains 'adenylate kinase 9':\n"); print(gff[grepl("adenylate kinase 9", product), .(gene_id, V1, V4, V5, product)])
ov_tab <- function(x, label) { g <- GRanges(x$chr, IRanges(x$start + 1L, x$end)); o <- findOverlaps(g, genes)
  data.table(set = label, chr = x$chr[queryHits(o)], start = x$start[queryHits(o)], meth.diff = x$meth.diff[queryHits(o)], gene_id = genes$gene_id[subjectHits(o)], product = genes$product[subjectHits(o)]) }
all <- rbind(ov_tab(p, "ploidy"), ov_tab(h, "heat"), ov_tab(i, "interaction"))
targets <- c("105326303", "105319299", "105329309", "105340590", gff[grepl("adenylate kinase 9", product), gene_id])
w <- dcast(all[gene_id %in% targets], gene_id + product + chr + start ~ set, value.var = "meth.diff")
for (cn in c("ploidy", "heat", "interaction")) if (!cn %in% names(w)) w[, (cn) := NA_real_]
w[, class := fifelse(is.na(ploidy) | is.na(heat), "in one list only", fifelse(ploidy > 0 & heat > 0, "both positive (jointly hyper)", fifelse(ploidy < 0 & heat < 0, "both negative (jointly hypo)", "opposite signs (contrasting)")))]
options(width = 230)
print(w[, .(gene_id, product = substr(product, 1, 34), chr, start, ploidy = round(ploidy, 1), heat = round(heat, 1), interaction = round(interaction, 1), class)])
cat("\nper gene: DML in ploidy list / heat list / both lists\n")
print(w[, .(ploidy_dml = sum(!is.na(ploidy)), heat_dml = sum(!is.na(heat)), shared = sum(!is.na(interaction)), classes = paste(unique(class[!is.na(interaction)]), collapse = "; ")), by = .(gene_id, product = substr(product, 1, 40))])

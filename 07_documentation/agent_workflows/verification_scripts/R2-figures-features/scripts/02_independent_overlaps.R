# Independent recomputation of step 07 tables, without GenomicRanges and without the
# pipeline helpers: plain interval arithmetic in 0-based half-open coordinates.
suppressPackageStartupMessages(library(data.table))
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
fdir <- file.path(repo, "03_input_files/genome_features")
t05 <- file.path(repo, "04_output/05_differential_methylation/tables")
t07 <- file.path(repo, "04_output/07_genome_features/tables")
files <- c(gene = "gene.gff", mRNA = "mRNA.gff", exon = "exon.gff", CDS = "CDS.gff", exon_UTR = "exonUTR.gff",
           intron = "intron.bed", nonCDS = "nonCDS.bed", lncRNA = "lncRNA.gff", upstream_flank = "upstream.gff",
           downstream_flank = "downstream.gff", flanks = "flanks.gff", intergenic = "intergenic.bed", TE = "rm.te.bed")
read_track0 <- function(f) {                       # -> merged 0-based half-open intervals per sequence
  p <- file.path(fdir, paste0("cgigas_uk_roslin_v1_", f, ".gz"))
  if (grepl("gff$", f)) { d <- fread(p, header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5)); setnames(d, c("chr", "s", "e")); d[, s := s - 1L] }
  else { d <- fread(p, header = FALSE, sep = "\t", quote = "", select = 1:3); setnames(d, c("chr", "s", "e")) }
  d <- d[e > s]
  setorder(d, chr, s, e)
  d[, grp := cumsum(c(TRUE, s[-1] > cummax(e)[-.N])), by = chr]      # touching intervals merged or not makes no difference for overlap
  d[, .(s = min(s), e = max(e)), by = .(chr, grp)][, grp := NULL][]
}
tr <- lapply(files, read_track0)
ov_one <- function(cpg, t) {                        # cpg: chr, start (0-based C). TRUE if [start, start+2) meets a merged interval
  out <- logical(nrow(cpg))
  idx <- split(seq_len(nrow(cpg)), cpg$chr)
  tt <- split(t, t$chr)
  for (ch in intersect(names(idx), names(tt))) {
    i <- idx[[ch]]; x <- tt[[ch]]; c0 <- cpg$start[i]
    k <- findInterval(c0 + 1L, x$s)                 # last interval with s <= c0 + 1, i.e. s < c0 + 2
    hit <- k > 0; hit[hit] <- x$e[k[hit]] > c0[hit]
    out[i] <- hit
  }
  out
}
res <- list()
for (mode in c("legacy", "corrected")) {
  u <- fread(file.path(t05, sprintf("united_cpg_%s.csv.gz", mode)))
  nl <- sum(grepl("^coverage", names(u)))
  pm <- sapply(seq_len(nl), function(j) 100 * u[[paste0("numCs", j)]] / (u[[paste0("numCs", j)]] + u[[paste0("numTs", j)]]))
  mp <- rowMeans(pm)
  cat(sprintf("\n[%s] CpGs %d; end - start: %s; strand: %s; coverage min %d; any NA pm: %s\n", mode, nrow(u), paste(unique(u$end - u$start), collapse = ","),
              paste(unique(u$strand), collapse = ","), min(as.matrix(u[, grep("^coverage", names(u)), with = FALSE])), anyNA(pm)))
  cat(sprintf("[%s] mean methylation < 10: %.2f%%; > 90: %.2f%%; eligible [10, 90]: %d (%.2f%%)\n", mode, 100 * mean(mp < 10), 100 * mean(mp > 90), sum(mp >= 10 & mp <= 90), 100 * mean(mp >= 10 & mp <= 90)))
  cpg <- u[, .(chr, start)]; elig <- mp >= 10 & mp <= 90
  dp <- fread(file.path(t05, sprintf("dml_ploidy_%s.csv", mode))); dh <- fread(file.path(t05, sprintf("dml_heat_%s.csv", mode))); ds <- fread(file.path(t05, sprintf("dml_shared_%s.csv", mode)))
  key <- function(d) paste(d$chr, d$start)
  stopifnot(identical(key(cpg[dp$united_row]), key(dp)), identical(key(cpg[dh$united_row]), key(dh)), identical(key(cpg[ds$united_row]), key(ds)))
  stopifnot(!anyDuplicated(key(cpg)))
  k <- key(cpg)
  sets <- list(`all analysed CpGs` = rep(TRUE, nrow(cpg)), `eligible CpGs` = elig, `ploidy DML` = k %in% key(dp), `heat DML` = k %in% key(dh),
               shared = k %in% key(ds), `shared, both hyper` = k %in% key(ds[class == "both_hyper"]), `shared, both hypo` = k %in% key(ds[class == "both_hypo"]),
               `shared, opposite` = k %in% key(ds[class == "opposite"]))
  cat(sprintf("[%s] DML not eligible: ploidy %d, heat %d, shared %d; unique DML %d\n", mode, sum(sets$`ploidy DML` & !elig), sum(sets$`heat DML` & !elig), sum(sets$shared & !elig), sum(sets$`ploidy DML` | sets$`heat DML`)))
  if (any((sets$`ploidy DML` | sets$`heat DML`) & !elig)) { w <- which((sets$`ploidy DML` | sets$`heat DML`) & !elig); print(data.table(cpg[w], mean_pm = mp[w], ploidy = sets$`ploidy DML`[w], heat = sets$`heat DML`[w])) }
  ov <- sapply(names(files), function(t) ov_one(cpg, tr[[t]]))
  mine <- rbindlist(lapply(names(sets), function(s) data.table(mode = mode, set = s, n_loci = sum(sets[[s]]), track = names(files), n_overlap = colSums(ov[sets[[s]], , drop = FALSE]))))
  theirs <- fread(file.path(t07, sprintf("feature_overlap_by_track_%s.csv", mode)))
  cmp <- merge(mine, theirs, by = c("mode", "set", "track"), suffixes = c("_mine", "_pipeline"))
  cat(sprintf("[%s] overlap-by-track rows compared: %d; n_loci differ: %d; n_overlap differ: %d\n", mode, nrow(cmp), cmp[n_loci_mine != n_loci_pipeline, .N], cmp[n_overlap_mine != n_overlap_pipeline, .N]))
  if (cmp[n_overlap_mine != n_overlap_pipeline, .N] > 0) print(cmp[n_overlap_mine != n_overlap_pipeline])
  eight <- c("exon_UTR", "CDS", "intron", "upstream_flank", "downstream_flank", "lncRNA", "TE", "intergenic")
  mine[, share := ifelse(track %in% eight, 100 * n_overlap / sum(n_overlap[track %in% eight]), NA_real_), by = set]
  cmp2 <- merge(mine, theirs, by = c("mode", "set", "track"))
  cat(sprintf("[%s] largest |legacy_share_pct difference|: %.3g; largest |pct_of_loci difference|: %.3g\n", mode, cmp2[, max(abs(share - legacy_share_pct), na.rm = TRUE)], cmp2[, max(abs(100 * n_overlap.x / n_loci.x - pct_of_loci))]))
  # exclusive class
  prio <- c("CDS", "lncRNA", "exon_UTR", "intron", "upstream_flank", "downstream_flank", "intergenic")
  first_hit <- max.col(ov[, prio], ties.method = "first"); cls <- ifelse(rowSums(ov[, prio]) == 0, "other", prio[first_hit])
  ex <- rbindlist(lapply(names(sets), function(s) { x <- cls[sets[[s]]]; te <- ov[sets[[s]], "TE"]
    data.table(mode = mode, set = s, class = c(prio, "other"), n = as.vector(table(factor(x, c(prio, "other")))),
               pct_in_TE = sapply(c(prio, "other"), function(cl) if (any(x == cl)) 100 * mean(te[x == cl]) else NA_real_)) }))
  ext <- fread(file.path(t07, sprintf("feature_class_exclusive_%s.csv", mode)))
  cmpx <- merge(ex, ext, by = c("mode", "set", "class"), suffixes = c("_mine", "_pipeline"))
  cat(sprintf("[%s] exclusive class rows compared: %d; n differ: %d; max |pct_in_TE diff|: %.3g; pct sums by set: %s\n", mode, nrow(cmpx), cmpx[n_mine != n_pipeline, .N],
              cmpx[, max(abs(pct_in_TE_mine - pct_in_TE_pipeline), na.rm = TRUE)], paste(round(ext[, sum(pct), by = set]$V1, 6), collapse = " ")))
  # what is "other"?
  oth <- cls == "other"
  cat(sprintf("[%s] class other: %d CpGs; on mito NC_001276.1: %d; in gene track: %d; in exon: %d; in TE: %d; in nonCDS: %d; in flanks: %d\n", mode, sum(oth), sum(oth & cpg$chr == "NC_001276.1"),
              sum(oth & ov[, "gene"]), sum(oth & ov[, "exon"]), sum(oth & ov[, "TE"]), sum(oth & ov[, "nonCDS"]), sum(oth & ov[, "flanks"])))
  cat(sprintf("[%s] analysed CpGs on mito: %d; DML on mito: ploidy %d heat %d\n", mode, sum(cpg$chr == "NC_001276.1"), sum(dp$chr == "NC_001276.1"), sum(dh$chr == "NC_001276.1")))
  cat(sprintf("[%s] CpGs in more than one of the 7 class tracks: %d (%.2f%%); CDS&intron %d; CDS&exon_UTR %d; lncRNA&CDS %d; lncRNA&intron %d; lncRNA&exon_UTR %d; lncRNA not in (exon_UTR or intron): %d of %d\n", mode,
              sum(rowSums(ov[, prio]) > 1), 100 * mean(rowSums(ov[, prio]) > 1), sum(ov[, "CDS"] & ov[, "intron"]), sum(ov[, "CDS"] & ov[, "exon_UTR"]), sum(ov[, "lncRNA"] & ov[, "CDS"]),
              sum(ov[, "lncRNA"] & ov[, "intron"]), sum(ov[, "lncRNA"] & ov[, "exon_UTR"]), sum(ov[, "lncRNA"] & !(ov[, "exon_UTR"] | ov[, "intron"] | ov[, "CDS"])), sum(ov[, "lncRNA"])))
  # Fisher tests
  tt <- fread(file.path(t07, sprintf("feature_association_tests_%s.csv", mode)))
  bg <- list(`all analysed CpGs` = rep(TRUE, nrow(cpg)), `eligible CpGs` = elig)
  my <- rbindlist(lapply(seq_len(nrow(tt)), function(r) { b <- tt$background[r]; s <- tt$set[r]; t <- tt$track[r]
    d <- sets[[s]]; o <- ov[, t]; comp <- bg[[b]] & !d
    a <- sum(d & o); bb <- sum(d & !o); cc <- sum(comp & o); dd <- sum(comp & !o)
    ft <- fisher.test(matrix(c(a, cc, bb, dd), 2))
    data.table(a = a, dml_total = a + bb, comp_in = cc, comp_total = cc + dd, or = unname(ft$estimate), p = ft$p.value, sample_or = (a / bb) / (cc / dd)) }))
  cat(sprintf("[%s] Fisher tests compared: %d; dml_in_track differ: %d; dml_total differ: %d; max rel OR diff: %.3g; max rel p diff: %.3g; max |pct_other diff|: %.3g; BH q max abs diff: %.3g\n", mode, nrow(tt),
              sum(my$a != tt$dml_in_track), sum(my$dml_total != tt$dml_total), max(abs(my$or / tt$odds_ratio - 1)), max(abs(my$p - tt$p_fisher) / pmax(tt$p_fisher, 1e-300)),
              max(abs(100 * my$comp_in / my$comp_total - tt$pct_other_cpg)), max(abs(p.adjust(my$p, "BH") - tt$q_BH))))
  res[[mode]] <- cbind(tt[, .(mode, background, set, track, dml_in_track, dml_total, odds_ratio, p_fisher, q_BH)], my[, .(comp_in, comp_total, sample_or)])
  # chi-squared
  cq <- fread(file.path(t07, sprintf("feature_class_chisq_%s.csv", mode)))
  for (s in c("ploidy DML", "heat DML", "shared")) {
    m <- table(sets[[s]], cls); m <- m[, colSums(m) > 0, drop = FALSE]
    ct <- suppressWarnings(chisq.test(m))
    cat(sprintf("[%s] chisq %s: mine X2 %.4f df %d p %.3g min expected %.3f | pipeline X2 %.4f df %d p %.3g min expected %.3f | classes with expected < 5: %d\n", mode, s, ct$statistic, ct$parameter, ct$p.value, min(ct$expected),
                cq[set == s, statistic], cq[set == s, df], cq[set == s, p], cq[set == s, min_expected], sum(ct$expected < 5)))
  }
  # regions, independent algorithm: single linkage on sorted starts per chr
  rg <- fread(file.path(t07, sprintf("dml_regions_%s.csv", mode)))
  for (nm in c("ploidy DML", "heat DML", "shared")) {
    d <- switch(nm, `ploidy DML` = dp, `heat DML` = dh, shared = ds)
    out <- d[order(chr, start), { g <- 1L; ids <- integer(.N); if (.N > 0) { ids[1] <- 1L; if (.N > 1) for (i in 2:.N) { if (start[i] - start[i - 1] > 1000) g <- g + 1L; ids[i] <- g } }; .(start = start, end = end, id = ids) }, by = chr]
    r <- out[, .(start = min(start), end = max(end), n_dml = .N), by = .(chr, id)]
    p <- rg[set == nm]
    m <- merge(r, p, by = c("chr", "start", "end"), all = TRUE)
    cat(sprintf("[%s] regions %s: mine %d, pipeline %d, unmatched %d, n_dml differ %d; span check failures %d; n_hyper+n_hypo != n_dml: %d; largest %d\n", mode, nm, nrow(r), nrow(p), m[is.na(n_dml.x) | is.na(n_dml.y), .N],
                m[!is.na(n_dml.x) & !is.na(n_dml.y) & n_dml.x != n_dml.y, .N], p[span_bp != end - start, .N], p[n_hyper + n_hypo != n_dml, .N], max(r$n_dml)))
    # min gap between consecutive regions on the same sequence
    setorder(r, chr, start); gap <- r[, .(gap = start[-1] - end[-.N]), by = chr]
    cat(sprintf("       smallest gap between neighbouring regions (start minus previous end): %s; pairs with start-to-start gap in (1000, 1002]: %d\n", min(gap$gap), out[, sum(diff(start) > 1000 & diff(start) <= 1002), by = chr][, sum(V1)]))
  }
  # linkage groups
  lg <- fread(file.path(t07, sprintf("dml_by_linkage_group_%s.csv", mode)))
  lgn <- sprintf("NC_0475%02d.1", 59:68); lab <- function(ch) ifelse(ch %in% lgn, paste0("LG", match(ch, lgn)), ifelse(ch == "NC_001276.1", "mitochondrial", "unplaced"))
  mine_lg <- data.table(sequence = lab(cpg$chr))[, .(n_cpg = .N), by = sequence]
  for (x in list(list(dp, "meth.diff", "p"), list(dh, "meth.diff", "h"), list(ds, "meth.diff_heat", "s"))) {
    a <- data.table(sequence = lab(x[[1]]$chr), v = x[[1]][[x[[2]]]])[, .(n = .N, hy = sum(v > 0), ho = sum(v < 0)), by = sequence]; setnames(a, c("n", "hy", "ho"), paste0(x[[3]], c("", "_hyper", "_hypo")))
    mine_lg <- merge(mine_lg, a, by = "sequence", all.x = TRUE)
  }
  for (j in names(mine_lg)[-1]) set(mine_lg, which(is.na(mine_lg[[j]])), j, 0L)
  m <- merge(mine_lg, lg, by = "sequence")
  cat(sprintf("[%s] linkage groups rows %d (pipeline %d): differing cells: n_cpg %d, ploidy %d/%d/%d, heat %d/%d/%d, shared %d/%d/%d\n", mode, nrow(m), nrow(lg), m[n_cpg != n_cpg_analysed, .N],
              m[p != ploidy_dml, .N], m[p_hyper != ploidy_dml_hyper, .N], m[p_hypo != ploidy_dml_hypo, .N], m[h != heat_dml, .N], m[h_hyper != heat_dml_hyper, .N], m[h_hypo != heat_dml_hypo, .N],
              m[s != shared, .N], m[s_hyper != shared_hyper, .N], m[s_hypo != shared_hypo, .N]))
  g <- fread(file.path(fdir, "cgigas_uk_roslin_v1_gene.gff.gz"), header = FALSE, sep = "\t", quote = "", select = 1)
  gl <- data.table(sequence = lab(g$V1))[, .(n_genes_mine = .N), by = sequence]
  m <- merge(lg, gl, by = "sequence", all.x = TRUE); m[is.na(n_genes_mine), n_genes_mine := 0L]
  cat(sprintf("[%s] genes per linkage group differ: %d; total genes in table %d of %d\n", mode, m[n_genes != n_genes_mine, .N], sum(lg$n_genes), nrow(g)))
  print(lg[, .(sequence, n_cpg_analysed, ploidy_dml, heat_dml, shared, n_genes, all3 = ploidy_dml + heat_dml + shared, pct_all3 = round(100 * (ploidy_dml + heat_dml + shared) / sum(ploidy_dml + heat_dml + shared), 2),
               pct_p = round(pct_of_ploidy_dml, 2), pct_h = round(pct_of_heat_dml, 2))])
  rm(u, pm, ov); gc(verbose = FALSE)
}
fwrite(rbindlist(res), "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R2-figures-features/out/02_fisher_recomputed.csv")

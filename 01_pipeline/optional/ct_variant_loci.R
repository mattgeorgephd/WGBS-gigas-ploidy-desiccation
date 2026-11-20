# ct_variant_loci.R
#
# Purpose   CpG sites at C-to-T sequence variants in the 22 libraries, found by
#           the rule of BS-Snper (Gao et al. 2015) as Venkataraman et al. 2022
#           applied it, computed in R from the sorted, deduplicated Bismark
#           alignments of the June 2021 run. BS-Snper itself is not run.
#           A CpG is listed when a call at the C of the CpG on the plus strand
#           (1-based position start + 1, start being the 0-based C position of
#           the coverage files) has REF C and ALT T or "T,G", whatever its
#           FILTER label, in the pooled call of all 22 libraries (the summed
#           per-library counts, which equals BS-Snper on a merged BAM) or in any
#           one library's call. BS-Snper defaults, except minimum coverage 5.
#           Every CpG site of the reference is screened once, so that one table
#           serves the 5x and the 10x analyses.
#
# Decisions of the author (07_documentation/reproducibility_plan.md, section 7,
# table "Decided"):
#   2026-11-20  R code only, BS-Snper is not run; exclusion is an option, off by
#               default (not wired into run_pipeline.R).
#   2026-11-20  Reads are counted with their full CIGAR (insertions, deletions
#               and clipping handled). This is the one stated difference from
#               2022: BS-Snper lays a read out with its last CIGAR operation
#               only. That behaviour is kept as an internal test switch
#               (cigar_mode = "bs_snper_last_op" of ct_count_window()), used only
#               by the regression test against the 2022 BS-Snper tables.
#   2026-11-20  Every CpG site of the genome is screened; plus-strand C only;
#               any FILTER; union of the pooled call and the 22 library calls.
#
# Run       WGBS_CACHE_DIR=/home/user/gannet_cache \
#             Rscript 01_pipeline/optional/ct_variant_loci.R
#           (from the repository root). Many hours: each alignment (3.3 to
#           6.8 GB) is downloaded, checked, counted and deleted in turn. The
#           run is resumable: start the same command again after an
#           interruption or a machine restart; finished libraries, finished
#           chromosomes and verified downloads are kept in the cache.
#           Options (testing only):
#             --libraries id1,id2   process only these libraries; no final tables
#             --keep-bam            do not delete an alignment after counting it
#
# Inputs    03_input_files/gannet_manifest_alignments.csv  (22 BAM and 22 index
#             rows: URL and size; gannet publishes no md5 for them)
#           03_input_files/metadata/wgbs_libraries.csv     (library ids, experiments)
#           Reference genome cgigas_uk_roslin_v1_genomic-mito.fa, fetched through
#             fetch_gannet_file() from 03_input_files/gannet_manifest.csv (md5 checked)
#           Alignments, downloaded into <gannet cache>/alignments/ one at a time
#             and deleted after use. Nothing is read from the outputs of other steps.
#
# Outputs   04_output/not_adopted/ct_variants/tables/
#             ct_variant_cpgs.csv.gz           one row per CpG with a C-to-T call
#                                              in any call set (see below)
#             ct_variant_calls.csv.gz          one row per (CpG, call set) with a
#                                              C-to-T call: 22 libraries, the pool
#                                              of all 22, and the pools of each
#                                              experiment (zr3534, zr3644)
#             ct_variant_assessability.csv.gz  one row per CpG of the genome:
#                                              pooled evidence-strand reads and
#                                              libraries able to give a call
#             ct_variant_run.csv               parameters, CIGAR mode, genome md5,
#                                              per-library size, md5, reads, times
#             ct_variant_log.txt               summary of the run
#           04_output/not_adopted/ct_variants/cache/  (git-ignored; per-library
#             counts at every CpG, enough to recompute any pooling or rule
#             parameter without the alignments; console log; run lock)
#
# Parameters  The rule parameters are defined below in one list, ct_variants,
#           with the names they will have in run_config.R (config$ct_variants);
#           when run_config.R defines config$ct_variants, that list is used
#           instead. The rule values are the 2022 values and are not toggles.
#
# Notes     Read classes as in BS-Snper: Watson = FLAG 0x10 and 0x80 both unset
#           or both set (FLAG 99, 147), Crick = otherwise (83, 163). Both mates
#           of a pair are counted, as BS-Snper does. Reads with FLAG 0x4 or MAPQ
#           below mapvalue are dropped; bases below minquali become N; a read
#           with more than max(1, round(L x errorate)) mismatches (bisulfite
#           conversions not counted; L = reference span) is dropped. At a
#           reference C, Watson T counts as C; Crick reads carry the genotype.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "Rsamtools", "GenomicAlignments", "Biostrings", "parallel"))

# ---- Parameters (to move to run_config.R as config$ct_variants) -------------------
ct_variants_default <- list(
  # BS-Snper rule, 2022 values (Venkataraman et al. 2022: defaults, mincover 5)
  minhetfreq = 0.1, minhomfreq = 0.85, minquali = 15L, mincover = 5L, minread2 = 2L,
  errorate = 0.02, mapvalue = 20L,
  # how the rule is applied, as in 2022 (fixed, not toggles)
  cigar = "full",            # author's decision of 2026-11-20 (the one difference from 2022)
  site = "C",                # the C of the CpG on the plus strand only
  filter = "any",            # PASS and Low records both list a CpG
  union = "pooled_all_and_each_library",
  # assessability: a CT call needs at least 3 Crick T reads (BS-Snper post-rule)
  min_crick_assessable = 3L,
  # computing
  window_bp = 500000L, max_reads_per_batch = 120000L,
  download_streams = 4L, download_range_mb = 256, download_retries = 10L,
  download_range_timeout_s = 1800L, user_agent = "wgbs-pipeline/1.0 (R; curl)"
)
ctv <- if (!is.null(config$ct_variants)) config$ct_variants else ct_variants_default
ct_code_version <- "ct_variant_loci 1.0 (2026-11-20)"

step <- "ct_variants"
step_rel <- file.path("not_adopted", "ct_variants")
cell_names <- c("w_A", "w_T", "w_C", "w_G", "c_A", "c_T", "c_C", "c_G")   # BS-Snper BSD order
q_names <- paste0("q", cell_names)
bgzf_eof <- as.raw(c(0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
                     0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00))

# =============================================================================
# Counting (one window of one alignment)
# =============================================================================

# Reads of one window, as scanBam returns them. A remote BamFile must be new for
# every query (a reused remote handle returns nothing after the first query).
ct_scan_window <- function(bf, chr, ws, we) {
  p <- ScanBamParam(which = GRanges(chr, IRanges(ws, we)), what = c("flag", "pos", "mapq", "cigar", "seq", "qual"))
  scanBam(bf, param = p)[[1]]
}

# Byte lookups (ASCII): base -> BSD cell, for Watson (T folded into C) and Crick.
ct_lut <- local({
  w <- integer(256); k <- integer(256)
  w[utf8ToInt("A") + 1] <- 1L; w[utf8ToInt("T") + 1] <- 3L; w[utf8ToInt("C") + 1] <- 3L; w[utf8ToInt("G") + 1] <- 4L
  k[utf8ToInt("A") + 1] <- 5L; k[utf8ToInt("T") + 1] <- 6L; k[utf8ToInt("C") + 1] <- 7L; k[utf8ToInt("G") + 1] <- 8L
  list(w = w, k = k)
})
rawA <- charToRaw("A"); rawC <- charToRaw("C"); rawG <- charToRaw("G"); rawT <- charToRaw("T"); rawN <- charToRaw("N")

# Counts at the CpG Cs of [ws, we] from the reads `r` (scanBam list). chr_seq: the
# chromosome as raw ASCII (upper case); cpg_pos: sorted 1-based C positions of the
# chromosome's CpGs. Returns the CpG indices (into cpg_pos) of the window, an
# integer matrix of counts (8 BSD cells) and of quality sums (8 cells), and read
# statistics over the reads that start in the window.
# cigar_mode "full": every M/=/X block of the CIGAR is laid out at its reference
# position; insertions and soft-clipped bases are not on the reference, deletions
# are reference positions without a base; L = reference span.
# cigar_mode "bs_snper_last_op" (internal test switch, regression test only):
# BS-Snper's layout, the read laid out from POS with its last CIGAR operation only
# (the first n bases for a final nM; nothing for a final I, D or S); L = n.
ct_count_window <- function(r, chr_seq, ws, we, cpg_pos, ctv, cigar_mode = "full") {
  stopifnot(cigar_mode %in% c("full", "bs_snper_last_op"))
  i0 <- findInterval(ws - 1L, cpg_pos) + 1L; i1 <- findInterval(we, cpg_pos)
  nw <- max(0L, i1 - i0 + 1L)
  cnt <- matrix(0L, nw, 8L, dimnames = list(NULL, cell_names)); qs <- cnt; colnames(qs) <- q_names
  flag <- r$flag; n <- length(flag)
  starts_here <- r$pos >= ws & r$pos <= we
  st <- c(records = sum(starts_here), unmapped = sum(starts_here & bitwAnd(flag, 4L) > 0L), mapq_pass = 0,
          mismatch_pass = 0, with_indel = 0, f99 = sum(starts_here & flag == 99L), f147 = sum(starts_here & flag == 147L),
          f83 = sum(starts_here & flag == 83L), f163 = sum(starts_here & flag == 163L))
  st["other_flag"] <- st["records"] - sum(st[c("f99", "f147", "f83", "f163")])
  keep <- which(bitwAnd(flag, 4L) == 0L & !is.na(r$mapq) & r$mapq >= ctv$mapvalue)
  if (nw == 0L || length(keep) == 0L) {
    return(list(idx = integer(0), cnt = cnt[0, , drop = FALSE], qs = qs[0, , drop = FALSE], stats = st))
  }
  batches <- split(keep, ceiling(seq_along(keep) / ctv$max_reads_per_batch))
  lim_q <- as.raw(ctv$minquali + 33L)
  for (sel in batches) {
    cig <- r$cigar[sel]; pos <- r$pos[sel]
    watson <- (bitwAnd(flag[sel], 16L) > 0L) == (bitwAnd(flag[sel], 128L) > 0L)
    seqs <- r$seq[sel]
    s_all <- charToRaw(as.character(unlist(seqs, use.names = FALSE)))
    q_all <- charToRaw(as.character(unlist(r$qual[sel], use.names = FALSE)))
    qoff <- c(0L, cumsum(width(seqs)))[seq_along(sel)]
    if (cigar_mode == "full") {
      qr <- cigarRangesAlongQuerySpace(cig, ops = c("M", "=", "X"))
      rr <- cigarRangesAlongReferenceSpace(cig, pos = pos, ops = c("M", "=", "X"))
      bid <- rep.int(seq_along(sel), lengths(qr))
      uq <- unlist(qr, use.names = FALSE); ur <- unlist(rr, use.names = FALSE)
      w <- width(uq)
      rid <- rep.int(bid, w)
      within <- sequence(w) - 1L
      qi <- rep.int(qoff[bid] + start(uq), w) + within
      rp <- rep.int(start(ur), w) + within
      L <- cigarWidthAlongReferenceSpace(cig)
      rm(qr, rr, uq, ur, bid)
    } else {
      lo <- substring(cig, nchar(cig))
      ll <- as.integer(sub("^.*?([0-9]+)[^0-9]$", "\\1", cig, perl = TRUE))
      w <- ifelse(lo == "M", ll, 0L)
      rid <- rep.int(seq_along(sel), w)
      within <- sequence(w) - 1L
      qi <- rep.int(qoff + 1L, w) + within
      rp <- rep.int(pos, w) + within
      L <- ll
    }
    st["with_indel"] <- st["with_indel"] + sum(pos >= ws & pos <= we & grepl("[ID]", cig))
    rm(within)
    rp[rp > length(chr_seq)] <- NA_integer_
    s <- s_all[qi]; q <- q_all[qi]; rf <- chr_seq[rp]
    rf[is.na(rp)] <- rawN
    rm(s_all, q_all, qi)
    masked <- s == rawN | rf == rawN | q < lim_q
    wv <- watson[rid]
    conv <- (wv & rf == rawC & s == rawT) | (!wv & rf == rawG & s == rawA)
    mm <- which(!masked & s != rf & !conv)
    rm(conv)
    nmm <- tabulate(rid[mm], nbins = length(sel))
    pass <- nmm <= pmax(1L, as.integer(L * ctv$errorate + 0.5))
    in_win <- pos >= ws & pos <= we
    st["mapq_pass"] <- st["mapq_pass"] + sum(in_win)
    st["mismatch_pass"] <- st["mismatch_pass"] + sum(in_win & pass)
    tg <- which(!masked & rf == rawC & rp >= ws & rp <= we)
    tg <- tg[pass[rid[tg]]]
    tg <- tg[chr_seq[rp[tg] + 1L] == rawG]
    if (length(tg) > 0L) {
      ci <- findInterval(rp[tg], cpg_pos)
      stopifnot(all(cpg_pos[ci] == rp[tg]))
      b <- as.integer(s[tg]) + 1L
      cell <- ifelse(wv[tg], ct_lut$w[b], ct_lut$k[b])
      key <- (ci - i0) * 8L + cell
      qv <- as.integer(q[tg]) - 33L
      cnt <- cnt + matrix(tabulate(key, nbins = nw * 8L), nw, 8L, byrow = TRUE)
      o <- order(key, method = "radix"); ks <- key[o]; cs <- cumsum(as.numeric(qv[o]))
      last <- which(c(ks[-1L] != ks[-length(ks)], TRUE))
      sums <- diff(c(0, cs[last]))
      qv_tab <- numeric(nw * 8L); qv_tab[ks[last]] <- sums
      qs <- qs + matrix(as.integer(qv_tab), nw, 8L, byrow = TRUE)
    }
    rm(s, q, rf, rp, rid, masked, wv, mm, tg)
  }
  any_cov <- which(rowSums(cnt) > 0L)
  list(idx = (i0:i1)[any_cov], cnt = cnt[any_cov, , drop = FALSE], qs = qs[any_cov, , drop = FALSE], stats = st)
}

# =============================================================================
# The BS-Snper decision at a reference C (printSnp candidate filter, Perl Bayes
# genotype with its priors and post-rules, and the record written by genotype())
# =============================================================================
# cnt, qs: matrices with the 8 BSD cells (counts and quality sums), reference C.
# Returns, for every row, candidate (TRUE/FALSE) and for candidates the genotype,
# QUAL, and for C-to-T genotypes (CT, TT, GT) alt, gt, dp, alt_reads, filter.
ct_call <- function(cnt, qs, ctv) {
  n <- nrow(cnt)
  g <- function(m, k) as.numeric(m[, k])
  wA <- g(cnt, 1); wT <- g(cnt, 2); wC <- g(cnt, 3); wG <- g(cnt, 4)
  cA <- g(cnt, 5); cT <- g(cnt, 6); cC <- g(cnt, 7); cG <- g(cnt, 8)
  m <- wA + wT + wC + wG + cA + cT + cC + cG
  rate <- function(v, d) d > 0 & (v / d) > ctv$minhetfreq
  cand <- rate(wA + cA, m) | rate(cT, cC + cT) | rate(wG + cG + cA, m)
  cand[is.na(cand)] <- FALSE
  out <- data.table(candidate = cand, genotype = NA_character_, qual = NA_integer_, alt = NA_character_,
                    gt = NA_character_, dp = NA_integer_, alt_reads = NA_integer_, filter = NA_character_)
  ix <- which(cand)
  if (length(ix) == 0L) return(out)
  s <- function(v) v[ix]
  wA <- s(wA); wT <- s(wT); wC <- s(wC); wG <- s(wG); cA <- s(cA); cT <- s(cT); cC <- s(cC); cG <- s(cG)
  mq <- function(k, nn) { v <- as.numeric(qs[ix, k]); ifelse(nn > 0, floor(v / nn + 0.5), 0) }
  qwA <- mq(1, wA); qwT <- mq(2, wT); qwC <- mq(3, wC); qwG <- mq(4, wG)
  qcA <- mq(5, cA); qcT <- mq(6, cT); qcC <- mq(7, cC); qcG <- mq(8, cG)
  e <- function(qv) as.numeric(sprintf("%.3e", 0.1^(qv / 10)))
  An <- wA; Tn <- cT; Cn <- wC + cC; Gn <- wG + cG
  nn <- lfactorial(An + Tn + Cn + Gn) - lfactorial(An) - lfactorial(Tn) - lfactorial(Cn) - lfactorial(Gn)
  eA <- pmin(e(qwA), e(qcA)); eT <- pmin(e(qwT), e(qcT)); eC <- pmin(e(qwC), e(qcC)); eG <- pmin(e(qwG), e(qcG))
  eCT <- e(qcT)
  hasA <- qwA > 0 | qcA > 0; hasT <- qcT > 0 | qwT > 0; hasC <- qcC > 0 | qwC > 0; hasG <- qwG > 0 | qcG > 0
  yes <- function(cond, v) ifelse(cond, v, NA_real_)
  L <- cbind(
    aa = yes(hasA, nn + An * log(1 - eA) + (Tn + Cn + Gn) * log(eA / 3)),
    tt = yes(hasT, nn + Tn * log(1 - eT) + (An + Cn + Gn) * log(eT / 3)),
    cc = yes(hasC, nn + Cn * log(1 - eC) + (An + Tn + Gn) * log(eC / 3)),
    gg = yes(hasG, nn + Gn * log(1 - eG) + (An + Cn + Tn) * log(eG / 3)),
    at = yes(hasA & (cT > 0 | wT > 0), nn + (An + Tn) * log((1 - (eA + eT) / 2) / 2) + (Cn + Gn) * log((eA + eT) / 4)),
    ac = yes(hasA & (cC > 0 | wC > 0), nn + (An + Cn) * log((1 - eA) / 2) + (Tn + Gn) * log(eA / 3)),
    ag = yes(hasA & (cG > 0 | wG > 0), nn + (An + Gn) * log((1 - eA) / 2) + (Cn + Tn) * log(eA / 3)),
    ct = yes(hasT & (wC > 0 | cC > 0), nn + (Cn + Tn) * log((1 - eT) / 2) + (An + Gn) * log(eCT / 3)),
    gt = yes(hasT & (wG > 0 | cG > 0), nn + (Tn + Gn) * log((1 - eT) / 2) + (An + Cn) * log(eT / 3)),
    cg = yes(hasC & (wG > 0 | cG > 0), nn + (Cn + Gn) * log((1 - eC) / 2) + (An + Tn) * log(eC / 2)))
  L[!is.na(L) & L == 0] <- NA_real_            # BS-Snper drops genotypes whose score is exactly 0
  prior_c <- c(aa = log(0.000083), tt = log(0.00033), cc = log(0.985), gg = log(0.000083), at = log(1.1) - 7 * log(10),
               ac = log(0.00017), ag = log(2.78) - 8 * log(10), ct = log(0.000667), gt = log(1.1) - 7 * log(10),
               cg = log(0.00017))
  L <- sweep(L, 2, prior_c[colnames(L)], "+")
  best <- rep(NA_integer_, length(ix)); bv <- rep(-Inf, length(ix))
  for (j in seq_len(ncol(L))) {                # first maximum wins, as which.max over this order
    up <- !is.na(L[, j]) & (is.na(best) | L[, j] > bv)
    best[up] <- j; bv[up] <- L[up, j]
  }
  nkept <- rowSums(!is.na(L))
  G <- ifelse(is.na(best), "NN", toupper(colnames(L)[best]))
  den <- rowSums(2.7^L, na.rm = TRUE)
  pq <- 1 - 2.7^bv / den
  qual <- ifelse(nkept > 1, ifelse(den == 0, 0, ifelse(pq == 0, 1000, suppressWarnings(as.integer(-10 * log10(pq))))), 0)
  # one genotype scored: QUAL from the number of reads of the homozygous allele (BS-Snper.pl 1543-1566)
  one <- which(nkept == 1)
  if (length(one)) {
    gb <- colnames(L)[best[one]]
    hom <- ifelse(gb == "aa", (wA + cA)[one], ifelse(gb == "tt", (cT + wT)[one], ifelse(gb == "cc", (wC + cC + wT)[one],
                  ifelse(gb == "gg", (wG + cG + cA)[one], NA_real_))))
    p1 <- ifelse(is.na(hom), 1, 1 - 1 / (1 + 0.5^hom))
    qual[one] <- ifelse(p1 == 0, 1000, suppressWarnings(as.integer(-10 * log10(p1))))
  }
  G[G == "CT" & cT == 0] <- "CC"
  G[G == "CT" & cT < 3] <- "NN"
  G[G == "AG" & wA == 0] <- "GG"
  G[G == "AG" & wA < 3] <- "NN"
  G[G == "TT" & cC > 1 & cC / (cC + cT) > 0.02] <- "NN"
  G[G == "AA" & wG > 1 & wG / (wA + wG) > 0.02] <- "NN"
  r3 <- function(x) as.numeric(sprintf("%.3f", x))
  alt <- gt <- filt <- rep(NA_character_, length(ix)); dp <- altn <- rep(NA_real_, length(ix))
  k <- G == "CT"                                # BS-Snper.pl "C>CT, only use crick"
  dp[k] <- (cC + cT)[k]; altn[k] <- cT[k]; alt[k] <- "T"; gt[k] <- "0/1"
  filt[k] <- ifelse(dp[k] >= ctv$mincover & qcT[k] >= ctv$minquali & altn[k] >= ctv$minread2 &
                      r3(altn[k] / dp[k]) >= ctv$minhetfreq, "PASS", "Low")
  k <- G == "TT"                                # "C>TT"
  dp[k] <- (cT + cC + wT + wC)[k]; altn[k] <- (cT + wT)[k]; alt[k] <- "T"; gt[k] <- "1/1"
  filt[k] <- ifelse(dp[k] >= ctv$mincover & qcT[k] >= ctv$minquali & altn[k] >= ctv$minread2 &
                      r3(altn[k] / dp[k]) >= ctv$minhomfreq, "PASS", "Low")
  k <- G == "GT"                                # "C>GT", ALT "T,G"
  varG <- wG + cG + wA; varT <- cT + wT; qG <- pmax(qwG, qcG)
  dp[k] <- (varG + varT)[k]; altn[k] <- varT[k]; alt[k] <- "T,G"; gt[k] <- "1/2"
  filt[k] <- ifelse(dp[k] >= ctv$mincover & qG[k] >= ctv$minquali & qcT[k] >= ctv$minquali & varG[k] >= ctv$minread2 &
                      varT[k] >= ctv$minread2 & r3(varG[k] / dp[k]) >= ctv$minhetfreq &
                      r3(varT[k] / dp[k]) >= ctv$minhetfreq, "PASS", "Low")
  set(out, ix, "genotype", G); set(out, ix, "qual", as.integer(qual)); set(out, ix, "alt", alt); set(out, ix, "gt", gt)
  set(out, ix, "dp", as.integer(dp)); set(out, ix, "alt_reads", as.integer(altn)); set(out, ix, "filter", filt)
  out
}

# C-to-T calls (REF C, ALT T or "T,G") of one call set, with the counts behind them.
ct_calls_table <- function(chr, cpg_pos, idx, cnt, qs, call_set, ctv) {
  x <- ct_call(cnt, qs, ctv)
  k <- which(x$genotype %in% c("CT", "TT", "GT"))
  if (length(k) == 0L) return(NULL)
  mean_q <- ifelse(cnt[k, , drop = FALSE] > 0, floor(qs[k, , drop = FALSE] / pmax(cnt[k, , drop = FALSE], 1) + 0.5), 0)
  colnames(mean_q) <- q_names
  cbind(data.table(chr = chr, start = cpg_pos[idx[k]] - 1L, call_set = call_set),
        x[k, .(genotype, alt, gt, filter, qual, dp, alt_reads)],
        as.data.table(cnt[k, , drop = FALSE]), as.data.table(mean_q))
}

# =============================================================================
# Downloads (one alignment at a time, ranged, checked, resumable)
# =============================================================================
ct_curl <- function(args) {
  out <- suppressWarnings(system2("curl", args, stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); list(status = if (is.null(st)) 0L else as.integer(st), out = out)
}

ct_head <- function(url, ctv) {
  for (k in seq_len(ctv$download_retries)) {
    h <- ct_curl(c("-f", "-sS", "-I", "-A", shQuote(ctv$user_agent), "--connect-timeout", "60", "--max-time", "300", shQuote(url)))
    cl <- grep("^content-length:", tolower(h$out), value = TRUE)
    if (h$status == 0L && length(cl) > 0) {
      return(list(bytes = as.numeric(sub("^content-length:\\s*", "", trimws(tail(cl, 1)))),
                  ranges = any(grepl("^accept-ranges:\\s*bytes", tolower(h$out))),
                  last_modified = sub("^[Ll]ast-[Mm]odified:\\s*", "", trimws(tail(grep("^last-modified:", h$out, ignore.case = TRUE, value = TRUE), 1)))))
    }
    wait <- min(600, 15 * 2^(k - 1)); log_msg("HEAD %s failed (curl status %d); retry in %.0f s", basename(url), h$status, wait)
    Sys.sleep(wait)
  }
  stop("HEAD request failed repeatedly: ", url, call. = FALSE)
}

ct_free_bytes <- function(dir) {
  x <- system2("df", c("-Pk", shQuote(dir)), stdout = TRUE)
  as.numeric(strsplit(trimws(x[length(x)]), "\\s+")[[1]][4]) * 1024
}

# One byte range into its own file: curl -f, HTTP 206, Content-Range a-b/total and
# the size checked; retried with backoff; a verified range leaves a marker file.
ct_fetch_range <- function(url, a, b, total, part, ctv) {
  ok_file <- paste0(part, ".ok")
  if (file.exists(ok_file) && file.exists(part) && file.size(part) == b - a + 1) return(list(ok = TRUE, attempts = 0L))
  tmp <- paste0(part, ".tmp"); hdr <- paste0(part, ".hdr")
  for (k in seq_len(ctv$download_retries)) {
    unlink(c(tmp, hdr))
    h <- ct_curl(c("-f", "-sS", "-A", shQuote(ctv$user_agent), "--connect-timeout", "60",
                   "--max-time", as.character(ctv$download_range_timeout_s), "--speed-limit", "20000", "--speed-time", "120",
                   "-r", sprintf("%.0f-%.0f", a, b), "-D", shQuote(hdr), "-o", shQuote(tmp), "-w", shQuote("%{http_code}"), shQuote(url)))
    code <- suppressWarnings(as.integer(tail(h$out, 1)))
    hl <- if (file.exists(hdr)) readLines(hdr, warn = FALSE) else character(0)
    cr <- trimws(sub("^[^:]*:", "", tail(grep("^content-range:", hl, ignore.case = TRUE, value = TRUE), 1)))
    want <- sprintf("bytes %.0f-%.0f/%.0f", a, b, total)
    good <- h$status == 0L && identical(code, 206L) && length(cr) == 1 && identical(cr, want) &&
      file.exists(tmp) && file.size(tmp) == b - a + 1
    if (good) {
      file.rename(tmp, part); unlink(hdr)
      writeLines(c(want, unname(tools::md5sum(part))), ok_file)
      return(list(ok = TRUE, attempts = k))
    }
    why <- sprintf("curl status %d, HTTP %s, Content-Range '%s', bytes %s", h$status, code, paste(cr, collapse = ""),
                   if (file.exists(tmp)) format(file.size(tmp)) else "none")
    wait <- min(600, 15 * 2^(k - 1)) * runif(1, 0.8, 1.2)
    message(sprintf("%s range %.0f-%.0f attempt %d failed (%s); retry in %.0f s", format(Sys.time(), "%H:%M:%S"), a, b, k, why, wait))
    Sys.sleep(wait)
  }
  unlink(c(tmp, hdr))
  list(ok = FALSE, attempts = ctv$download_retries)
}

# Download one alignment (BAM and index) into dest_dir and check it: size against
# HEAD Content-Length and the manifest, BGZF end-of-file block, md5 recorded, then
# the whole file decoded once (records counted) and compared with the index.
# A verified alignment leaves <bam>.verified.rds; parts survive a restart.
ct_fetch_alignment <- function(bam_row, bai_row, dest_dir, ctv) {
  bam <- file.path(dest_dir, bam_row$file); bai <- file.path(dest_dir, bai_row$file)
  marker <- paste0(bam, ".verified.rds")
  if (file.exists(marker) && file.exists(bam) && file.exists(bai) && file.size(bam) == bam_row$bytes) {
    v <- readRDS(marker); log_msg("%s: verified alignment found in the cache", bam_row$library_id); return(v)
  }
  t0 <- Sys.time()
  hd <- ct_head(bam_row$url, ctv)
  if (!isTRUE(all.equal(hd$bytes, as.numeric(bam_row$bytes))) || !hd$ranges) {
    stop(sprintf("%s: server reports %.0f bytes (ranges %s), manifest %.0f: the file has changed or ranges are refused",
                 bam_row$file, hd$bytes, hd$ranges, bam_row$bytes), call. = FALSE)
  }
  total <- hd$bytes
  parts_dir <- paste0(bam, ".parts"); dir.create(parts_dir, showWarnings = FALSE)
  rb <- ctv$download_range_mb * 2^20
  a <- seq(0, total - 1, by = rb); b <- pmin(a + rb - 1, total - 1)
  parts <- file.path(parts_dir, sprintf("part_%04d.bin", seq_along(a)))
  have <- sum(file.exists(parts) & file.exists(paste0(parts, ".ok")))
  need <- total * 2 - have * rb + 1e9
  free <- ct_free_bytes(dest_dir)
  if (free < need) stop(sprintf("Not enough disk space for %s: %.1f GB free, %.1f GB needed", bam_row$file, free / 1e9, need / 1e9), call. = FALSE)
  log_msg("%s: downloading %.2f GB in %d ranges, %d streams (%d ranges already in the cache)", bam_row$library_id,
          total / 1e9, length(a), ctv$download_streams, have)
  res <- mclapply(seq_along(a), function(i) ct_fetch_range(bam_row$url, a[i], b[i], total, parts[i], ctv),
                  mc.cores = ctv$download_streams, mc.preschedule = FALSE)
  ok <- vapply(res, function(z) is.list(z) && isTRUE(z$ok), logical(1))
  if (!all(ok)) stop(sprintf("%s: %d of %d ranges failed after retries; rerun to resume", bam_row$file, sum(!ok), length(ok)), call. = FALSE)
  t_dl <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # concatenate the parts in order, check, then delete them
  tmp <- paste0(bam, ".tmp"); con_out <- file(tmp, "wb")
  for (p in parts) {
    con_in <- file(p, "rb")
    repeat { x <- readBin(con_in, "raw", 64 * 2^20); if (length(x) == 0L) break; writeBin(x, con_out) }
    close(con_in)
  }
  close(con_out)
  if (file.size(tmp) != total) stop(sprintf("%s: joined size %.0f, expected %.0f", bam_row$file, file.size(tmp), total), call. = FALSE)
  con <- file(tmp, "rb"); seek(con, total - 28); eof <- readBin(con, "raw", 28); close(con)
  if (!identical(eof, bgzf_eof)) stop(bam_row$file, ": the BGZF end-of-file block is missing", call. = FALSE)
  file.rename(tmp, bam)
  unlink(parts_dir, recursive = TRUE)
  # index (small, one request) after the BAM, so that it is not older than the data
  for (k in seq_len(ctv$download_retries)) {
    h <- ct_curl(c("-f", "-sS", "-A", shQuote(ctv$user_agent), "--connect-timeout", "60", "--max-time", "600",
                   "-o", shQuote(paste0(bai, ".tmp")), shQuote(bai_row$url)))
    if (h$status == 0L && file.exists(paste0(bai, ".tmp")) && file.size(paste0(bai, ".tmp")) == bai_row$bytes) break
    Sys.sleep(min(600, 15 * 2^(k - 1)))
  }
  if (!file.exists(paste0(bai, ".tmp")) || file.size(paste0(bai, ".tmp")) != bai_row$bytes) stop("Index download failed: ", bai_row$url, call. = FALSE)
  file.rename(paste0(bai, ".tmp"), bai)
  md5 <- unname(tools::md5sum(bam))
  t1 <- Sys.time()
  records <- countBam(bam)$records                      # decodes every BGZF block and record once
  ix <- idxstatsBam(bam, index = bai)
  indexed <- sum(as.numeric(ix$mapped)) + sum(as.numeric(ix$unmapped))
  if (records != indexed) stop(sprintf("%s: %.0f records decoded, %.0f in the index", bam_row$file, records, indexed), call. = FALSE)
  v <- list(library_id = bam_row$library_id, bytes = total, md5 = md5, last_modified = hd$last_modified, records = records,
            mapped = sum(as.numeric(ix$mapped)), download_s = round(t_dl), check_s = round(as.numeric(difftime(Sys.time(), t1, units = "secs"))),
            ranges = length(a), range_retries = sum(vapply(res, function(z) max(0L, z$attempts - 1L), integer(1))))
  saveRDS(v, marker)
  log_msg("%s: %.2f GB in %.0f s (%.1f MB/s), md5 %s, %.0f records decoded = index", bam_row$library_id, total / 1e9,
          t_dl, total / 1e6 / max(1, t_dl), md5, records)
  v
}

# =============================================================================
# One library, genome-wide (resumable per chromosome)
# =============================================================================
ct_count_settings <- function(ctv) {
  list(code = ct_code_version, mapvalue = ctv$mapvalue, minquali = ctv$minquali, errorate = ctv$errorate, cigar = ctv$cigar)
}

ct_chr_cpg <- function(fa, chr, len) {
  dna <- getSeq(fa, GRanges(chr, IRanges(1L, len)))[[1]]
  list(seq = charToRaw(as.character(dna)), cpg = start(matchPattern("CG", dna)))
}

ct_count_library <- function(lib, bam, bai, fasta, seqs, ctv, lib_dir) {
  dir.create(file.path(lib_dir, "counts"), recursive = TRUE, showWarnings = FALSE)
  fa <- FaFile(fasta); open(fa); on.exit(close(fa), add = TRUE)
  # Not opened: an opened BamFile returns no reads after its first scanBam query
  # (seen with Rsamtools 2.18 on local and remote files), so each query opens it.
  bf <- BamFile(bam, index = bai)
  ix <- idxstatsBam(bam, index = bai)
  in_index <- setNames(as.numeric(ix$mapped) + as.numeric(ix$unmapped), as.character(ix$seqnames))
  for (j in seq_len(nrow(seqs))) {
    chr <- seqs$chr[j]; f <- file.path(lib_dir, "counts", paste0(chr, ".rds"))
    if (file.exists(f)) next
    t0 <- Sys.time()
    cs <- ct_chr_cpg(fa, chr, seqs$len[j])
    wst <- seq(1L, seqs$len[j], by = ctv$window_bp)
    res <- vector("list", length(wst))
    for (k in seq_along(wst)) {
      ws <- wst[k]; we <- min(seqs$len[j], ws + ctv$window_bp - 1L)
      res[[k]] <- ct_count_window(ct_scan_window(bf, chr, ws, we), cs$seq, ws, we, cs$cpg, ctv, ctv$cigar)
    }
    res <- res[!vapply(res, is.null, logical(1))]
    x <- list(chr = chr, n_cpg = length(cs$cpg),
              idx = unlist(lapply(res, `[[`, "idx")),
              cnt = do.call(rbind, c(list(matrix(0L, 0, 8, dimnames = list(NULL, cell_names))), lapply(res, `[[`, "cnt"))),
              qs = do.call(rbind, c(list(matrix(0L, 0, 8, dimnames = list(NULL, q_names))), lapply(res, `[[`, "qs"))),
              stats = Reduce(`+`, lapply(res, `[[`, "stats")),
              seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    # every record of the chromosome must have been seen once (by its start position)
    if (x$stats[["records"]] != in_index[[chr]]) {
      stop(sprintf("%s %s: %.0f records read, %.0f in the index; the alignment was not read completely",
                   lib, chr, x$stats[["records"]], in_index[[chr]]), call. = FALSE)
    }
    saveRDS(x, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    if (seqs$len[j] > 5e6) log_msg("%s %s: %d CpGs, %d with reads, %.0f s", lib, chr, x$n_cpg, length(x$idx), x$seconds)
  }
}

ct_read_chr <- function(lib_dir, chr) readRDS(file.path(lib_dir, "counts", paste0(chr, ".rds")))

# =============================================================================
# Final tables (all 22 libraries counted)
# =============================================================================
ct_finalize <- function(libs, seqs, ctv, cache, out_tables, fasta) {
  log_msg("final tables: per-library calls, pooled calls (all 22 libraries; zr3534; zr3644), assessability")
  experiments <- split(libs$library_id, libs$zymo_project)
  calls <- list(); assess <- list(); per_chr <- list()
  fa <- FaFile(fasta); open(fa); on.exit(close(fa), add = TRUE)
  for (j in seq_len(nrow(seqs))) {
    chr <- seqs$chr[j]
    x1 <- ct_read_chr(file.path(cache, libs$library_id[1]), chr)
    n <- x1$n_cpg
    cs <- ct_chr_cpg(fa, chr, seqs$len[j]); stopifnot(length(cs$cpg) == n)
    pool <- lapply(c(all = "all", experiments), function(z) list(cnt = matrix(0, n, 8), qs = matrix(0, n, 8)))
    n_ass <- lapply(experiments, function(z) integer(n))
    for (l in libs$library_id) {
      x <- ct_read_chr(file.path(cache, l), chr); stopifnot(x$n_cpg == n)
      e <- libs$zymo_project[libs$library_id == l]
      for (pn in c("all", e)) {
        pool[[pn]]$cnt[x$idx, ] <- pool[[pn]]$cnt[x$idx, ] + x$cnt
        pool[[pn]]$qs[x$idx, ] <- pool[[pn]]$qs[x$idx, ] + x$qs
      }
      ok <- (x$cnt[, "c_C"] + x$cnt[, "c_T"]) >= ctv$min_crick_assessable
      n_ass[[e]][x$idx[ok]] <- n_ass[[e]][x$idx[ok]] + 1L
      calls[[length(calls) + 1]] <- ct_calls_table(chr, cs$cpg, x$idx, x$cnt, x$qs, l, ctv)
    }
    for (pn in names(pool)) {
      p <- pool[[pn]]; has <- which(rowSums(p$cnt) > 0)
      colnames(p$cnt) <- cell_names; colnames(p$qs) <- q_names
      calls[[length(calls) + 1]] <- ct_calls_table(chr, cs$cpg, has, p$cnt[has, , drop = FALSE], p$qs[has, , drop = FALSE],
                                                   if (pn == "all") "pooled_all" else paste0("pooled_", pn), ctv)
    }
    crick <- as.integer(pool$all$cnt[, 7] + pool$all$cnt[, 6])
    a <- data.table(chr = chr, gap = c(cs$cpg[1] - 1L, diff(cs$cpg)), crick_pooled = crick)
    for (e in names(n_ass)) set(a, j = paste0("n_assessable_", e), value = n_ass[[e]])
    assess[[j]] <- a
    per_chr[[j]] <- data.table(chr = chr, n_cpg = n, first_start = cs$cpg[1] - 1L, last_start = cs$cpg[n] - 1L)
    rm(pool, n_ass, cs)
  }
  calls <- rbindlist(calls)
  assess <- rbindlist(assess)
  setorder(calls, chr, start, call_set)
  # one row per CpG with a call in any call set
  union_sets <- c("pooled_all", libs$library_id)
  cp <- calls[, .(listed = any(call_set %in% union_sets),
                  pooled_all = c(genotype[call_set == "pooled_all"], "")[1],
                  pooled_all_filter = c(filter[call_set == "pooled_all"], "")[1],
                  pooled_zr3534 = c(genotype[call_set == "pooled_zr3534"], "")[1],
                  pooled_zr3644 = c(genotype[call_set == "pooled_zr3644"], "")[1],
                  n_libraries_called = sum(call_set %in% libs$library_id),
                  libraries_called = paste(call_set[call_set %in% libs$library_id], collapse = ";"),
                  any_pass = any(filter[call_set %in% union_sets] == "PASS")), by = .(chr, start)]
  # assessability of each listed CpG, joined by position
  assess[, start := cumsum(gap) , by = chr]
  cp <- merge(cp, assess[, c("chr", "start", "crick_pooled", grep("^n_assessable_", names(assess), value = TRUE)), with = FALSE],
              by = c("chr", "start"), all.x = TRUE, sort = FALSE)
  cp[, pos1 := start + 1L]
  setcolorder(cp, c("chr", "start", "pos1"))
  chr_order <- seqs$chr
  cp <- cp[order(match(chr, chr_order), start)]; calls <- calls[order(match(chr, chr_order), start, call_set)]
  assess[, start := NULL]
  write_table(cp, file.path(out_tables, "ct_variant_cpgs.csv.gz"))
  write_table(calls, file.path(out_tables, "ct_variant_calls.csv.gz"))
  write_table(assess, file.path(out_tables, "ct_variant_assessability.csv.gz"))
  list(cpgs = cp, calls = calls, assess = assess, per_chr = rbindlist(per_chr))
}

# =============================================================================
# Main
# =============================================================================
ct_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  t_start <- start_step("ct_variant_loci (optional)", "CpG sites at C-to-T variants, BS-Snper rule in R, 22 libraries")
  only <- NULL; keep_bam <- FALSE
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--libraries" && i < length(args)) { only <- strsplit(args[i + 1], ",")[[1]]; i <- i + 2 }
    else if (args[i] == "--keep-bam") { keep_bam <- TRUE; i <- i + 1 }
    else stop("Cannot read the argument: ", args[i], " (options: --libraries id1,id2  --keep-bam)", call. = FALSE)
  }
  if (!identical(ctv$cigar, "full") || !identical(ctv$site, "C") || !identical(ctv$filter, "any")) {
    stop("ct_variants: cigar must be \"full\", site \"C\" and filter \"any\" (the author's decisions of 2026-11-20)", call. = FALSE)
  }
  out_tables <- output_dir(step_rel, "tables")
  cache <- output_dir(step_rel, "cache")
  # run lock: a second run on the same cache stops; a lock left by a dead process is taken over
  lock <- file.path(cache, "run.lock")
  if (file.exists(lock)) {
    pid <- suppressWarnings(as.integer(readLines(lock, warn = FALSE)[1]))
    cmd <- if (!is.na(pid) && file.exists(sprintf("/proc/%d/cmdline", pid))) paste(readLines(sprintf("/proc/%d/cmdline", pid), warn = FALSE), collapse = " ") else ""
    if (!is.na(pid) && pid != Sys.getpid() && grepl("ct_variant_loci", cmd)) stop("Another run is active (PID ", pid, "); lock ", lock, call. = FALSE)
  }
  writeLines(c(Sys.getpid(), format(Sys.time())), lock)
  on.exit(unlink(lock), add = TRUE)

  libs <- as.data.table(read_libraries())[, .(library_id, zymo_project)]
  man <- fread(input_path("gannet_manifest_alignments.csv"))
  bam_rows <- man[grepl("\\.bam$", file)]; bai_rows <- man[grepl("\\.bai$", file)]
  if (!setequal(bam_rows$library_id, libs$library_id) || !setequal(bai_rows$library_id, libs$library_id) || nrow(bam_rows) != 22) {
    stop("gannet_manifest_alignments.csv must hold one BAM and one index row for each of the 22 libraries", call. = FALSE)
  }
  if (!is.null(only)) {
    bad <- setdiff(only, libs$library_id); if (length(bad)) stop("Unknown library: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  run_libs <- if (is.null(only)) libs$library_id else only

  fasta <- fetch_gannet_file("cgigas_uk_roslin_v1_genomic-mito.fa", config, verify_md5 = TRUE)
  genome_md5 <- read_gannet_manifest()[read_gannet_manifest()$file == basename(fasta), "md5"]
  if (!file.exists(paste0(fasta, ".fai"))) indexFa(fasta)
  si <- seqinfo(FaFile(fasta))
  seqs <- data.table(chr = seqnames(si), len = as.integer(seqlengths(si)))
  aln_dir <- gannet_cache_dir(config, "alignments")
  settings <- ct_count_settings(ctv)
  log_msg("CIGAR mode %s; rule minhetfreq %.2f minhomfreq %.2f minquali %d mincover %d minread2 %d errorate %.2f mapvalue %d",
          ctv$cigar, ctv$minhetfreq, ctv$minhomfreq, ctv$minquali, ctv$mincover, ctv$minread2, ctv$errorate, ctv$mapvalue)
  log_msg("libraries to process: %d; cache %s; alignments %s", length(run_libs), cache, aln_dir)

  for (l in run_libs) {
    lib_dir <- file.path(cache, l); dir.create(lib_dir, showWarnings = FALSE)
    done <- file.path(lib_dir, "done.rds")
    sf <- file.path(lib_dir, "count_settings.rds")
    if (file.exists(sf) && !identical(readRDS(sf), settings)) {
      log_msg("%s: cached counts were made with other settings; counting again", l)
      unlink(file.path(lib_dir, c("counts", "done.rds")), recursive = TRUE)
    }
    if (file.exists(done)) { log_msg("%s: finished earlier (cache)", l); next }
    saveRDS(settings, sf)
    t0 <- Sys.time()
    fi <- ct_fetch_alignment(bam_rows[library_id == l], bai_rows[library_id == l], aln_dir, ctv)
    bam <- file.path(aln_dir, bam_rows[library_id == l, file]); bai <- file.path(aln_dir, bai_rows[library_id == l, file])
    t1 <- Sys.time()
    ct_count_library(l, bam, bai, fasta, seqs, ctv, lib_dir)
    t2 <- Sys.time()
    st <- Reduce(`+`, lapply(seqs$chr, function(chr) ct_read_chr(lib_dir, chr)$stats))
    ncall <- 0L; ncall_pass <- 0L; nsites <- 0L
    for (chr in seqs$chr) {
      x <- ct_read_chr(lib_dir, chr); nsites <- nsites + length(x$idx)
      cl <- ct_calls_table(chr, seq_len(x$n_cpg), x$idx, x$cnt, x$qs, l, ctv)
      if (!is.null(cl)) { ncall <- ncall + nrow(cl); ncall_pass <- ncall_pass + sum(cl$filter == "PASS") }
    }
    summ <- c(fi, as.list(st), list(cpg_with_reads = nsites, ct_calls = ncall, ct_calls_pass = ncall_pass,
                                    fetch_s = round(as.numeric(difftime(t1, t0, units = "secs"))),
                                    count_s = round(as.numeric(difftime(t2, t1, units = "secs"))),
                                    max_rss_mb = ct_peak_rss_mb(), finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    saveRDS(summ, done)
    log_msg("%s: counted in %.1f min; %.0f records, %.0f used (MAPQ and mismatch filters); %d CpGs with reads; %d C-to-T calls (%d PASS); peak memory %.0f MB",
            l, summ$count_s / 60, st["records"], st["mismatch_pass"], nsites, ncall, ncall_pass, summ$max_rss_mb)
    if (!keep_bam) { unlink(c(bam, bai, paste0(bam, ".verified.rds"))); log_msg("%s: alignment deleted", l) }
  }

  finished <- libs$library_id[file.exists(file.path(cache, libs$library_id, "done.rds"))]
  if (!is.null(only) || length(finished) < 22) {
    log_msg("%d of 22 libraries finished; the final tables are written when all 22 are", length(finished))
    return(invisible(NULL))
  }
  fin <- ct_finalize(libs, seqs, ctv, cache, out_tables, fasta)
  ct_write_run_record(libs, ctv, cache, out_tables, fin, genome_md5, seqs)
  finish_step(file.path(step_rel), t_start)
}

ct_peak_rss_mb <- function() {
  s <- readLines("/proc/self/status", warn = FALSE)
  v <- grep("^VmHWM:", s, value = TRUE)
  if (length(v)) as.numeric(gsub("[^0-9]", "", v)) / 1024 else NA_real_
}

ct_write_run_record <- function(libs, ctv, cache, out_tables, fin, genome_md5, seqs) {
  lib_rows <- rbindlist(lapply(libs$library_id, function(l) {
    s <- readRDS(file.path(cache, l, "done.rds"))
    data.table(library_id = l, experiment = libs$zymo_project[libs$library_id == l], bam_bytes = s$bytes, bam_md5 = s$md5,
               bam_last_modified = s$last_modified, records = s$records, mapped = s$mapped, flag_99 = s$f99, flag_147 = s$f147,
               flag_83 = s$f83, flag_163 = s$f163, other_flags = s$other_flag, mapq_pass = s$mapq_pass,
               mismatch_pass = s$mismatch_pass, with_indel = s$with_indel, cpg_with_reads = s$cpg_with_reads,
               ct_calls = s$ct_calls, ct_calls_pass = s$ct_calls_pass, download_s = s$download_s, check_s = s$check_s,
               count_s = s$count_s, range_retries = s$range_retries, peak_memory_mb = round(s$max_rss_mb), finished = s$finished)
  }))
  par <- data.table(item = c("code_version", "cigar_mode", "site", "filter", "union", names(ctv)[1:7], "min_crick_assessable",
                             "genome", "genome_md5", "cpg_sites_screened", "bs_snper_reference"),
                    value = c(ct_code_version, ctv$cigar, ctv$site, ctv$filter, ctv$union,
                              vapply(ctv[1:7], function(v) format(v), ""), ctv$min_crick_assessable,
                              "cgigas_uk_roslin_v1_genomic-mito.fa", genome_md5, sum(fin$per_chr$n_cpg),
                              "BS-Snper 1.1 (Gao et al. 2015), rule as used by Venkataraman et al. 2022"))
  # long format: section parameter rows first, then one block of rows per library
  long_lib <- melt(lib_rows[, lapply(.SD, as.character)], id.vars = "library_id", variable.name = "item", value.name = "value")
  run <- rbind(data.table(section = "parameter", library_id = "", item = par$item, value = par$value),
               data.table(section = "library", library_id = long_lib$library_id, item = as.character(long_lib$item), value = long_lib$value))
  write_table(run, file.path(out_tables, "ct_variant_run.csv"))
  # summary log
  cp <- fin$cpgs; cl <- fin$calls; a <- fin$assess
  ne <- grep("^n_assessable_", names(a), value = TRUE)
  n_ass <- rowSums(as.matrix(a[, ne, with = FALSE]))
  sets <- c("pooled_all", "pooled_zr3534", "pooled_zr3644", libs$library_id)
  lines <- c(sprintf("ct_variant_loci: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             sprintf("Rule: BS-Snper as used by Venkataraman et al. 2022 (C of the CpG, plus strand; REF C, ALT T or T,G; any FILTER; pooled call of all 22 libraries plus each library). CIGAR: %s.", ctv$cigar),
             sprintf("CpG sites screened: %d", nrow(a)),
             sprintf("CpGs listed (pooled_all or any library): %d (%.2f%%)", sum(cp$listed), 100 * sum(cp$listed) / nrow(a)),
             sprintf("  with at least one PASS record: %d; only Low records: %d", sum(cp$listed & cp$any_pass), sum(cp$listed & !cp$any_pass)),
             sprintf("  from the pooled call: %d; only from library calls: %d; only from the pooled call: %d",
                     sum(cp$pooled_all != ""), sum(cp$listed & cp$pooled_all == ""), sum(cp$pooled_all != "" & cp$n_libraries_called == 0)),
             sprintf("CpGs called only by an experiment pool (not listed under the 2022 rule): %d", sum(!cp$listed)),
             sprintf("Assessability (evidence = Crick C + T reads at the C after the BS-Snper read filters; a library can give a CT call with at least %d):", ctv$min_crick_assessable),
             sprintf("  CpGs with fewer than %d pooled evidence reads: %d (%.2f%%)", ctv$min_crick_assessable, sum(a$crick_pooled < ctv$min_crick_assessable), 100 * mean(a$crick_pooled < ctv$min_crick_assessable)),
             sprintf("  CpGs where no library has %d evidence reads: %d (%.2f%%); all 22 libraries have them: %d (%.2f%%)", ctv$min_crick_assessable,
                     sum(n_ass == 0), 100 * mean(n_ass == 0), sum(n_ass == 22), 100 * mean(n_ass == 22)),
             "Calls by call set (C-to-T records at CpG Cs; PASS; genotypes CT/TT/GT):",
             vapply(sets, function(s) { x <- cl[call_set == s]
               sprintf("  %-14s %8d  PASS %8d  CT %8d  TT %8d  GT %6d", s, nrow(x), sum(x$filter == "PASS"), sum(x$genotype == "CT"),
                       sum(x$genotype == "TT"), sum(x$genotype == "GT")) }, ""),
             "Per library: see ct_variant_run.csv (bytes, md5, records, reads used, times).",
             sprintf("Output sizes (MB): cpgs %.1f, calls %.1f, assessability %.1f",
                     file.size(file.path(out_tables, "ct_variant_cpgs.csv.gz")) / 1e6,
                     file.size(file.path(out_tables, "ct_variant_calls.csv.gz")) / 1e6,
                     file.size(file.path(out_tables, "ct_variant_assessability.csv.gz")) / 1e6))
  writeLines(lines, file.path(out_tables, "ct_variant_log.txt"))
  cat(lines, sep = "\n")
}

if (!identical(Sys.getenv("CT_VARIANTS_FUNCTIONS_ONLY"), "1")) ct_main()

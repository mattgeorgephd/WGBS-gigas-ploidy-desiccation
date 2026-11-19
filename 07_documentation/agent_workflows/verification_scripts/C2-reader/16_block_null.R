# 16_block_null.R (scratch)
# Does it matter statistically that the surviving ploidy foreground is one accession block?
# Test: use as "genes of interest" every universe transcript in the block XM_034462509.1 to XM_034479694.1,
# whether or not it holds a DML, and run the legacy test (BP). Compare with random sets of the same size.
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C2-reader"
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
.libPaths(c("/home/user/gannet_cache/legacy_r_lib", .libPaths()))
suppressPackageStartupMessages(library(topGO)); stopifnot(as.character(packageVersion("GO.db")) == "3.13.0")
geneID2GO <- suppressWarnings(topGO::readMappings(file.path(repo, "03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz")))
# one entry per transcript here, so that margins are counted in transcripts (the legacy run counted map lines)
geneID2GO <- geneID2GO[!duplicated(names(geneID2GO))]
geneNames <- names(geneID2GO)
blk <- geneNames >= "XM_034462509.1" & geneNames <= "XM_034479694.1"
cat(sprintf("universe transcripts %d; in the block %d (%.1f%%)\n", length(geneNames), sum(blk), 100 * mean(blk)))
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
mk <- function(sel) { g <- factor(as.integer(sel), levels = c(0, 1)); names(g) <- geneNames; g }
td <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = mk(blk), annot = annFUN.gene2GO, gene2GO = geneID2GO))
pv <- function(sel) { td2 <- updateGenes(td, mk(sel)); score(suppressMessages(getSigGroups(td2, test.stat))) }
p_blk <- pv(blk)
cat(sprintf("foreground = every universe transcript of the block: BP terms tested %d; p < 0.01: %d (%.1f%%); p < 1e-6: %d\n",
            length(p_blk), sum(p_blk < 0.01), 100 * mean(p_blk < 0.01), sum(p_blk < 1e-6)))
set.seed(20261119)
rnd <- t(sapply(1:10, function(i) { p <- pv(seq_along(geneNames) %in% sample(length(geneNames), sum(blk))); c(sum(p < 0.01), sum(p < 1e-6)) }))
cat(sprintf("10 random sets of the same size: terms with p < 0.01: median %g (range %d to %d); p < 1e-6: median %g (max %d)\n",
            median(rnd[, 1]), min(rnd[, 1]), max(rnd[, 1]), median(rnd[, 2]), max(rnd[, 2])))
# overlap with the legacy ploidy result
gl <- readRDS(file.path(S, "gentable_ploidy_BP_legacy.rds")); pl <- as.numeric(sub("^< *", "", gl$classic)); names(pl) <- gl$GO.ID
ge <- readRDS(file.path(S, "gentable_ploidy_BP_exact.rds")); pe <- as.numeric(sub("^< *", "", ge$classic)); names(pe) <- ge$GO.ID
sig_l <- names(pl)[pl < 0.01]; sig_e <- names(pe)[pe < 0.01]
cat(sprintf("legacy ploidy BP terms with p < 0.01: %d; of these, block-associated (p < 0.01 for the block itself): %d (%.0f%%)\n",
            length(sig_l), sum(p_blk[sig_l] < 0.01), 100 * mean(p_blk[sig_l] < 0.01)))
cat(sprintf("complete-read ploidy BP terms with p < 0.01: %d; of these, block-associated: %d (%.0f%%)\n",
            length(sig_e), sum(p_blk[sig_e] < 0.01), 100 * mean(p_blk[sig_e] < 0.01)))
only_l <- setdiff(sig_l, sig_e)
cat(sprintf("terms with p < 0.01 under the legacy reader only: %d; block-associated: %d (%.0f%%)\n", length(only_l), sum(p_blk[only_l] < 0.01), 100 * mean(p_blk[only_l] < 0.01)))
cat(sprintf("legacy ploidy BP top ten: block p-values: %s\n", paste(signif(p_blk[gl$GO.ID[1:10]], 2), collapse = " ")))

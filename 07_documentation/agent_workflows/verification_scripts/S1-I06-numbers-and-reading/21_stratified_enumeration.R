# A valid test of the defender's argument "compare with relabellings that are balanced in the high-level libraries":
# the five high-level libraries of I-49 keep their real ploidy labels; the labels of the other 17 libraries are
# reassigned within cohort in every possible way (21 x 210 = 4,410 labellings, the observed one included).
# Statistic: CpGs whose pooled difference (all 22 libraries) exceeds 20 points. Legacy mode.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
R <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S1-I06-numbers-and-reading"
lib <- fread(file.path(R, "03_input_files/metadata/wgbs_libraries.csv"))
u <- fread(file.path(R, "04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz"))
cs <- as.matrix(u[, paste0("numCs", 1:22), with = FALSE]) + 0; ns <- cs + as.matrix(u[, paste0("numTs", 1:22), with = FALSE]); rm(u); invisible(gc(verbose = FALSE))
# only CpGs that can exceed 20 points matter for speed: keep CpGs whose range of library values is above 20 points
pm <- cs / ns; rng <- apply(pm, 1, function(x) max(x) - min(x)); keepc <- rng > 0.2; cat("CpGs whose library values span more than 20 points:", sum(keepc), "of", length(keepc), "\n")
cs <- cs[keepc, ]; ns <- ns[keepc, ]; rm(pm)
trip <- as.integer(lib$ploidy == "triploid"); high <- lib$library_id %in% c("zr3534_3", "zr3534_4", "zr3534_5", "zr3644_14", "zr3644_15")
hn <- which(lib$cohort == "heat_experiment" & !high); cn <- which(lib$cohort == "control_experiment" & !high)
cat("non-high heat libraries:", length(hn), "with", sum(trip[hn]), "triploid; non-high control libraries:", length(cn), "with", sum(trip[cn]), "triploid\n")
ch <- combn(length(hn), sum(trip[hn])); cc <- combn(length(cn), sum(trip[cn]))
L <- matrix(0L, 22, ncol(ch) * ncol(cc)); k <- 0L
for (a in seq_len(ncol(ch))) for (b in seq_len(ncol(cc))) { k <- k + 1L; z <- trip * high; z[hn[ch[, a]]] <- 1L; z[cn[cc[, b]]] <- 1L; L[, k] <- z }
stopifnot(all(colSums(L) == 11)); io <- which(colSums(L != trip) == 0); stopifnot(length(io) == 1)
rs_c <- rowSums(cs); rs_n <- rowSums(ns); n <- integer(ncol(L)); t0 <- Sys.time()
for (from in seq(1L, ncol(L), by = 60L)) { cols <- from:min(from + 59L, ncol(L)); c1 <- cs %*% L[, cols, drop = FALSE]; n1 <- ns %*% L[, cols, drop = FALSE]
  n[cols] <- as.integer(colSums(abs(100 * (c1 / n1 - (rs_c - c1) / (rs_n - n1))) > 20)) }
cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
obs <- n[io]; cat("observed:", obs, "(must be 1038)\n")
cat(sprintf("all %d labellings with the five high-level libraries fixed: median %g, 5th to 95th percentile %g to %g, min %g, max %g; labellings with at least the observed count: %d (p = %.4f)\n",
            length(n), median(n), quantile(n, .05), quantile(n, .95), min(n), max(n), sum(n >= obs), mean(n >= obs)))
# by what is reassigned: heat non-high labels real, control reassigned; and the reverse
heat_real <- colSums(L[hn, , drop = FALSE] != trip[hn]) == 0; ctrl_real <- colSums(L[cn, , drop = FALSE] != trip[cn]) == 0
cat(sprintf("heat cohort at its real labels, control non-high reassigned (%d): median %g; at least observed %d\n", sum(heat_real), median(n[heat_real]), sum(n[heat_real] >= obs)))
cat(sprintf("control cohort at its real labels, heat non-high reassigned (%d): median %g; at least observed %d; counts: %s\n", sum(ctrl_real), median(n[ctrl_real]), sum(n[ctrl_real] >= obs),
            paste(sort(n[ctrl_real]), collapse = " ")))
fwrite(data.table(labelling = seq_along(n), n_above = n, observed = seq_along(n) == io, heat_real = heat_real, ctrl_real = ctrl_real), file.path(S, "21_stratified_enumeration.csv"))

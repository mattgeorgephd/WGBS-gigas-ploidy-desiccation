suppressPackageStartupMessages({library(data.table); library(methylKit)})
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml/t05"
fit <- readRDS(file.path(scratch, "equiv_res.rds"))$legacy$fit
lib <- fread("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv"); setorder(lib, legacy_sample_index)
trt <- factor(lib$treatment, levels = c("control", "heat"))
u <- fit$united
a1 <- getData(suppressMessages(calculateDiffMeth(u, covariates = data.frame(heat = trt), mc.cores = 1)))
a2 <- getData(suppressMessages(calculateDiffMeth(u, covariates = data.frame(heat = trt), mc.cores = 2)))
b1 <- getData(suppressMessages(calculateDiffMeth(u, covariates = data.frame(HEAT = c(rep("yes", 10), rep("no", 12))), mc.cores = 1)))
cat("factor covariate, 1 core vs 2 cores: identical p", identical(a1$pvalue, a2$pvalue), "; identical q", identical(a1$qvalue, a2$qvalue), "; identical diff", identical(a1$meth.diff, a2$meth.diff), "\n")
cat("factor covariate (pipeline) vs character yes/no (legacy), 1 core: identical p", identical(a1$pvalue, b1$pvalue), "; max abs p diff", max(abs(a1$pvalue - b1$pvalue)),
    "; max rel p diff where p < 0.01:", max((abs(a1$pvalue - b1$pvalue) / b1$pvalue)[b1$pvalue < 0.01]), "; max rel q diff where q < 0.05:", max((abs(a1$qvalue - b1$qvalue) / b1$qvalue)[b1$qvalue < 0.05]), "\n")
w <- which.max(abs(a1$pvalue - b1$pvalue)); cat("largest difference at p =", b1$pvalue[w], "vs", a1$pvalue[w], "\n")

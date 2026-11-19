suppressPackageStartupMessages({library(data.table); library(methylKit)})
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/R1-coverage-dml/t05"
res <- readRDS(file.path(scratch, "equiv_res.rds"))
fit <- res$legacy$fit
lib <- fread("/home/user/WGBS-gigas-ploidy-desiccation/03_input_files/metadata/wgbs_libraries.csv"); setorder(lib, legacy_sample_index)
lib$treatment <- factor(lib$treatment, levels = c("control", "heat"))
u <- fit$united[1:300, ]
w <- NULL
r <- withCallingHandlers(suppressMessages(calculateDiffMeth(u, covariates = data.frame(heat = lib$treatment), mc.cores = 1)),
                         warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
print(table(w))

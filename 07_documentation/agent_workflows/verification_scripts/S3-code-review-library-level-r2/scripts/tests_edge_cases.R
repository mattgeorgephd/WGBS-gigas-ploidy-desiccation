# tests_edge_cases.R
# The functions of 06_diagnostics/library_methylation_level.R, taken from the parsed file (the file is only read),
# run on small synthetic inputs to look for failure scenarios. Nothing of the repository is written.
suppressPackageStartupMessages(library(data.table)); setDTthreads(1)
repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
S <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/S3-code-review-library-level-r2"
ex <- parse(file.path(repo, "06_diagnostics/library_methylation_level.R"), keep.source = FALSE)
defs <- c("above", "at_threshold", "above_largest_gap", "pooled", "difference", "relabel_reference", "round_cols")
env <- new.env()
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% defs) eval(e, env)
stopifnot(all(defs %in% ls(env)))
assign("cutoff", 20, env); assign("tie_tolerance", 1e-9, env)
say <- function(...) cat(..., "\n", sep = "")

say("== 1. above_largest_gap ==")
g <- env$above_largest_gap(c(9.46, 9.9, 10.0, 10.1, 10.2))        # the largest gap is below the second smallest value
say("largest gap at the bottom, values 9.46 9.9 10.0 10.1 10.2: high = ", paste(g$high, collapse = " "), " (", sum(g$high), " of 5 libraries marked high)")
g <- env$above_largest_gap(c(1, 2, 3, 3.5))
say("two equal largest gaps, values 1 2 3 3.5: high = ", paste(g$high, collapse = " "), "; gap ", g$gap, "; second_gap ", g$second_gap, " (equal, no warning)")
g <- env$above_largest_gap(c(1, 1, 1))
say("all values equal: high = ", paste(g$high, collapse = " "), "; gap ", g$gap)
g <- tryCatch(env$above_largest_gap(c(1, NA, 3, 10)), error = function(e) paste("error:", conditionMessage(e)))
say("a missing value (1 NA 3 10): high = ", paste(g$high, collapse = " "), " (sort() drops the NA silently)")

say("\n== 2. relabel_reference against exact integer arithmetic on random counts ==")
exact_counts <- function(m, strata, n_lib) {
  libs <- unlist(lapply(strata, function(s) c(s$a, s$b)))
  per <- lapply(strata, function(s) { x <- c(s$a, s$b); if (length(x) == 1) list(if (length(s$a)) x else integer(0)) else combn(x, length(s$a), simplify = FALSE) })
  grid <- expand.grid(lapply(per, seq_along))
  sapply(seq_len(nrow(grid)), function(i) {
    a <- unlist(lapply(seq_along(per), function(k) per[[k]][[grid[i, k]]])); b <- setdiff(libs, a)
    ca <- rowSums(m$C[, a, drop = FALSE]); na <- rowSums(m$N[, a, drop = FALSE]); cb <- rowSums(m$C[, b, drop = FALSE]); nb <- rowSums(m$N[, b, drop = FALSE])
    sum(5 * abs(cb * na - ca * nb) > na * nb)
  })
}
set.seed(1)
n_lib <- 9L; assign("n_lib", n_lib, env)
ncpg <- 20000
N <- matrix(sample(10:40, ncpg * n_lib, replace = TRUE), ncpg, n_lib); p <- matrix(rep(sample(c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1), ncpg, replace = TRUE), n_lib), ncpg, n_lib)
C <- matrix(rbinom(ncpg * n_lib, N, pmin(1, pmax(0, p + rnorm(ncpg * n_lib, 0, 0.15)))), ncpg, n_lib); storage.mode(C) <- "integer"; storage.mode(N) <- "integer"
# rows built to be exact ties: every library of group a at 60%, every library of group b at 80% (span exactly 20)
C[1:50, 1:4] <- 6L; N[1:50, 1:4] <- 10L; C[1:50, 5:8] <- 8L; N[1:50, 5:8] <- 10L
m <- list(C = C, N = N)
for (case in list(list(name = "4 against 4", strata = list(list(a = 1:4, b = 5:8))),
                  list(name = "2 against 5", strata = list(list(a = 1:2, b = 3:7))),
                  list(name = "two strata, 2+1 against 2+3", strata = list(list(a = 1:2, b = 3:4), list(a = 5L, b = 6:8))))) {
  r <- env$relabel_reference(m, case$strata)
  e <- exact_counts(m, case$strata, n_lib)
  say(case$name, ": ", length(r$counts), " groupings; counts equal to exact arithmetic in all: ", identical(as.numeric(sort(r$counts)), as.numeric(sort(e))),
      "; observed ", r$observed, "; mirror present: ", any(r$is_mirror), "; CpGs kept by the screen: ", r$n_cpg_screened, " of ", ncpg)
}
# what a rule without margin would do with the 50 built ties, with and without the screen
a <- 1:4; b <- 5:8
d_mat <- 100 * (rowSums(C[, b]) / rowSums(N[, b]) - rowSums(C[, a]) / rowSums(N[, a])); d_script <- env$difference(m, a, b)
pm <- 100 * C[, 1:8] / N[, 1:8]; span <- apply(pm, 1, max) - apply(pm, 1, min)
say("the 50 built ties (60% against 80%): value in the matrix form ", sprintf("%.17g", d_mat[1]), ", in the form of difference() ", sprintf("%.17g", d_script[1]), ", span ", sprintf("%.17g", span[1]))
say("  counted by '> 20' in the matrix form: ", sum(abs(d_mat[1:50]) > 20), "; kept by the screen 'span > 20': ", sum(span[1:50] > 20),
    " -> without the margin the screen would remove CpGs that the unscreened rule counts")
say("  counted by above() (with the margin): ", sum(env$above(d_mat[1:50])), "; at_threshold(): ", sum(env$at_threshold(d_mat[1:50])))

say("\n== 3. relabel_reference with a stratum of one library (combn(x, m) reads a single number x as 1:x) ==")
strata1 <- list(list(a = 1:2, b = 3:4), list(a = 7L, b = integer(0)))
r <- tryCatch(env$relabel_reference(m, strata1), error = function(e) paste("error:", conditionMessage(e)))
if (is.list(r)) say("strata {1,2 | 3,4} and {7 | none}: groupings enumerated ", length(r$counts), " (6 exist); group a of the first few: ",
                    paste(sapply(r$group_a[1:8], paste, collapse = ","), collapse = " ; ")) else say(r)
strata0 <- list(list(a = integer(0), b = 1:3))
r <- tryCatch(env$relabel_reference(m, strata0), error = function(e) paste("error:", conditionMessage(e)))
if (is.list(r)) say("empty group a: groupings ", length(r$counts), "; counts ", paste(r$counts, collapse = ","), "; observed ", r$observed) else say("empty group a: ", r)

say("\n== 4. round_cols ==")
x <- data.table(median_of_counts = c(4018.5, 123456.5, 1234567), p05 = c(2970.75, 10297.25, 98765.45), p = c(2 / 252, 1.23456789e-10, 0.05), n_int = c(13028L, 1234567L, 3L))
y <- env$round_cols(copy(x))
print(cbind(x, y)); say("integer column untouched: ", identical(x$n_int, y$n_int))
com <- file.path(repo, "06_diagnostics/library_methylation_level")
ref <- fread(file.path(com, "relabelling_reference.csv"))
say("committed relabelling_reference.csv: largest reference value ", max(ref$reference_p95), "; values with more than 6 significant digits before rounding are impossible to see after it; ",
    "half-integers or quantiles of counts of 100,000 or more would be altered; here the largest count is ", max(ref$n_above_threshold_observed))

say("\n== 5. integer overflow in the expected overlap ==")
o <- tryCatch(withCallingHandlers(50000L * 50000L / 65480L, warning = function(w) { say("warning: ", conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) NA)
say("sum(logical) * sum(logical) with two counts of 50,000: ", o)
agr <- fread(file.path(com, "contrast_agreement.csv"))
say("largest product of the two counts in the committed table: ", max(as.numeric(agr$n_above_in_heat_cohort) * agr$n_above_in_control_cohort), " (limit 2147483647)")

say("\n== 6. cut() with quantile breaks ==")
s <- c(rep(0, 40), rnorm(60))
r <- tryCatch(cut(s, quantile(s, 0:10 / 10), include.lowest = TRUE, labels = FALSE), error = function(e) paste("error:", conditionMessage(e)))
say("slopes with 40% exact zeros: ", if (is.character(r)) r else "no error")
s <- c(1:20) + 0
r <- cut(s, quantile(s, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
say("20 distinct values: sizes of the tenths ", paste(table(r), collapse = " "), "; a value equal to an inner break goes to the lower tenth")

say("\n== 7. part 8 on synthetic lists: one CpG in total, an empty list file, a position missing from the united table ==")
loop <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("for")) && identical(e[[2]], as.name("mode")), as.list(ex))[[1]]
body <- as.list(loop[[4]])[-1]
i_files <- which(vapply(body, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("files")), logical(1)))
part8 <- body[i_files:(i_files + 2)]                                   # files <- ; present <- ; if (!any(present)) ... else ...
run_part8 <- function(list_rows, united_keys, values) {
  tmp <- file.path(S, "out", "part8_test"); unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
  e <- new.env(parent = globalenv())
  e$corrected_lists <- names(list_rows); e$sens_dir <- tmp; e$mode <- "legacy"
  for (l in names(list_rows)) if (!is.null(list_rows[[l]])) fwrite(list_rows[[l]], file.path(tmp, sprintf("dml_%s_legacy.csv", l)))
  e$grp <- list("diploid control" = 1:2, "diploid heat" = 3:4, "triploid control" = 5:6, "triploid heat" = 7:8); e$n_lib <- 8L
  e$lib <- data.frame(library_id = paste0("L", 1:8)); e$m <- list(key = united_keys); e$pm <- values
  e$value_rows <- list(); e$apart_rows <- list()
  res <- tryCatch({ for (st in part8) eval(st, e); list(values = rbindlist(e$value_rows), apart = rbindlist(e$apart_rows)) }, error = function(err) paste("error:", conditionMessage(err)))
  unlink(tmp, recursive = TRUE); res
}
keys <- paste("chr1", c(100L, 200L, 100000L)); vals <- rbind(c(0, 0, 0, 0, 0, 0, 50, 60), c(10, 20, 30, 40, 50, 60, 70, 80), c(0, 0, 0, 0, 0, 0, 0, 0.4))
r <- run_part8(list(A = data.table(chr = "chr1", start = 100L, qvalue = 0.001)), keys, vals)
if (is.list(r)) { say("one CpG in one list: "); print(r$values); print(r$apart) } else say("one CpG in one list: ", r)
r <- run_part8(list(A = data.table(chr = "chr1", start = 100L, qvalue = 0.001), B = data.table(chr = character(0), start = integer(0), qvalue = numeric(0))), keys, vals)
if (is.list(r)) { say("a second list whose file holds a header and no row: "); print(r$apart) } else say("a second list whose file holds a header and no row: ", r)
r <- run_part8(list(A = data.table(chr = character(0), start = integer(0), qvalue = numeric(0))), keys, vals)
if (is.list(r)) { say("the only list is empty: values rows ", nrow(r$values), "; apart table: "); print(r$apart) } else say("the only list is empty: ", r)
r <- run_part8(list(A = data.table(chr = "chr1", start = 300L, qvalue = 0.001)), keys, vals)
say("a called position that is not in the united table: ", if (is.list(r)) "no error" else r)
r <- run_part8(list(A = data.table(chr = "chr1", start = 100000, qvalue = 0.001)), keys, vals)     # start written as a double: fwrite writes 1e+05
say("start = 100000 stored as a double in the list file (fwrite writes 1e+05): ", if (is.list(r)) paste("no error; rows", nrow(r$values)) else r)
r <- run_part8(list(A = data.table(chr = "chr1", start = c(100L, 200L, 100000L), qvalue = 0.001)), keys, vals)
if (is.list(r)) { say("three CpGs; the third has group 'triploid heat' at 0 and 0.4 against 0 elsewhere (ranges touch at 0): "); print(r$values[, 1:9]); print(r$apart) } else say(r)

say("\n== 7b. the only list is empty: what is written ==")
r <- run_part8(list(A = data.table(chr = character(0), start = integer(0), qvalue = numeric(0))), keys, vals)
if (is.list(r)) { print(r$values); str(as.list(r$values[1])) } else say(r)
say("vapply with a zero-length FUN.VALUE returns: "); str(vapply(c("a", "b"), function(g) logical(0), logical(0)))

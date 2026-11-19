repo <- "/home/user/WGBS-gigas-ploidy-desiccation"
source(file.path(repo, "02_R_functions/sample_sheet.R"))
tmpd <- tempfile(); dir.create(file.path(tmpd, "metadata"), recursive = TRUE)
lib0 <- read.csv(file.path(repo, "03_input_files/metadata/wgbs_libraries.csv"), stringsAsFactors = FALSE)
input_path <- function(...) file.path(tmpd, ...)
try_case <- function(label, f) {
  x <- lib0; x <- f(x); write.csv(x, file.path(tmpd, "metadata/wgbs_libraries.csv"), row.names = FALSE, na = "")
  r <- try(read_libraries(), silent = TRUE)
  if (inherits(r, "try-error")) { cat(sprintf("%-55s -> stops: %s\n", label, trimws(conditionMessage(attr(r, "condition"))))); return(invisible()) }
  s <- sample_sets(r)
  cat(sprintf("%-55s -> accepted; NA ploidy %d, NA treatment %d; set sizes all/diploid/triploid/control/heat = %s; treat_heat = %s\n", label, sum(is.na(r$ploidy)), sum(is.na(r$treatment)),
              paste(lengths(s), collapse = "/"), paste(as.integer(r$treatment == "heat"), collapse = "")))
}
try_case("as committed", identity)
try_case("treatment relabelled 'heat treated'", function(x) { x$treatment[x$treatment == "heat"] <- "heat treated"; x })
try_case("one ploidy typed 'Triploid'", function(x) { x$ploidy[22] <- "Triploid"; x })
try_case("duplicated legacy_sample_index (two rows 5)", function(x) { x$legacy_sample_index[6] <- 5L; x })
try_case("rows shuffled", function(x) x[sample(22), ])
try_case("23 rows", function(x) rbind(x, transform(x[1, ], library_id = "zr3534_11", legacy_sample_index = 23L)))

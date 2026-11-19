# legacy_reader_mechanism.R
#
# Purpose   How the legacy enrichment script lost most of the ploidy annotation
#           table and a quarter of the heat table when it read them, and what
#           the loss did to the test (issue I-40, with the universe gap of
#           I-39). The legacy script read three tables (one line per DML,
#           transcript and GO term) with read.table() and its default quote
#           characters, the double quote and the apostrophe. GO term names hold
#           apostrophes ("mRNA 3'-UTR binding"). An apostrophe without a
#           partner on its line opens a quoted string, and every line up to the
#           next such apostrophe goes into one field of one row. Step 09
#           repeats that read in order to reproduce the draft. This script
#           documents it for the three legacy tables (ploidy, heat, shared):
#             (a) lines, transcripts and genes in the file and in the rows of
#                 four reads (default quote characters, no quote character,
#                 the double quote alone, the apostrophe alone), with the
#                 warnings R gives for each read;
#             (b) the lines with an unpaired apostrophe, the merged rows, and
#                 two interventions: with the apostrophes replaced on those
#                 lines alone the default read returns every line, and with
#                 the apostrophes replaced on all other lines it returns the
#                 rows it returned before;
#             (c) what is lost: the transcripts and genes that are in no
#                 returned row, and the share of the skipped lines that are
#                 theirs;
#             (d) where the survivors are: the longest run of surviving
#                 transcripts in accession order, and surviving and lost
#                 transcripts and genes by linkage group;
#             (e) the overlap with the universe gap of I-39: transcripts and
#                 genes admitted by the universe alone, by the reader alone,
#                 and by both;
#             (f) survival of a gene by the number of its transcripts in the
#                 table, and of a transcript by whether it carries a term with
#                 an unpaired apostrophe;
#             (g) a position-only control for the BP test of the ploidy table.
#                 Nearly all transcripts of interest of that test lie in one
#                 accession run. To see what position alone does to the test,
#                 every universe transcript of the run is taken as the list of
#                 interest, with no DML information, and tested as step 09
#                 tests (transcripts as genes, the legacy universe with its
#                 repeated names, classic Fisher test of topGO, GO.db 3.13.0).
#                 The number of BP terms with p below the cutoff is compared
#                 with that of random lists of the same size from the
#                 universe. That list is many times as long as the legacy
#                 list, so random lists of the size of the legacy list are
#                 also drawn, from the run and from the whole universe. Every
#                 list is tested twice: with the totals of the test counted in
#                 lines of the universe map, as in step 09 (issue I-43), and
#                 counted in transcripts.
#
# Run       Rscript 06_diagnostics/legacy_reader_mechanism.R      (about 3 minutes)
#
# Reads     03_input_files/annotation/legacy/DML-Cov10-20-BOTH-{PLOIDY,HEAT,COMBO}.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz
#           03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz
#           04_output/08_gene_annotation/tables/transcript_table.csv.gz
#           04_output/09_go_enrichment_legacy/tables/legacy_reader_loss.csv, legacy_funnel.csv, legacy_fisher_ploidy_BP.csv
#           GO.db 3.13.0 (side library in the cache, as step 09)
# Writes    06_diagnostics/legacy_reader_mechanism.csv            long format: table, quantity, value, definition
#           06_diagnostics/legacy_reader_by_linkage_group.csv     one row per table and sequence
#           06_diagnostics/legacy_reader_position_control.csv     one row per kind of list and unit of the totals
#
# Notes     - The three tables are the legacy files, which the legacy script
#             read. Step 09 reads copies that step 08 regenerates. The script
#             stops unless its row, transcript and gene counts equal those of
#             legacy_reader_loss.csv and legacy_funnel.csv, and its p-values
#             for the legacy list those of legacy_fisher_ploidy_BP.csv.
#           - Which line starts which row is predicted from the quote rule of
#             scan(), which read.table() calls: in a file with a separator a
#             quote character opens a quoted string wherever in a field it
#             stands, only the same character closes it, and two of them in a
#             row stand for the character itself (R sources, src/main/scan.c,
#             fillBuffer()). The script stops unless each of the four reads
#             returns exactly the predicted rows.
#           - A transcript survives when it is named in the transcript column
#             of a row of the default read, a gene when one of its transcripts
#             survives. A skipped line is a line that does not start a row.
#           - Accession order is the order of the number in the transcript
#             name (XM_<number>.<version>). It does not depend on the collation
#             of the R session.
#           - The linkage groups are also reported as two halves in the order
#             of config$linkage_groups, because the annotated mRNAs in the
#             accession range of the ploidy run lie on the first five and on
#             unplaced scaffolds (quantity mrnas_in_run_range_by_sequence).
#           - Position control: the two fixed lists are tested by topGO itself.
#             The random lists are tested with the same hypergeometric tail
#             computed from the term members of the topGO object, which the
#             script shows to be identical to the p-values of topGO for the
#             two fixed lists, in both units. 200 lists of each of the three
#             random kinds (600 draws), set.seed(20261119); the same lists are
#             used for both units.
#           - Warning texts are part of the output. They are those of the R
#             version named in the first row of legacy_reader_mechanism.csv.
#           - Legacy mode only: the tables and the test belong to the legacy
#             analysis.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
go_db_version <- use_legacy_go_db(config)            # must come before topGO is loaded
load_packages(c("data.table", "topGO"))
setDTthreads(1)                                      # the tables are small
Sys.setenv(LANGUAGE = "en")                          # warnings are quoted in the output: keep them in English

n_draws <- 200L                                      # random lists of each kind in the position control
seed <- 20261119L
tables <- c(ploidy = "PLOIDY", heat = "HEAT", shared = "COMBO")   # the name of each table in its legacy file name
quote_settings <- list(                              # the four reads of part (a): the quote characters of each
  default_quotes     = c("\"", "'"),
  no_quote_character = character(0),
  double_quote_only  = "\"",
  apostrophe_only    = "'")
replacement <- "^"                                   # stands in for an apostrophe in the interventions
out_file <- function(name) here("06_diagnostics", name)
step09 <- function(name) output_path("09_go_enrichment_legacy", "tables", name)
log_msg("%s; GO.db %s (ontology of %s); topGO %s", R.version.string, go_db_version, go_source_date(),
        as.character(packageVersion("topGO")))
# The prediction below takes the default quote characters of read.table() to be these two.
stopifnot(identical(formals(utils::read.table)$quote, paste(quote_settings$default_quotes, collapse = "")))

# ==========================================================================================
# Functions
# ==========================================================================================

# Evaluate an expression; return its value and the distinct warnings it gave.
with_warnings <- function(expr) {
  messages <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(messages))
}

# read.table() as the legacy script called it, with the quote characters as the
# one argument that varies. `from` is a file path or a connection.
read_with_quotes <- function(from, quote_chars) {
  as.data.table(utils::read.table(from, header = FALSE, sep = "\t", col.names = dml_annotation_columns(),
                                  quote = paste(quote_chars, collapse = "")))
}
# The same read of lines held in memory (for the interventions).
read_lines_with_quotes <- function(text_lines, quote_chars) {
  con <- textConnection(text_lines)
  on.exit(close(con))
  read_with_quotes(con, quote_chars)
}

# State at the end of one line under the quote rule of scan() (see Notes).
#   state     "" outside a quoted string, else the quote character that is open
#   quote_at  positions of the quote characters of the line
#   quote_is  the quote characters themselves
state_after_line <- function(state, quote_at, quote_is) {
  i <- 1L
  n <- length(quote_is)
  while (i <= n) {
    if (state == "") {
      state <- quote_is[i]                           # opens a quoted string
    } else if (quote_is[i] == state) {
      doubled <- i < n && quote_is[i + 1L] == state && quote_at[i + 1L] == quote_at[i] + 1L
      if (doubled) i <- i + 1L else state <- ""      # two in a row are a literal character; one closes the string
    }                                                # the other quote character inside a string is a literal character
    i <- i + 1L
  }
  state
}

# The rows scan() makes of the lines of a file: a row ends at the end of a line
# that is not inside a quoted string, and at the end of the file. Returns the
# first and the last line of every row, whether the file ends inside a quoted
# string, and for every line the state it would end in if it were read alone.
predict_rows <- function(text_lines, quote_chars) {
  n <- length(text_lines)
  state_after <- rep("", n)
  alone <- rep("", n)
  if (length(quote_chars) > 0) {
    pattern <- paste0("[", paste(quote_chars, collapse = ""), "]")
    with_quote <- grep(pattern, text_lines)
    quote_at <- gregexpr(pattern, text_lines[with_quote])
    quote_is <- regmatches(text_lines[with_quote], quote_at)
    state <- ""
    after <- character(length(with_quote))
    for (j in seq_along(with_quote)) {
      state <- state_after_line(state, quote_at[[j]], quote_is[[j]])
      after[j] <- state
      alone[with_quote[j]] <- state_after_line("", quote_at[[j]], quote_is[[j]])
    }
    # A line without a quote character ends in the state of the last line that had one.
    state_after <- c("", after)[findInterval(seq_len(n), with_quote) + 1L]
  }
  ends_row <- state_after == ""
  ends_inside_string <- !ends_row[n]
  ends_row[n] <- TRUE
  last_line <- which(ends_row)
  list(first_line = c(1L, head(last_line, -1L) + 1L), last_line = last_line,
       ends_inside_string = ends_inside_string, alone = alone)
}

# Stop unless read.table() returned the predicted rows: the first seven fields
# of a row are those of its first line, the last two those of its last line,
# and the rows that hold a line break are the rows predicted to hold more than
# one line. `complete` is the table with every line as a row.
check_predicted_rows <- function(returned, predicted, complete, label) {
  first <- predicted$first_line
  last <- predicted$last_line
  ok <- nrow(returned) == length(first)
  if (ok) {
    closed <- rep(TRUE, length(last))                # a row that the end of the file cuts short has no last fields
    holds_break <- last > first
    if (predicted$ends_inside_string) {
      closed[length(last)] <- FALSE
      holds_break[length(last)] <- TRUE              # the line break that ends the file is inside the string
    }
    ok <- identical(returned$transcript, complete$transcript[first]) &&
      identical(returned$geneID, complete$geneID[first]) &&
      identical(returned$start, complete$start[first]) &&
      identical(returned$GOID, complete$GOID[first]) &&
      identical(returned$GOSlim[closed], complete$GOSlim[last][closed]) &&
      identical(returned$GOcat[closed], complete$GOcat[last][closed]) &&
      identical(grepl("\n", returned$GOterm, fixed = TRUE), holds_break)
  }
  if (!ok) stop("read.table() did not return the rows predicted from the quote rule: ", label, call. = FALSE)
  invisible(TRUE)
}

# The number in a transcript name, XM_<nine digits>.<version>.
accession_number <- function(id) {
  if (!all(grepl("^XM_[0-9]{9}\\.[0-9]+$", id))) stop("A transcript name is not XM_<nine digits>.<version>", call. = FALSE)
  as.integer(substr(id, 4L, 12L))
}

# One row of a step 09 table, or stop.
one_row <- function(d, what) {
  if (nrow(d) != 1) stop("Expected exactly one row of ", what, call. = FALSE)
  d
}

count_text <- function(x) sprintf("%d", as.integer(x))
percent_text <- function(part, whole) if (whole == 0) NA_character_ else sprintf("%.1f", 100 * part / whole)

# ==========================================================================================
# 1. Definitions of the quantities of legacy_reader_mechanism.csv
# ==========================================================================================
lg_first <- config$linkage_groups[1:5]
lg_last  <- config$linkage_groups[6:10]
lg_first_text <- sprintf("the first five linkage groups (%s to %s)", lg_first[1], lg_first[5])
lg_last_text  <- sprintf("the last five linkage groups (%s to %s)", lg_last[1], lg_last[5])
read_text <- c(
  default_quotes     = "the default quote characters of read.table(), the double quote and the apostrophe, as in the legacy script",
  no_quote_character = "no quote character (quote = \"\")",
  double_quote_only  = "the double quote as the only quote character",
  apostrophe_only    = "the apostrophe as the only quote character")
definitions <- c(
  r_version = "R version that made the reads and gave the warnings reported in this table",
  universe_map_lines = "Lines of the legacy universe map (one line per blastx alignment of a transcript)",
  universe_transcripts = "Distinct transcripts of the legacy universe",
  universe_transcripts_first_five_linkage_groups = paste("Universe transcripts on", lg_first_text),
  percent_universe_first_five_linkage_groups = paste("Percent of the universe transcripts on", lg_first_text),
  universe_transcripts_last_five_linkage_groups = paste("Universe transcripts on", lg_last_text),
  percent_universe_last_five_linkage_groups = paste("Percent of the universe transcripts on", lg_last_text),
  universe_transcripts_unplaced_scaffolds = "Universe transcripts on unplaced scaffolds",
  lines_in_file = "Lines of the annotation table (one line per DML, transcript and GO term)",
  transcripts_in_file = "Distinct transcripts named in the table",
  genes_in_file = "Distinct genes named in the table",
  percent_rows_default_quotes = "Rows of the default read as a percent of the lines of the table",
  lines_with_apostrophe = "Lines that hold at least one apostrophe (all apostrophes are in GO term names)",
  lines_with_unpaired_apostrophe = paste("Lines with an odd number of apostrophes outside a double-quoted field.",
                                         "Read alone with the default quote characters, such a line ends inside a quoted string"),
  merged_rows_default_quotes = "Rows of the default read that hold more than one line of the table (a line break inside a field)",
  lines_in_largest_merged_row = "Lines of the table in the largest merged row of the default read",
  table_ends_inside_quoted_string = paste("Whether the last quoted string of the default read is still open at the end of the table",
                                          "(yes or no); R then warns of an end of file within a quoted string"),
  skipped_lines = "Lines that do not start a row of the default read: lines of the table minus rows returned",
  rows_after_replacing_unpaired_apostrophes = paste("Rows of the default read after every apostrophe on the lines with an unpaired",
                                                    "apostrophe, and on no other line, was replaced by a caret"),
  repaired_read_returns_every_line = paste("Whether that read returns one row per line, with the transcript and the GO ID of the line,",
                                           "and gives no warning (yes or no)"),
  rows_after_replacing_all_other_apostrophes = paste("Rows of the default read after every apostrophe on all other lines was replaced",
                                                     "by a caret and the lines with an unpaired apostrophe were left as they are"),
  control_read_returns_the_same_rows = "Whether that read returns rows that start on the same lines as those of the default read (yes or no)",
  transcripts_lost = "Transcripts of the table that are in no row of the default read",
  genes_lost = "Genes of the table that are in no row of the default read",
  skipped_lines_of_lost_transcripts = "Skipped lines that belong to transcripts in no row of the default read",
  percent_skipped_lines_of_lost_transcripts = "Those lines as a percent of the skipped lines",
  skipped_lines_of_lost_genes = "Skipped lines that belong to genes in no row of the default read",
  percent_skipped_lines_of_lost_genes = "Those lines as a percent of the skipped lines",
  runs_of_surviving_transcripts = "Runs of surviving transcripts when the transcripts of the table are put in accession order",
  longest_surviving_run_transcripts = "Transcripts in the longest such run",
  longest_surviving_run_first = "First transcript of the longest run",
  longest_surviving_run_last = "Last transcript of the longest run",
  percent_surviving_transcripts_in_longest_run = "Transcripts of the longest run as a percent of the surviving transcripts",
  universe_transcripts_in_run_range = "Universe transcripts whose accession lies between the first and the last transcript of the longest run",
  mrnas_in_run_range_by_sequence = paste("Annotated XM mRNAs (transcript table of step 08) whose accession lies in that range,",
                                         "counted by sequence; sequences with none are left out"),
  tested_transcripts_first_five_linkage_groups = paste("Transcripts that survive and are in the legacy universe, on", lg_first_text),
  percent_tested_transcripts_first_five_linkage_groups = "Those transcripts as a percent of all that survive and are in the universe",
  tested_transcripts_last_five_linkage_groups = paste("Transcripts that survive and are in the legacy universe, on", lg_last_text),
  percent_tested_transcripts_last_five_linkage_groups = "Those transcripts as a percent of all that survive and are in the universe",
  tested_transcripts_unplaced_scaffolds = "Transcripts that survive and are in the legacy universe, on unplaced scaffolds",
  transcripts_admitted_by_universe_alone = "Transcripts of the table that are in the legacy universe: admitted if the universe were the only filter",
  transcripts_admitted_by_reader_alone = "Transcripts of the table that survive the default read: admitted if the reader were the only filter",
  transcripts_admitted_by_both = paste("Transcripts that survive the default read and are in the legacy universe:",
                                       "the transcripts of interest of the legacy test"),
  genes_admitted_by_universe_alone = "Genes with at least one transcript of the table in the legacy universe",
  genes_admitted_by_reader_alone = "Genes with at least one transcript that survives the default read",
  genes_admitted_by_both = "Genes with at least one transcript that survives the default read and is in the legacy universe",
  percent_surviving_transcripts_in_universe = "Percent of the surviving transcripts that are in the legacy universe",
  percent_lost_transcripts_in_universe = "Percent of the lost transcripts that are in the legacy universe (empty when no transcript is lost)",
  genes_with_one_transcript = "Genes with one transcript in the table",
  genes_with_one_transcript_surviving = "Of those genes, the ones that survive the default read",
  percent_genes_with_one_transcript_surviving = "Percent of the genes with one transcript in the table that survive",
  genes_with_two_to_five_transcripts = "Genes with two to five transcripts in the table",
  genes_with_two_to_five_transcripts_surviving = "Of those genes, the ones that survive the default read",
  percent_genes_with_two_to_five_transcripts_surviving = "Percent of the genes with two to five transcripts in the table that survive",
  genes_with_more_than_five_transcripts = "Genes with more than five transcripts in the table",
  genes_with_more_than_five_transcripts_surviving = "Of those genes, the ones that survive the default read",
  percent_genes_with_more_than_five_transcripts_surviving = "Percent of the genes with more than five transcripts in the table that survive",
  transcripts_with_unpaired_apostrophe_term = "Transcripts with at least one line that holds an unpaired apostrophe",
  transcripts_with_unpaired_apostrophe_term_surviving = "Of those transcripts, the ones that survive the default read",
  percent_other_transcripts_surviving = "Percent of all other transcripts of the table that survive the default read")
for (r in names(quote_settings)) {
  definitions[paste0("rows_", r)] <- sprintf("Rows returned by read.table(header = FALSE, sep = tab) with %s", read_text[[r]])
  definitions[paste0("transcripts_", r)] <- "Distinct transcripts in the transcript column of that read"
  definitions[paste0("genes_", r)] <- "Distinct genes in the gene column of that read"
  definitions[paste0("warnings_", r)] <- "Warnings R gave for that read, distinct messages joined by ' | '; 'none' when R gave none"
}
stopifnot(!anyDuplicated(names(definitions)))

mechanism <- list()
add <- function(table_name, quantity, value) {
  if (!quantity %in% names(definitions)) stop("No definition for quantity ", quantity, call. = FALSE)
  if (length(value) != 1) stop("Quantity ", quantity, " has not one value", call. = FALSE)
  mechanism[[length(mechanism) + 1L]] <<- data.table(table = table_name, quantity = quantity, value = as.character(value),
                                                    definition = definitions[[quantity]])
  invisible(NULL)
}
yes_no <- function(x) if (isTRUE(x)) "yes" else "no"
add("session", "r_version", R.version.string)

# ==========================================================================================
# 2. The legacy universe and the transcript table
# ==========================================================================================
gene2go <- topGO::readMappings(input_path("annotation", "geneid2go-union1x-BOTH.tab.gz"))
map_names <- names(gene2go)                          # one name per line of the map; names repeat, as in the legacy run
universe <- unique(map_names)
universe_accession <- accession_number(universe)
# Every line of a transcript carries the same GO IDs, so one entry per transcript is well defined.
go_text <- vapply(gene2go, paste, character(1), collapse = ",")
stopifnot(uniqueN(data.table(map_names, go_text)) == length(universe), !anyDuplicated(universe_accession))

tx_table <- fread(output_path("08_gene_annotation", "tables", "transcript_table.csv.gz"), colClasses = c(gene_id = "character"))
stopifnot(!anyDuplicated(tx_table$transcript_id), !any(tx_table$chr == config$mito_sequence),
          setequal(universe, tx_table[in_legacy_universe == TRUE, transcript_id]),
          tx_table[, uniqueN(chr), by = gene_id][, all(V1 == 1L)])
sequence_levels <- c(config$linkage_groups, "unplaced scaffolds")
sequence_of <- function(chr_names) ifelse(chr_names %in% config$linkage_groups, chr_names, "unplaced scaffolds")
count_by_sequence <- function(sequence_names) as.integer(table(factor(sequence_names, levels = sequence_levels)))
xm <- tx_table[grepl("^XM_", transcript_id), .(transcript_id, chr)]          # the curated NM mRNAs have their own numbers
xm[, `:=`(accession = accession_number(transcript_id), sequence = sequence_of(chr))]
universe_sequence <- sequence_of(tx_table$chr[match(universe, tx_table$transcript_id)])
universe_by_sequence <- count_by_sequence(universe_sequence)

add("universe", "universe_map_lines", count_text(length(map_names)))
add("universe", "universe_transcripts", count_text(length(universe)))
add("universe", "universe_transcripts_first_five_linkage_groups", count_text(sum(universe_sequence %in% lg_first)))
add("universe", "percent_universe_first_five_linkage_groups", percent_text(sum(universe_sequence %in% lg_first), length(universe)))
add("universe", "universe_transcripts_last_five_linkage_groups", count_text(sum(universe_sequence %in% lg_last)))
add("universe", "percent_universe_last_five_linkage_groups", percent_text(sum(universe_sequence %in% lg_last), length(universe)))
add("universe", "universe_transcripts_unplaced_scaffolds", count_text(sum(universe_sequence == "unplaced scaffolds")))

# ==========================================================================================
# 3. The three tables: reads, mechanism, loss (parts a to f)
# ==========================================================================================
step09_loss <- fread(step09("legacy_reader_loss.csv"))
step09_funnel <- fread(step09("legacy_funnel.csv"))
by_sequence <- list()
survivors <- list()                                  # per table: the transcripts of the default read and the longest run
for (tbl in names(tables)) {
  path <- input_path("annotation", "legacy",
                     sprintf("DML-Cov10-20-BOTH-%s.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz", tables[[tbl]]))
  file_lines <- readLines(path)
  n_file_lines <- length(file_lines)
  complete <- read_dml_annotation(path, "exact")     # every line a row
  fields <- tstrsplit(paste0(file_lines, "\t"), "\t", fixed = TRUE)    # the added tab keeps a last field that is empty
  # Input as the description above assumes it: ten fields on every line, no comment character, no carriage return,
  # the replacement character absent, quote characters in GO term names only, lines sorted by transcript accession.
  stopifnot(n_file_lines > 0, length(fields) == length(dml_annotation_columns()), !anyNA(fields[[10]]),
            nrow(complete) == n_file_lines, identical(complete$transcript, fields[[1]]), identical(complete$GOterm, fields[[8]]),
            !any(grepl("#", file_lines, fixed = TRUE)), !any(grepl("\r", file_lines, fixed = TRUE)),
            !any(grepl(replacement, file_lines, fixed = TRUE)),
            !any(grepl("[\"']", unlist(fields[-8]))),
            !is.unsorted(accession_number(complete$transcript)))

  # ---- (a) four reads --------------------------------------------------------------------
  add(tbl, "lines_in_file", count_text(n_file_lines))
  add(tbl, "transcripts_in_file", count_text(uniqueN(complete$transcript)))
  add(tbl, "genes_in_file", count_text(uniqueN(complete$geneID)))
  reads <- list(); predicted <- list(); read_warnings <- list()
  for (r in names(quote_settings)) {
    read <- with_warnings(if (r == "default_quotes") read_dml_annotation(path, "legacy") else read_with_quotes(path, quote_settings[[r]]))
    reads[[r]] <- read$value
    read_warnings[[r]] <- if (length(read$warnings) > 0) paste(read$warnings, collapse = " | ") else "none"
    predicted[[r]] <- predict_rows(file_lines, quote_settings[[r]])
    check_predicted_rows(reads[[r]], predicted[[r]], complete, paste(tbl, r))
    add(tbl, paste0("rows_", r), count_text(nrow(reads[[r]])))
    add(tbl, paste0("transcripts_", r), count_text(uniqueN(reads[[r]]$transcript)))
    add(tbl, paste0("genes_", r), count_text(uniqueN(reads[[r]]$geneID)))
    add(tbl, paste0("warnings_", r), read_warnings[[r]])
  }
  default_read <- reads$default_quotes
  rows_default <- predicted$default_quotes
  add(tbl, "percent_rows_default_quotes", percent_text(nrow(default_read), n_file_lines))
  # Step 09 made the same two reads of the copies that step 08 regenerates.
  s9_complete <- one_row(step09_loss[set == tbl & grepl("^every line", reader)], "legacy_reader_loss.csv")
  s9_default  <- one_row(step09_loss[set == tbl & grepl("^read.table defaults", reader)], "legacy_reader_loss.csv")
  stopifnot(s9_complete$rows == n_file_lines, s9_complete$transcripts == uniqueN(complete$transcript),
            s9_complete$genes == uniqueN(complete$geneID), s9_default$rows == nrow(default_read),
            s9_default$transcripts == uniqueN(default_read$transcript), s9_default$genes == uniqueN(default_read$geneID))

  # ---- (b) unpaired apostrophes, merged rows, interventions ----------------------------------
  # A line with an unpaired apostrophe is one that, read alone, ends inside a string that an apostrophe opened.
  # In these tables that is the same as an odd number of apostrophes outside double-quoted fields, and no line
  # read alone ends inside a string that a double quote opened.
  has_apostrophe <- grepl("'", file_lines, fixed = TRUE)
  unpaired <- rows_default$alone == "'"
  outside_double_quotes <- gsub("\"[^\"]*\"", "", file_lines)
  apostrophes_outside <- nchar(outside_double_quotes) - nchar(gsub("'", "", outside_double_quotes, fixed = TRUE))
  stopifnot(identical(unpaired, apostrophes_outside %% 2L == 1L), !any(rows_default$alone == "\""))
  lines_per_row <- rows_default$last_line - rows_default$first_line + 1L
  add(tbl, "lines_with_apostrophe", count_text(sum(has_apostrophe)))
  add(tbl, "lines_with_unpaired_apostrophe", count_text(sum(unpaired)))
  add(tbl, "merged_rows_default_quotes", count_text(sum(grepl("\n", default_read$GOterm, fixed = TRUE))))
  add(tbl, "lines_in_largest_merged_row", count_text(max(lines_per_row)))
  add(tbl, "table_ends_inside_quoted_string", yes_no(rows_default$ends_inside_string))
  add(tbl, "skipped_lines", count_text(n_file_lines - nrow(default_read)))

  # The interventions are read from memory. That route returns what the read of the file returns.
  stopifnot(identical(as.data.frame(suppressWarnings(read_lines_with_quotes(file_lines, quote_settings$default_quotes))),
                      as.data.frame(default_read)))
  repaired <- file_lines
  repaired[unpaired] <- gsub("'", replacement, file_lines[unpaired], fixed = TRUE)
  repaired_read <- with_warnings(read_lines_with_quotes(repaired, quote_settings$default_quotes))
  add(tbl, "rows_after_replacing_unpaired_apostrophes", count_text(nrow(repaired_read$value)))
  add(tbl, "repaired_read_returns_every_line",
      yes_no(nrow(repaired_read$value) == n_file_lines && length(repaired_read$warnings) == 0 &&
               identical(repaired_read$value$transcript, complete$transcript) && identical(repaired_read$value$GOID, complete$GOID)))
  others <- has_apostrophe & !unpaired
  control <- file_lines
  control[others] <- gsub("'", replacement, file_lines[others], fixed = TRUE)
  control_read <- suppressWarnings(read_lines_with_quotes(control, quote_settings$default_quotes))
  add(tbl, "rows_after_replacing_all_other_apostrophes", count_text(nrow(control_read)))
  add(tbl, "control_read_returns_the_same_rows",
      yes_no(nrow(control_read) == nrow(default_read) && identical(control_read$transcript, default_read$transcript) &&
               identical(control_read$GOID, default_read$GOID) && identical(control_read$start, default_read$start)))

  # ---- (c) what is lost ------------------------------------------------------------------------
  complete[, `:=`(starts_row = FALSE, unpaired_line = unpaired)]
  complete[rows_default$first_line, starts_row := TRUE]
  per_tx <- complete[, .(gene_id = as.character(geneID[1L]), genes_named = uniqueN(geneID), chr = chr[1L],
                         sequences_named = uniqueN(chr), n_lines = .N, n_rows = sum(starts_row),
                         has_unpaired_term = any(unpaired_line)), by = transcript]
  in_table <- match(per_tx$transcript, tx_table$transcript_id)
  stopifnot(all(per_tx$genes_named == 1L), all(per_tx$sequences_named == 1L), !anyNA(in_table),
            identical(tx_table$gene_id[in_table], per_tx$gene_id), identical(tx_table$chr[in_table], per_tx$chr))
  per_tx[, `:=`(survives = n_rows > 0L, in_universe = transcript %chin% universe,
                accession = accession_number(transcript), sequence = sequence_of(chr))]
  stopifnot(setequal(per_tx[survives == TRUE, transcript], default_read$transcript), !anyDuplicated(per_tx$accession))
  per_gene <- per_tx[, .(sequence = sequence[1L], n_transcripts = .N, n_lines = sum(n_lines), survives = any(survives),
                         in_universe = any(in_universe), tested = any(survives & in_universe)), by = gene_id]
  stopifnot(setequal(per_gene[survives == TRUE, gene_id], as.character(default_read$geneID)))
  n_skipped <- n_file_lines - nrow(default_read)
  add(tbl, "transcripts_lost", count_text(per_tx[survives == FALSE, .N]))
  add(tbl, "genes_lost", count_text(per_gene[survives == FALSE, .N]))
  add(tbl, "skipped_lines_of_lost_transcripts", count_text(per_tx[survives == FALSE, sum(n_lines)]))
  add(tbl, "percent_skipped_lines_of_lost_transcripts", percent_text(per_tx[survives == FALSE, sum(n_lines)], n_skipped))
  add(tbl, "skipped_lines_of_lost_genes", count_text(per_gene[survives == FALSE, sum(n_lines)]))
  add(tbl, "percent_skipped_lines_of_lost_genes", percent_text(per_gene[survives == FALSE, sum(n_lines)], n_skipped))

  # ---- (d) where the survivors are -------------------------------------------------------------
  setorder(per_tx, accession)
  runs <- rle(per_tx$survives)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1L
  surviving_runs <- which(runs$values)
  longest <- surviving_runs[which.max(runs$lengths[surviving_runs])]
  run_first <- per_tx$transcript[run_start[longest]]
  run_last <- per_tx$transcript[run_end[longest]]
  run_range <- c(per_tx$accession[run_start[longest]], per_tx$accession[run_end[longest]])
  universe_in_range <- universe_accession >= run_range[1] & universe_accession <= run_range[2]
  mrnas_in_range <- count_by_sequence(xm[accession >= run_range[1] & accession <= run_range[2], sequence])
  n_surviving <- per_tx[survives == TRUE, .N]
  add(tbl, "runs_of_surviving_transcripts", count_text(length(surviving_runs)))
  add(tbl, "longest_surviving_run_transcripts", count_text(runs$lengths[longest]))
  add(tbl, "longest_surviving_run_first", run_first)
  add(tbl, "longest_surviving_run_last", run_last)
  add(tbl, "percent_surviving_transcripts_in_longest_run", percent_text(runs$lengths[longest], n_surviving))
  add(tbl, "universe_transcripts_in_run_range", count_text(sum(universe_in_range)))
  add(tbl, "mrnas_in_run_range_by_sequence",
      paste(sprintf("%s: %d", sequence_levels[mrnas_in_range > 0], mrnas_in_range[mrnas_in_range > 0]), collapse = "; "))
  tested_by_sequence <- count_by_sequence(per_tx[survives == TRUE & in_universe == TRUE, sequence])
  n_tested <- sum(tested_by_sequence)
  add(tbl, "tested_transcripts_first_five_linkage_groups", count_text(sum(tested_by_sequence[sequence_levels %in% lg_first])))
  add(tbl, "percent_tested_transcripts_first_five_linkage_groups",
      percent_text(sum(tested_by_sequence[sequence_levels %in% lg_first]), n_tested))
  add(tbl, "tested_transcripts_last_five_linkage_groups", count_text(sum(tested_by_sequence[sequence_levels %in% lg_last])))
  add(tbl, "percent_tested_transcripts_last_five_linkage_groups",
      percent_text(sum(tested_by_sequence[sequence_levels %in% lg_last]), n_tested))
  add(tbl, "tested_transcripts_unplaced_scaffolds", count_text(tested_by_sequence[sequence_levels == "unplaced scaffolds"]))

  seq_tab <- data.table(
    table = tbl, sequence = sequence_levels,
    transcripts_in_table = count_by_sequence(per_tx$sequence),
    transcripts_surviving = count_by_sequence(per_tx[survives == TRUE, sequence]),
    transcripts_lost = count_by_sequence(per_tx[survives == FALSE, sequence]),
    genes_in_table = count_by_sequence(per_gene$sequence),
    genes_surviving = count_by_sequence(per_gene[survives == TRUE, sequence]),
    genes_lost = count_by_sequence(per_gene[survives == FALSE, sequence]),
    transcripts_in_universe = count_by_sequence(per_tx[in_universe == TRUE, sequence]),
    transcripts_tested = tested_by_sequence,
    genes_in_universe = count_by_sequence(per_gene[in_universe == TRUE, sequence]),
    genes_tested = count_by_sequence(per_gene[tested == TRUE, sequence]),
    universe_transcripts = universe_by_sequence,
    universe_transcripts_in_run_range = count_by_sequence(universe_sequence[universe_in_range]),
    annotated_mrnas_in_run_range = mrnas_in_range)
  count_columns <- setdiff(names(seq_tab), c("table", "sequence"))
  seq_tab <- rbind(seq_tab, cbind(data.table(table = tbl, sequence = "all"), seq_tab[, lapply(.SD, sum), .SDcols = count_columns]))
  seq_tab[, `:=`(percent_transcripts_surviving = round(100 * transcripts_surviving / transcripts_in_table, 1),
                 percent_genes_surviving = round(100 * genes_surviving / genes_in_table, 1),
                 percent_of_tested_transcripts = round(100 * transcripts_tested / n_tested, 1),
                 percent_of_universe_transcripts = round(100 * universe_transcripts / length(universe), 1))]
  by_sequence[[tbl]] <- seq_tab

  # ---- (e) overlap with the universe gap -------------------------------------------------------
  add(tbl, "transcripts_admitted_by_universe_alone", count_text(per_tx[in_universe == TRUE, .N]))
  add(tbl, "transcripts_admitted_by_reader_alone", count_text(n_surviving))
  add(tbl, "transcripts_admitted_by_both", count_text(n_tested))
  add(tbl, "genes_admitted_by_universe_alone", count_text(per_gene[in_universe == TRUE, .N]))
  add(tbl, "genes_admitted_by_reader_alone", count_text(per_gene[survives == TRUE, .N]))
  add(tbl, "genes_admitted_by_both", count_text(per_gene[tested == TRUE, .N]))
  add(tbl, "percent_surviving_transcripts_in_universe", percent_text(n_tested, n_surviving))
  add(tbl, "percent_lost_transcripts_in_universe",
      percent_text(per_tx[survives == FALSE & in_universe == TRUE, .N], per_tx[survives == FALSE, .N]))
  s9_reader <- one_row(step09_funnel[set == tbl & grepl("^3 ", stage)], "legacy_funnel.csv")
  s9_tested <- one_row(step09_funnel[set == tbl & grepl("^4 ", stage)], "legacy_funnel.csv")
  stopifnot(s9_reader$genes == per_gene[survives == TRUE, .N], s9_reader$transcripts == n_surviving,
            s9_tested$genes == per_gene[tested == TRUE, .N], s9_tested$transcripts == n_tested)

  # ---- (f) survival by number of transcripts and by annotation ---------------------------------
  size_classes <- list(one = c(1, 1), two_to_five = c(2, 5), more_than_five = c(6, Inf))
  for (k in names(size_classes)) {
    in_class <- per_gene[n_transcripts >= size_classes[[k]][1] & n_transcripts <= size_classes[[k]][2]]
    suffix <- if (k == "one") "transcript" else "transcripts"
    add(tbl, sprintf("genes_with_%s_%s", k, suffix), count_text(nrow(in_class)))
    add(tbl, sprintf("genes_with_%s_%s_surviving", k, suffix), count_text(in_class[survives == TRUE, .N]))
    add(tbl, sprintf("percent_genes_with_%s_%s_surviving", k, suffix), percent_text(in_class[survives == TRUE, .N], nrow(in_class)))
  }
  add(tbl, "transcripts_with_unpaired_apostrophe_term", count_text(per_tx[has_unpaired_term == TRUE, .N]))
  add(tbl, "transcripts_with_unpaired_apostrophe_term_surviving", count_text(per_tx[has_unpaired_term == TRUE & survives == TRUE, .N]))
  add(tbl, "percent_other_transcripts_surviving",
      percent_text(per_tx[has_unpaired_term == FALSE & survives == TRUE, .N], per_tx[has_unpaired_term == FALSE, .N]))

  survivors[[tbl]] <- list(transcripts = per_tx[survives == TRUE, transcript], run_range = run_range,
                           run_first = run_first, run_last = run_last)
  log_msg("%-6s %6d lines; default read %6d rows, %4d of %4d transcripts, %4d of %4d genes; %2d lines with an unpaired apostrophe, %2d merged rows; warnings: %s",
          tbl, n_file_lines, nrow(default_read), n_surviving, nrow(per_tx), per_gene[survives == TRUE, .N], nrow(per_gene),
          sum(unpaired), sum(lines_per_row > 1L), read_warnings$default_quotes)
}
mechanism <- rbindlist(mechanism)
by_sequence <- rbindlist(by_sequence)
stopifnot(!anyDuplicated(mechanism[, .(table, quantity)]))
fwrite(mechanism, out_file("legacy_reader_mechanism.csv"))
fwrite(by_sequence, out_file("legacy_reader_by_linkage_group.csv"))

# ==========================================================================================
# 4. Position-only control for the BP test of the ploidy table (part g)
# ==========================================================================================
cutoff <- config$go_pvalue
legacy_list <- intersect(survivors$ploidy$transcripts, universe)                 # the transcripts of interest of the legacy test
run_pool <- universe[universe_accession >= survivors$ploidy$run_range[1] & universe_accession <= survivors$ploidy$run_range[2]]
log_msg("position control: legacy list %d transcripts; accession run %s to %s holds %d of the %d universe transcripts",
        length(legacy_list), survivors$ploidy$run_first, survivors$ploidy$run_last, length(run_pool), length(universe))

# ---- topGO as step 09 calls it, on the map as it is and on one entry per transcript ---------
bp_fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
interest_factor <- function(names_in_order, of_interest) {
  f <- factor(as.integer(names_in_order %in% of_interest), levels = 0:1)
  names(f) <- names_in_order
  f
}
topgo_object <- function(go_map, of_interest) {
  suppressMessages(new("topGOdata", ontology = "BP", allGenes = interest_factor(names(go_map), of_interest),
                       annot = annFUN.gene2GO, gene2GO = go_map))
}
topgo_p <- function(object, go_map, of_interest) {
  object <- updateGenes(object, interest_factor(names(go_map), of_interest))
  score(suppressMessages(getSigGroups(object, bp_fisher)))
}
maps <- list(lines = gene2go, transcripts = gene2go[!duplicated(map_names)])
fixed_lists <- list(legacy_list = legacy_list, accession_run = run_pool)
topgo_results <- list()
objects <- list()
for (u in names(maps)) {
  objects[[u]] <- topgo_object(maps[[u]], legacy_list)
  for (l in names(fixed_lists)) topgo_results[[paste(u, l)]] <- topgo_p(objects[[u]], maps[[u]], fixed_lists[[l]])
  log_msg("topGO, totals in %s: %d BP terms; p < %g for %d (legacy list) and %d (accession run)", u, length(topgo_results[[paste(u, "legacy_list")]]),
          cutoff, sum(topgo_results[[paste(u, "legacy_list")]] < cutoff), sum(topgo_results[[paste(u, "accession_run")]] < cutoff))
}

# ---- The same test computed directly ---------------------------------------------------------
# topGO counts the transcripts of a term and those of interest among them once per name, and the two totals of
# the test once per entry of the map that has a BP term (issue I-43). With one entry per transcript the totals
# are counts of transcripts too. The p-value is the upper tail of the hypergeometric distribution, as in
# fisher.test(alternative = "greater"); a term with no transcript of interest has p = 1.
term_members <- genesInTerm(objects$lines)           # term -> transcripts that carry it or a descendant of it
bp_terms <- names(term_members)
entries <- genes(objects$lines)                      # map lines of the transcripts that have a BP term
test_names <- unique(entries)
entries_of_name <- tabulate(match(entries, test_names), nbins = length(test_names))
member_term <- rep(seq_along(term_members), lengths(term_members))
member_name <- match(unlist(term_members, use.names = FALSE), test_names)
term_size <- unname(lengths(term_members))
stopifnot(!anyNA(member_name), setequal(test_names, genes(objects$transcripts)),
          identical(sort(bp_terms), sort(names(genesInTerm(objects$transcripts)))))
direct_test <- function(of_interest, totals_in) {
  in_list <- test_names %chin% of_interest
  hits <- tabulate(member_term[in_list[member_name]], nbins = length(bp_terms))
  n_interest <- if (totals_in == "lines") sum(entries_of_name[in_list]) else sum(in_list)
  n_total <- if (totals_in == "lines") length(entries) else length(test_names)
  p <- rep(1, length(bp_terms))
  scored <- hits > 0
  p[scored] <- stats::phyper(hits[scored] - 1, n_interest, n_total - n_interest, term_size[scored], lower.tail = FALSE)
  names(p) <- bp_terms
  list(p = p, transcripts_in_test = sum(in_list), of_interest_as_counted = n_interest, total_as_counted = n_total)
}
for (u in names(maps)) for (l in names(fixed_lists)) {
  from_topgo <- topgo_results[[paste(u, l)]]
  if (!identical(direct_test(fixed_lists[[l]], u)$p[names(from_topgo)], from_topgo)) {
    stop("The direct test does not give the p-values of topGO: ", u, ", ", l, call. = FALSE)
  }
}

# ---- The legacy result of step 09 and its top terms ------------------------------------------
step09_bp <- fread(step09("legacy_fisher_ploidy_BP.csv"))
p_legacy <- topgo_results[["lines legacy_list"]]
top_terms <- step09_bp[rank <= config$go_legacy_topn, go_id]
stopifnot(setequal(step09_bp$go_id, names(p_legacy)),
          isTRUE(all.equal(unname(p_legacy[step09_bp$go_id]), step09_bp$p_fisher, tolerance = 1e-9)),
          sum(p_legacy < cutoff) == sum(step09_bp$below_p_cutoff),
          length(top_terms) == config$go_legacy_topn,
          setequal(names(sort(p_legacy))[seq_len(config$go_legacy_topn)], top_terms),
          sort(p_legacy)[config$go_legacy_topn] < sort(p_legacy)[config$go_legacy_topn + 1L])

# ---- Lists: two fixed, three kinds of random list --------------------------------------------
set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
draw_from <- function(pool, size) pool[sample.int(length(pool), size)]
kinds <- data.table(
  list_of_interest = c("legacy_list", "accession_run", "random_universe_size_of_run", "random_run_size_of_legacy_list",
                       "random_universe_size_of_legacy_list"),
  description = c(
    "Transcripts of interest of the legacy test: returned by the default read of the ploidy table and in the legacy universe",
    "Every universe transcript whose accession lies in the longest surviving run of the ploidy table; no DML information",
    "Random universe transcripts, as many as the accession run holds",
    "Random universe transcripts of the accession run, as many as the legacy list holds",
    "Random universe transcripts, as many as the legacy list holds"),
  uses_dml_information = c("yes", "no", "no", "no", "no"))
all_lists <- c(
  list(legacy_list = list(legacy_list), accession_run = list(run_pool)),
  list(random_universe_size_of_run = replicate(n_draws, draw_from(universe, length(run_pool)), simplify = FALSE)),
  list(random_run_size_of_legacy_list = replicate(n_draws, draw_from(run_pool, length(legacy_list)), simplify = FALSE)),
  list(random_universe_size_of_legacy_list = replicate(n_draws, draw_from(universe, length(legacy_list)), simplify = FALSE)))
stopifnot(identical(names(all_lists), kinds$list_of_interest))

per_list <- rbindlist(lapply(names(all_lists), function(kind_name) {
  rbindlist(lapply(seq_along(all_lists[[kind_name]]), function(i) {
    rbindlist(lapply(names(maps), function(u) {
      res <- direct_test(all_lists[[kind_name]][[i]], u)
      data.table(list_of_interest = kind_name, list_number = i, unit = u,
                 transcripts_in_list = length(all_lists[[kind_name]][[i]]),
                 transcripts_in_test = res$transcripts_in_test, of_interest_as_counted = res$of_interest_as_counted,
                 total_as_counted = res$total_as_counted,
                 terms_below = sum(res$p < cutoff), top_terms_below = sum(res$p[top_terms] < cutoff),
                 top_terms_not_below = if (all(res$p[top_terms] < cutoff)) "none" else paste(top_terms[res$p[top_terms] >= cutoff], collapse = "; "))
    }))
  }))
}))
unit_text <- c(lines = "lines of the universe map, as in step 09", transcripts = "transcripts, one entry per transcript")
control <- per_list[, .(
  lists = .N, transcripts_in_list = transcripts_in_list[1L],
  transcripts_with_bp_term_median = median(as.numeric(transcripts_in_test)),
  of_interest_as_counted_median = median(as.numeric(of_interest_as_counted)), total_as_counted = total_as_counted[1L],
  bp_terms = length(bp_terms), p_cutoff = cutoff,
  terms_below_cutoff_median = median(as.numeric(terms_below)), terms_below_cutoff_min = min(terms_below),
  terms_below_cutoff_q05 = unname(quantile(terms_below, 0.05)), terms_below_cutoff_q95 = unname(quantile(terms_below, 0.95)),
  terms_below_cutoff_max = max(terms_below),
  legacy_top_terms = length(top_terms),
  legacy_top_terms_below_cutoff_median = median(as.numeric(top_terms_below)), legacy_top_terms_below_cutoff_min = min(top_terms_below),
  legacy_top_terms_below_cutoff_max = max(top_terms_below),
  lists_with_no_legacy_top_term_below_cutoff = sum(top_terms_below == 0L),
  legacy_top_terms_not_below_cutoff = if (.N == 1L) top_terms_not_below[1L] else NA_character_),
  by = .(unit, list_of_interest)]
stopifnot(per_list[, uniqueN(transcripts_in_list) == 1L && uniqueN(total_as_counted) == 1L, by = .(unit, list_of_interest)][, all(V1)])
control <- merge(kinds, control, by = "list_of_interest", sort = FALSE)
control[, totals_counted_in := unit_text[unit]]
control[, unit_order := match(unit, names(maps))]
control[, list_order := match(list_of_interest, kinds$list_of_interest)]
setorder(control, unit_order, list_order)
control[, c("unit", "unit_order", "list_order") := NULL]
setcolorder(control, c("totals_counted_in", "list_of_interest", "description", "uses_dml_information"))
fwrite(control, out_file("legacy_reader_position_control.csv"))

# ==========================================================================================
# 5. Summary for the log
# ==========================================================================================
options(width = 250)
show <- function(quantities) {
  wide <- dcast(mechanism[quantity %in% quantities], quantity ~ table, value.var = "value")
  print(wide[match(quantities, quantity), c("quantity", intersect(names(tables), names(wide))), with = FALSE], row.names = FALSE)
}
cat("\nRows returned by each read, and warnings of the default read\n")
show(c("lines_in_file", paste0("rows_", names(quote_settings)), "warnings_default_quotes"))
cat("\nMechanism and interventions\n")
show(c("lines_with_unpaired_apostrophe", "merged_rows_default_quotes", "rows_after_replacing_unpaired_apostrophes",
       "repaired_read_returns_every_line", "rows_after_replacing_all_other_apostrophes", "control_read_returns_the_same_rows"))
cat("\nWhat is lost, and where the survivors are\n")
show(c("transcripts_in_file", "transcripts_default_quotes", "transcripts_lost", "genes_in_file", "genes_default_quotes", "genes_lost",
       "percent_skipped_lines_of_lost_transcripts", "longest_surviving_run_transcripts", "longest_surviving_run_first",
       "longest_surviving_run_last", "transcripts_admitted_by_both", "tested_transcripts_first_five_linkage_groups",
       "percent_tested_transcripts_first_five_linkage_groups", "genes_admitted_by_universe_alone", "genes_admitted_by_reader_alone",
       "genes_admitted_by_both", "percent_genes_with_one_transcript_surviving", "percent_genes_with_more_than_five_transcripts_surviving"))
cat(sprintf("\nUniverse transcripts on %s: %s%%\n", lg_first_text,
            mechanism[table == "universe" & quantity == "percent_universe_first_five_linkage_groups", value]))
cat(sprintf("\nPosition-only control, BP test of the ploidy table: terms with p < %g (%d random lists of each kind)\n", cutoff, n_draws))
print(control[, .(totals = sub(",.*$", "", totals_counted_in), list_of_interest, lists, transcripts_in_list,
                  median = terms_below_cutoff_median, min = terms_below_cutoff_min, max = terms_below_cutoff_max,
                  top_terms_median = legacy_top_terms_below_cutoff_median, top_terms_max = legacy_top_terms_below_cutoff_max,
                  lists_with_no_top_term = lists_with_no_legacy_top_term_below_cutoff)], row.names = FALSE)
log_msg("wrote legacy_reader_mechanism.csv (%d rows), legacy_reader_by_linkage_group.csv (%d rows), legacy_reader_position_control.csv (%d rows)",
        nrow(mechanism), nrow(by_sequence), nrow(control))

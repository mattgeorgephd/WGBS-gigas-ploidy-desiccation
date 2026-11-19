# universe_file_cut.R
#
# Purpose   Show where the gene universe of the legacy GO enrichment lost 15,666
#           transcripts (issue I-39). The universe was read off the first column
#           of a 20 GB table on gannet,
#             data/gene_annotations/BOTH/union_1x-BOTH.GeneIDs.geneOverlap.transcriptIDs.GOAnnot
#           which a notebook cell wrote with `join ... | sort | uniq > file`; the
#           saved output of that cell is "uniq: write error". The table is
#           sorted by transcript accession. At one place a line stops after five
#           fields and the text of a line that belongs far later in the sort
#           order follows on the same line. Everything between the two
#           accessions is missing from the file, and so from the universe.
#           This script fetches a few hundred bytes around that place with an
#           HTTP range request (nothing is downloaded in full) and reports
#             - the cut line and the offset of the cut;
#             - the last transcript before the cut and the first after it;
#             - the first and last annotated transcript absent from the
#               universe, from the pipeline's own table, for comparison.
#
# Run       Rscript 06_diagnostics/universe_file_cut.R        (a few seconds; needs gannet)
#
# Reads     the table on gannet (two requests: headers, and 4 KB around the cut)
#           04_output/08_gene_annotation/tables/transcript_table.csv.gz
# Writes    06_diagnostics/universe_file_cut.csv              one row: the facts
#           06_diagnostics/universe_file_cut_context.txt      the lines around the cut
#
# Notes     - The offset was found by an independent verifying agent with a
#             binary search over the sorted file (range requests of 4 KB). This
#             script checks that offset; it does not search. If the file on
#             gannet is replaced, the check fails with a message.
#           - One cut was located. The file was not scanned end to end; a
#             further cut that removes no whole transcript would not change the
#             universe.
#           - A cut at a multiple of 4,096 bytes that is followed by the start of
#             a later line is what a failed block write in the middle of a
#             buffered stream leaves behind. Why the write failed is not known.

suppressPackageStartupMessages(library(here))
source(here("run_config.R"))
invisible(lapply(list.files(here("02_R_functions"), pattern = "\\.R$", full.names = TRUE), source))
load_packages(c("data.table", "curl"))

url <- paste0("https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/",
              "data/gene_annotations/BOTH/union_1x-BOTH.GeneIDs.geneOverlap.transcriptIDs.GOAnnot")
expected_size <- 20149928521          # bytes, as on gannet on 2026-11-19
expected_cut  <- 5283561472           # byte offset at which the later line starts
n_fields <- 10                        # fields of a complete line

# ---- File size ---------------------------------------------------------------------------
head_req <- curl::curl_fetch_memory(url, handle = curl::new_handle(nobody = TRUE))
headers <- curl::parse_headers_list(head_req$headers)
size <- as.numeric(headers[["content-length"]])
if (!isTRUE(all.equal(size, expected_size))) {
  stop(sprintf("The table on gannet has %s bytes, expected %s: it is not the file that was examined.",
               format(size, big.mark = ","), format(expected_size, big.mark = ",")))
}

# ---- Bytes around the cut ----------------------------------------------------------------
window <- 2048
from <- expected_cut - window; to <- expected_cut + window - 1
chunk <- curl::curl_fetch_memory(url, handle = curl::new_handle(range = sprintf("%.0f-%.0f", from, to)))
stopifnot(chunk$status_code == 206, length(chunk$content) == 2 * window)
text <- rawToChar(chunk$content)
line_starts <- c(1L, gregexpr("\n", text, fixed = TRUE)[[1]] + 1L)       # positions in the chunk, 1-based
lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
offsets <- from + line_starts[seq_along(lines)] - 1                      # byte offset of each line start in the file
fields <- lengths(strsplit(lines, "\t", fixed = TRUE))
inner <- seq(2, length(lines) - 1)                                       # the first and last line of the chunk are partial
bad <- inner[fields[inner] != n_fields]
if (length(bad) != 1) stop("Expected exactly one malformed line in the window, found ", length(bad))

cut_line <- lines[bad]
parts <- strsplit(cut_line, "\t", fixed = TRUE)[[1]]
# The later line starts where a second transcript accession appears.
is_accession <- grepl("^[XN][MR]_[0-9]+\\.[0-9]+$", parts)
stopifnot(sum(is_accession) == 2, is_accession[1])
second <- which(is_accession)[2]
kept_part <- paste(parts[seq_len(second - 1)], collapse = "\t")
cut_offset <- offsets[bad] + nchar(kept_part, type = "bytes") + 1        # +1 for the tab after the kept fields
before <- parts[1]; after <- parts[second]

# ---- The gap in the universe, from the pipeline's table ------------------------------------
tx <- fread(output_path("08_gene_annotation", "tables", "transcript_table.csv.gz"))
annotated <- tx[in_legacy_go_map == TRUE][order(transcript_id)]
absent <- annotated[in_legacy_universe == FALSE]
runs <- rle(annotated$in_legacy_universe)
ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
longest <- which(!runs$values & runs$lengths == max(runs$lengths[!runs$values]))
gap_first <- annotated$transcript_id[starts[longest]]; gap_last <- annotated$transcript_id[ends[longest]]
pos <- function(id) match(id, annotated$transcript_id)

facts <- data.table(
  file = basename(url), file_bytes = size,
  cut_offset_bytes = cut_offset, cut_offset_is_multiple_of_4096 = cut_offset %% 4096 == 0,
  cut_offset_as_expected = isTRUE(all.equal(cut_offset, expected_cut)),
  bytes_kept_of_cut_line = nchar(kept_part, type = "bytes") + 1, fields_kept_of_cut_line = second - 1,
  last_transcript_before_cut = before, first_transcript_after_cut = after,
  annotated_mRNAs_absent_from_universe = nrow(absent),
  longest_run_absent = max(runs$lengths[!runs$values]), run_first = gap_first, run_last = gap_last,
  transcript_before_cut_is_just_before_run = isTRUE(pos(before) == pos(gap_first) - 1),
  transcript_after_cut_is_just_after_run = isTRUE(pos(after) == pos(gap_last) + 1),
  checked_on = as.character(Sys.Date()))
fwrite(facts, here("06_diagnostics", "universe_file_cut.csv"))

show <- seq(max(2, bad - 3), min(length(lines) - 1, bad + 3))
context <- c(sprintf("File: %s", url),
             sprintf("Lines around byte %s (the offset of each line start, then the line; tabs shown as ' | '):", format(expected_cut, big.mark = ",")),
             sprintf("%s%14.0f  %s", ifelse(show == bad, ">>", "  "), offsets[show], gsub("\t", " | ", lines[show], fixed = TRUE)))
writeLines(context, here("06_diagnostics", "universe_file_cut_context.txt"))

cat(sprintf("File size %s bytes.\n", format(size, big.mark = ",")))
cat(sprintf("A line of %s stops after %d fields; a line of %s begins at byte %s (%s x 4,096%s).\n",
            before, second - 1, after, format(cut_offset, big.mark = ","), format(cut_offset %/% 4096, big.mark = ","),
            if (cut_offset %% 4096 == 0) "" else sprintf(" + %d", cut_offset %% 4096)))
cat(sprintf("Longest run of annotated mRNAs absent from the universe: %s, %s to %s.\n",
            format(facts$longest_run_absent, big.mark = ","), gap_first, gap_last))
cat(sprintf("Transcript before the cut is the one just before that run: %s; transcript after the cut is the one just after it: %s.\n",
            facts$transcript_before_cut_is_just_before_run, facts$transcript_after_cut_is_just_after_run))

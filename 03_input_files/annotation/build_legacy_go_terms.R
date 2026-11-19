# build_legacy_go_terms.R
#
# Purpose   Recover the part of the lost GO-to-GOslim table that the legacy
#           annotation used, and write it as legacy_go_goslim_terms.tsv.gz.
#
#           The legacy notebook joined each transcript's GO IDs to
#           `GO-GOslim.sorted` (GO ID, term name, GOslim bin, ontology letter), a
#           2017 file on a server that no longer responds. No copy has been
#           found. The three legacy annotation tables in legacy/ carry those four
#           columns on every line, so the rows of the lost table can be read back
#           from them for every GO ID that occurs in a transcript of a gene that
#           holds a legacy DML. That is enough to regenerate the three tables
#           exactly (06_diagnostics/check_legacy_outputs.R confirms it), and it
#           is not enough for anything else: GO IDs of other transcripts are not
#           recoverable this way.
#
# Run       Rscript 03_input_files/annotation/build_legacy_go_terms.R
#           Run once; the output is tracked. Rerun only if the files in legacy/
#           change.
#
# Reads     03_input_files/annotation/legacy/DML-Cov10-20-BOTH-{PLOIDY,HEAT,COMBO}.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz
# Writes    03_input_files/annotation/legacy_go_goslim_terms.tsv.gz

suppressPackageStartupMessages({library(here); library(data.table)})

cols <- c("transcript", "geneID", "chr", "start", "end", "meth.diff", "GOID", "GOterm", "GOSlim", "GOcat")
tabs <- lapply(c("PLOIDY", "HEAT", "COMBO"), function(s) {
  fread(here("03_input_files", "annotation", "legacy",
             sprintf("DML-Cov10-20-BOTH-%s.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz", s)),
        header = FALSE, sep = "\t", quote = "", col.names = cols)
})
terms <- unique(rbindlist(tabs)[, .(GOID, GOterm, GOSlim, GOcat)])
setorder(terms, GOID, GOSlim)

# One GO ID never carries two ontology letters; a GO ID can sit in several GOslim bins.
stopifnot(terms[, uniqueN(GOcat), by = GOID][, all(V1 == 1)])

out <- here("03_input_files", "annotation", "legacy_go_goslim_terms.tsv.gz")
fwrite(terms, out, sep = "\t", quote = FALSE, compress = "gzip")
cat(sprintf("%d rows for %d GO IDs written to %s\n", nrow(terms), uniqueN(terms$GOID), out))
print(terms[, .(rows = .N, go_ids = uniqueN(GOID)), by = GOcat])

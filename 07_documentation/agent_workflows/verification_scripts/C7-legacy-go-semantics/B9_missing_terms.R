suppressPackageStartupMessages(library(data.table))
scratch <- "/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/C7-legacy-go-semantics"
P <- readRDS(file.path(scratch, "A_port.rds")); E <- lapply(P$DMLGOterms_exact, as.data.table); out <- P$out
full <- fread(cmd = "zcat /home/user/WGBS-gigas-ploidy-desiccation/03_input_files/annotation/geneid2go.tab.gz", header = FALSE, sep = "\t", col.names = c("id", "go"), quote = "")
full <- unique(full)
long <- full[, .(go = trimws(unlist(strsplit(go, ",", fixed = TRUE)))), by = id][grepl("^GO:[0-9]{7}$", go)]
long <- unique(long)
tab <- unique(rbindlist(E)[, .(transcript, GOID)])
tx <- unique(tab$transcript)
m <- long[id %chin% tx]
cat("transcripts in the three DML annotation tables:", length(tx), "; of these named in the complete map:", uniqueN(m$id), "\n")
cat("distinct GO IDs on their map lines:", uniqueN(m$go), "; distinct GO IDs in the tables:", uniqueN(tab$GOID),
    "; map GO IDs that never occur in the tables:", length(setdiff(unique(m$go), unique(tab$GOID))),
    "; table GO IDs not on the map lines:", length(setdiff(unique(tab$GOID), unique(m$go))), "\n")
m[, in_table := paste(id, go) %chin% paste(tab$transcript, tab$GOID)]
cat("(transcript, GO ID) pairs on the map lines:", nrow(m), "; present as table rows:", sum(m$in_table), sprintf("(%.1f%%)", 100 * mean(m$in_table)), "\n")
cat("is presence a property of the GO ID alone? GO IDs with some pairs present and some absent:", m[, .(a = any(in_table), b = any(!in_table)), by = go][a & b, .N], "\n")
top <- unique(unlist(lapply(out, function(r) r$top10)))
cat("distinct top-ten terms over the nine tests:", length(top), "; never in any DML annotation table:", sum(!top %chin% tab$GOID),
    "; of those, named directly on a map line:", sum(!top %chin% tab$GOID & top %chin% long$go), "\n")
per <- sapply(out, function(r) sum(!r$top10 %chin% tab$GOID)); print(per); cat("sum over the nine tests:", sum(per), "of 90\n")

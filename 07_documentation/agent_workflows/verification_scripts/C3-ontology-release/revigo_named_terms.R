suppressPackageStartupMessages(library(data.table))
scr <- getwd(); source(file.path(scr, "parse_archive.R"))
R <- "/home/user/WGBS-gigas-ploidy-desiccation/08_archive"
l <- readLines(file.path(R, "figures/revigo_figure/Revigo-PLOIDY-BP.R")); l <- l[grepl('c\\("GO:[0-9]{7}"', l)]
rv <- data.table(go_id = sub('.*c\\("(GO:[0-9]{7})","([^"]*)".*', "\\1", l), description = sub('.*c\\("(GO:[0-9]{7})","([^"]*)".*', "\\2", l))
fish <- as.data.table(parse_archive(file.path(R, "GOterms/BOTH-PLOIDY/PLOIDY-BP-FisherTestResults.csv")))
fish[, p := classic_numeric(classic)]
dml_terms <- fread(file.path(R, "GOterms/BOTH-PLOIDY/ploidy-BP-GOterms.csv"))
rv[, in_dml_term_list := go_id %in% dml_terms$GOID]
rv <- merge(rv, fish[, .(go_id, annotated, significant, p, rank)], by = "go_id", all.x = TRUE)
cat(sprintf("REVIGO export: %d terms; in DML-associated GO list (%d ids): %d; in Fisher table: %d; with Fisher p < 0.01: %d; p >= 0.01: %d; p >= 0.05: %d; p >= 0.5: %d\n",
            nrow(rv), nrow(dml_terms), sum(rv$in_dml_term_list), sum(!is.na(rv$p)), sum(rv$p < 0.01, na.rm = TRUE), sum(rv$p >= 0.01, na.rm = TRUE),
            sum(rv$p >= 0.05, na.rm = TRUE), sum(rv$p >= 0.5, na.rm = TRUE)))
cat("REVIGO ids not in the DML-associated list:\n"); print(rv[in_dml_term_list == FALSE])
cat(sprintf("Fisher p<0.01 terms (%d): in DML-associated list %d; in REVIGO export %d\n", sum(fish$p < 0.01), sum(fish[p < 0.01]$go_id %in% dml_terms$GOID), sum(fish[p < 0.01]$go_id %in% rv$go_id)))
named <- c("DNA replication", "homologous chromosome pairing at meiosis", "chromosome segregation", "ATP generation from ADP",
           "proton transmembrane transport", "calcium ion transport")
cat("\nTerms named in Results 3.4 for ploidy-DML (Figure 6):\n")
print(rv[description %in% named, .(go_id, description, in_dml_term_list, annotated, significant, p, rank)])
print(fish[term %in% short_term(named), .(go_id, term, annotated, significant, classic, rank)])

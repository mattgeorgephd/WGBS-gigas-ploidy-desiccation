db <- news()
db <- as.data.frame(unclass(db)[c("Version", "Date", "Category", "Text")], stringsAsFactors = FALSE)
v <- numeric_version(sub("[^0-9.].*$", "", db$Version))
db <- db[v > "4.1.0" | (v == "4.1.0"), ]                      # 4.1.0 and later
pat <- "scan|readLines|read\\.|pushBack|text[- ]mode|line ending|CRLF|carriage|embedded nul|quot(e|ing)|apostrophe|comment\\.char|fill *=|blank\\.lines|skipNul|connection"
hit <- db[grepl(pat, db$Text, ignore.case = TRUE), ]
hit <- hit[order(numeric_version(sub("[^0-9.].*$", "", hit$Version))), ]
cat("entries 4.1.0 and later:", nrow(db), "; matching the broad pattern:", nrow(hit), "\n")
for (i in seq_len(nrow(hit))) cat(sprintf("\n[R %s | %s]\n%s\n", hit$Version[i], hit$Category[i], substr(gsub("\\s+", " ", hit$Text[i]), 1, 420)))

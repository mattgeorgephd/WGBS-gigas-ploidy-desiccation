# Search R's own NEWS database for entries that mention the table readers.
db <- news()                      # R NEWS for the 4.x series, as shipped with this R
db <- as.data.frame(unclass(db)[c("Version", "Date", "Category", "Text")], stringsAsFactors = FALSE)
cat("versions covered:", paste(unique(db$Version)[c(1, length(unique(db$Version)))], collapse = " ... "), "; entries:", nrow(db), "\n")
pat <- "read\\.table|read\\.csv|read\\.delim|\\bscan\\(|count\\.fields|readTableHead|type\\.convert|quote *=|quoted string|embedded (quote|newline)|EOF within"
hit <- db[grepl(pat, db$Text), ]
ver <- numeric_version(sub("[^0-9.].*$", "", hit$Version))
hit <- hit[order(ver), ]
for (i in seq_len(nrow(hit))) cat(sprintf("\n[R %s | %s | %s]\n%s\n", hit$Version[i], hit$Date[i], hit$Category[i], gsub("\\s+", " ", hit$Text[i])))

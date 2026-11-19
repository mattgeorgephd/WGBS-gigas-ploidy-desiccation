#!/bin/bash
# Europe PMC literature search with generic terms. Usage: epmc_search.sh <label> <query (already URL-encoded)> [pageSize]
W=/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/L2-mitochondria-and-library-preparation
label="$1"; q="$2"; n="${3:-25}"
d="$W/downloads/search_$label"; mkdir -p "$d"
curl -sS -A "literature-check/1.0" --max-time 60 -o "$d/search.json" -w "search $label HTTP %{http_code} size %{size_download}\n" "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=$q&format=json&resultType=lite&pageSize=$n&sort=CITED%20desc"

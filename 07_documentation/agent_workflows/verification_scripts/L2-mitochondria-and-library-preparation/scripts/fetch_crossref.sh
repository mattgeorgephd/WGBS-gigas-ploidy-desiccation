#!/bin/bash
# Fetch Crossref and Europe PMC records by DOI. Usage: fetch_crossref.sh <label> <doi>
# Each record goes into its own folder under downloads/. Only the DOI is sent.
W=/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/L2-mitochondria-and-library-preparation
label="$1"; doi="$2"
d="$W/downloads/crossref_$label"; mkdir -p "$d"
curl -sS -A "literature-check/1.0" --max-time 60 -o "$d/crossref.json" -w "crossref $label HTTP %{http_code}\n" "https://api.crossref.org/works/$doi"
e="$W/downloads/epmc_$label"; mkdir -p "$e"
curl -sS -A "literature-check/1.0" --max-time 60 -o "$e/epmc.json" -w "epmc $label HTTP %{http_code}\n" "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22$doi%22&format=json&resultType=core"

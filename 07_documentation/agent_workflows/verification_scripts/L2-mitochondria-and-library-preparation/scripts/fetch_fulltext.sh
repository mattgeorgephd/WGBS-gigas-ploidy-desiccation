#!/bin/bash
# Fetch open-access full text XML from Europe PMC. Usage: fetch_fulltext.sh <label> <PMCID>
W=/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf/L2-mitochondria-and-library-preparation
label="$1"; pmcid="$2"
d="$W/downloads/fulltext_$label"; mkdir -p "$d"
curl -sS -A "literature-check/1.0" --max-time 90 -o "$d/$pmcid.xml" -w "fulltext $label $pmcid HTTP %{http_code} size %{size_download}\n" "https://www.ebi.ac.uk/europepmc/webservices/rest/$pmcid/fullTextXML"

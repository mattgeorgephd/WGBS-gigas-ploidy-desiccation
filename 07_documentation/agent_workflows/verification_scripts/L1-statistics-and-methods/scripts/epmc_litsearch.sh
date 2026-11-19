#!/bin/bash
# Europe PMC literature search with generic terms; one request at a time
W="$1"; key="$2"; query="$3"
out="$W/downloads/epmc_litsearch/$key.json"
code=$(curl -sS --max-time 60 -A "literature-check/1.0" -G -o "$out" -w "%{http_code}" "https://www.ebi.ac.uk/europepmc/webservices/rest/search" --data-urlencode "query=$query" --data-urlencode "format=json" --data-urlencode "resultType=lite" --data-urlencode "pageSize=15")
echo "$key $code"

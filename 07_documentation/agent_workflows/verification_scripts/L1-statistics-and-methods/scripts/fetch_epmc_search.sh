#!/bin/bash
# Europe PMC record by DOI (core result: abstract, PMCID, access flags); one request at a time with a time limit
W="$1"
while IFS=$'\t' read -r key doi; do
  out="$W/downloads/epmc_search/$key.json"
  if [ ! -s "$out" ]; then
    code=$(curl -sS --max-time 40 -A "literature-check/1.0" -G -o "$out" -w "%{http_code}" "https://www.ebi.ac.uk/europepmc/webservices/rest/search" --data-urlencode "query=DOI:\"$doi\"" --data-urlencode "format=json" --data-urlencode "resultType=core")
    echo "$key $doi $code"
    sleep 1
  else
    echo "$key $doi cached"
  fi
done < "$W/scripts/dois.tsv"

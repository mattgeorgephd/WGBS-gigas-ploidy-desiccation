#!/bin/bash
# fetch Crossref records, one at a time, generic client name
W="$1"
while IFS=$'\t' read -r key doi; do
  out="$W/downloads/crossref/$key.json"
  if [ ! -s "$out" ]; then
    code=$(curl -sS -A "literature-check/1.0" -o "$out" -w "%{http_code}" "https://api.crossref.org/works/$doi")
    echo "$key $doi $code"
    sleep 1
  else
    echo "$key $doi cached"
  fi
done < "$W/scripts/dois.tsv"

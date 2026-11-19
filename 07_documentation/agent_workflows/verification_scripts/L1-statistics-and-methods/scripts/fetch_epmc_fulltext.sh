#!/bin/bash
# Europe PMC full text XML of articles with a PMCID; one request at a time with a time limit
W="$1"; LIST="$2"
while IFS=$'\t' read -r key pmcid; do
  out="$W/downloads/epmc_fulltext/${key}_${pmcid}.xml"
  if [ ! -s "$out" ]; then
    code=$(curl -sS --max-time 150 -A "literature-check/1.0" -o "$out.part" -w "%{http_code}" "https://www.ebi.ac.uk/europepmc/webservices/rest/$pmcid/fullTextXML")
    if [ "$code" = "200" ]; then mv "$out.part" "$out"; else mv "$out.part" "$out.http$code"; fi
    echo "$key $pmcid $code $(stat -c %s "$out" 2>/dev/null)"
    sleep 1
  else
    echo "$key $pmcid cached $(stat -c %s "$out")"
  fi
done < "$LIST"
echo DONE

Check S5-literature: index of scripts and saved material (read-only on the repository; nothing sent to outside services but DOIs, titles and generic search terms).
scripts/lookup_doi.py            Crossref and Europe PMC records by DOI (metadata, abstracts) -> downloads/crossref, downloads/epmc, text/lookup_batch*.txt
scripts/epmc_fulltext.py         open-access full text from Europe PMC -> downloads/epmc_fulltext/*.xml, text/<label>_<PMCID>.txt
scripts/epmc_search.py, crossref_search.py   literature searches (generic terms) -> text/search_*.txt, text/crsearch_*.txt
scripts/html_to_text.py          saved article pages to text (PMC page of Park and Wu 2016, Leek et al. 2010; DSS user's guide; Springer abstract page)
scripts/methylkit_help.R         help page of calculateDiffMeth of the installed methylKit 1.28.0 -> text/methylkit_help.txt; text/methylkit_vignette_excerpt.txt
scripts/measure_fig2_bars.py     bar heights of Figure 2 of Wreczycka et al. 2017 from the embedded 300 dpi image -> text/fig2_bar_measurements.txt
scripts/mito_cpg_check.R         the 363 mitochondrial CpGs of united_cpg_legacy.csv.gz by library and by CpG -> mito/mito_cpg_check.out, mito/*.csv
downloads/mdc/16722oa.pdf        Wreczycka et al. 2017, publisher's version from the repository of the Max Delbrueck Center (https://edoc.mdc-berlin.de/16722/); text in text/wreczycka2017_mdc16722_raw.txt
downloads/fao/cd1852en.pdf       Vignier, Adams and Lovatelli 2024 (FAO technical paper 698); text in text/fao_cd1852en.txt
locus/                           united-table rows of NC_047565.1:35030701 and its neighbour
Not reachable: bioRxiv full text (HTTP 429, three attempts), ScienceDirect (HTTP 403), PMC page of Korthauer et al. (browser check), Europe PMC full text of author manuscripts (HTTP 500).

Second reading of the literature, reader L1-statistics-and-methods, 2026-11-19.
Read-only on the repository. Requests to outside services carried DOIs, titles of published works, names of public code
repositories and generic search terms only, with the client name "literature-check/1.0".

WHAT WAS READ, AND WHERE (all files below are under this folder)

Bibliographic records
  downloads/crossref/*.json            Crossref, https://api.crossref.org/works/<doi>, for the 17 DOIs of scripts/dois.tsv and the two
                                       records of McCullagh and Nelder, Generalized Linear Models (10.1007/978-1-4899-3242-6, 10.1201/9780203753736)
  downloads/epmc_search/*.json         Europe PMC record and abstract by DOI (webservices/rest/search, resultType=core)
  downloads/biorxiv_api/wreczycka_details.json   bioRxiv API, versions of preprint 10.1101/109512 (v1 2017-02-17, v2 2017-08-09)
  downloads/phipson2010_arxiv/api.xml  arXiv record 1603.05766 (journal reference: Stat. Appl. Genet. Molec. Biol. 9 (2010), Issue 1, Article 39)

Full texts retrieved by this reader
  Wreczycka et al. 2017, article       downloads/mdc_pdf/16722oa.pdf  from https://edoc.mdc-berlin.de/id/eprint/16722/1/16722oa.pdf
                                       (publisher's version, 11 pages; same SHA-256 as the first reader's copy); text/wreczycka2017_raw.txt, _layout.txt
  Wreczycka et al., preprint v2        downloads/biorxiv_fulltext/109512v2.full.html  from https://www.biorxiv.org/content/10.1101/109512v2.full
                                       text/wreczycka2017_preprint_v2.txt
  Wreczycka et al., code repository    downloads/github_bimsb_clone/repo  (https://github.com/BIMSBbioinfo/Strategies_for_analyzing_BS-seq,
                                       commit 48f3fb1 of 2017-08-10; the article names it in its Conclusions). Read, not run:
                                       diff_meth_tools/compareUsingSimulated.Rmd and .md, functions/dataSim2.R,
                                       compareUsing_CTCF_RRBS_pairsCelllines.Rmd and .md, diff_meth_figs/diff_meth-sim-vis-1.png and -vis-2.png
  Akalin et al. 2012                   downloads/epmc_fulltext/akalin2012_PMC3491415.xml   (Europe PMC fullTextXML)   text/akalin2012_PMC3491415.txt
  Feng, Conneely and Wu 2014           downloads/epmc_fulltext/feng2014_PMC4005660.xml                               text/feng2014_PMC4005660.txt
  Lazic 2010                           downloads/epmc_fulltext/lazic2010_PMC2817684.xml                              text/lazic2010_PMC2817684.txt
  Lea, Tung and Zhou 2015              downloads/epmc_fulltext/lea2015_PMC4657956.xml                                text/lea2015_PMC4657956.txt
  Liu et al. 2012                      downloads/epmc_fulltext/liu2012_PMC3491382.xml                                text/liu2012_PMC3491382.txt
  Silliman et al. 2023                 downloads/epmc_fulltext/silliman2023_PMC10468963.xml                          text/silliman2023_PMC10468963.txt
  Venkataraman et al. 2022             downloads/epmc_fulltext/venkataraman2022_PMC9351233.xml                       text/venkataraman2022_PMC9351233.txt
  Park and Wu 2016                     downloads/pmc_html/park2016_PMC12157722.html  (https://pmc.ncbi.nlm.nih.gov/articles/PMC12157722/)  text/park2016_PMC12157722.txt
  Lea et al. 2017                      downloads/ncbi_efetch/lea2017_PMC5656403.xml  (NCBI E-utilities efetch, db=pmc; author manuscript)  text/lea2017_PMC5656403.txt
  Leek et al. 2010                     downloads/ncbi_efetch/leek2010_PMC3880143.xml (NCBI E-utilities efetch, db=pmc; author manuscript)  text/leek2010_PMC3880143.txt
  Gelman and Stern 2006                downloads/gelman2006/signif4.pdf  (author's copy of the published article,
                                       https://sites.stat.columbia.edu/gelman/research/published/signif4.pdf)  text/gelman2006_signif4_raw.txt
  Nieuwenhuis et al. 2011              downloads/nieuwenhuis2011_author/NieuwenhuisEtAl2011.pdf  (author's copy of the published article,
                                       https://www.ejwagenmakers.com/2011/NieuwenhuisEtAl2011.pdf)  text/nieuwenhuis2011_author_copy.txt
  Phipson and Smyth 2010               downloads/phipson2010_arxiv/1603.05766.pdf  (https://export.arxiv.org/pdf/1603.05766, the authors' copy)  text/phipson2010_arxiv1603.05766.txt
  DSS user's guide 2.60.0              downloads/bioc_dss/DSS_vignette.html  (https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html)  text/DSS_vignette.txt
  dmrseq vignette 1.32.0               downloads/bioc_dmrseq/dmrseq_vignette.html (https://bioconductor.org/packages/release/bioc/vignettes/dmrseq/inst/doc/dmrseq.html)  text/dmrseq_vignette.txt
  methylKit 1.28.0 (installed)         text/methylkit_calculateDiffMeth_help.txt; vignette source /usr/local/lib/R/site-library/methylKit/doc/methylKit.Rmd
  stats::anova.glm, R 4.3.3            text/stats_anova_glm_help_R433.txt
  Bismark documentation                downloads/bismark_clone/repo/docs/src/content/docs/usage/*.md (https://github.com/FelixKrueger/Bismark, commit 4eb52d2)
  BS-Snper README and source           downloads/bssnper_clone/repo (https://github.com/hellbelly/BS-Snper, commit d4a2446 of 2021-04-15)

Abstract only
  Hurlbert 1984                        abstract in the Crossref record; the publisher's page answered HTTP 403
  Korthauer et al. 2019                abstract in the Europe PMC record; Europe PMC full text HTTP 500, NCBI efetch refused by the publisher,
                                       PMC page behind a reCAPTCHA check, bioRxiv HTTP 429

Uncited work found in the search and read in full (Europe PMC fullTextXML)
  Huh et al. 2019, Brief Bioinform 20:33, doi:10.1093/bib/bbx077            text/huh2019_PMC6357555.txt
  Piao et al. 2021, Int J Environ Res Public Health 18:7975, doi:10.3390/ijerph18157975   text/piao2021_PMC8345583.txt
  Venkataraman et al. 2024, Environ Epigenet 10:dvae018, doi:10.1093/eep/dvae018         text/venkataraman2024_PMC11556341.txt
  Tran et al. 2018, Genes 9:75, doi:10.3390/genes9020075 (not relevant to the statements)  text/tran2018_PMC5852571.txt
  searches: downloads/epmc_litsearch/s1.json to s6.json (queries in scripts/epmc_litsearch.sh calls, generic terms)

Refused or failed (not worked around)
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5656403/, PMC3880143, PMC6587918   reCAPTCHA page (after the first PMC page had been served)
  Europe PMC fullTextXML of PMC5656403, PMC3880143, PMC6587918, PMC12157722   HTTP 500
  NCBI efetch of PMC12157722 and PMC6587918                                   "The publisher of this article does not allow downloading of the full text in XML form"
  https://www.biorxiv.org/content/10.1101/183210v3.full                        HTTP 429
  https://esajournals.onlinelibrary.wiley.com/doi/10.2307/1942661              HTTP 403
  https://www.nature.com/articles/nn.2886.pdf                                  preview page of subscription content (abstract only)
  https://edoc.mdc-berlin.de/16393/ (repository record of the preprint)        HTTP 410
  https://api.github.com/repos/BIMSBbioinfo/...                                HTTP 403 (API not enabled for the session); the public repository was read through a shallow git clone

MEASUREMENTS OF FIGURE 2 OF WRECZYCKA ET AL. 2017
  scripts/measure_fig2.py          bars of the image embedded in the PDF (downloads/mdc_pdf_images/p5-000.png)      -> text/fig2_measurement_L1.txt
  scripts/measure_repo_fig.py      bars of the uncompressed rendering in the article's code repository (vis-2)       -> text/fig2_measurement_repo_png_L1.txt
  scripts/measure_repo_counts.py   bars of the counts figure of the repository (TP, FN, FP, TN of 2,500 each; vis-1) -> text/fig_counts_measurement_repo_png_L1.txt
  Result: unchanged CpGs called, of 2,500, at effect sizes 5, 10, 15, 20, 25:
    default test                     714, 730, 755, 803, 764   (28.6% to 32.1%)
    chi-squared with the correction  123, 115, 118, 162, 143   (4.6% to 6.5%)
    F with the correction             19,  21,  20,  32,  20   (0.8% to 1.3%)
    DSS                               about 2, 0, 4, 16, 25     (at most 1.0%)
    limma                             about 0, 4, 3, 8, 7       (at most 0.3%)
  Changed CpGs called (true positives): default 1124, 1499, 1752, 2069, 2208; chi-squared with the correction 257, 515, 718, 1069, 1273;
    F with the correction 69, 138, 176, 246, 243; DSS about 4, 30, 179, 495, 971; limma about 5, 51, 126, 247, 315.
  The calls were made with getMethylDiff(difference = 5, qvalue = 0.01) (compareUsingSimulated.Rmd, lines 102 and 216).

OTHER SCRIPTS
  scripts/fetch_crossref.sh, fetch_epmc_search.sh, fetch_epmc_fulltext.sh, epmc_litsearch.sh   downloads
  scripts/show_crossref.py, show_epmc.py, show_abstract.py, show_search.py                    print records
  scripts/jats_to_text.py, html_to_text.py                                                    conversion to text (standard library only)
  scripts/kwic.py                                                                             sentences that match a pattern
  scripts/find_quote.py                                                                       check that a quotation occurs verbatim in a text
  scripts/fig2_profile.py                                                                     pixel profile through one bar

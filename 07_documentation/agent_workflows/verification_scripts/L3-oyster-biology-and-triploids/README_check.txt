Second reading of the literature statements of issues.md, reader L3-oyster-biology-and-triploids (2026-11-19).
Scope: I-49 item 2 (*The animal.*), I-49 "Why it matters" (point on published alternatives), I-49 item 7 (one sentence),
I-49 "What would settle it" (genotypes from bisulfite reads), I-49 "Needs from authors" (three questions), I-08 (qualification
of 2026-11-19), I-01 "Needs from authors" (sentence citing Hurlbert). Read-only on the repository. Outside services received
DOIs, titles of published works, one BioProject accession named in a published article, and generic search terms only;
client name "literature-check/1.0"; no e-mail address. One request at a time.
The log was edited by the session at 19:11 (corrections of readers L1 and L2); the passages of this scope were re-read after
that edit and the quoted clauses checked again (scripts/check_quotes.py, scripts/quotes.txt: each occurs once).

SCRIPTS (scripts/)
  fetch_records.py        Crossref and Europe PMC records by DOI -> downloads/records_<label>/
  show_records.py         prints bibliographic fields and abstracts of the saved records -> text/records_batch1..5.txt
  epmc_fulltext.py        open-access full text XML from Europe PMC -> downloads/epmc_fulltext_<label>/, text/<label>_<PMCID>.txt
  html_to_text.py         saved publisher pages to text (parser of the standard library; scripts and styles dropped)
  openalex_locations.py   open-access locations of a work by DOI -> downloads/openalex_<label>/, text/openalex_locations.txt
  search.py               Crossref and Europe PMC searches with generic terms -> downloads/search_*/, text/search_*.txt
  check_quotes.py, quotes.txt   occurrence of each quoted clause in issues.md

WORKS OF THE SCOPE: WHAT WAS READ AND WHERE
  Johnson, Sirovy and Kelly 2022, J Anim Ecol 91:1135, doi:10.1111/1365-2656.13645
      FULL TEXT, published version (13 pages) from the NSF Public Access Repository: https://par.nsf.gov/servlets/purl/10324761
      (record https://par.nsf.gov/biblio/10324761-...). downloads/nsfpar_johnson2022/10324761.pdf. The text layer of the PDF is
      garbled (font encoding), so the pages were read as images; crops of Table 2 and three passages in pages/.
      Refused: Wiley (HTTP 403, browser check), NOAA repository (HTTP 403, "Access Denied").
  Johnson and Kelly 2020, Evol Appl 13:945, doi:10.1111/eva.12912
      FULL TEXT, Europe PMC PMC7232765 -> text/johnson2020_PMC7232765.txt
  Jiang, Li, Yu and Kong 2016, Mar Biotechnol 18:124, doi:10.1007/s10126-015-9674-4
      ABSTRACT (Europe PMC MED 26585587; https://link.springer.com/article/10.1007/s10126-015-9674-4 is a subscription preview).
      Tissue not stated in the abstract. text/springer_jiang2016.txt
  Olson and Roberts 2014, bioRxiv, doi:10.1101/012831
      ABSTRACT (Crossref, Europe PMC PPR35466); versions from https://api.biorxiv.org/details/biorxiv/10.1101/012831
      (v1 2014-12-16, v2 2015-02-18, v3 2015-03-13, not published in a journal). Full text refused: HTTP 429 (error 1015).
  Teng, Xu, Liu, Yu, Kong and Li 2025, Mar Biotechnol 27:31, doi:10.1007/s10126-024-10409-6
      ABSTRACT and data availability statement of the article (Springer page, subscription preview) -> text/springer_teng2025.txt
      FULL TEXT of the preprint of the same study: Teng et al. 2024, Research Square, doi:10.21203/rs.3.rs-3866716/v1,
      https://www.researchsquare.com/article/rs-3866716/v1 -> downloads/researchsquare_teng2024/, text/teng2024_researchsquare_v1.txt
      BioProject PRJNA801419 (named in the article): NCBI E-utilities summary -> downloads/ncbi_bioproject_PRJNA801419/
      ("The whole genome bisulfite sequencing of the Pacific oyster somatic tissues").
      Earlier preprint of the same group: Li et al. 2022, bioRxiv, doi:10.1101/2022.02.22.481396 (abstract only).
  Silliman, Spencer, White and Roberts 2023, Genome Biol Evol 15:evad013, doi:10.1093/gbe/evad013
      FULL TEXT, Europe PMC PMC10468963 -> text/silliman2023_PMC10468963.txt
  Sun, Yu and Li 2022, Mar Biotechnol 24:18, doi:10.1007/s10126-021-10083-y
      ABSTRACT (Europe PMC MED 35041105; Springer page is a subscription preview) -> text/springer_sun2022mb.txt
  Sun, Li and Yu 2022, Gene 820:146260, doi:10.1016/j.gene.2022.146260
      ABSTRACT (Europe PMC MED 35121028). ScienceDirect refuses automated access (HTTP 403).
  Valk, Engelsma, Megens, Kamermans, Murk and Nijland 2025, BMC Genomics 26:567, doi:10.1186/s12864-025-11736-1
      FULL TEXT, Europe PMC PMC12144841 -> text/valk2025_PMC12144841.txt
  Hurlbert 1984, Ecol Monogr 54:187, doi:10.2307/1942661
      ABSTRACT (Crossref, OpenAlex). Wiley refuses automated access (HTTP 403); no open copy listed by OpenAlex.
  Lazic 2010, BMC Neurosci 11:5, doi:10.1186/1471-2202-11-5
      FULL TEXT, Europe PMC PMC2817684 -> text/lazic2010_PMC2817684.txt
  Lazic and Essioux 2013, BMC Neurosci 14:37, doi:10.1186/1471-2202-14-37
      FULL TEXT, Europe PMC PMC3661356 -> text/lazic2013_PMC3661356.txt
  Lea, Tung and Zhou 2015, PLoS Genet 11:e1005650, doi:10.1371/journal.pgen.1005650
      FULL TEXT, Europe PMC PMC4657956 -> text/lea2015_PMC4657956.txt
  Lea, Vilgalys, Durst and Tung 2017, Nat Ecol Evol 1:1074, doi:10.1038/s41559-017-0229-0
      FULL TEXT (author manuscript), Europe PMC PMC5656403 -> text/lea2017_PMC5656403.txt
  Guo, DeBrosse and Allen 1996, Aquaculture 142:149, doi:10.1016/0044-8486(95)01243-5
      METADATA ONLY (Crossref, OpenAlex; no abstract in either). ScienceDirect HTTP 403. Known here through the FAO manual.
  Vignier, Adams and Lovatelli 2024, FAO Fisheries and Aquaculture Technical Paper 698, doi:10.4060/cd1852en
      FULL TEXT, 68 pages, retrieved from https://openknowledge.fao.org/bitstreams/e8732dc8-6111-4fda-a1ec-882caad1960c/download
      (same file as the first reader's copy, SHA-256 8674c625...be4d) -> downloads/fao_pdf/, text/fao_cd1852en_layout.txt
  Miller, Elliott, Vaillancourt, Kube and Koutoulis 2014, Aquaculture 433:318, doi:10.1016/j.aquaculture.2014.06.025
      ABSTRACT (Europe PMC AGR IND605323399; UTAS repository record on figshare, file under embargo). ScienceDirect HTTP 403.
  Miller, Elliott, Vaillancourt, Koutoulis and Henshall 2016, Aquac Res 47:1288, doi:10.1111/are.12588
      ABSTRACT (Europe PMC AGR IND605153131; UTAS repository record on figshare, file under embargo). Wiley HTTP 403.

FURTHER WORKS READ WHILE LOOKING FOR CONTRADICTING OR LIMITING EVIDENCE
  Wang, Li, Wang, Que, Zhang and Li 2021, Heredity 126:10, doi:10.1038/s41437-020-0351-7
      FULL TEXT on the publisher's page (free to read): https://www.nature.com/articles/s41437-020-0351-7 -> text/wang2021heredity_nature.txt
  Wang, Li, Wang, Zhang and Li 2021, Environ Res 197:111058, doi:10.1016/j.envres.2021.111058        ABSTRACT (Europe PMC)
  Jiang, Li, Yu and Kong 2017, Genes Genom 39:997, doi:10.1007/s13258-017-0563-2                      ABSTRACT (Springer page) -> text/springer_jiang2017.txt
  Covelo-Soto, Leunda, Perez-Figueroa and Moran 2015, Anim Genet 46:280, doi:10.1111/age.12287        ABSTRACT (Europe PMC)
  Salem, Al-Tobasei, Ali and Kenney 2022, Genes 13:1151, doi:10.3390/genes13071151                    FULL TEXT, PMC9319582
  Beato, Sanchez-Baizan, Felip and Piferrer 2025, BMC Genomics 26:1004, doi:10.1186/s12864-025-12220-6  FULL TEXT, PMC12593797
  Zhou et al. 2016, Int J Mol Sci 17:1299, doi:10.3390/ijms17081299                                   FULL TEXT, PMC5000696
  Han et al. 2021, Sci Rep 11:7564, doi:10.1038/s41598-021-87278-9                                    FULL TEXT, PMC8027599
  Ren et al. 2022, Epigenetics Chromatin 15:19, doi:10.1186/s13072-022-00455-6 (allotriploid hybrids) FULL TEXT, PMC9123727 (skimmed, not used)
  Xiao et al. 2013, PLoS ONE 8:e56409, doi:10.1371/journal.pone.0056409 (allotetraploid hybrids)      ABSTRACT (not used)
  Gawra et al. 2023, Sci Adv 9:eadh8990, doi:10.1126/sciadv.adh8990                                   FULL TEXT, PMC10491289 (skimmed, not used)
  Gavery and Roberts 2017, PeerJ 5:e4147 (review)                                                     FULL TEXT, PMC5723431 (searched for the tissue of Jiang et al. 2016; not stated)

SEARCHES (text/search_*.txt)
  Europe PMC: works citing Jiang et al. 2016 (CITES:26585587_MED); titles on ploidy and methylation in animals; abstracts on
  ploidy and methylation in molluscs; abstracts on triploidy and methylation in fish and other animals; title of Johnson, Sirovy and Kelly.
  Crossref: four generic queries on triploid or polyploid bivalves and DNA methylation; title of Johnson, Sirovy and Kelly.
  NSF Public Access Repository: title of Johnson, Sirovy and Kelly.

ACCESS REFUSED (not circumvented)
  onlinelibrary.wiley.com (besjournals, esajournals, pdfdirect): HTTP 403, browser check, three works, one attempt each.
  www.sciencedirect.com: HTTP 403, one attempt (Miller et al. 2014); not tried for Guo et al. 1996 and Sun, Li and Yu 2022.
  www.biorxiv.org: HTTP 429 (error 1015), one attempt (Olson and Roberts).
  repository.library.noaa.gov: HTTP 403 ("Access Denied"), one attempt.
  Europe PMC full text XML of PMC7852555 (Wang et al. 2021, free but not open access): HTTP 500; read on the publisher's page instead.
  figshare.utas.edu.au records 22917161 and 22917155 (Miller et al.): files under embargo.

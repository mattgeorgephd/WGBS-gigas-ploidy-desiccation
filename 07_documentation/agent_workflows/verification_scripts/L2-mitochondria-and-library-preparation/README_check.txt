Second reading L2-mitochondria-and-library-preparation (2026-11-19)
Scope: issues.md, I-01 update (first point, from "Published work supports the inference." to the note on how the sources were read),
I-49 item 2 (mitochondrial clause, "In the control cohort ..." sentence, *Library preparation.* part, CHH as measure of non-conversion),
and the row of the corrections table on the mitochondrial genome.
Read-only on the repository. Requests to outside services carried only DOIs, PMCIDs, titles and generic search terms; client name "literature-check/1.0".

SCRIPTS (scripts/)
  fetch_crossref.sh     Crossref record (https://api.crossref.org/works/<doi>) and Europe PMC core record by DOI -> downloads/crossref_<label>, downloads/epmc_<label>
  show_record.py        prints the bibliographic fields of both records (run with python3 -I)
  fetch_fulltext.sh     Europe PMC open-access full text (https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML) -> downloads/fulltext_<label>
  jats_to_text.py       JATS XML to plain text by regular expressions, no XML parser (python3 -I) -> text/<label>_<PMCID>.txt
  epmc_search.sh, show_search.py   Europe PMC searches with generic terms -> downloads/search_<label>
  show_openalex.py      open-access locations of one work from OpenAlex
  count_quotes.py, quotes.txt   checks that each quoted clause of the log occurs exactly once in issues.md

WORKS OF THE SCOPE (DOI checked at Crossref; all nine DOIs resolve to the cited work, journal, volume and first page or article number)
  Hong et al. 2013, Mol Cell Biol 33:2683, doi:10.1128/MCB.00220-13           ABSTRACT ONLY (Europe PMC record). Full text not retrievable: Europe PMC fullTextXML HTTP 500 (not in the open-access subset);
                                                                                NCBI efetch db=pmc returns front matter only ("The publisher of this article does not allow downloading of the full text in XML form");
                                                                                https://pmc.ncbi.nlm.nih.gov/articles/PMC3700126/ returns a reCAPTCHA challenge; https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC3700126&blobtype=pdf HTTP 403;
                                                                                https://mcb.asm.org/content/mcb/33/14/2683.full.pdf HTTP 403; https://doi.org/10.1128/MCB.00220-13 -> tandfonline.com HTTP 403. Not retried, not circumvented.
  Liu et al. 2016, Sci Rep 6:23421, doi:10.1038/srep23421                      full text, PMC4800444 -> text/liu2016_PMC4800444.txt
  Mechta et al. 2017, Front Genet 8:166, doi:10.3389/fgene.2017.00166          full text, PMC5671948 -> text/mechta2017_PMC5671948.txt
  Owa et al. 2018, PLoS ONE 13:e0192722, doi:10.1371/journal.pone.0192722      full text, PMC5805350 -> text/owa2018_PMC5805350.txt
  Leroux et al. 2024, Sci Rep 14:31544, doi:10.1038/s41598-024-83368-6         full text, PMC11682101 -> text/leroux2024_PMC11682101.txt (the log marks it as abstract only)
  Ji et al. 2014, Front Genet 5:341, doi:10.3389/fgene.2014.00341              full text, PMC4204604 -> text/ji2014_PMC4204604.txt (the log marks it as abstract only)
  Olova et al. 2018, Genome Biol 19:33, doi:10.1186/s13059-018-1408-2          full text, PMC5856372 -> text/olova2018_PMC5856372.txt
    its correction, Genome Biol 20:43, doi:10.1186/s13059-019-1656-9           full text, PMC6385436 -> text/olova2019corr_PMC6385436.txt (replaces Additional file 3 only)
  Lea et al. 2017, Nat Ecol Evol 1:1074, doi:10.1038/s41559-017-0229-0         full text (author manuscript), PMC5656403 -> text/lea2017_PMC5656403.txt
  Wang et al. 2014, BMC Genomics 15:1119, doi:10.1186/1471-2164-15-1119        full text, PMC4378273 -> text/wang2014_PMC4378273.txt
  Supplementary files of these works were not opened.

FURTHER WORKS READ (not cited in the scoped passages)
  Olson and Roberts 2014, Front Physiol 5:224, doi:10.3389/fphys.2014.00224    full text, PMC4060414 (C. gigas male gametes: mitochondrial genome unmethylated)
  Zhou et al. 2019, Sci Rep 9:10383, doi:10.1038/s41598-019-46875-5            full text retrieved, PMC6637168; read: results on library preparation, mitochondrial paragraph, discussion, methods
  Morrison et al. 2021, Epigenetics Chromatin 14:28, doi:10.1186/s13072-021-00401-y   full text retrieved, PMC8214260; read: design, passages on conversion controls and replicates
  Yin, Gong and Mao 2026, Epigenetics Chromatin 19:44, doi:10.1186/s13072-026-00692-z  full text, PMC13688742
  Guitton, Nido and Tzoulis 2022, Nucleic Acids Res 50:9190, doi:10.1093/nar/gkac701   full text, PMC9458446
  Guitton et al. 2022, Epigenetics 17:906, doi:10.1080/15592294.2022.2045754   full text retrieved, PMC9423827; read: abstract, introduction, results on alignment and on the methylation landscape, discussion
  Patil et al. 2019, Nucleic Acids Res 47:10072, doi:10.1093/nar/gkz762        full text retrieved, PMC6821263; read: abstract, methods on pre-processing and noise filtering, first three results sections
  Dou et al. 2019, Genome Res 29:1622, doi:10.1101/gr.234021.117               full text retrieved, PMC6771398; read: abstract, introduction, first two results sections, passages on controls and artifacts
  Foox et al. 2021, Genome Biol 22:332, doi:10.1186/s13059-021-02529-2         full text retrieved, PMC8650396; read: passages on replicates and limitations
  Abstracts only (in the retrieved full text): Bicci et al. 2021 (Nucleic Acids Res 49:12757, doi:10.1093/nar/gkab1179), Matsuda et al. 2018 (Sci Rep 8:5801, doi:10.1038/s41598-018-24251-z),
    Goldsmith et al. 2021 (Sci Rep 11:8032, doi:10.1038/s41598-021-87457-8), Morris et al. 2020 (Environ Epigenet 6:dvaa001, doi:10.1093/eep/dvaa001), Lueth et al. 2021 (Front Aging Neurosci 13:713084, doi:10.3389/fnagi.2021.713084)

SEARCHES (Europe PMC, https://www.ebi.ac.uk/europepmc/webservices/rest/search)
  search_mito_meth_title, search_mito_meth_2   mitochondrial DNA methylation, bisulfite, non-CpG, nanopore, artifact, controversy (2018 to 2026)
  search_wgbs_libprep                           WGBS library preparation, benchmark, bias, batch effect
  search_pico, search_batch_wgbs, search_batch_abs, search_epiqc   Pico Methyl-Seq; batch effects in bisulfite sequencing; EpiQC
  OpenAlex (https://api.openalex.org/works/doi:10.1128/MCB.00220-13) for open copies of Hong et al.

NOT DONE
  The numbers of the log computed from the project data were not recomputed (out of scope by instruction).
  issues.md changed while this check ran (other sessions); the scoped passages were read again at the end and are unchanged (line numbers then: 122, 1069, 1073, 1216).

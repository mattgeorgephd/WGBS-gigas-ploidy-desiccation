# annotation

Inputs for functional annotation of genes that contain differentially methylated loci.

## Files

| File | Contents | Origin |
|---|---|---|
| `20210605-cgigas-roslin-mito-blastx.outfmt6.gz` | blastx of the Roslin v1 transcripts (plus mitochondrial genes) against UniProt Swiss-Prot, BLAST tabular format 6, 93,983 rows. Query IDs have the form `lcl|<sequence>_mrna_<transcript>_<n>`; subject IDs the form `sp|<accession>|<entry name>`. | gannet project directory, `data/`; produced for Venkataraman et al. (2022) |
| `geneid2go-union1x-BOTH.tab.gz` | The transcript-to-GO map used by the legacy topGO analysis as its gene universe: one line per blastx hit with the transcript ID and a comma-separated list of GO IDs, 48,263 lines for 27,469 transcripts. 15,666 annotated mRNAs that are adjacent when sorted by accession (XM_011455141.3 to XM_034461018.1) are missing from it, and 49 others (issue I-39). Used by step 09 only. | gannet project directory, `analyses/GOterm-annotation/BOTH/`; md5 of the uncompressed file `87bbfa137a501e4d37ec8404abda08ca` |
| `geneid2go.tab.gz` | The complete transcript-to-GO map from the same notebook: 92,197 lines for 47,266 transcript names, 45,346 of them with GO terms. The source of the gene-to-GO map of steps 08 and 10. | gannet project directory, `data/geneid2go.tab` (June 2021); md5 of the uncompressed file `841dc4861b95d71dcc28ecd087fb321f` |
| `legacy_go_goslim_terms.tsv.gz` | GO ID, term name, GOslim bin, and ontology letter for the 4,564 GO IDs that occur in the legacy annotation tables: the recoverable part of the lost table `GO-GOslim.sorted`. Used only to regenerate the legacy tables. | Built by `build_legacy_go_terms.R` from the files in `legacy/` |
| `build_legacy_go_terms.R` | Builds the file above. Run once. | this repository |
| `legacy/DML-Cov10-20-BOTH-{PLOIDY,HEAT,COMBO}.GeneIDs.geneOverlap.transcriptIDs.GOAnnot.gz` | The three legacy annotation tables (one line per DML, transcript, and GO term; ten columns, no header) that the legacy enrichment script read. Reference for the checks in `06_diagnostics/`, and source of `legacy_go_goslim_terms.tsv.gz`. COMBO is the set the draft calls interaction-DML. | gannet project directory, `analyses/GOterm-annotation/BOTH/` |
| `legacy/DML-Cov10-20-BOTH-{PLOIDY,HEAT,COMBO}.geneOverlap.transcriptIDs.gz` | Legacy intermediate: one line per DML and transcript. Reference for the checks. | same |
| `legacy/cgigas_uk_roslin_v1_mRNA.transcriptID.geneID.gz` | Legacy intermediate: gene ID and transcript ID for the 63,341 mRNAs. Reference for the checks. | gannet project directory, `analyses/GOterm-annotation/` |

## Things to know

- The GO terms in the maps came from a 2017 UniProt-to-GO table (`uniprot-SP-GO.sorted`) on a server that no longer responds; no copy of that table has been found in the project locations. The maps are therefore the frozen annotation source: GO terms as UniProt held them in July 2017, attached to blastx hits of June 2021. See issues I-33 and I-42 in [`issues.md`](../../issues.md).
- The maps are keyed by transcript, one line per blastx alignment, so a transcript appears on as many lines as its best hit has aligned segments; all lines of a transcript are identical. The legacy analysis treated each transcript as a gene; see issue I-03.
- In the legacy maps, transcript IDs were rebuilt with an `XM_` prefix regardless of their true prefix. Of the 45,346 names with GO terms, 43,184 are `XM_` mRNAs, 123 are curated `NM_` mRNAs written as `XM_` (and therefore lost in the legacy analysis), and 2,038 are non-coding `XR_` transcripts written as `XM_`. The pipeline recovers the true IDs from the blastx query IDs (`true_transcript_ids()` in `02_R_functions/annotation.R`) and uses the mRNAs only.
- The files in `legacy/` are outputs of the legacy analysis, not inputs of the corrected one. They are kept here, under the GitHub-first rule, because the checks need them and because the lost GOslim table can only be recovered from them.
- Transcript-to-gene and transcript-to-product relations are read from `../genome_features/cgigas_uk_roslin_v1_mRNA.gff.gz`.

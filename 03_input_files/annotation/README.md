# annotation

Inputs for functional annotation of genes that contain differentially methylated loci.

## Files

| File | Contents | Origin |
|---|---|---|
| `20210605-cgigas-roslin-mito-blastx.outfmt6.gz` | blastx of the Roslin v1 transcripts (plus mitochondrial genes) against UniProt Swiss-Prot, BLAST tabular format 6, 93,983 rows. Query IDs have the form `lcl|<sequence>_mrna_<transcript>_<n>`; subject IDs the form `sp|<accession>|<entry name>`. | gannet project directory, `data/`; produced for Venkataraman et al. (2022) |
| `geneid2go-union1x-BOTH.tab.gz` | The transcript-to-GO map used by the legacy topGO analysis: one line per blastx hit with the transcript ID and a comma-separated list of GO IDs, 48,263 lines for 27,469 transcripts. | gannet project directory, `analyses/GOterm-annotation/BOTH/` |

## Things to know

- The GO terms in the map came from a 2017 UniProt-to-GO table (`uniprot-SP-GO.sorted`) on a server that no longer responds; no copy of that table has been found in the project locations. The map itself is therefore the frozen annotation source. See issue I-33 in [`issues.md`](../../issues.md).
- The map is keyed by transcript, and 5,874 transcripts appear on more than one line. The legacy analysis treated each transcript as a gene; see issue I-03.
- In the legacy map, transcript IDs were rebuilt with an `XM_` prefix regardless of their true prefix, so the 157 curated `NM_` transcripts in the blastx file are mislabelled there. The pipeline parses the true transcript ID from the blastx query ID.
- Transcript-to-gene and transcript-to-product relations are read from `../genome_features/cgigas_uk_roslin_v1_mRNA.gff.gz`.

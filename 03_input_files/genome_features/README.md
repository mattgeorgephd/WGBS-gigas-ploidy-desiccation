# genome_features

Genome feature tracks for the *Crassostrea gigas* Roslin v1 assembly (RefSeq GCF_902806645.1) with the mitochondrial genome (NC_001276.1). These are the tracks the legacy analysis intersected with CpG and DML positions. They were generated for Venkataraman et al. (2022, BMC Genomics 23:556) and are archived on OSF, doi:10.17605/OSF.IO/YGCTB, folder `genome-feature-files/`; the notebook that built them is `code/08-Generating-Genome-Feature-Tracks.ipynb` in the same archive.

All files are gzipped copies of the originals. `md5_uncompressed.txt` holds the md5 of each file before compression. The sizes of all thirteen match the OSF copies, and the md5 of the CDS and exon tracks, which were downloaded from OSF, match the OSF records.

GFF files use 1-based closed coordinates. BED files use 0-based half-open coordinates. The pipeline converts both to the same convention when it reads them (`02_R_functions/`).

## Files

| File | Feature | Source of the copy |
|---|---|---|
| `cgigas_uk_roslin_v1_gene.gff.gz` | Genes (38,208: 31,371 protein coding, 5,904 lncRNA, others) | this repository (legacy `old_WGBS/WGBS/GFF/`) |
| `cgigas_uk_roslin_v1_mRNA.gff.gz` | mRNA transcripts (63,341) | this repository |
| `cgigas_uk_roslin_v1_exon.gff.gz` | Exons | OSF |
| `cgigas_uk_roslin_v1_CDS.gff.gz` | Coding sequence | OSF |
| `cgigas_uk_roslin_v1_exonUTR.gff.gz` | Untranslated regions of exons | this repository |
| `cgigas_uk_roslin_v1_intron.bed.gz` | Introns | this repository |
| `cgigas_uk_roslin_v1_nonCDS.bed.gz` | Complement of coding sequence | this repository |
| `cgigas_uk_roslin_v1_lncRNA.gff.gz` | Long non-coding RNA | this repository |
| `cgigas_uk_roslin_v1_upstream.gff.gz` | 1 kb upstream flanks | this repository |
| `cgigas_uk_roslin_v1_downstream.gff.gz` | 1 kb downstream flanks | this repository |
| `cgigas_uk_roslin_v1_flanks.gff.gz` | Upstream and downstream flanks together | this repository |
| `cgigas_uk_roslin_v1_intergenic.bed.gz` | Intergenic regions | this repository |
| `cgigas_uk_roslin_v1_rm.te.bed.gz` | Transposable elements (RepeatMasker) | this repository |

Flank width, checked in the tracks themselves: 1,000 bp for 78% of records, shorter where a neighbouring feature or the end of a sequence intervenes, and zero length for a few records, which the pipeline drops on import.

The tracks overlap one another: a transposable element inside an intron is in both tracks, and exon contains CDS and UTR. See issue I-18 in [`issues.md`](../../issues.md) for what that means for percentages.

The exon and CDS tracks were previously fetched by the legacy notebooks from `owl.fish.washington.edu`, which no longer responds.

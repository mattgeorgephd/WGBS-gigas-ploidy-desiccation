# 03_input_files

Every input the pipeline reads. Nothing in this folder is produced by the pipeline, and nothing in `04_output/` is needed to rebuild it.

For what the pipeline does with these files, see the [root README](../README.md) and [`01_pipeline/README.md`](../01_pipeline/README.md).

## The rule: GitHub first

A file the pipeline needs lives in this folder if it fits in GitHub. The working ceiling is 50 MB per file, after gzip where the format allows it. R reads the gzipped files directly; none needs to be unpacked by hand.

Only files above that ceiling stay on gannet. They are listed in `gannet_manifest.csv` with URL, size, and md5. `01_pipeline/01_fetch_inputs.R` downloads them into `gannet_cache/` (git-ignored) and verifies the checksums. To put the cache somewhere else, for example on a larger drive, set the environment variable `WGBS_CACHE_DIR` before running.

## Contents

| Path | What it is | Origin | Read by |
|---|---|---|---|
| `metadata/` | Library and animal sample sheets | lab records | all steps |
| `elisa/` | Raw plate absorbance and plate map for the global 5-mC ELISA (N = 32) | plate reader exports, 2019-02-06 | step 03 |
| `qpcr/` | Raw Cq values, nine per-gene replicate tables and one consolidated table | qPCR instrument exports | step 11 |
| `bismark_reports/` | Bismark alignment, deduplication, and summary reports for the 22 libraries, all from the June 2021 alignment | gannet `/panopea/061021-big/` | steps 02, 03 |
| `genome_features/` | Thirteen feature tracks for the Roslin v1 assembly, gzipped | Venkataraman et al. 2022, OSF doi:10.17605/OSF.IO/YGCTB | steps 07, 08 |
| `annotation/` | blastx hits against Swiss-Prot and the transcript-to-GO map | gannet project directory | steps 08, 09 |
| `gannet_manifest.csv` | The inputs that stay on gannet: 22 coverage files, 9.13 GB in total | gannet project directory, `cov_files/` | steps 01, 04, 05 |
| `gannet_cache/` | Download cache; created on first run; not tracked | | |

Each subfolder has its own README with a file inventory.

## The coverage files

The 22 files named `<library>_R1.CpG_report.merged_CpG_evidence.cov` are the starting point of the methylation analysis. Each is the output of Bismark `coverage2cytosine --merge_CpG --zero_based`: one row per CpG dinucleotide with reads from both strands combined, and six tab-separated columns with no header:

| Column | Content |
|---|---|
| 1 | Sequence name (RefSeq accession; NC_047559.1 to NC_047568.1 are the ten linkage groups, NC_001276.1 the mitochondrial genome, NW_ unplaced scaffolds) |
| 2 | Start, zero-based position of the C |
| 3 | End, start + 2 |
| 4 | Percent methylation |
| 5 | Count of methylated reads |
| 6 | Count of unmethylated reads |

They are about 440 MB each, which is why they are not in GitHub. How they were produced is documented in [`07_documentation/upstream_processing/`](../07_documentation/upstream_processing/README.md).

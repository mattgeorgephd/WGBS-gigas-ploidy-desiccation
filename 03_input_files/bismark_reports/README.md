# bismark_reports

Bismark reports for the 22 libraries, all from the June 2021 alignment run (Bismark 0.22.3), the run that produced the coverage files used throughout the analysis. Copied unchanged from gannet `/panopea/061021-big/` on 2026-11-18. The job script of that run is in `07_documentation/upstream_processing/`.

## Files

| File | Contents |
|---|---|
| `061021-big_bismark_summary_report.txt` | Summary table for all 78 libraries of the run (this study uses 22 of them): read pairs, aligned pairs, duplicates removed, and methylated and unmethylated cytosine counts by context after deduplication. Step 02 reads this file. |
| `<library>_R1_bismark_bt2_PE_report.txt` | Per-library alignment report: input files, Bismark version, bowtie2 options, mapping efficiency, and methylation counts before deduplication. |
| `<library>_R1_bismark_bt2_pe.deduplication_report.txt` | Per-library deduplication report. |

## Things to know

- Two sets of methylation counts exist for each library: before deduplication (alignment report) and after (summary report). The pipeline uses the counts after deduplication, which are the reads the coverage files are built from. The legacy Figure 1B used counts before deduplication, taken from a mix of two alignment runs; see issue I-13 in [`issues.md`](../../issues.md).
- The per-library reports that were in the repository before the reorganization are in `08_archive/data/PE_reports/`.

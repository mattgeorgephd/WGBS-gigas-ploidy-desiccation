# upstream_processing

Records of the steps that ran before anything in this repository: read alignment and methylation extraction on the University of Washington Mox cluster. Nothing here is executed by the R pipeline. The files are kept so that the Methods can be checked against what was run.

All files were copied unchanged from gannet on 2026-11-18.

## Files

| File | Source on gannet | What it is |
|---|---|---|
| `061021-big_job.sh` | `/panopea/061021-big/job.sh` | Job script of the June to July 2021 alignment of all *C. gigas* WGBS libraries together (78 libraries, including the 22 used in this study). Bismark 0.22.3, bowtie2 2.3.4.1, samtools 1.9. This run produced the coverage files used for every differential methylation analysis in the manuscript. |
| `061021-big_job2.sh` | `/panopea/061021-big/job2.sh` | Second submission of the same job with the alignment block commented out and a different methylation extractor call (`--comprehensive --merge_non_CpG`). The directory also holds two SLURM logs and a `timeout/` folder, which suggests the first job reached its time limit and this one completed the downstream steps. |
| `061021-big_bismark_summary_report.txt` | `/panopea/061021-big/bismark_summary_report.txt` | Bismark summary for all 78 libraries of that run: read pairs, aligned, duplicates removed, methylated and unmethylated counts by context. |
| `030521-ronrosM_job.sh` | `/panopea/030521-ronrosM/job.sh` | Job script of the March 2021 alignment of the ten heat-treated libraries alone (zr3534). Bismark 0.21.0. Superseded by the June run. Identical in content to the legacy `old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh`, which is named for fastp but contains no trimming command. |
| `030521-ronrosM_bismark_summary_report.txt` | `/panopea/030521-ronrosM/bismark_summary_report.txt` | Bismark summary of the March run. |

## What these establish

- **Alignment settings, both runs:** `bismark --non_directional -score_min L,0,-0.6`, paired end, against the Roslin v1 assembly plus mitochondrial genome (`roslin_M`); then `deduplicate_bismark --paired`; then `bismark_methylation_extractor`; then `coverage2cytosine --merge_CpG --zero_based`, which produces the `*_R1.CpG_report.merged_CpG_evidence.cov` files (one row per CpG, both strands combined, zero-based start).
- **The same reads went into both runs.** The ten zr3534 libraries have identical read-pair counts in the two summaries, so the June run used the same trimmed reads as the March run under shorter file names.
- **The project's coverage files are the June outputs.** The files in the project's `cov_files/` directory on gannet have the same sizes as the June outputs, and the two that were checksummed (zr3534_7 and zr3644_14) are identical.
- **Trimming is not documented here.** The March inputs are named `*.fastp-trim.20201202.fq.gz`; no trimming command, version, or parameter set exists in this repository or in these directories. See issue I-28 in [`issues.md`](../../issues.md).
- **The heat-treated and control libraries differ as sequencing batches** in depth and duplication rate. See issue I-01.

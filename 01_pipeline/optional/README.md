# 01_pipeline/optional

Steps that are not part of the paper run. `run_pipeline.R` does not start them; each is run by hand and writes under `04_output/not_adopted/`. The optional analyses that use them are switches of `run_pipeline.R`:

| Analysis | Adopted in the paper? | Switch of `run_pipeline.R` | Setting in `run_config.R` (list `optional`) | Id or token | Outputs | What it changes | Issue |
|---|---|---|---|---|---|---|---|
| The paper | yes | none: `Rscript run_pipeline.R` | every setting at its default | `legacy` | `04_output/<step>/`, `05_manuscript/` | at least 10 reads per CpG in each of the 22 libraries, no upper coverage filter, first line of each coverage file dropped, logistic regression without overdispersion correction, DML at q < 0.01 and a difference of more than 20 points; steps 01 to 09 and 11 to 13 | |
| Coverage filter of the Methods | no | `--coverage-filter` | `coverage_filter_of_methods = TRUE` | `corrected` | `04_output/not_adopted/corrected/` | CpGs at or above the 99.9th percentile of coverage in a library removed; every line of each coverage file read | I-05, I-38 |
| 5x standard | no | `--coverage 5` | `coverage_standard = 5L` | `cov5` | `04_output/not_adopted/cov5/` | 5 reads per CpG in each library (2,297,236 CpGs tested instead of 480,579) | Venkataraman et al. 2022 |
| Overdispersion correction | no | `--overdispersion` | `overdispersion_chisq = TRUE` | `mnchisq` | `04_output/not_adopted/mnchisq/` | methylKit `overdispersion = "MN"`, `test = "Chisq"`; no other test is offered | I-06 |
| C-to-T variant exclusion | no | `--exclude-ct-variants` | `exclude_ct_variants = TRUE` | `ctexcl` | `04_output/not_adopted/ctexcl/` | CpGs listed by `01_pipeline/optional/ct_variant_loci.R` removed from the CpG sets of step 04 and, after the test, from the tested CpGs and the DML (fits, p-values and q-values unchanged); the run refuses to start without the list | I-08 |
| Any combination of the four | no | the switches together | the settings together | the tokens joined by `-` in the order above | `04_output/not_adopted/<id>/`, for example `corrected-cov5-ctexcl/` | each of its parts | |
| Gene-level GO enrichment | no | `--gene-level-go`, alone (on the paper's DML) or with any of the above | `gene_level_go = TRUE` | | `04_output/not_adopted/<id>/10_go_enrichment_gene_level/` and `.../11_semantic_similarity/` (the paper's DML: `<id>` = `legacy`) | step 10 and the gene-level part of step 11 | I-41 |
| Diagnostics | no | `--diagnostics` (all) or `--diagnostics name1,name2` | `diagnostics = "all"` or names | | `06_diagnostics/` | the checks and controls of the audit, listed in `06_diagnostics/diagnostics.csv` | |
| Side by side | no | made by every optional run | | | `04_output/not_adopted/comparisons/`, `05_manuscript/not_adopted/` | the tables that list analyses, from every analysis folder on disk, the paper first; the manuscript values of every analysis ("not run" where absent) | |

## `ct_variant_loci.R`: CpG sites at C-to-T variants

### What it does

It lists the CpG sites of the reference genome where the 22 libraries carry a C-to-T sequence variant at the C of the CpG. At such a site a T in the bisulfite reads is not evidence of an unmethylated C, so Venkataraman et al. 2022 (BMC Genomics 23:556), the study on which the draft's methods are modelled, removed these sites from the CpGs they characterised and from their candidate DML. This step reproduces their rule in R, so that the analyses can exclude the same kind of sites behind the option `--exclude-ct-variants` of `run_pipeline.R` (token `ctexcl`; setting `exclude_ct_variants` in the list `optional` of `run_config.R`), off by default.

The rule is that of BS-Snper (Gao et al. 2015, version 1.1) as the 2022 study used it, ported from BS-Snper's C and Perl source (`07_documentation/agent_workflows/scope_inputs/variant_rule_spec.md`):

- Reads: dropped if unmapped or MAPQ below 20; bases with quality below 15 become N; a read is dropped if it has more mismatches than max(1, round(0.02 x its reference span)), not counting bisulfite conversions. Both mates of a pair are counted, as BS-Snper does.
- Strands: Watson reads (FLAG 99, 147) and Crick reads (FLAG 83, 163). At a reference C only the Crick reads show the genotype; a Watson T counts as C.
- Per position: the candidate filter, the Bayesian genotype with BS-Snper's priors and post-rules, and the FILTER label (PASS or Low) of the record BS-Snper writes.
- A CpG is **listed** when a call at the C of the CpG on the plus strand (1-based position `start + 1`, where `start` is the 0-based C position of the coverage files) has REF C and ALT T or "T,G" (genotypes CT, TT, GT), whatever its FILTER label, either in the pooled call of all 22 libraries or in the call of any one library. The pooled call uses the summed per-library counts, which equals running BS-Snper on a merged BAM, because every BS-Snper filter acts on single reads.
- Parameters: BS-Snper defaults except minimum coverage 5 (`minhetfreq` 0.1, `minhomfreq` 0.85, `minquali` 15, `mincover` 5, `minread2` 2, `errorate` 0.02, `mapvalue` 20). `mincover` only sets the PASS/Low label, so it does not change which sites are listed.

### Decisions of the author

Recorded in `07_documentation/reproducibility_plan.md`, section 7, table "Decided":

| Date | Decision |
|---|---|
| 2026-11-20 | Variants are called by R code only, applying the rule of BS-Snper to the sorted, deduplicated alignments; BS-Snper itself is not run. Exclusion of the listed sites is an option, off by default, so the default run still reproduces the draft. |
| 2026-11-20 | Reads are counted with their full CIGAR (insertions, deletions and clipping handled), not with BS-Snper's last CIGAR operation. This is the one stated difference from 2022 (below). |
| 2026-11-20 | Every CpG site of the genome is screened once, so that one table serves the 5x and the 10x analyses. |
| 2026-11-20 | As in 2022: the C of the CpG on the plus strand only; every call whatever its FILTER label; the union of one pooled call over all 22 libraries and the 22 per-library calls. |
| 2026-11-20 | Placement and outputs as in `07_documentation/agent_workflows/scope_inputs/settings_build_spec.md`, section 6.14. |

### The one difference from 2022: full CIGAR

BS-Snper keeps only the last operation of each read's CIGAR string and lays the read out from its start with that operation alone, so a read with an insertion or deletion is partly dropped and partly shifted by the indel length (31% of the alignments carry an indel). This step lays every aligned block out at its reference position. The behaviour of BS-Snper is kept inside the script as an internal test switch (`cigar_mode = "bs_snper_last_op"` of `ct_count_window()`), used only for the regression test; the run itself refuses any CIGAR setting other than `"full"`.

On the test window of the 2022 alignments (NC_047566.1:24,250,001-24,500,000, 4,510 CpGs, eight animals), BS-Snper's union of nine calls is 394 sites; with the full CIGAR this step lists 421 (389 shared, 32 more, 5 fewer).

### Validation

- **Regression against BS-Snper (T1).** With the last-operation switch on, the port reproduces every C-to-T record at CpGs that BS-Snper wrote for the 2022 study in the test window: 876 records of the eight animals and 160 records of the merged BAM (1,036 records), with the same positions, all eight strand counts, genotype, ALT, FILTER and QUAL. The union is 394 against 394. The merged-BAM records are reproduced from the summed per-animal counts.
- **Agreement with the scoping pilot (T2).** For library zr3644_14, at the pilot's 6,170 test loci, the eight counts and eight quality sums equal the pilot's full-CIGAR counts at every locus, and the 216 C-to-T calls are the same, with the same genotype and FILTER.
- **Download and read checks.** Every byte range must come back as HTTP 206 with the expected Content-Range and size; the joined file must have the size of the HEAD Content-Length and of the manifest and end with the BGZF end-of-file block; it is then decoded once in full and its record count must equal the count in its index. While counting, the records read on each chromosome must equal the index's count for that chromosome, so a query that silently returns nothing (seen with a reused Rsamtools handle) stops the run instead of passing as empty.

### How to run

From the repository root, with the download cache on a drive with at least 15 GB free:

```sh
WGBS_CACHE_DIR=/home/user/gannet_cache \
  nohup setsid Rscript 01_pipeline/optional/ct_variant_loci.R \
  >> 04_output/not_adopted/ct_variants/cache/console.log 2>&1 < /dev/null &
```

The run is resumable. After an interruption or a machine restart, start the same command again: verified byte ranges, verified alignments, finished chromosomes and finished libraries are kept and skipped. A second run started while one is active stops at once (lock file `cache/run.lock`; a lock left by a dead process is taken over).

Progress: `tail 04_output/not_adopted/ct_variants/cache/console.log` (one line per large chromosome and per library) and `ls 04_output/not_adopted/ct_variants/cache/*/done.rds` (finished libraries). The final tables are written when all 22 libraries are finished.

Testing options: `--libraries id1,id2` processes only these libraries and writes no final table; `--keep-bam` keeps the alignment after counting.

The rule parameters are in `run_config.R` as the list `config$ct_variants`, with the values of the list `ct_variants_default` at the top of the script (the 2022 values; not toggles, decision 10.8 of the settings specification). The script uses `config$ct_variants` when it is defined, as it is since 2026-11-20, and its own list otherwise; the two hold the same values. (The run of the screen started on 2026-11-20 at 06:15 UTC read its settings before the list moved; the values are the same.)

### The option that uses the list: `--exclude-ct-variants`

With the option on (alone, `Rscript run_pipeline.R --exclude-ct-variants`, or with any other switch), the CpGs with `listed` TRUE in `ct_variant_cpgs.csv.gz` are removed:

- in step 04, from every library before any count (after its first line and its upper cutoff are taken), so the CpG counts, distributions, correlations and PCA of the analysis are over the other CpGs;
- in step 05, after the test: the fit is that of the same analysis without the option (`methylkit_<fit id>.rds`, no refit), and the listed CpGs leave the tested CpGs, the test tables and the DML lists; p-values, q-values and differences of the other CpGs are those of the fit on all CpGs, as in Venkataraman et al. 2022. `ct_exclusion_summary_<id>.csv` gives per contrast the tested CpGs and DML before and after and the numbers removed, and the kept CpGs and kept DML that the screen could not assess (fewer than `min_crick_assessable` = 3 Crick evidence reads at the C, summed over the 22 libraries, so that no call set could have listed them; from `ct_variant_assessability.csv.gz`, placed on the CpGs of the reference and checked against `ct_variant_cpg_index.csv`). `ct_excluded_cpgs_<id>.csv.gz` lists the removed tested CpGs with their p, q and difference per contrast and whether each was a DML;
- the later steps (06 to 09, and 10 and 11 with `--gene-level-go`) read the reduced tables of the analysis.

A run with the option refuses to start, with the command above, while the list is absent.

### Inputs

| Input | Source |
|---|---|
| 22 alignments and their indexes, `<library>_R1_bismark_bt2_pe.deduplicated.sorted.bam(.bai)` | gannet `/panopea/061021-big/`, listed with URL and size in `03_input_files/gannet_manifest_alignments.csv` (gannet publishes no md5 for them; the md5 of each download is recorded in `ct_variant_run.csv`). Downloaded one at a time into `<gannet cache>/alignments/` and deleted after counting |
| Reference genome `cgigas_uk_roslin_v1_genomic-mito.fa` | fetched by the step itself through `fetch_gannet_file()` (`03_input_files/gannet_manifest.csv`, md5 checked) |
| Library ids and experiments | `03_input_files/metadata/wgbs_libraries.csv` |

Nothing is read from the outputs of other steps.

### Outputs

In `04_output/not_adopted/ct_variants/tables/`:

| File | Content | Expected size |
|---|---|---|
| `ct_variant_cpgs.csv.gz` | One row per CpG with a C-to-T call in any call set: `chr`, `start` (0-based C, as in the coverage files), `pos1`, `listed` (the 2022 rule: pooled call of all 22 or any library), the genotype and FILTER of the pooled call, the genotypes of the two experiment pools, `n_libraries_called`, `libraries_called`, `any_pass`, and the assessability of the site (next row) | 5 to 15 MB (estimate: 0.3 to 0.9 million sites) |
| `ct_variant_calls.csv.gz` | One row per (CpG, call set) with a C-to-T call: call sets `pooled_all`, `pooled_zr3534`, `pooled_zr3644` and the 22 libraries; genotype, ALT, GT, FILTER, QUAL, depth, ALT reads, the eight strand counts and their mean qualities (BS-Snper's BSD and BSQ). The union of 2022 is rebuilt from it with `call_set %in% c("pooled_all", <library ids>)` | 15 to 25 MB (about 11 bytes per row; 40,656 rows for zr3644_14 alone) |
| `ct_variant_assessability.csv.gz` | One row per CpG site of the genome (13,245,836): `crick_pooled` (Crick C plus T reads at the C after the read filters, summed over the 22 libraries) and `n_assessable_zr3534`, `n_assessable_zr3644` (libraries of each experiment with at least 3 such reads, the minimum for a CT call). No coordinate columns: see below | about 23 MB |
| `ct_variant_cpg_index.csv` | One row per reference sequence: `chr`, `n_cpg`, `first_row`, `first_start`, `last_start`, and `md5_starts` (md5 of the 0-based starts written one per line), to rebuild and check the coordinates of the assessability rows | 17 kB |
| `ct_variant_run.csv` | `section, library_id, item, value`: the parameters, CIGAR mode, code version and genome md5; per library the alignment's size, md5 and date, records, reads by FLAG, reads passing MAPQ and the mismatch filter, reads with an indel, CpGs with reads, calls, download, check and counting times, retries and peak memory | 20 kB |
| `ct_variant_log.txt` | Summary: sites screened and listed, PASS and Low shares, pooled-only and library-only shares, assessability, calls per call set, output sizes | 3 kB |

`04_output/not_adopted/ct_variants/cache/` (git-ignored by the `.gitignore` of that folder) keeps, per library, the eight strand counts and quality sums at every CpG with reads (about 30 to 60 MB per library). From it any pooling (for example per experiment) or other rule parameters can be recomputed without the alignments. It also holds the console log and the run lock.

**Why the assessability table has no coordinates.** A site that no call lists is clean only if enough evidence reads were seen there; later steps must be able to report the CpGs of their analysed set (2,297,236 at 5x, 480,579 at 10x in legacy mode) that could not be assessed. That needs a value for every CpG. With `chr` and `start` (even written as gaps between neighbouring CpGs) the table is about 45 MB, too close to the 50 MB limit for committed files; without them it is about 23 MB. The rows are the CpGs of the reference in a fixed order (sequences in the order of the FASTA index, positions increasing, every plus-strand "CG"), which `ct_assessability_coordinates(fasta)` in the script rebuilds in about a minute; `ct_variant_cpg_index.csv` checks the rebuild per sequence. Listed sites carry the same three values, by coordinate, in `ct_variant_cpgs.csv.gz`.

### Cost

Measured on 4 cores and 15 GB of memory, one R process:

| Part | Time |
|---|---|
| Download (4 ranged streams, gannet at 20 to 28 MB/s) | 150 to 200 s per 3.3 to 5.7 GB alignment; up to about 6 minutes for the largest (6.8 GB) |
| Full decode and checks | about 1 minute per 3.3 GB |
| Counting, genome-wide | 9.7 minutes for zr3644_14 (40 million records); about 20 minutes for 77 million |
| Final tables | about 6 minutes |
| Whole run, 22 libraries (109 GB of alignments) | about 7 to 9 hours if gannet keeps its speed; slower downloads add time |

Peak memory: about 1.9 GB while counting the smallest library, about 2.6 GB for the final tables. Disk: twice the alignment being downloaded (parts and joined file, at most 13.5 GB), plus about 1 GB of cache.

### Caveats (from the scoping and its review)

- Plus-strand C only, as in 2022: a CpG whose G carries a variant (a C-to-T change on the minus strand) is not listed.
- Every call lists a site, whatever its FILTER, as in 2022; a single Crick T read with no other reads gives a TT (Low) call.
- Both mates of overlapping pairs are counted (as BS-Snper does), and the mismatch filter counts the variant base itself, so variant reads pass the read filter less often than reference reads.
- The model is diploid; at the same depth a single-copy allele in a triploid is called less often than a heterozygous allele in a diploid, so masking sensitivity depends on ploidy, and the listed set grows with the number of libraries and with depth.

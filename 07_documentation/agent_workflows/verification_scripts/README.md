# Scripts of the verifying agents

One folder per verification task of the workflow `audit-verification` (claims C1 to C8), per code review (R1 to R6), and per check of the workflow `check-sensitivity-text` that wrote scripts (S1 to S5). Each folder holds the scripts the agent wrote and its small result tables, copied from the scratch directory of the audit session.

**What these are.** A record of how each number in a **Verification** paragraph of [`issues.md`](../../../issues.md) was obtained. The agents were instructed to refute a statement of the issue log with code that shares nothing with the pipeline.

**What these are not.** They are not part of the pipeline and will not run from a clone: paths point at the session's scratch directory, and large intermediates (downloads from gannet, copies of the repository, package sources) were not copied. Numbers that a revision of the manuscript would rely on are being moved into `06_diagnostics/` scripts that do run from the repository.

| Folder | Statement checked | Reported in |
|---|---|---|
| `C1-universe` | The gene universe of the legacy GO enrichment is missing a block of 15,666 transcripts | I-39 |
| `C2-reader` | `read.table` skipped most of the ploidy annotation table | I-40, I-43 |
| `C3-ontology-release` | The legacy tables are reproduced only with GO.db 3.13.0 | I-37, I-42, I-43, I-22 |
| `C4a-eligible-theory` | The eligible-CpG background and the gene-level enrichment (reasoning and calibration) | I-41, I-06 |
| `C4b-eligible-reimplementation` | The gene-level enrichment, rebuilt without topGO or goseq | I-41 |
| `C5-feature-odds` | The odds ratios of DML by genome feature | I-41, I-14 |
| `C6-coverage-filter-locus` | The coverage filter, the first-line loss, and the locus on NC_047563.1 | I-05, I-09, I-38 |
| `C7-legacy-go-semantics` | What Figure 5, Figure 9 and Table 1 count | I-02, I-07, I-20, I-21, I-42, I-45 |
| `C8-counts-and-correlations` | CpG counts, the genome total, the correlation range | I-16, I-44 |
| `S1-I06-numbers-and-reading` | Issue I-06 as first written after the sensitivity run: every number recomputed from the cached fits; methylKit's corrected test rebuilt in base R and applied to relabelled libraries (p-value shares, list lengths, corrected calls); the dispersion reference; where the correction is not applied | I-06, I-01, I-02, I-07 |
| `S2-I49-evidence-and-alternatives` | Issue I-49 as first written: numbers recomputed, and eight alternative explanations tested (sequencing depth, read depth and pair similarity, the level and pair similarity, other rules for the high libraries, the slope analysis, the draft's own hypothesis, bisulfite conversion, the ELISA) | I-49, I-01, I-06, I-08 |
| `S3-code-review-library-level-r2` | The code of `06_diagnostics/library_methylation_level.R`: all 13 tables recomputed with other code (long-format sums, bit-pattern enumeration of the regroupings, exact integer arithmetic for every threshold decision); two timed runs of a copy of the script (byte-identical tables); edge cases of its functions; the CpGs at exactly 20 points; exact Spearman p-values; the second largest gap as the cut for the high libraries. The copy of the script and the 13 tables of each run were not kept | I-49; `06_diagnostics/README.md`; the code pass |
| `S4-updates-and-cross-file` | The updates of I-01, I-02, I-04, I-07 and I-08 after the sensitivity run, and every other file that repeats them: numbers recomputed from the united table and the cached fits; the interaction contrast without the five high libraries against all 4,410 regroupings; the acute heat shock by ploidy; how much of a pooled list lies above 20 points inside each cohort under relabelling; sets of five heat libraries outside the range of the other 17; mitochondrial CpGs (methylation and relative coverage); Figure 1B on counts before deduplication; a spike-in of a common ploidy difference | I-01, I-02, I-04, I-07, I-08, I-29, I-49 |
| `R1` to `R6` | Code reviews of the pipeline steps. The reviewers wrote no scripts worth keeping; their findings and the second agents' verdicts are in `../results/audit-verification.json`, and the confirmed ones are sorted by area in `../scripts/code-pass-inputs/` | the code pass (not yet run) |

The agents' full reports are in [`../results/audit-verification.json`](../results/audit-verification.json).

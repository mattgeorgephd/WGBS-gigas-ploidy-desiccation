# Scripts of the verifying agents

One folder per verification task of the workflow `audit-verification` (claims C1 to C8) and per code review (R1 to R6). Each folder holds the scripts the agent wrote and its small result tables, copied from the scratch directory of the audit session.

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
| `C7-legacy-go-semantics` | What Figure 5, Figure 9 and Table 1 count | I-02, I-07, I-20, I-21, I-45 |
| `C8-counts-and-correlations` | CpG counts, the genome total, the correlation range | I-16, I-44 |
| `R1` to `R6` | Code reviews of the pipeline steps (in progress) | to be added |

The agents' full reports are in [`../results/audit-verification.json`](../results/audit-verification.json).

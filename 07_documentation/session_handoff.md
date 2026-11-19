# Session handoff

How to continue this work from the repository alone, for a new session or a person. The working session rewrites the "State" and "In flight" sections whenever it pushes.

Last updated: 2026-11-19 04:15 UTC.

## If the session ends or its machine is replaced

The message "This session is out of disk, memory or time. Start a new session to get a fresh machine." means the current machine cannot continue. Nothing that was pushed is lost.

1. Start a new session in the same cloud environment (`WGBS-oyster-ploidy`) with this repository attached.
2. Give it this instruction: *Check out branch `reproducibility-audit`, read `07_documentation/session_handoff.md`, then `issues.md` and `07_documentation/reproducibility_plan.md`, and continue from "Next actions". Commit and push to `reproducibility-audit`.*
3. The new machine sets itself up: the environment's setup script installs R 4.3.3 and the packages (copy at `07_documentation/cloud_environment_setup.sh`; if R is missing, run that script).
4. Restore the inputs that are not in git (9.8 GB, about ten minutes): `Rscript 01_pipeline/01_fetch_inputs.R`. Optional: set `WGBS_CACHE_DIR` to a directory outside the repository, and `WGBS_CORES=4`.
5. Confirm that the machine reproduces the repository before changing anything: `Rscript run_pipeline.R --steps 6:10`, then `git status` (only PDF files and `session_info.txt` may show as changed), then `Rscript 06_diagnostics/check_legacy_outputs.R` (329 of 331 comparisons match).

What a new session cannot get from git: the download cache (step 4 restores it); a text copy of the manuscript (export the Google Doc linked in the root README); the agents' full transcripts (their returned results are in `07_documentation/agent_workflows/results/`).

## Rules of the project

- R only in the pipeline. Numbered folders, a README in each. One script per step in `01_pipeline/`, helpers in `02_R_functions/`, settings in `run_config.R`.
- GitHub first: any file that fits (50 MB after gzip) lives in the repository; gannet holds only what does not.
- The treatment is "heat" or "heat-treated". "Desiccation" appears only in legacy names. Legacy labels: `ronit` = heat-treated animals, `yaamini` = control animals, `both` = the 22-library comparison.
- Legacy mode reproduces the draft exactly; corrected analyses sit beside it. A number goes into `issues.md` or a README only after it was recomputed.
- Work on branch `reproducibility-audit`. No pull request unless the author asks. Commit and push after every finished unit of work; label snapshots of unreviewed work as such.
- Questions that only the authors can answer are collected in `07_documentation/reproducibility_plan.md`, section 7. Do not guess the answers.

## State

| Item | State |
|---|---|
| Audit | `issues.md`: 42 issues, six of them bearing on the main conclusions. |
| Steps 01 to 10 | Done. Legacy outputs reproduced: 329 of 331 checks in `06_diagnostics/legacy_check_results.csv`; the two differences come from the legacy files. |
| Step 11 (semantic similarity) and step 12 (qPCR) | Written by agents; committed as unreviewed snapshots; independent review and fixes in progress. `01_pipeline/README.md` still lists them as planned. |
| Step 13 (manuscript values) | The script works. `03_input_files/manuscript/manuscript_values_map.csv` holds three starter rows; the full inventory is being produced (see "In flight"). |
| Sensitivity analyses | `06_diagnostics/dml_sensitivity.R` written by an agent and tested on 20,000 CpGs only; unreviewed. The full model fits have not been run (about two hours on four cores: `Rscript 06_diagnostics/dml_sensitivity.R`). |
| Verification of the audit | In progress (see "In flight"). Five of nine claim reports received, none yet carried into `issues.md`; all are in `agent_workflows/results/audit-verification.json`. Summary below. |

### Verification reports received and still to be carried into `issues.md`

- **C1, incomplete GO universe (I-39): confirmed.** The cause is now shown in the file on gannet (a line cut at a byte offset that is a multiple of 4,096, followed by a line from the far side of the gap). Six small corrections to the issue text. The control-only universe has the same defect.
- **C2, file-reading fault (I-40): numbers confirmed.** Four corrections to the issue text; the shared set lost no transcripts, so its tests are unaffected by this fault. New defect to log: the legacy Fisher tests counted their totals in map lines (48,263) and their term counts in transcripts (27,469).
- **C3, ontology release (I-37, I-28): reproduction exact, one inference withdrawn.** The tables pin GO.db 3.13.0 only; topGO 2.34.0, 2.44.0 and 2.54.0 give identical files. The version strings in the Methods are verbatim from Venkataraman et al. 2022. To log: three "-new" tables in the archive come from a later run on GO.db 3.16.0; and the REVIGO figures (Figures 6 to 8) were drawn from all GO terms attached to DML-bearing transcripts, not from terms with p < 0.01, so most plotted and named terms are not enriched even by the draft's own test. This changes I-22, the header of step 09, and the premise of step 11.
- **C4a, the eligible-universe argument (I-41): central result holds, one conclusion must be withdrawn.** "No GO term enriched" was reproduced independently at the configured settings; the broadest term ("organelle") crosses q = 0.05 under some settings, so the statement needs its sensitivity. To withdraw: "DML are about twice as likely as other eligible CpGs to lie in gene bodies". Lists made by shuffling ploidy labels show the same odds ratio, so it is a property of the selection rule. The eligible background is fairer than all CpGs but not fair.
- **C5, feature odds ratios (I-41 table, I-14): every number reproduced independently.** Corrections: an odds ratio of 2.5 was worded as "twice as likely" (94.5% against 87.2% of loci); ploidy DML also stay nominally enriched in lncRNA and exon UTR, both through the one locus of I-09; the four chromosome percentages of Results 3.3 are not reproduced (I-25); the header of step 07 omits one output file. Matched on methylation level and between-animal variability, the gene-body odds ratio is 1.13 (ploidy) and 1.01 (heat), in agreement with C4a.
- **Step 12 (qPCR), build and review reports received.** The legacy analysis is reproduced (11 of 11 ANOVA tables). No effect survives correction across genes. The mitochondrial abundance statements of the Discussion are not supported; the HSP90 statements are partly supported, and the HSP90 heat effect depends on normalization to actin, which differs between groups. The review found one medium and seven low items; the fix stage is running.
- **Step 11 (semantic similarity), build report received.** It found independently what C3 found: Figures 6 to 8 were drawn from all GO terms annotated to DML-bearing transcripts. Of the 21 terms that Results 3.4 names as enriched, 2 have p < 0.01 in the draft's own test. It also reports that the archived final image of Figure 8 shows the ploidy plot, and that rrvgo's threshold means the opposite of the comment in `run_config.R` (REVIGO's 0.7 corresponds to 0.3 in rrvgo); both to be checked. Review is running.
- **New finding, checked a second time by the working session:** with ploidy labels shuffled within cohort (1,000 permutations), a median of 751 CpGs pass the 20-point rule (5th to 95th percentile 514 to 1,316) against 1,040 observed; 12% of permutations reach the observed count (14% with the coverage filter). The q-value filter removes none of the CpGs above 20 points. The number of ploidy DML is therefore not distinguishable from chance differences among animals. This belongs with I-06 and bears on every ploidy result. Script: scratch only so far; to be added to `06_diagnostics/dml_sensitivity.R` as a permutation analysis.

## In flight when this file was last updated

Four agent runs (scripts and partial results in `07_documentation/agent_workflows/`). If the machine is lost, the finished agents' results are in `results/`; the rest must be run again.

| Run | Finished | Still to come | If lost |
|---|---|---|---|
| `audit-verification` | claims C1, C2, C3, C4a, C5 | claims C4b, C6, C7, C8 (eligible-universe result by independent re-implementation; coverage filter and the 2.2 kb locus; legacy GO definitions and Table 1; CpG counts and correlations), six code reviews, six verifications | Run the script again; or check by hand the claims that matter most to the authors (C4, C5, C6). |
| `build-steps-11-12` | build of both steps; review of step 12 | review of step 11; fix for both | The snapshot of both scripts is in git. Review each (run twice, compare outputs, recompute key numbers), then mark done in `01_pipeline/README.md`. |
| `manuscript-values`, first batch (Methods, Results 3.1 to 3.3, figure legends) | extraction and mapping for all four parts | verification for four | Extracted statements and finished mappings are in `results/manuscript-values-batch1.json`. Run the script for the unfinished parts. |
| `build-dml-sensitivity` | none | build, review, fix | The snapshot of the script is in git. Review it, then run it in full. |

## Next actions, in order

1. Carry the verification reports into the repository: correct `issues.md` where a report found an error, fix confirmed code defects, rerun the steps they touch, rerun `06_diagnostics/check_legacy_outputs.R`.
2. Finish steps 11 and 12: apply confirmed review findings; add their rows to `01_pipeline/README.md` and `04_output/README.md`; update `issues.md` (I-22 for step 11, I-29 for step 12).
3. Manuscript values: run `manuscript-values.js` for the remaining parts (`results34_35`, `discussion`), which need the outputs of steps 11 and 12. Assemble all verified rows into `03_input_files/manuscript/manuscript_values_map.csv` (write the CSV with a CSV library so that quotes inside expressions are escaped), run step 13, and read every row that "differs".
4. Sensitivity analyses: after review, run `06_diagnostics/dml_sensitivity.R` in full, and write the results into `issues.md` (I-02 interaction, I-04 acute heat-shock libraries, I-06 overdispersion, I-07 ploidy within controls).
5. Issue I-15: the genomic distribution of CpGs at 1x in Results 3.1 has no pipeline output yet. Extend step 04 to write the overlap of its CpG sets with the feature tracks (the draft's numbers appear to come from the heat-only 10x set), then rerun step 04 (8 minutes, 8 GB).
6. `renv.lock` for the audit environment; a full run from a clean clone; final pass over the READMEs; remove obsolete items from the "still to be written" list in `06_diagnostics/README.md`.
7. Put the open decisions to the authors (`reproducibility_plan.md`, section 7: thirteen items).

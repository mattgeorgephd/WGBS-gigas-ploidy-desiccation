# Session handoff

How to continue this work from the repository alone, for a new session or a person. The working session rewrites the "State" and "In flight" sections whenever it pushes.

Last updated: 2026-11-19 03:25 UTC.

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
| Verification of the audit | In progress (see "In flight"). Reports received so far, not yet carried into `issues.md`: claim C1 (incomplete GO universe, I-39) confirmed, with the cause now shown directly in the file on gannet and six small corrections; claim C2 (file-reading fault, I-40) confirmed in its numbers, with four corrections to the issue text and one new defect to log (the legacy Fisher tests counted their totals in map lines and their term counts in transcripts). Both reports are in `agent_workflows/results/audit-verification.json`. |

## In flight when this file was last updated

Four agent runs (scripts and partial results in `07_documentation/agent_workflows/`). If the machine is lost, the finished agents' results are in `results/`; the rest must be run again.

| Run | Finished | Still to come | If lost |
|---|---|---|---|
| `audit-verification` | claims C1, C2 | claims C3 to C8 (ontology release; eligible-universe argument, by theory and by independent re-implementation; feature odds ratios; coverage filter and the 2.2 kb locus; legacy GO definitions and Table 1; CpG counts and correlations), six code reviews, six verifications | Run the script again; or check by hand the claims that matter most to the authors (C4, C5, C6). |
| `build-steps-11-12` | none | build, review, fix for both steps | The snapshot of both scripts is in git. Review each (run twice, compare outputs, recompute key numbers), then mark done in `01_pipeline/README.md`. |
| `manuscript-values`, first batch (Methods, Results 3.1 to 3.3, figure legends) | extraction for all four parts; mapping for Results 3.1 | mapping for three parts; verification for four | Extracted statements and finished mappings are in `results/manuscript-values-batch1.json`. Run the script for the unfinished parts. |
| `build-dml-sensitivity` | none | build, review, fix | The snapshot of the script is in git. Review it, then run it in full. |

## Next actions, in order

1. Carry the verification reports into the repository: correct `issues.md` where a report found an error, fix confirmed code defects, rerun the steps they touch, rerun `06_diagnostics/check_legacy_outputs.R`.
2. Finish steps 11 and 12: apply confirmed review findings; add their rows to `01_pipeline/README.md` and `04_output/README.md`; update `issues.md` (I-22 for step 11, I-29 for step 12).
3. Manuscript values: run `manuscript-values.js` for the remaining parts (`results34_35`, `discussion`), which need the outputs of steps 11 and 12. Assemble all verified rows into `03_input_files/manuscript/manuscript_values_map.csv` (write the CSV with a CSV library so that quotes inside expressions are escaped), run step 13, and read every row that "differs".
4. Sensitivity analyses: after review, run `06_diagnostics/dml_sensitivity.R` in full, and write the results into `issues.md` (I-02 interaction, I-04 acute heat-shock libraries, I-06 overdispersion, I-07 ploidy within controls).
5. Issue I-15: the genomic distribution of CpGs at 1x in Results 3.1 has no pipeline output yet. Extend step 04 to write the overlap of its CpG sets with the feature tracks (the draft's numbers appear to come from the heat-only 10x set), then rerun step 04 (8 minutes, 8 GB).
6. `renv.lock` for the audit environment; a full run from a clean clone; final pass over the READMEs; remove obsolete items from the "still to be written" list in `06_diagnostics/README.md`.
7. Put the open decisions to the authors (`reproducibility_plan.md`, section 7: thirteen items).

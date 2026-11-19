# agent_workflows

Records of the multi-agent runs used in this audit: the scripts that defined them and the results the agents returned. They are kept for two reasons. A machine reset loses everything that is not in git, and these runs take hours. And the verification of the audit's claims should itself be open to inspection.

Nothing here is part of the analysis pipeline, and no pipeline step reads these files.

## scripts/

Workflow scripts for the Workflow tool of Claude Code (JavaScript). Each defines the agents, their briefs, and how results pass between them. They are as launched, apart from one sentence reworded in `build-dml-sensitivity.js`.

| Script | What it does |
|---|---|
| `audit-verification.js` | Nine agents each try to refute one major claim of `issues.md` from the data. Six reviewers examine areas of the pipeline code, one of them by running steps 06 to 10 from a clean clone. Each reviewer's findings are re-examined by a second agent. Read-only on the repository. |
| `build-steps-11-12.js` | Builds pipeline steps 11 (semantic similarity) and 12 (qPCR): one agent writes the step, a second reviews it independently, a third fixes what the review confirmed. |
| `manuscript-values.js` | For a list of manuscript sections given as `args`: one agent lists every quantitative statement, a second ties each to an R expression over the pipeline outputs and tests it, a third recomputes each value independently and looks for missed statements. The rows it returns become `03_input_files/manuscript/manuscript_values_map.csv`. |
| `build-dml-sensitivity.js` | Builds `06_diagnostics/dml_sensitivity.R` (overdispersion, tested interaction, without the four acute heat-shock libraries, within-cohort contrasts) on a subset of CpGs, with independent review and fixes. The full model fits are run separately. |
| `eval_expressions.R` | Evaluates candidate expressions for the manuscript values map the way step 13 does; used by the agents of `manuscript-values.js`. |

To run one again in a new session, edit the constants at the top of the script: `SCRATCH` (any scratch directory) and `MANUSCRIPT` (a plain-text copy of the manuscript draft with one paragraph per line, exported from the Google Doc linked in the root README; `manuscript-values.js` refers to its line numbers, so the section line ranges in that script must be checked against the copy). `manuscript-values.js` also expects a catalog of output tables at `SCRATCH/output_catalog.txt` (one line per table of `04_output/*/tables` with its column names) and `eval_expressions.R` in `SCRATCH`.

## results/

One JSON file per run, holding what each finished agent returned (`agent`, `phase`, `result`) and the names of agents that had not finished when the file was written. The files are rewritten as agents finish. They are raw reports: a claim marked "confirmed" or a finding marked "high" is one agent's conclusion until it has been checked and carried into `issues.md`, the code, or the manuscript values map. Once a run has been integrated, its file stays as the record of what was reported.

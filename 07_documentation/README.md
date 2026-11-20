# 07_documentation

Project documentation that is not code and not data: the plan for the reproducible pipeline, the audit behind it, and records of processing that happened before this repository.

For the list of open problems and decisions, see [`issues.md`](../issues.md) at the repository root.

## Files

| Path | Contents |
|---|---|
| `session_handoff.md` | How to continue the work in a new session: recovery steps, project rules, current state, work in flight, next actions. Rewritten whenever the working session pushes. |
| `reproducibility_plan.md` | The plan: goal, rules (R only, GitHub first for data), target folder layout, the thirteen pipeline steps, phases, validation targets, and the decisions needed from the authors. |
| `toggles_migration/` | One-time scripts of the settings layout of 2026-11-20: `rewrite_maps.py` (the two manuscript maps: `expression_corrected` became the template `expression_optional`; the items map got `optional` and `comparisons` rows), `split_tables.R` (the tables that listed the legacy and the corrected mode were split into the paper's copy, the copy of `04_output/not_adopted/corrected/` and the side-by-side copy of `not_adopted/comparisons/`), and `settings_table.md` (the table of the analyses and switches that the READMEs repeat). Not run by the pipeline. Deviation from the build spec and from the R-only rule of the repository: `rewrite_maps.py` is Python (run with `python3 -I`), because the maps must keep their CRLF line endings and quoting byte for byte, which the csv module writes and its round-trip check verifies; the review fixes of the values map of 2026-11-20 (nine templates of `expression_optional`) were made the same way. |
| `agent_workflows/` | Scripts and returned results of the multi-agent runs used to verify the audit and to build parts of the pipeline. |
| `audit/gannet_manifest_2026-11-18.tsv` | Every file in the project directory on gannet as crawled on 2026-11-18 (843 files, 217 GB): path, modification time, size. |
| `audit/elisa_qpcr_audit.md` | Detailed audit of the ELISA, the percent CpG methylation summary, and the qPCR material: computational path, recomputed statistics, per-sample values, script inventory. |
| `upstream_processing/` | Alignment job scripts copied from gannet, with a README on what they establish about how the coverage files were produced. |
| `protocols/` | Kit manual for the MethylFlash global 5-mC ELISA (Epigentek P-1030), with margin notes from the 2018 pilot. |
| `sequencing_reports/` | Zymo Research service reports for the WGBS libraries of project zr3534. |
| `cloud_environment_setup.sh` | Setup script for a Linux container (Ubuntu 24.04): installs R 4.3.3 and the Bioconductor 3.18 packages the pipeline uses. Stored as the setup script of the `WGBS-oyster-ploidy` cloud environment. |

## Conventions

- The treatment is called "heat" or "heat-treated" in everything written from 2026-11-18 on. Legacy files call it "desiccation".
- Legacy dataset labels: `ronit` = heat-treated animals only, `yaamini` = control animals only, `both` = the 22-library comparison in the manuscript.

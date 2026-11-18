# 06_diagnostics

Checks that the pipeline reproduces what it should. These scripts are not part of the analysis; they compare its outputs against three references and report every difference.

## The three comparisons

1. **Regenerated against legacy outputs.** DML lists, per-chromosome counts, feature overlap counts, and topGO tables in `08_archive/` (and the few legacy tables that exist only on gannet) are compared row by row with the outputs of the pipeline run in legacy mode.
2. **Regenerated against the manuscript.** Each number quoted in the text is compared with the regenerated value; the result is `05_manuscript/manuscript_values.csv`.
3. **Legacy mode against corrected mode.** Where the plan defines a corrected analysis (the upper coverage filter, the interaction test, the unit and background of the GO test), the two sets of results are tabulated side by side.

## Targets for legacy mode

The values that must be matched exactly are listed in `07_documentation/reproducibility_plan.md`, section 6. One is already confirmed outside the pipeline: 480,579 CpGs with at least 10x coverage in all 22 libraries, reproduced from the coverage files on gannet.

## Planned files

| File | Purpose |
|---|---|
| `check_legacy_outputs.R` | Comparison 1 |
| `check_modes.R` | Comparison 3 |
| `legacy_reference/` | Small legacy tables that exist only on gannet, copied here so the checks run from the repository |

No scripts exist here yet; each is added with the pipeline step it checks.

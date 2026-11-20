# Action plan: remove the F test of the DML model, keep the correction with the chi-squared test only

Scope: the instruction of M. George of 2026-11-20 ("remove all mention of f-test"; the overdispersion option is `overdispersion = "MN"`, `test = "Chisq"` only, as in Venkataraman et al. 2022). The records under `07_documentation/agent_workflows/` stay as they are (author, 00:45 UTC; plan line 164). Nothing in this plan has been applied. Every repository path is relative to `/home/user/WGBS-gigas-ploidy-desiccation` (branch `reproducibility-audit`).

Inputs merged: the term sweep (`../F1-terms/hits.csv`, 2,305 lines) and the dependency trace (`../F2-dependencies/statements.csv`, 67 statements; `analyses.csv`; `derived/`). The corrected and extended statement table of this pass is `statements_checked.csv` (83 rows) in this folder.

State of the files read: HEAD 471764a, then 671feca at the end of this pass (that commit touched only 04_output/04_cpg_coverage/ and agent_workflows/results/scope-author-changes.json). `git diff --stat 435f971 HEAD` on issues.md, README.md, 06_diagnostics/, the plan, the handoff, 03_input_files/manuscript/ and 05_manuscript/ gives no change, so the line numbers of both inventories (taken at 435f971) still hold. md5 of issues.md edc445fea497749b86f65c921c022dba, as the trace recorded. The working tree differs from HEAD only in three files under `04_output/04_cpg_coverage/`. The hardening pass is rerunning the pipeline now (step 05 was running at the time of reading), so line numbers and outputs can still move. Check them again before editing.

## 0. Cross-check of the two inventories (task step 1)

What was checked, and how:
- Every hit of the sweep in the two F classes ("named explicitly" and "test not named") in the 8 text and code files was matched against the line ranges of the trace (`code/crosscheck.py`). 40 lines were not covered. Each was opened and is now covered in `statements_checked.csv` or explained below.
- Every statement of the trace was located by its text in its file (`code/locate.py`). 61 of 67 were found on the stated line. The other 6 were checked by eye: the line was right and only the paraphrase failed to match (S04, S12, S15, S18, S33, S35).
- Lines confirmed by opening the files (well over 40):
  - issues.md: 24, 39, 122, 123, 166-170, 289, 291, 295, 299, 300, 302, 304, 306, 307, 311, 313-319, 321, 323, 327-337, 339-341, 346-348, 350-353, 356-358, 361, 363, 365, 367, 369, 395, 418, 1111, 1191, 1197, 1199, 1202, 1211, 1213, 1215, 1225, 1226
  - README.md: 13
  - 06_diagnostics/README.md: 11, 15, 44-49, 80-145, 155-161
  - plan: 100, 114, 141, 153, 164, 166, 169-171
  - handoff: 23, 40, 53, 71, 95, 97
  - values map: rows 283, 337, 347, 462, 466, 476, 496, 521, 527, 529, 533, 537
  - 06_diagnostics/dml_sensitivity.R: 1-185, 366-475, 530-600, 985-1115
  - 06_diagnostics/library_methylation_level.R: 48-76, 240-252, 400-470, 530-540
- The coverage was also checked in reverse, from trace statements to sweep hits (inline Python on both CSVs). Three statements of the trace have no matching F-class line in the sweep: S06 (issues.md 170, missed), S51 (plan 100, classed "unrelated") and S62 (values-map rows, classed "unrelated"). In each, the sweep's class is wrong.

Wrong line numbers in the trace, now fixed in `statements_checked.csv`:

| Statement | Trace says | Correct | Why it matters |
|---|---|---|---|
| S40-S50 (06_diagnostics/README.md) | 48, 49, 51, 54, 76, 77, 91-101, 105, 106, 108, 126 | 83, 84, 86, 89, 111, 112, 126-136, 140, 141, 143, 161 | Every line is off by 35. Lines 48-49 of that file are the ELISA and %mCpG two-way ANOVA rows ("interaction F"). An edit made at the trace's line numbers would hit an F statistic of the draft. |
| S09 (issues.md) | 300 | 299-300 | The benchmark percentages for the F test are on 299 |
| S67 (library_methylation_level.R) | 50-62 | 54-62 | Lines 50-53 belong to part 7 |

Statements missing from the trace (added as X01-X16 in `statements_checked.csv`):
- issues.md:
  - 122 (I-01: "138 of the CpGs with q < 0.01 in the corrected heat contrast are mitochondrial")
  - 170 (also missed by the sweep)
  - 299 (benchmark)
  - 307 (I-06 update: "what a corrected list means (item 7)")
  - 311 (item 2 heading; it holds)
  - 1213 and 1225 (corrections-table rows on the benchmark)
- Values map: rows R3.5-049 (line 347), D4.9-002 (521) and D3.5x-001 (537). They carry the same "none when variation among animals is modelled" as the four rows of S62 and rest on the A2 F fit. Neither inventory lists them.
- Handoff: 40, 53, 71. Plan: 164. dml_sensitivity.R: 388-391 (the chi-squared rows are defined by reference to the F ids) and 1005 (a comment).

Sweep hits that need no action: issues.md 289, 304, 657, 1215 and 1226 (chi-squared wording that holds); 06_diagnostics/README.md 11 and 124 (generic; the count of 22 fits stays); dml_sensitivity.R 95-96, 376-377, 472, 551 and 567 (generic code that passes `spec$test`); library_methylation_level.R 536; plan 114.

Spot checks of the trace's numbers, independent of its code:
- `code/spotcheck_derivation.R` derives chi-squared p-values from the cached F fits. It reproduces 10 DML (A4 controls), 1,093 (A4 heat-treated) and 14 (A1 ploidy). At dispersion not above 1 it gives 2, 172 and 1. Run with ulimit -v 1.5 GB.
- Inline Python on the stored relabelling tables gives the following.
  - Pooled ploidy, from `07_relabel_level_pooled_a.csv` and `_b.csv`:
    - chi-squared: 86, 69 and 61 of 300 relabellings reach the observed shares; mean 7.19%, 2.16% and 0.47%; 5th to 95th percentile below 0.05 is 4.6% to 11.7%.
    - F: 78, 95 and 119 of 300.
  - Heat share: 12 of 300 reach 13.2%.
  - Heat cohort (`08_relabel_level_heat_experiment.csv`): 0 of 125. Controls (`08_relabel_level_control_experiment.csv`): 39, 28 and 11 of 150.
- The shares and statistics in `06_diagnostics/dml_sensitivity/dml_counts.csv` agree with the trace (cat of the file).

## (a) Code changes

### 06_diagnostics/dml_sensitivity.R (1,115 lines)

Registry. Replace each F entry with its chi-squared twin and keep the 22 fits (16 legacy, 6 corrected):

| Lines | Now | Change to |
|---|---|---|
| 384-387 | `A1_ploidy_mn_f`, `A1_heat_mn_f` ("MN", "F", both_modes) | delete |
| 388-391 | `A1_ploidy_mn_chisq`, `A1_heat_mn_chisq` (legacy only); question "As A1_ploidy_mn_f with the chi-squared test, the less stringent form of the correction" | keep both ids; set `modes = both_modes`; question: "Which ploidy [heat] DML remain when the variation among animals is estimated (overdispersion "MN", chi-squared test, the setting of Venkataraman et al. 2022)? Issue I-06." |
| 394-395 | `A2_interaction_mn_f` ("MN", "F", both_modes) | `A2_interaction_mn_chisq`, "MN", "Chisq", both_modes; "The same interaction test with the variation among animals estimated. Issues I-02 and I-06." |
| 402-403 | `A4_ploidy_in_control_mn_f` ("MN", "F") | `A4_ploidy_in_control_mn_chisq`, "MN", "Chisq"; "As A4_ploidy_in_control with the variation among animals estimated." |
| 406-407 | `A4_ploidy_in_heat_mn_f` | `A4_ploidy_in_heat_mn_chisq`, "MN", "Chisq" (same wording) |

The cache signature (lines 548-551) holds `analysis_id`, mode, input md5, libraries, indicator, covariates, overdispersion and test. It does not hold `modes` or the question. So the cached legacy fits of `A1_*_mn_chisq` stay valid when their modes widen, and only 6 fits are new (section b).

Branches:
- 373 (registry comment): `test  "Chisq" (the only test offered; write it explicitly, methylKit's default is another test)`.
- 425-427: test validation becomes `spec$test != "Chisq"` → bad("only the chi-squared test is offered (author, 2026-11-20)"). The check at 426-427 (none needs Chisq) then folds into it.
- 452 (figure label): `if (spec$overdispersion == "MN") "; MN, chi-squared test" else ""`. The uniqueness check at 465 still passes: one MN row per family and contrast.
- 457-459 (`reference_distribution`): for MN, "chi-squared with 1 degree of freedom, deviance divided by the dispersion estimate (estimates not above 1 set to 1)"; otherwise "chi-squared with 1 degree of freedom".
- 997-999: `statistic_needed_for_a_single_call = stats::qchisq(p_single, 1, lower.tail = FALSE)`, without the branch. The column then holds 31.4 (legacy) or 31.3 (corrected) for every row. Keep the column (the README explains it) or drop it; that is the code pass's choice.
- 1005 comment (optional): "Analyses that differ only in the overdispersion setting share their effect size."
- 567 already passes `test = spec$test`. Keep `"Chisq"` explicit in the registry, because `formals(methylKit::calculateDiffMeth)$test` lists "F" first (trace, verified).

Header text:
- 13-17 (A1, A2): "A1 overdispersion (issue I-06): the same two contrasts with the variation among animals estimated ("MN") and the chi-squared test, the setting of Venkataraman et al. 2022. A2 a tested ploidy by heat interaction (issue I-02), without and with "MN", cross-tabulated with ...".
- 31-32: "Every analysis runs in mode "legacy"; A0, A1 and A2 also run in mode "corrected"."
- 100-105: end the "MN" note at "... referred to chi-squared with 1 degree of freedom. The script passes test = "Chisq" explicitly, because calculateDiffMeth defaults to another test when the correction is on." Delete "or to F with 1 and the residual degrees of freedom (test "F"). The F distribution is used only where the dispersion is above 1."
- 107-122 ("Which of the two tests") becomes a note "The test". Keep the benchmark's rating of the default test (lowest specificity, best overall accuracy) and of the chi-squared test with the correction ("a good compromise between stringent and relaxed methods"), and the three samples per group. Delete the F clause (113-114) and the stats::anova.glm sentence (117-118). Keep: "Under relabelled ploidy labels it does not hold its level on these libraries (issues.md, I-06, item 4), and where the dispersion estimate is not above 1 the correction is not applied (item 7)." Whether the description of the benchmark may name the F test at all is an author decision (open points).
- 170-179, note (3): delete the sentence "For the F test this statistic rises steeply ... not the absence of differences" (lines 176-179). Keep the definition.

### 06_diagnostics/library_methylation_level.R (541 lines), part 8

| Line | Change |
|---|---|
| 54-57 | "8. The CpGs that the sensitivity analyses call with the overdispersion correction (the lists dml_<analysis>_mn_chisq_<mode>.csv): percent methylation ..." Delete "and the F test" and "the two lists made with the chi-squared test are not read". |
| 73 | `06_diagnostics/dml_sensitivity/dml_<analysis>_mn_chisq_<mode>.csv   (part 8; skipped ...)` |
| 247 | `corrected_lists <- c("A1_ploidy_mn_chisq", "A1_heat_mn_chisq", "A2_interaction_mn_chisq", "A4_ploidy_in_control_mn_chisq", "A4_ploidy_in_heat_mn_chisq")` |
| 415 | message "no dml_*_mn_chisq_<mode>.csv under ..." |
| 536 | no change |

The code of part 8 (411-445) works on positions only and needs no other change. The output grows from 145 + 9 to about 1,418 + 398 CpGs (trace, derived). The file names `corrected_dml_*` can stay; renaming them would break the citations in issues.md I-08 and I-49.

### Other readers of the F outputs

None. `git grep` for `_mn_f|dml_counts|permutation_reference|interaction_vs_shared|corrected_dml_` over the repository (excluding agent_workflows and 08_archive) finds only the two scripts above. `check_legacy_outputs.R` reads step 05's own `dml_counts_by_threshold`. The values map reads only uncorrected rows of the sensitivity tables: `dml_counts.csv` row `A4_ploidy_in_control` and `interaction_vs_shared_corrected.csv` row `A2_interaction` (lines 157, 271, 283, 324, 346, 433, 447, 520, 525 and others). The regenerated tables must keep those ids and values.

### Interface with the toggle design (lane 3)

When step 05 gets its overdispersion option, it should take one switch (correction on or off) with no choice of test, and pass `overdispersion = "MN", test = "Chisq"` explicitly. Nothing in 01_pipeline/, 02_R_functions/, run_config.R or run_pipeline.R mentions the F test today (sweep, verified).

## (b) Fits to run

Six new methylKit fits. The minutes are those of the F fits they replace (field `minutes` of the cached fits, from the trace), at `WGBS_CORES=3` on a shared machine. The chi-squared test costs the same: the same glm and the same dispersion estimate.

| New fit | Mode | Replaces | Minutes |
|---|---|---|---|
| A1_ploidy_mn_chisq | corrected | A1_ploidy_mn_f corrected | 10.2 |
| A1_heat_mn_chisq | corrected | A1_heat_mn_f corrected | 9.6 |
| A2_interaction_mn_chisq | legacy | A2_interaction_mn_f legacy | 8.7 |
| A2_interaction_mn_chisq | corrected | A2_interaction_mn_f corrected | 15.0 |
| A4_ploidy_in_control_mn_chisq | legacy | A4_ploidy_in_control_mn_f | 6.4 |
| A4_ploidy_in_heat_mn_chisq | legacy | A4_ploidy_in_heat_mn_f | 5.7 |
| Total | | | about 56 |

The cheapest way to get them:
1. Wait until the hardening pass has finished and is integrated. Then confirm that `md5sum 04_output/05_differential_methylation/tables/united_cpg_{legacy,corrected}.csv.gz` still gives `8e7f5d67bc3fb7306ddef29a3894110f` and `ebec40d3e333f787a091f0ae529da9b6`. These are the `united_md5` of the cache signatures, read from `fit_A4_ploidy_in_heat_mn_f_legacy.rds` and `fit_A2_interaction_mn_f_corrected.rds`. If either changed, all 22 fits are refitted, about 3 hours at 3 cores (script header).
2. Edit the registry and run test mode once: `WGBS_SENSITIVITY_SUBSET=20000 Rscript 06_diagnostics/dml_sensitivity.R`, 10 to 15 minutes on 2 cores.
3. Full run: `WGBS_CORES=3 Rscript 06_diagnostics/dml_sensitivity.R`. The 16 unchanged fits come from the cache. Expect about 56 minutes of fitting plus about 10 for the tables and the permutation reference, so about 66 minutes. It must not overlap with the other workflow's 10 GB jobs.
4. Check the new fits against the scratch derivation. It predicts every count exactly: `../F2-dependencies/derived/chisq_counts_derived.csv` and `chisq_dml_lists_derived.csv`. Compare n_dml, n_significant, n_in_primary and the shares, and require identical DML sets. The validated pairs agreed to 2e-8 relative in p.
5. `Rscript 06_diagnostics/library_methylation_level.R`, about 2 minutes and 1.5 GB.
6. After the values-map edits, rerun step 13 only (`01_pipeline/13_manuscript_values.R`) to rewrite `05_manuscript/manuscript_values.csv`.

A stop-gap that is not a substitute: `../F2-dependencies/code/derive_chisq_from_cache.R` gives all the numbers from the cache in about 5 minutes. The repository should hold real fits, so use it only as the check in step 4. It needs the cached F fits, so delete them only after step 4.

New relabelling computations: only if the author chooses "recompute" over "withdraw" (open points). The F-only statements that cite them cannot be rewritten in full without them. Each is new code outside the records, 1 to 4 minutes per the trace:
- calls under relabelling at each list's p threshold (scratch 15_ and 15b_ with `test = "Chisq"`): I-06 item 7, I-07, the README tail
- the interaction share under relabelling (`int_mn_chisq` was never run): I-02, item 2
- the split by the five high libraries (18_) and by the two high control libraries (19_): item 2, I-07
- the corrected p of the DML under relabelling (11_): I-01, item 7 "low power"

The scope rule of 2026-11-19 (handoff line 23: "Do not start sensitivity analyses, new diagnostics ...") weighs against adding them. The instruction of 2026-11-20 (handoff 95 (iii): "recompute it with the chi-squared test or withdraw it") allows either.

## (c) Output files

Delete (after check (b)4):
- Tracked, `06_diagnostics/dml_sensitivity/`: `dml_A1_ploidy_mn_f_legacy.csv`, `dml_A1_ploidy_mn_f_corrected.csv`, `dml_A1_heat_mn_f_legacy.csv`, `dml_A1_heat_mn_f_corrected.csv`, `dml_A2_interaction_mn_f_legacy.csv`, `dml_A2_interaction_mn_f_corrected.csv`, `dml_A4_ploidy_in_control_mn_f_legacy.csv`, `dml_A4_ploidy_in_heat_mn_f_legacy.csv` (8 files). The script's stale-file note (line 1060) lists them until they are gone.
- Git-ignored: `06_diagnostics/dml_sensitivity/cache/fit_*_mn_f_*.rds` (8 files, 49 MB per the sweep), plus `test_subset/dml_*_mn_f_*.csv` (8) and `test_subset/cache/fit_*_mn_f_*.rds` (8). `find 06_diagnostics -iname '*_mn_f_*'` lists all 32 today.

Regenerated by the full run of `dml_sensitivity.R`:
- `analyses.csv`: the test, reference_distribution and figure_label columns.
- `dml_counts.csv`: 8 F rows replaced by `_mn_chisq` rows, giving 22 rows.
- `interaction_vs_shared_legacy.csv` and `_corrected.csv`: 6 rows each, `A2_interaction_mn_f` becomes `A2_interaction_mn_chisq`.
- `permutation_reference.csv`: only the `analyses` column changes, for example "A0_ploidy; A1_ploidy_mn_chisq". The counts need no model.
- `figure_dml_sensitivity.png` and `.pdf`.
- `equivalence_with_step05.csv`, `dispersion_summary.csv` and `session_info.txt`, rewritten with the same content apart from session info.
- 6 new lists: `dml_A1_{ploidy,heat}_mn_chisq_corrected.csv`, `dml_A2_interaction_mn_chisq_{legacy,corrected}.csv`, `dml_A4_ploidy_in_{control,heat}_mn_chisq_legacy.csv`.

Regenerated by `library_methylation_level.R`: `corrected_dml_library_values.csv` and `corrected_dml_group_apart.csv`. The other 11 tables of that folder should come out byte-identical; check with `git diff --stat`.

Then step 13: `05_manuscript/manuscript_values.csv` (rows D4.3-016, D4.9-008, D4.9-010, D4.9-014 at lines 476, 527, 529, 533).

## (d) Text changes

Actions used below:
- **R**: replace with chi-squared numbers that exist now, in a repository table.
- **RF**: replace after the new fit. The number is predicted exactly by the scratch derivation.
- **RS**: replace after the script writes a new column. No repository table holds the number now.
- **W**: withdraw.
- **N**: needs a new relabelling computation (withdraw unless the author chooses to recompute).
- **A**: author decision.

Numbers in quotes are proposed wording. `claims.csv` in this folder holds 35 of them with their sources (section h).

### issues.md

| Line | Where | Action | Proposed wording / note |
|---|---|---|---|
| 24 | Updates, fourth pass (dated) | R + add dated entry | Reword "... added that methylKit's corrected test leaves the test uncorrected at CpGs whose dispersion estimate is not above 1, so that its lists are not lists of established differences" to "... added that methylKit's correction is not applied at CpGs whose dispersion estimate is not above 1 (I-06, item 7)". Add a dated entry: "2026-11-2x: the sensitivity analysis fits the correction with the chi-squared test only (author, 2026-11-20); I-01, I-02, I-06, I-07, I-08 and I-49 rewritten on its results." |
| 39 | Summary row I-06 | R/RF | "... methylKit's corrected test (chi-squared, the setting of Venkataraman et al. 2022) leaves 14 ploidy and 185 heat CpGs (157 with the coverage filter) and does not hold its level with relabelled libraries". This replaces "2 ploidy and 16 heat CpGs (6 with the coverage filter), and relabelled libraries get such calls too". |
| 122 | I-01, mitochondrial paragraph | RF/RS | "it is why 275 of the 1,124 CpGs with q < 0.01 in the corrected heat contrast are mitochondrial (I-06, item 6)". The script must write q < 0.01 by coverage-filter status (already on the code-pass list, handoff 40). |
| 123 | I-01, cohort contrast | R + RF + N | "(13% of the eligible CpGs below 0.05, where relabelled ploidy labels give 7.2% on average; 12 of 300 relabellings of the ploidy label ... reach that share) and 185 CpGs pass the DML thresholds, 157 with the coverage filter (I-06)". The sentence "Of the draft's 2,351 heat DML, 82% have a corrected p below 0.05, a share that 7 of 150 lists ... reach" becomes "88% (2,062) have a chi-squared corrected p below 0.05" plus N for the relabelling clause, or W. "which CpGs differ cannot be named" now rests on item 4 (A). |
| 167 | I-02, With the correction | RF | "*With the correction* it calls 242 CpGs, 231 with the coverage filter, none of them among the 200 (none of the 159); at 5 of the 242 methylKit applied no correction because the dispersion estimate is not above 1 (I-06, item 7)". The last clause is RS. |
| 168 | I-02, ploidy difference differs | RF + N/W | "With the correction, 23% of the eligible CpGs still have p < 0.05". The clause "a share that 1 of 300 relabellings ... reaches" is W, or N (`int_mn_chisq`). |
| 170 | I-02, Consequence for the options | RF | "... uncorrected, it lists 3.5% of all tested CpGs; corrected, it lists 242 CpGs (231 with the coverage filter), none of them among the draft's 200." The conclusion (no rescue of the set) holds. |
| 291 | I-06, Proposed resolution | R | "(Superseded ...: methylKit's corrected test is not a calibrated primary analysis on these libraries; items 4 and 7.)" |
| 295 | Rating | R | Add "(with the correction 14 of the 1,038 ploidy DML and 185 of the 2,351 heat DML remain; item 4)". |
| 299-300 | Correction paragraph, benchmark | A | Keeping the description of the cited article as it is means accepting the F test in a quotation. If it is withdrawn: on 299, "... against 4.6% to 6.5% for the chi-squared test with the correction and at most 1.0% for DSS and limma"; the bar counts lose "19 to 32 for the F test with the correction". On 300, drop the quotation "F-test with overdispersion has similar results to DSS" and shorten "limma, DSS and the F test with the correction can be applied ..." to "limma and DSS ...". That alters a quotation and needs a marked ellipsis. |
| 302 | Correction paragraph, second part | A + R | "(note "Which of the two tests") is confirmed on its four points" becomes "(note "The test") ... on three points" if the F point is withdrawn. "Items 2, 4 and 7 below measure the three tests" becomes "the two tests". The quoted help text "automatically switches to the F-test. The Chisq-test can be manually chosen" is A (quotation of the documentation). |
| 306 | Update opening | none | 22 fits, 16 legacy and 6 corrected, stays true |
| 307 | Update, second paragraph | R | "It changed several readings, among them that of item 7." |
| 311 | Item 2 heading | keep | Holds: 86 of 300 relabellings reach the chi-squared ploidy share. The line's last sentence ("A calibrated test ... gives 5%, 1% and 0.1%") should be followed by "The chi-squared form of the correction does not give these shares even where groups do not differ (item 4); read its shares against relabelled libraries." |
| 313-319 | Item 2 table | R/RF | Column header "With the correction (chi-squared test)". Rows: Ploidy "7.3%, 2.2%, 0.46%" (R); Heat "13%, 5.5%, 1.7%" (R); Interaction "23%, 12%, 4.1%" (RF); Controls "10%, 3.3%, 0.81%" (RF); Heat-treated "25%, 13%, 5.2%" (RF). |
| 321 | Item 2 paragraph | R + N/W | "With the ploidy labels reassigned within cohort (300 relabellings, each at 3,000 eligible CpGs), the uncorrected test gives 33%, 21% and 12% on average, the same as with the real labels, and the corrected test gives 7.2%, 2.2% and 0.47% on average (4.6% to 12% below 0.05 in nine relabellings of ten). With the correction, the ploidy contrast of the draft has the p-value distribution of relabelled libraries: 86 of the 300 relabellings give at least the share below 0.05 that the real labels give at the same CpGs (69 at 0.01, 61 at 0.001). The ploidy contrast among the controls is of the same kind: 39 of 150 regroupings of the 12 controls reach its share below 0.05 (28 at 0.01, 11 at 0.001). The ploidy contrast inside the heat cohort keeps an excess that relabelling does not give: none of the 125 other groupings of the ten heat-treated libraries reaches 25%. The cohort contrast cannot be relabelled; 12 of the 300 relabellings of the ploidy label reach its 13%, so ..." Three sentences are W or N: "7.6% ... 4.5%" (18_); "mean 6.6% over 70 ... 4.4% over 80" (19_); "1 of the 300 relabellings reaches the 20% of the interaction". Source of the numbers given: the stored tables `07_relabel_level_pooled_{a,b}.csv` and `08_relabel_level_{heat,control}_experiment.csv` (columns `mn_chisq1-3`) of S1. |
| 323 | Item 3, last sentence | R | "With the correction, 168 and 1,124 CpGs have q < 0.01." |
| 327-337 | Item 4 table | R/RF | One correction column, "With the correction (chi-squared test)": 14 (14); 185 (157); 242 (231); not run; not run; 10; 1,093; not run; not run. |
| 340 | Item 4, one-sided lists | RF | "101 of the 125 corrected calls" becomes "933 of the 1,093 corrected calls". |
| 341 | Item 4, chi-squared paragraph | W + R | Keep: "The chi-squared form of the correction does not hold its level with 19 residual degrees of freedom, because it refers a statistic with an estimated dispersion to the distribution that a known dispersion would give: with relabelled libraries it gives 7.2%, 2.2% and 0.47% of p-values below 0.05, 0.01 and 0.001 (`10_combine_relabel.log`). Its lists are therefore not lists with 1% false discoveries." Withdraw "an F variable ... (2.1e-5 against 2.1e-8)" and "Its 14 and 185 CpGs are not a middle course between the two other columns". |
| 346 | Item 6 | RF/RS | "Ploidy: 14 CpGs, the same in both modes. Heat: 185 CpGs, 157 in corrected mode. Of the 28 that drop, 2 are at CpGs that the filter removes and 26 lose q < 0.01 (q 0.0067 to 0.0098 in legacy mode), because 391 of the 1,124 CpGs with q < 0.01 in legacy mode are high-coverage CpGs that the filter removes (275 of them mitochondrial), and a q-value depends on how many other CpGs are significant. With the Benjamini-Hochberg adjustment in place of methylKit's SLIM the lists have 14 and 13 ploidy CpGs and 175 and 150 heat CpGs." Withdraw "The corrected heat list is therefore not stable at its threshold." The filter-status and Benjamini-Hochberg numbers are RS (script). |
| 347 | Item 6, ploidy pair | RF | "NC_047562.1:40,685,928 has q = 1.2e-8 in both modes, and NC_047565.1:35,030,701 differs by 20.2 points. Both are called inside the heat cohort (q = 1.9e-6 and 1.0e-9; 51.8 and 40.9 points) and neither among the controls (q = 0.44 and 0.79; 10.0 and 0.0 points): both are differences of the heat cohort alone. 10 of the 14 corrected ploidy CpGs are also called inside the heat cohort and none among the controls. The 10 CpGs called among the controls include 3 of the draft's 1,038; none of the 242 interaction CpGs is among its 200." Withdraw the four named control and interaction CpGs and "leaves the list under the Benjamini-Hochberg adjustment". |
| 348 | Item 6, NC_047565.1:35,030,701 | keep | Holds: in both chi-squared lists, dispersion 0.54, so the p-value is the uncorrected one. Re-cite the regenerated `corrected_dml_library_values.csv`. |
| 350 | Item 7, first paragraph | W + RF/RS | "With the chi-squared test a single CpG needs a statistic above 31.4 (p below 2.1e-8) to be called on its own, in every analysis (`dml_counts.csv`, column `statistic_needed_for_a_single_call`). The correction is not applied where the dispersion estimate is not above 1: methylKit then sets it to 1, and the corrected p-value equals the uncorrected one (`analyses.csv`, column `reference_distribution`). The estimate is not above 1 at 1.0% of the eligible CpGs with 22 libraries, 3.3% with the 12 controls and 6.1% with the 10 heat-treated libraries; such CpGs are 1 of the 14 ploidy CpGs, 2 of the 185 heat CpGs, 5 of the 242 interaction CpGs, 2 of the 10 among the controls and 172 of the 1,093 among the heat-treated." The counts are RS. Withdraw "F statistic above 84 ... above 474", "This holds only where ..." and "those CpGs supply most of the corrected lists ...". |
| 351 | Item 7, relabelling | N or W | Every number is F-specific: 0.14, 2.2 and 5.4 CpGs below 2.1e-8 per relabelling; medians 0, 2 and 10; "8 of 200"; "95 of 150"; "4 to 41, largest 69"; "80%, all, 99.9%". Proposed if withdrawn: "The corrected lists are not lists with 1% false discoveries, because the corrected test does not hold its level on these libraries (item 4)." |
| 352 | Item 7, third paragraph | R | "A corrected list is not a list of established differences, and a short one is not evidence that two groups are alike. These data cannot name the CpGs at which groups differ. Item 2, read against relabelled libraries, is the better guide to whether groups differ at all." |
| 353 | Item 7, low power | RF + N/W | "Of the draft's 1,038 ploidy DML, 862 (83%) have a corrected p below 0.05 and 116 (11%) below 0.001" (48, 4.6%, are below 1e-4). The relabelling comparison ("72% and 4.7% ... 27 and 48") is N or W, and so is the heading sentence "Low power does not explain the result". "A looser threshold gives no list" becomes "A looser threshold gives a list: with q < 0.05 the corrected test calls 37 ploidy CpGs above 20 points, with q < 0.1 it calls 52." The last two numbers are RS: the script must write counts by q threshold, or cite the cached fit. |
| 356 | What follows, Ploidy | RF | "of the 1,093 CpGs called there with the correction, 179 are among the 1,038". "its CpGs cannot be named" now rests on item 4 (A to confirm). |
| 357 | What follows, Heat | RF/RS | "... which CpGs differ cannot be named: 185 CpGs pass with the correction, 157 with the coverage filter, 175 and 150 with another adjustment for multiple testing, from a test that does not hold its level on these libraries (item 4)." |
| 358 | What follows, heat cohort | RF + N/W | "The two list sizes, 12,405 uncorrected and 1,093 corrected, do not count established loci: relabelled libraries give a median of 4,019.5 CpGs above 20 points (`permutation_reference.csv`), the corrected test does not hold its level (item 4), and 172 of the 1,093 were tested without the correction because their dispersion estimate is not above 1." "about ten corrected calls (item 7)" is N or W. |
| 361 | Options 1 | R | "... methylKit's corrected test is not by itself a calibrated primary analysis: its chi-squared form does not hold its level on these libraries (item 4), and it leaves the test uncorrected wherever the dispersion estimate is not above 1 (item 7)." |
| 363 | Options 3 | keep | "does not have the weakness of item 7" still holds, because the mechanism is the same for the chi-squared test. |
| 369 | Decision of the author 2026-11-20 | A | Proposed: delete "; the F test is not offered as an option" (the preceding "with the chi-squared test only" carries it). The rest holds: 14 and 185, all in the draft's lists, per dml_counts.csv `n_in_primary`. |
| 395 | I-07 update | RF + N/W | "With the overdispersion correction 10 CpGs remain, 3 of them among the draft's 1,038, and the share of corrected p-values below 0.05 (10%) is reached by 39 of the 150 regroupings of the controls (I-06, item 2)." "as many as regrouped control libraries give" (15b_) and "regroupings that keep the two high control libraries in one group" (19_) are N or W. |
| 418 | I-08 update | RF (part 8) | Rewrite: "At a minority of them the test applied no correction because the dispersion estimate is not above 1 (172 of the 1,093 CpGs of the heat-cohort list; NC_047565.1:35,030,701)". Withdraw "so these are CpGs at which the libraries of each group agree closely, not CpGs that passed a test of the variation among animals". "the only CpG in both corrected lists" becomes "one of the four CpGs in both corrected lists of the draft's contrasts, and the one with the smallest p-value". "Of the 1,093 CpGs called for ploidy inside the heat cohort, 36 have all five triploids outside the range of the other 17 libraries, and 16 have all five diploids outside it ... at 24 of the 52 the group is separated from the nearest other library by less than 3.5 points, and at 10 (eight with the triploids apart, two with the diploids apart) by 11 to 19 points; at 14 of the 16 diploid CpGs the diploids are the higher group". The 98, 34 and 14 counts over all 480,579 CpGs do not depend on the test; keep them. |
| 1111 | I-49 item 7 | RF (part 8) | "Of the 1,093 CpGs that separate the ploidy groups of the heat cohort with the overdispersion correction (I-06), 36 have all five heat-cohort triploids outside the range of the other 17 libraries and 16 have all five heat-cohort diploids outside it; the other 1,041 show no group apart." "at 12 of the 17" becomes "at 24 of the 52". "the CpG with the smallest p-value in the two corrected lists" holds. |
| 1191, 1197, 1199, 1202, 1211, 1213, 1225 | Corrections table (dated records) | A | 1197 holds as worded. 1199 names the F threshold and rests on "most corrected calls ... no correction": withdraw the row, or reduce it to the option-1 correction. 1191, 1202 and 1211 carry F numbers (17 of 125; 124 of 125): keep them as dated records with a note "superseded: F test removed 2026-11-20", or update them (52 of 1,093; 172 of 1,093). 1213 and 1225 are the benchmark (as 299-300). |

### README.md

| Line | Action | Proposed wording |
|---|---|---|
| 13, point (8) | R + W | "... and methylKit's corrected test (chi-squared, the setting of the paper the Methods follow) leaves 14 ploidy CpGs and 185 heat CpGs, from a test that does not hold its level on these libraries." Withdraw "calls of a kind that relabelled libraries also produce". |

### 06_diagnostics/README.md (correct line numbers; the trace's are off by 35)

| Line | Action | Proposed wording |
|---|---|---|
| 83 | R | "\| A1 \| Overdispersion correction (`overdispersion = "MN"`, `test = "Chisq"`, the setting of Venkataraman et al. 2022) \| I-06 \|" |
| 84, 86 | keep | "without and with the correction" holds |
| 89 | R | "All analyses run in `legacy` mode. A0, A1 and A2 also run in `corrected` mode (I-05)." |
| 111 | W + R | Replace the two F sentences with: "In the analyses with the correction (`_mn_chisq`) the shares cannot be read against these nominal values: with relabelled ploidy labels the chi-squared form gives 7.2%, 2.2% and 0.47% on average (`issues.md`, I-06, items 2 and 4); read them against relabelled libraries." |
| 112 | W + R | Keep the definition. "The column gives the chi-squared statistic at that p-value (31.4 in legacy mode, 31.3 in corrected mode); with the correction it is the deviance divided by the dispersion estimate that must exceed it. The correction applies only where the dispersion estimate of a CpG is above 1; where it is not, methylKit sets the estimate to 1, so the "corrected" p-value is the uncorrected one (`analyses.csv`, column `reference_distribution`). A corrected list is not a list with 1% false discoveries, because the corrected test does not hold its level (item 3)." Withdraw "For the F test it rises steeply ...", the two "A corrected list that is short ..." sentences, "and a statistic of 31.4 is enough" and "Most CpGs of the corrected lists are of this kind, and relabelled libraries get such calls too". |
| 124 | keep | 22 fits, 16 and 6 |
| 126-136 | R/RF | Header "DML with the correction (chi-squared test)". Rows: 14 (14); 185 (157); 242 (231); not run; not run; 10; 1,093; not run; not run. Share column: "33% and 7.3%", "43% and 13%", "50% and 23%", "37% and 10%", "50% and 25%". |
| 140 | keep | "none with it (... none of the 159 in corrected mode)" holds with chi-squared (0 and 0); re-source to A2_interaction_mn_chisq |
| 141 | W + R | "What an independent check added (...): the chi-squared form of the correction does not hold its level under relabelling (7.2%, 2.2% and 0.47% of p-values below 0.05, 0.01 and 0.001)." Withdraw the F calibration clause and the "median of 0, 2 and 10 corrected calls" tail (N if recomputed). |
| 143 | RS | "(391 of the 1,124 corrected heat CpGs of legacy mode are removed by the filter, 275 of them mitochondrial), and the lists under the Benjamini-Hochberg adjustment (175 and 150 heat CpGs)" |
| 161 | RF (part 8) | "... at the CpGs that `dml_sensitivity.R` calls with the overdispersion correction (chi-squared test; 1,418 CpGs in legacy mode, 398 in corrected mode) ..." Delete the last sentence ("The CpGs called with the correction and the chi-squared test ... only where the F test also calls them (2 and 16)"). |

### 07_documentation/reproducibility_plan.md

| Line | Action | Proposed wording |
|---|---|---|
| 100 | RF | "... lists 16,971 CpGs without the overdispersion correction and 242 with it (231 with the coverage filter), none of them among the draft's 200 ..." |
| 141 | RF | "with the correction it calls 242, none of them among the 200" |
| 153 | R + W | "methylKit's corrected test (chi-squared, as in Venkataraman et al. 2022) leaves 14 ploidy CpGs and 185 heat CpGs (157 with the coverage filter), but it does not hold its level on these libraries, and it applies no correction where its dispersion estimate is not above 1." Withdraw "and relabelled libraries get such calls too". |
| 164, 166, 169 | A | Decision records that quote the instruction. If reworded: on 166, "**Keep one model with the overdispersion correction in the sensitivity analysis, the chi-squared one, and remove the other test of the DML model and every mention of it** ..."; on 169, delete "The F test is not coded as an option of the pipeline."; on 164 (3), "... that name the removed test ...". |

### 07_documentation/session_handoff.md

| Line | Action | Note |
|---|---|---|
| 40 | update when done | The Benjamini-Hochberg item stays on the code-pass list. Numbers change (175 and 150). |
| 53 | R (after the header rewrite) | "confirmed on its four points" becomes "three points" if the F point of the header note is withdrawn |
| 71 | A | Record of this scoping workflow ("(2) the F test of the DML model: ...") |
| 95 | keep until done, then reword | (iii) "done on <date>: the DML model with the correction is fitted with the chi-squared test only" |
| 97 | W | Replace the parenthesis "(the F test results of that diagnostic, 2 and 16, stay where they are as part of the earlier sensitivity analysis)" with "(the diagnostic fits the correction with the chi-squared test only; author, 2026-11-20)". It now contradicts plan 166 and handoff 95. |

## (e) Rows of 03_input_files/manuscript/manuscript_values_map.csv

Change the note column only. Expressions are unaffected: they read `table1_group_means.csv` or `united_cpg`. The p-values were checked against the cached fits by the trace: F 0.043, 0.23, 0.18, 0.019 and 0.15 equal the cached F fits; the binomial p-values 1.3e-16, 0.0038, 0.00008, 0.00001 and 0.018 equal the uncorrected fits.

| Row | Line | Now | New (chi-squared, derived; RF) |
|---|---|---|---|
| D4.3-016 | 476 | "p = 0.043 once the variation among libraries is allowed for (F test; I-06)" | "p = 0.016 once the variation among libraries is allowed for (overdispersion correction, chi-squared test; I-06)". Source: A4_ploidy_in_heat at NC_047564.1:8,425,110, p_chisq 0.01623 |
| D4.9-008 | 527 | "p = 0.23 ... (F test; I-06)" | "p = 0.22 ... (overdispersion correction, chi-squared test; I-06)"; A2_interaction at CAD, 0.2174 |
| D4.9-010 | 529 | "p = 0.18 ... (F test; I-06)" | "p = 0.16 ..."; A2_interaction at nesprin-1, 0.1596 |
| D4.9-014 | 533 | "p = 0.019 and 0.15 ... (F test; I-06)" | "p = 0.010 and 0.13 ..."; A2_interaction at AK9 5,569,229 and 5,569,232: 0.01003 and 0.1306 |

"In a scratch test on the read counts of step 05" should then cite a repository source. Either the script writes the p-values at the named CpGs (a small `named_loci.csv` in `dml_sensitivity/`; a new output, so A under the scope rule), or the note says "from the fits of `06_diagnostics/dml_sensitivity.R` (cached, not written to a table)".

Keep as worded: rows R3.4-029 (283), R3.5-039 (337), R3.5-049 (347), D4.3-002 (462), D4.6-002 (496), D4.9-002 (521) and D3.5x-001 (537). "none when variation among animals is modelled" holds with the chi-squared test: 0 of the 200 and 0 of the 159 (derived). Their source becomes `interaction_vs_shared_<mode>.csv` row `A2_interaction_mn_chisq`. Then rerun step 13 so that `05_manuscript/manuscript_values.csv` follows.

## (f) Must NOT be touched

F statistics of the draft's analyses of variance and other F statistics. These are results of the paper (191 lines in 17 files, sweep class "analysis of variance"):
- 01_pipeline/03_global_methylation.R: 113, 115, 118-119, 121-122, 130, 136
- 01_pipeline/12_qpcr.R: 10, 33, 40, 483, 521, 526, 543, 558-559, 563, 567-568, 572-573, 586, 588, 590-591, 605, 614, 619, 686, 688, 1013, 1026, 1137, 1176
- 01_pipeline/README.md: 22. 03_input_files/qpcr/README.md: 11
- 03_input_files/manuscript/manuscript_values_map.csv and 05_manuscript/manuscript_values.csv: 106-117, 120-125, 492, 494, 510, 536 (R3.1-016 to R3.1-035 and others)
- 03_input_files/qpcr/legacy/dct_analysis_NA-to-45_COX1convert.md: 45 lines from 67 to 652
- 04_output/03_global_methylation/tables/anova_by_scale.csv (columns F_sequential, F_type2); 04_output/12_qpcr/tables/qpcr_anova_legacy.csv (F_value); qpcr_cox1.csv 1, 7-13; qpcr_manuscript_statements.csv 7-9, 11
- 06_diagnostics/README.md: 48-49 (ELISA and %mCpG two-way ANOVA, "interaction F"). The trace's wrong numbers for S40-S41 point here.
- 06_diagnostics/check_legacy_outputs.R: 45, 52-53. legacy_check_results.csv: 4-7
- 07_documentation/audit/elisa_qpcr_audit.md: 34 lines (9 to 535)
- 07_documentation/reproducibility_plan.md: 131, 133
- issues.md: 106, 682, 1051, 1078, 1187, 1204 (for example "F(1,13) = 6.79" in 1187)

Look-alikes:
- GO term name "ATP synthase complex, catalytic core F(1)" in `04_output/09_go_enrichment_legacy/tables/legacy_fisher_{heat,ploidy,shared}_CC.csv` and the six `04_output/10_go_enrichment_gene_level/tables/go_gene_level_*.csv` (30 lines in 9 files).
- GO ontology code "F" in 01_pipeline/08_gene_annotation.R:192, 09_go_enrichment_legacy.R:122, 148, 150 and 11_semantic_similarity.R:187. Gene names such as "member F5" in the dml_gene_overlap tables.
- The chi-squared test of feature class, in 01_pipeline/07_genome_features.R (29, 66, 78, 124, 128, 243-255), 01_pipeline/README.md:17, the feature_class_chisq tables, 04_output/README.md:35 and the values map row 223.
- "corrected test" of the censored regression in 01_pipeline/12_qpcr.R:118-124.
- The uncorrected sensitivity analyses A0, A2_interaction, A3, A4_ploidy_in_control, A4_ploidy_in_heat and A5, their rows, lists and cached fits. The values map reads A4_ploidy_in_control and A2_interaction (rows of lines 157, 271, 283, 324, 346, 433, 447, 520, 525 and others); keep these ids unchanged.
- The two cached chi-squared legacy fits `fit_A1_{ploidy,heat}_mn_chisq_legacy.rds`.
- `dispersion_summary.csv`, the dispersion estimate text (Pearson chi-squared over df).
- "F-score" in the benchmark text (issues.md 299): a search for "F-" must not catch it.
- Values-map ids such as F6-001 and F9-001 (figure rows).
- Spreadsheet ranges such as "A1:K63" in elisa_qpcr_audit.md.
- The words "A1/A2/A4" in GO and Fisher tables (177 lines in 31 files).

Also untouched: 08_archive/ and all of 07_documentation/agent_workflows/ (author's answer).

## (g) Historical records under 07_documentation/agent_workflows/

The author answered at 00:45 UTC on 2026-11-20 to leave them as they are (plan line 164; handoff 95). The three options are listed here for completeness, without choosing.

Facts:
- The folder holds 650 tracked files and 21 MB: 623 under verification_scripts/, 15 under scripts/, 11 under results/, and README.md.
- The sweep found 301 lines in 55 files that name the F test or an `_mn_f` id (`../F1-terms/agent_records_F_mentions.csv`). By folder: S1-I06 32 files, S3 7, results/ 6, S4 5, scripts/ 3, L1 2.
- The sweep's count is a lower bound. Column names such as `mn_f1` in `08_relabel_level_heat_experiment.csv`, `08_relabel_level_control_experiment.csv` and `19_control_regroupings_by_high_libraries.R` are not in its list.

| Option | What it removes | What is lost |
|---|---|---|
| Leave (the author's answer) | nothing | Nothing. The check of section h must exclude the folder. The F mentions stay in 55+ files. |
| Delete the files concerned (about 57 files) | the 55 listed files and the 3 above | The source of the chi-squared numbers that the rewritten text rests on: the relabelling tables `07_relabel_level_pooled_{a,b}.csv` and `08_relabel_level_*.csv`, which hold both tests in one file, and `10_combine_relabel.log` (7.2%, 2.2%, 0.47%). Also 11 files cited by name in issues.md and the handoff: 10_combine_relabel.log, 04_item6.log, 08_relabel_level_within_cohort.log, 14_correction_not_applied.log, 15_corrected_lists_heat.log, 15b_corrected_lists_pooled_control.log, 20_phi_at_corrected_calls.log, S4 03_group_apart.log (the test-independent 98 and 34 CpGs of I-08 and I-49), results/audit-verification.json, scripts/scope-author-changes.js and results/scope-author-changes.json (the record of this scoping workflow). The scripts that produced test-independent numbers (06_sample_and_parity.R, 12b_binomial_null.R, 13_what_the_text_omits.R) also go. Editing the files instead of deleting them would break the record. |
| Delete the folder | all 650 files | Every recomputation behind the statements marked as checked once by an agent. The live documents cite about 55 record files by name (basename match, an upper bound): S1 16, S2 12, S4 7, scripts/code-pass-inputs 3, and others. issues.md carries 26 references to the folder and the handoff 21. The code-pass inputs (`scripts/code-pass-inputs/extra_items.txt`, groups A and F) and the hardening workflow's own results (`results/pipeline-hardening.json`) also go. Every citation would have to be withdrawn or the numbers regenerated by repository scripts first. |

## (h) Checks after the change

1. No mention left. Run `code/check_no_f_mention.sh /home/user/WGBS-gigas-ploidy-desiccation [keep_list]`. It is read-only and does three things:
   - It runs `git grep -n -I -i -E` with the pattern
     `\bF[- ]?tests?\b|\bF-test|test *= *"F"|"MN", *"F"|"Chisq", *"F"|test == "F"|_mn_f\b|_mn_f_|\bmn_f[0-9]?\b|qf\(|\bF statistic|\bF with 1 and|\bF variable|\bF threshold|\bF distribution|\bcorrected F\b|statistic above (84|474|89|250)`
     over the repository, excluding `07_documentation/agent_workflows` and `08_archive`. It removes lines that are in the ANOVA files of (f) and carry an ANOVA marker (anova, aov(, F_value, F_sequential, F_type2, F(n,m), interaction F, two-way, qPCR, ELISA, Tobit, censored). It removes lines containing a string of the optional keep list, for the passages the author decides to keep (benchmark, decision records).
   - It runs `find 06_diagnostics -iname '*_mn_f_*'`, which covers the cache and test_subset.
   - It runs `pdftotext` on every tracked PDF outside the records and greps for "F test" and `_mn_f`.

   It exits 0 only if the remaining count, the file count and the PDF count are all 0. Baseline today (`baseline_check.txt`): 267 lines with a term, 10 of them ANOVA look-alikes, 257 left; 32 files named `*_mn_f_*`; 5 PDF lines (figure_dml_sensitivity.pdf).
2. No F-only statement left, named or not. `git grep -n -F -f code/f_only_phrases.txt` over issues.md, README.md, 06_diagnostics/README.md, the plan, the handoff, the values map and the two scripts must return nothing outside the keep list. The list holds 55 fixed strings of F results, such as "16 heat CpGs", "not stable at its threshold", "A looser threshold gives no list", "5.6%, 1.2%" and "124 of the 125". It hits 47 lines today (the line 466 of the values map that an earlier version caught, "0.36 with the coverage filter", is excluded). Read it together with `statements_checked.csv`: every row must be closed.
3. Every number equals its source. Run `python3 -I code/check_numbers.py claims.csv --repo <repo>`. For each claim it requires the anchor on exactly one line of the file and the written value on that line, and it recomputes the value from `06_diagnostics/dml_sensitivity/dml_counts.csv` or from the stored relabelling tables. Adjust the anchors to the final wording and add a claim for each new number. Run before the rewrite with `--scratch ../F2-dependencies/derived` (values only): all 35 claims pass today against the derivation.
4. Fits equal the derivation: identical DML sets for the 6 new fits and for the two existing chi-squared ones, against `chisq_dml_lists_derived.csv`.
5. Nothing else moved: `git diff --stat` shows no change outside the files named in (a), (c), (d) and (e). The 11 untouched tables of `library_methylation_level/` and the uncorrected rows of `dml_counts.csv` are byte- or value-identical. `06_diagnostics/check_legacy_outputs.R` and step 13 report what they reported before, apart from the four notes.

## 3. Conclusions of issues.md that change in substance

These change but can be kept, rewritten on the chi-squared fit:
1. The corrected lists are 14 ploidy and 185 heat (157 with the filter), not 2 and 16 (6). The summary row, README (8), plan 153, I-01 and I-06 What follows all change.
2. The corrected heat list is largely stable under the coverage filter: 157 of 185 kept, Benjamini-Hochberg 175 and 150. "Not stable at its threshold" is reversed (item 6).
3. The ploidy pair is not at the edge: q = 1.2e-8, not 0.0099. Both CpGs are called inside the heat cohort; the text says neither is. "Differences of the heat cohort alone" holds and is strengthened (item 6).
4. The corrected interaction is a list of 242 CpGs (231), not 2 (I-02, plan 100 and 141, the diagnostics README). "Does not rescue the overlap set" still holds: none of the 200 or 159.
5. The heat-cohort list is 1,093 (179 of them in the 1,038), not 125 (7). The controls list is 10 (3), not 2 (0) (I-06 What follows, I-07).
6. CpGs with no correction applied are a minority of the chi-squared lists (1 of 14, 2 of 185, 5 of 242, 2 of 10, 172 of 1,093), not most. The ground for "a corrected list is not a list of established differences" must move to item 4: the chi-squared form gives 7.2%, 2.2% and 0.47% under relabelling (I-06 item 7, I-08, Options 1, header line 24).
7. In the I-08 and I-49 group-apart reading: 36 and 16 of the 1,093, not 11 and 6 of 125. 24 of the 52 are separated by less than 3.5 points, not 12 of 17. 14 of the 16 diploid CpGs have the diploids higher, not all six. NC_047565.1:35,030,701 is one of four CpGs in both corrected lists, not the only one.
8. A looser threshold gives a list: 37 at q < 0.05 and 52 at q < 0.1, not 3 and 6.
9. Item 2 must be read against relabelled libraries, not against nominal shares, because the chi-squared form does not give 5% where groups do not differ. The readings that follow still hold: pooled ploidy like relabelling (86 of 300), controls like regroupings (39 of 150), heat cohort beyond every regrouping (0 of 125), cohort contrast 12 of 300.
10. I-01 mitochondrial count: 275 of the 1,124 CpGs with q < 0.01, not 138 of 287.
11. "Which CpGs differ cannot be named" (I-01, I-06 What follows): the conclusion can stay. It now rests on the failure of the corrected test to hold its level, not on the shortness of the lists. The author should confirm the reading for lists of 185 and 1,093.

These cannot be kept at all (no chi-squared counterpart, or the counterpart reverses them). Withdraw them, or recompute where marked N:
- The F-statistic thresholds (above 84 with 22 libraries, above 474 with 10) and the reading that single CpGs can hardly be called from this number of animals (item 7, diagnostics README item 4, header note 3). With the chi-squared test the threshold is 31.4 everywhere.
- The F against chi-squared tail comparison (2.1e-5 against 2.1e-8) and "not a middle course" (item 4).
- "The corrected test holds its level on average under relabelling (5.6%, 1.2%, 0.14%; 3.5% to 9.3% in nine of ten)" (item 2, diagnostics README 111 and 141). The chi-squared form does the opposite.
- "Most corrected calls are at CpGs where no correction was applied" (item 7, corrections row 1199, I-08).
- The relabelling tail counts of item 7: 0.14, 2.2 and 5.4 per relabelling below 2.1e-8; medians 0, 2 and 10; "about ten calls from an arbitrary grouping"; "two CpGs among the controls are what regrouped controls give" (I-07); and "1 of 300 relabellings reaches the 20% of the interaction" (I-02, item 2). N.
- The splits by the high libraries (7.6% against 4.5%; 6.6% over 70 against 4.4% over 80) (item 2, I-07). N.
- "Low power does not explain the result for ploidy" (75% against 72% under relabelling) (item 7). The observed share is now 83%, and the relabelling share is not computed. N.
- I-01: "82% of the heat DML have corrected p < 0.05, reached by 7 of 150 relabelled lists". The observed share is now 88%; the comparison is N.
- The four named control and interaction CpGs of item 6.
- The literature passages on the F test of the benchmark, and the decision records that name the F test, are withdrawn only if the author decides so.

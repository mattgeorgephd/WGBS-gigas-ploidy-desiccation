# Report on the overnight work of 2026-11-21 (for the author)

This covers the author's instructions of about 07:00 UTC on 2026-11-21. The work ran from 06:58 to 14:45 UTC, and every step was committed and pushed to `reproducibility-audit` as it finished. The table with commit ids is `07_documentation/overnight_plan_2026-11-21.md`. Nothing here is adopted: every analysis except `legacy` is optional, is off by default, and writes only under `not_adopted/`.

## 1. What was asked and what was done

| Asked | Done | Where |
|---|---|---|
| Run `--coverage-filter --coverage 5 --exclude-ct-variants` | Done (5x fit `corrected-cov5`, 88 min) | `04_output/not_adopted/corrected-cov5-ctexcl/`; issues.md, I-51 (update) |
| Fix the `run_config.R` comment calling 20 points a "minimum" | Done ("more than this", strict as in `getMethylDiff`) | 8a7b64f |
| Run every variation at 50 points as well as 20 | Done for the ten analyses run at 20 points (switch `--difference 50`, token `d50`, reviewed twice before use) | issues.md, I-52; `04_output/not_adopted/README.md`, "The option `d50`" |
| Run the analysis closest to Venkataraman et al. 2022 | Done at 20 and 50 points (5x fit `corrected-cov5-mnchisq`, 89 min) | `corrected-cov5-mnchisq-ctexcl/`, `corrected-cov5-mnchisq-ctexcl-d50/` |
| A table comparing all the variations, committed | Done: 20 analyses | `04_output/not_adopted/comparisons/analysis_comparison.md` (and `.csv`), copied into `04_output/not_adopted/README.md`; written by `01_pipeline/optional/compare_analyses.R` |
| Regular check-ins; save work so a restart loses nothing | Hourly scheduled check-ins; one commit per analysis; 5x fits saved part by part. The container did not restart. | |
| Next steps: final default run, `renv.lock`, clean-clone run, README pass | Done (see section 3) | |
| Checks for errors once everything had finished | Done (section 4) | |

Also run, on the session's judgment: combination C2 (`mnchisq-ctexcl`, with the gene-level GO) and its 50-point version. The handoff listed C2 among the runs still to do, and no document had dropped it.

## 2. Results

The full table is `04_output/not_adopted/comparisons/analysis_comparison.md`. The main points:

- **At 50 points the lists nearly vanish.** With the paper's settings otherwise, 1 ploidy and 4 heat DML remain, all in the paper's lists, and no locus is in both lists. Across the ten 50-point analyses the counts run from 1 to 6 ploidy and 2 to 92 heat DML, with no locus in both lists in any of them (I-52). At 10x, 1,037 of the paper's 1,038 ploidy DML and 2,347 of its 2,351 heat DML have differences between 20 and 50 points.
- **The analysis closest to Venkataraman et al. 2022** has 96 ploidy and 875 heat DML at 20 points (19 in both lists), and 1 and 47 at 50 points. One difference from 2022 remains: the variant rule counts reads with their full CIGAR (your decision of 2026-11-20; I-08).
- **The overdispersion correction is the setting that changes most at 20 points.** At 5x it shortens the lists about 72-fold (ploidy) and 18-fold (heat). At 50 points it changes nothing at 10x: the 1 and 4 DML pass the corrected test.
- **The upper coverage filter changes the 5x lists little** (ploidy 7,046 to 6,902, heat 15,530 to 15,376).
- **The C-to-T exclusion removes 23.8% of the ploidy and 22.8% of the heat DML** of the analysis closest to 2022, against the 19.7% that 2022 reports. The two studies are not directly comparable.

Whether any of these lists are real differences is not shown by the runs. The corrected test does not hold its level with relabelled libraries at 10x (I-06), and that has not been checked at 5x.

## 3. Repository

- **Final default run** (`Rscript run_pipeline.R`): `tests/check_paper_outputs.R` passes, with 189 paper files byte-identical (d99825c).
- **Run from a clean clone:** a fresh clone of d99825c, with no fit cache, ran the whole pipeline. Step 05 refitted the paper's model in 19.3 minutes, and T1 passes: 189 paper files byte-identical and the input inventory identical. It used the download cache of the gannet inputs on this machine, so the download from gannet itself was not repeated.
- **`renv.lock`:** 201 packages, R 4.3.3, Bioconductor 3.18. Restore instructions are in README, "Running the pipeline". A restore from the lockfile on a clean machine has not been tested.
- **README pass:** the status paragraph is shortened, the counts are current (52 issues), and the table rows of the settings are identical in every README. The handoff and the plan are updated, with decisions 18 (I-51) and 19 (I-52) added to section 7.

## 4. Checks for errors

| Check | Result |
|---|---|
| `tests/check_settings.R` | 62 of 62 |
| `tests/check_overdispersion.R` | pass (14 and 185 DML identical to the diagnostic lists) |
| `tests/check_paper_recipe_legacy.R` | 6 of 6 |
| `tests/check_paper_outputs.R` | pass, after the final default run and in the clean clone |
| `tests/check_cached_fit.R`, `check_step13_not_run.R`, `check_ct_exclusion.R`, `check_fit_by_sequence.R` (each in a fresh scratch copy) | all pass: cached fit 11 of 11; step 13 not run 17 of 17; C-to-T exclusion 23 of 23; fit by sequence `identical()` to the fit in one piece for `legacy` and `mnchisq`, with resume, stale-part and fingerprint checks (9 of 9 each); 13:46 to 14:42 UTC |
| Paper snapshots around every optional run | no paper output changed in any run |
| Independent review of the new issues.md text (I-51 update, I-52) | 14 findings, all applied (bd0be97). Among them: two wrong loci named in I-52, item 1; "every setting of Venkataraman et al. 2022", which ignored the full-CIGAR difference; step 13 counts; run times. Venkataraman et al. 2022 settings confirmed against the article text. |
| Final check of the outputs of every analysis and of the documentation (independent agent) | No wrong number in any output table (counts, lists, bed files, subsets, step 13 states, about 120 values quoted in the documentation). 11 documentation findings, applied (b0b6c35): stale statuses, out-of-date handoff, run times and file sizes quoted as estimates, inconsistent copies of the settings table, one em dash. One code finding: `analyses_included.csv` lagged after a run ending with step 09, fixed in `compare_analyses.R`, which now rebuilds it. |

Errors of the session found during the night, all corrected: I-52 named the wrong two heat loci for `ctexcl-d50`, and the overnight plan table gave three start times 2 to 17 minutes late (corrected from the logs).

## 5. For you to decide or know

1. **I-51 and I-52 (plan, section 7, items 18 and 19).** Which coverage standard and which difference threshold the paper uses, and whether to report the other as a sensitivity analysis.
2. **Eleven step 13 statements apply the threshold to group means** (R3.5-029, R3.5-031, R3.5-033, R3.5-047, D4.2-004, D4.3-005, D4.3-016, D4.9-008, D4.9-010, D4.9-012, D4.9-014). Under `d50` they are evaluated at 50 points. Should they keep 20, the threshold of the draft's text?
3. **Done on 2026-11-21 (evening), at your request: the remaining six combinations at both thresholds** (section 6). All 16 combinations of the four settings now exist at 20 and 50 points.
4. **Step 04 runs again for every `d50` id** (about 9 minutes and 9.2 GB each), although its outputs equal those of the analysis without `d50`. Reusing them would need a change of step 04.
5. **The two step 13 errors of `cov5`** (R3.5-046, D4.9-012) remain. Selecting the start 42,498 in their two templates would let them be evaluated (I-51, item 6). Not changed, because the values map needs your approval.
6. **Step 09 does not rewrite `comparisons/analyses_included.csv`**, so the file lags after a run that ends with step 09. `compare_analyses.R` now rebuilds it. The pipeline step itself was not changed.
7. **One tracked legacy file is over 50 MB** (`08_archive/old_WGBS/WGBS/DML/ALL-CpG-Cov10.csv`, 85 MB; I-34, plan decision 10). The license is also still open (I-34).
8. **The `GO_bias_by_CpG_count` figure** is listed as "not produced" in step 13's log for analyses run without `--gene-level-go`. This is expected: step 10 makes it.

## 6. Follow-up of 2026-11-21 evening: the remaining six combinations

Asked: "Run the remaining six combinations at both thresholds". Twelve runs, 19:20 UTC to 00:00 UTC on 2026-11-22, one at a time, each with a paper snapshot before and after (no paper file changed) and a commit of its own. The fits of four were in the cache; `cov5-mnchisq` needed a new 5x fit by sequence. A container restart at about 21:32 UTC stopped that fit after 6 of its 13 ploidy parts; it resumed from the saved import and parts and finished at 22:55 UTC. The table of all 32 analyses is `04_output/not_adopted/comparisons/analysis_comparison.md`.

| Analysis | CpGs tested | Ploidy DML | Heat DML | In both lists | At 50 points: ploidy, heat |
|---|---|---|---|---|---|
| `corrected-ctexcl` | 427,349 | 635 | 1,596 | 100 | 1, 2 |
| `corrected-mnchisq-ctexcl` | 427,349 | 12 | 118 | 4 | 1, 2 |
| `corrected-cov5` | 2,265,759 | 9,791 | 21,200 | 1,711 | 6, 92 |
| `corrected-cov5-mnchisq` | 2,265,759 | 126 | 1,133 | 27 | 1, 62 |
| `cov5-mnchisq` | 2,297,236 | 126 | 1,191 | 27 | 1, 63 |
| `cov5-mnchisq-ctexcl` | 2,153,017 | 96 | 923 | 19 | 1, 48 |

What the full set adds:
- Across all 16 combinations at 50 points, the ploidy list has 1 to 6 DML, always including NC_047561.1:6,740,817, and the heat list 2 to 92. No combination has a locus in both lists (I-52).
- The upper filter changes nothing at 50 points: with and without it, the 5x lists are the same 6 and 92 DML (4 and 63 with the exclusion). At 20 points it removes few: 9,942 to 9,791 ploidy DML at 5x, 1,191 to 1,133 heat DML with the correction.
- With the correction, the 5x ploidy list is 126 with or without the filter, and 96 after the exclusion.
- Step 13 also reports its two known errors for `corrected-cov5` (R3.5-046, D4.9-012), with the same cause as for `cov5` (I-51, item 6).

These results are recorded in issues.md (I-08, I-51 and I-52), in `04_output/not_adopted/README.md` and in the handoff.

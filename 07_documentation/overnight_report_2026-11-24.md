# Report on the runs of 2026-11-22 to 2026-11-24 (for the author)

This covers your instructions of 2026-11-22: the q < 0.05 variations (about 04:30 UTC), and, before you went to bed (about 04:45 UTC), the 25-point columns, the control-only ploidy contrast, the tested interaction and the full factorial ("Everything, ~2 days"). The work ran from 04:27 UTC on 2026-11-22 to about 15:30 UTC on 2026-11-24. Every analysis was committed and pushed to `reproducibility-audit` as it finished. The table with commit ids is `07_documentation/overnight_plan_2026-11-21.md`, section "Night of 2026-11-22", rows 16 to 24.

Nothing here is adopted. Every analysis except `legacy` is optional, is off by default, and writes only under `not_adopted/`. No paper file changed in any run: snapshots were taken before and after each run, and `git diff` from the merge of the options to the last run shows no paper output changed.

## 1. What was asked and what was done

| Asked | Done | Where |
|---|---|---|
| Run the variations at q < 0.05 (all 32: the 16 combinations at 20 and 50 points) | Switch `--qvalue 0.05` (token `q05`), reviewed independently, 94 of 94 settings checks; all 32 runs | issues.md, I-52, "The q-value cutoff"; `04_output/not_adopted/README.md`, "The option `q05`" |
| Counts at 10, 15 and 25 points (25 is the methylKit default) as columns of the comparison table | Done; the script stops if a threshold sweep disagrees with the DML counts | `01_pipeline/optional/compare_analyses.R`; second table of `analysis_comparison.md` |
| Ploidy contrast on the 12 control libraries (I-07), in the full factorial | Switch `--controls-only-ploidy` (token `ctrl`); 64 analyses | issues.md, I-07, updates of 2026-11-23 and 2026-11-24 |
| Tested ploidy by heat interaction (I-02), in the full factorial | Switch `--interaction` (token `inter`); 64 analyses | issues.md, I-02, updates of 2026-11-23 and 2026-11-24 |
| A table comparing all the variations | 192 analyses | `04_output/not_adopted/comparisons/analysis_comparison.md` (and `.csv`) |
| Regular check-ins; save progress in case of a restart | Hourly scheduled check-ins; one commit per analysis; 5x fits saved part by part; done flags let the batch driver resume | |

The 192 analyses are the 64 combinations of six settings that do not change the contrasts (`corrected`, `cov5`, `mnchisq`, `ctexcl`, `d50`, `q05`, each on or off; the paper is all six off), each with the paper's contrasts, with `ctrl` and with `inter`. `ctrl` and `inter` needed 16 new fits. For each, only the new contrast was fitted and the base fit was taken from the cache: 6.5 to 13.6 minutes for step 05 at 10x and 27 to 69 minutes at 5x. The peak memory of the machine was 10.4 GB. The other 112 runs reuse those fits.

Built and tested before use:
- On simulated data, the control-only fit gives the p- and q-values of the diagnostic A4 exactly (`tests/check_contrast_variants.R`).
- On the real data, `ctrl` reproduces A4's 2,768 ploidy DML and `inter` reproduces A2's 16,971 interaction DML. `tests/check_contrast_variants_real.R` passes for `ctrl`, `mnchisq-ctrl`, `inter`, `mnchisq-inter`, `corrected-inter` and `corrected-mnchisq-inter` (rerun today).
- An independent review found no correctness defect. Its two findings were applied.

## 2. Results

**q < 0.05 (I-52).**
- *Without the overdispersion correction the cutoff hardly matters*, because the 20-point threshold binds. At 10x it adds one heat DML; at 5x the lists grow by 1% to 5%.
- *With the correction the lists are 2.1 to 2.6 times as long.* For example, `mnchisq` goes from 14 to 37 ploidy and from 185 to 410 heat DML.
- *At 50 points* nothing material changes, and no analysis has a locus in both lists.

**Control-only ploidy contrast (I-07).**
- *Without the correction, at 20 points,* the control-only list is longer than the all-sample list: 1,985 to 3,326 DML at 10x and 14,741 to 27,725 at 5x, 1.9 to 3.8 times as long.
  - 77% to 84% of it is higher in triploids. That is the direction of the difference in overall methylation level between the control groups (I-49).
  - At the paper's settings a list of this length is what regroupings of the control libraries give (update of 2026-11-19); this was not checked for the other settings.
- *With the correction* it has 5 to 53 DML at 10x and 106 to 398 at 5x.
- *It almost never shares a locus with the heat list.* Under the correction there is no shared locus in 14 of the 16 analyses at 20 points. Each of the other two shares one locus, at q < 0.05, and that locus is not among the all-sample shared loci.
- *The all-sample shared loci are not ploidy differences among the controls.* The loci that the all-sample corrected lists share (4 at 10x, 27 at 5x) have control-only ploidy q-values of 0.25 to 0.92 (update of 2026-11-23).
- **Reading:** the shared set ("interaction-DML") depends on the heat-treated libraries, which are also a separate experiment and sequencing project (I-01).

**Tested interaction (I-02).**
- *Without the correction* the interaction list holds 3% to 5% of all tested CpGs, and 86% to 91% of it has one sign: heat-cohort diploids and control-cohort triploids are the more methylated groups. A quarter to more than half of the shared loci are among them (26% to 56%).
- *With the correction* it has 194 to 2,006 DML at 10x and 1,427 to 9,044 at 5x. None of the shared loci is among them in any of the 16 analyses.
- *At 50 points* it keeps 61 to 7,876 DML. The interaction difference is a difference of two differences, so a threshold chosen for the main effects does not carry over to it (I-52, update of 2026-11-24).

**What follows.** The conclusions of I-02 and I-07 hold at every combination of settings. Under the correction, the shared set comes from the all-sample ploidy contrast, it does not reappear among the controls, and the tested interaction does not contain it. Nothing in the runs shows a ploidy-specific response to heat that the confounding of I-01 and the five high-level libraries of I-49 would not also produce.

**Checked afterwards at 5x, at your request (2026-11-24, `06_diagnostics/dml_corrected_test_5x.R`).** The relabelling check of the 10x diagnostic, made for the corrected 5x lists of `cov5-mnchisq-ctrl` and `cov5-mnchisq-inter` (the paper's settings with 5x coverage):
- *Control-only ploidy:* the real list has 120 DML; 150 regroupings of the control libraries give a median of 135.5, and 92 of them reach 120. No evidence of a ploidy difference at single CpGs among the controls, as at 10x.
- *Interaction:* the real list has 1,718 DML; 200 relabellings give a median of 136.5 and at most 1,451, so none reaches it. The relabellings that come closest keep most of the real grouping of the heat-treated libraries. The ploidy difference is not the same in the heat cohort as among the controls, but the heat cohort is also a separate experiment (I-01), so this is a ploidy by cohort interaction, not evidence of a ploidy-specific response to heat.
- *The corrected test exceeds its level at 5x too:* over all relabellings, 8.6% (controls) and 6.7% (interaction) of p-values fall below 0.05, against 5%.
**Checked on 2026-11-25, at your request:**
- *Every corrected control-only and interaction analysis* (64; `06_diagnostics/dml_corrected_test_settings.R`). All 32 control-only lists are of the length that regroupings of the control libraries give (at least 48 of 150 regroupings reach each). All 32 interaction lists are longer than relabelled lists (no relabelling reaches 30 of them; one of 200 reaches the other two), and in every one the closest relabellings keep the real grouping of the heat cohort.
- *The five high-level libraries and the 5x interaction list* (`06_diagnostics/interaction_high_libraries_5x.R`). Without them 118 of the 1,718 DML stay called; without any other set of five of the same make-up, a median of 701 (134 to 1,189). Genome-wide without them the corrected test calls 526, against a median of 2,314 without the other sets. What remains is the pairing of D11 and D12 against the heat-cohort triploids.
- *Reading.* The control-only contrast gives no evidence of a ploidy difference at single CpGs at any setting. The interaction is a difference between the heat cohort's ploidy groups and the controls', carried mostly by the five libraries and otherwise by D11 and D12, in a cohort that is also a separate experiment (I-01). It is not evidence of a ploidy-specific response to heat. (issues.md, I-02, I-06, I-07 and I-49, updates of 2026-11-25.)
- *Not checked by these scripts:* the uncorrected lists at settings other than the paper's, and the corrected ploidy and heat lists with all 22 libraries outside the paper's settings.

## 3. Incidents

1. **Out-of-memory kill caused by the session (2026-11-23, about 04:46 UTC).**
   - A runner was started outside the session's task tracking and then stopped, but its step-04 R process survived.
   - The next run (`inter-d50-q05`) was killed for lack of memory, with two step-04 copies at once.
   - The orphan was killed and the run repeated. The runner now refuses to start while any pipeline R process is running.
   - Nothing wrong was committed.
2. **Container restart (2026-11-24, about 00:15 UTC).** The run in progress was lost two minutes in and was rerun. The caches, done flags and the scratch directory survived.
3. **Disk.** With 9.9 GB free before the last 15 runs, 434 git-ignored per-CpG tables (`diffmeth_*`, `united_cpg_*`) of finished `d50` and `q05` analyses were deleted. None was tracked, read by a test or cited in issues.md, and the tables of the base `ctrl` and `inter` analyses that the real-data test reads were kept. Free space went from 9.9 to 27 GB. Any of them is regenerated in half a minute with the fit in the cache: `WGBS_ANALYSIS=<id> Rscript 01_pipeline/05_differential_methylation.R`.

## 4. Checks

| Check | Result |
|---|---|
| `tests/check_settings.R` | 134 of 134 (rerun today) |
| `tests/check_contrast_variants_real.R` for the six base analyses it covers | all pass (rerun today) |
| Paper snapshots around every one of the 160 runs (32 `q05`, 128 `ctrl` and `inter`) | no paper output changed |
| Files changed outside `not_adopted/` and the documentation, from the merge of the options to the last run | one: `compare_analyses.R` (the planned 25-point columns) |
| Comparison table against the analysis folders | 192 rows, one per analysis folder; none missing |
| Independent check of the new issues.md, README and plan text against `analysis_comparison.csv` and the per-analysis lists (all 4 x 192 list lengths equal the table) | 6 findings, all applied: one wrong ratio (3.9 for 3.8), a range quoted with the wrong end points, a denominator (16 for 32 analyses in the 50-point column), three statements wider than their evidence (the relabelling checks were made at the paper's settings only, not at 10x in general); and four statements of earlier dates corrected where they stand (issues.md, corrections table, rows of 2026-11-24) |

## 5. For you to decide or know

1. **Adoption.** None of the 192 analyses is adopted. The decisions of I-02 (rename the shared set; whether to report the tested interaction), I-07 (the words "baseline", "under control conditions", "controls only") and I-06 (the model) are still yours. The runs give the counts for any choice.
2. **A threshold for the interaction**, if it is reported. The 20 or 50 points of the main effects do not transfer to a difference of differences.
3. **The fit caches exist only on this machine** (your decision of 2026-11-22: no external storage). If the container is reclaimed, each analysis is rebuilt by rerunning it, and its fit is made again first. Step 05 of the 16 `ctrl` and `inter` fits alone took about 7.5 hours (386 minutes for the eight at 5x), and each 5x base fit 88 to 107 minutes. The paper's own fit takes 19 minutes.
4. **The draft pull request was not created.** GitHub refused it with 403 "Resource not accessible by integration". Either reconnect GitHub for Claude and approve the Claude app on the repository, then ask again, or open it yourself from https://github.com/mattgeorgephd/WGBS-gigas-ploidy-desiccation/compare/master...reproducibility-audit. The description written for it is in the session's scratch directory (`pr_body.md`), which does not survive a new machine.
5. **Still open from earlier:**
   - the empty Google Drive folder: whether to delete it;
   - a license (I-34);
   - the tracked legacy file over 50 MB (I-34);
   - items 2 and 4 to 6 of the report of 2026-11-21.

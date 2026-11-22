# Overnight work plan, 2026-11-21 (author asleep; approved by the author)

The author's instructions (2026-11-21, about 07:00 UTC):
- Run the combination of all three options, `--coverage-filter --coverage 5 --exclude-ct-variants`.
- Fix the `run_config.R` comment (done, 8a7b64f).
- Run every variation of the analyses with a 50-point difference cutoff (Venkataraman et al. 2022) as well as the paper's 20 points.
- Also run the analysis closest to Venkataraman et al. 2022, `--coverage-filter --coverage 5 --overdispersion --exclude-ct-variants`.
- Make a table comparing all the variations, and commit it.
- Check in regularly and save work along the way, so a restart loses nothing.
- When everything has finished, run checks to find errors, so the author can review.
- Then go on with the next steps: the final default run, `renv.lock`, a run from a clean clone, and the README pass.

Long jobs run as tracked background tasks, one 5x fit at a time; the machine has 15 GB and a 5x fit peaks near 9.5 GB. A stopped 5x fit is restarted with the same command; it resumes from its saved parts. Before any restart, run `pgrep -af "05_differential"` and kill any leftover process.

| # | Work | State |
|---|---|---|
| 1 | `corrected-cov5-ctexcl` steps 4:5 (5x fit `corrected-cov5`) | done 08:42 UTC (878d4b6): ploidy 6,902, heat 15,376 DML; 88 min |
| 2 | 50-point switch (`--difference 50`, token `d50`), built by an agent, then an independent review | built f93496c, review fixes 797bfda, re-review fixes 033411b (done). At 10x: 1 ploidy and 4 heat DML. A d50 run reruns step 04 (9.2 GB), so not beside a 5x fit |
| 3 | `corrected-cov5-mnchisq-ctexcl` steps 4:5 (5x fit `corrected-cov5-mnchisq`) | done 10:19 UTC (5a0ea9f): ploidy 96, heat 875 DML (before C-to-T exclusion 126 and 1,133); 89 min, peak 9.04 GB |
| 4 | Steps 6:9,13 of 1 and 3 | done 10:31 UTC (cdf9742, c73b4a5, fc6349b); paper files unchanged both times; step 13: only the 2 known cov5 errors |
| 5 | d50 version of every analysis: legacy, corrected, mnchisq, corrected-mnchisq, ctexcl, cov5, cov5-ctexcl, corrected-cov5-ctexcl, corrected-cov5-mnchisq-ctexcl | done 12:21 UTC: all nine committed (c302b2a, a5ae0f1, 6573719, a7415ae, a1b6b6f, a6517c5, d84d8ff, b6d827b, 81d98e5); 10 to 14 min each; paper files unchanged each time; step 13: 3 to 10 statements not evaluable per d50 analysis (documented state), only the 2 known cov5 errors |
| 6 | Comparison table of all analyses at 20 and 50 points (script, CSV, README table) | done 12:25 UTC (6b98683): 18 analyses; table in 04_output/not_adopted/README.md |
| 7 | Results into issues.md, then an independent review | done (c0dfb9b; review applied bd0be97: 14 findings, among them two wrong loci in I-52, item 1, and 'every setting of Venkataraman et al. 2022' narrowed for the full CIGAR) |
| 8 | Final default run with tests/check_paper_outputs.R | done 12:41 UTC (d99825c): T1 passes, 189 paper files byte-identical |
| 9 | `renv.lock`, run from a clean clone, final README pass (task 17) | renv.lock (8a8566e), README (4a97807, 6dcbbe9); clean-clone run of `Rscript run_pipeline.R` (fresh clone of d99825c, no fit cache): done 13:19 UTC, step 05 refitted in 19.3 min, T1 passes (189 paper files byte-identical) |
| 10 | Final error checks: every test in tests/, the outputs of every analysis, documentation consistency; a report for the author | check_settings 62/62, check_overdispersion, check_paper_recipe_legacy 6/6 pass; outputs and documentation check done (no wrong number; 11 documentation findings applied, b0b6c35); C2 (mnchisq-ctexcl, planned in the handoff and missed) done 13:33 (e34e7ac; 12 ploidy, 142 heat DML) and mnchisq-ctexcl-d50 done 13:44 (42fd560; 1, 2); comparison table now 20 analyses; scratch-copy tests all pass (check_cached_fit 11/11, check_step13_not_run 17/17, check_ct_exclusion 23/23, check_fit_by_sequence identical for legacy and mnchisq), 14:42 UTC. Report for the author: 07_documentation/overnight_report_2026-11-21.md. DONE |

Each run has `tests/snapshot_paper.R` before and after it. Progress is logged in the "State" column and in commits.

## Follow-up of 2026-11-21 (author's instruction: "Run the remaining six combinations at both thresholds")

The six combinations of the four settings not run before, each at 20 and 50 points (12 runs), one at a time, each with a paper snapshot before and after and a commit:

| # | Work | State |
|---|---|---|
| 11 | Batch 1 (10x, cached fits): `corrected-ctexcl`, `corrected-ctexcl-d50`, `corrected-mnchisq-ctexcl`, `corrected-mnchisq-ctexcl-d50` | done 20:12 UTC (193b307, ff4c3b8, 4552c84, a2411bb): corrected-ctexcl 635 ploidy, 1,596 heat DML; corrected-mnchisq-ctexcl 12, 118; both d50 versions 1, 2; paper unchanged each time |
| 12 | Batch 2 (5x, cached fits): `corrected-cov5`, `corrected-cov5-d50`, `corrected-cov5-mnchisq`, `corrected-cov5-mnchisq-d50` | done 21:03 UTC (53266aa, c2101c8, 689bca5, 0fd6514): corrected-cov5 9,791 ploidy, 21,200 heat DML; corrected-cov5-mnchisq 126, 1,133; d50: 6, 92 and 1, 62; paper unchanged each time |
| 13 | New 5x fit `cov5-mnchisq` by sequence: `--coverage 5 --overdispersion --steps 4:5` (about 1.7 h; resumes from saved parts if stopped) | step 04 done (2756933); step 05 stopped by a container restart at about 21:32 UTC after 6 of 13 ploidy parts; resumed from the saved import and parts at 21:34 UTC; done 22:55 UTC (af65ec7): 126 ploidy, 1,191 heat DML; paper unchanged |
| 14 | Batch 4: `cov5-mnchisq`, `cov5-mnchisq-d50`, `cov5-mnchisq-ctexcl`, `cov5-mnchisq-ctexcl-d50` | done 23:59 UTC (f549e30, c267aaf, ac3fcc1, 796cccf); paper unchanged each time |
| 15 | Comparison table (32 analyses), issues.md (I-51, I-52), READMEs, report; independent check | table and documents done 2026-11-22 (761eb37); independent check: no wrong number, 7 wording findings applied |

## Night of 2026-11-22: the q < 0.05 variations (author asleep; approved by the author)

The author's instructions (2026-11-22, about 04:30 UTC): "add more run variations based on the q < 0.05"; scope chosen: all 32 (the 16 combinations of corrected, cov5, mnchisq, ctexcl, each at 20 and 50 points, with q < 0.05); "Conduct regular check-ins and save progress in case of a restart." The literature evidence was added to I-06 and I-52 first (1e4b6fc).

| # | Work | State |
|---|---|---|
| 16 | Build the switch `--qvalue 0.05` (token `q05`, last in the token order; setting `dml_qvalue_standard`; stamp column `dml_qvalue`; fit reused, as for d50), by an agent in the working tree, not committed by it | done (93 checks; paper run T1 passes; q05 1,038 ploidy, 2,352 heat DML) |
| 17 | Independent review of the build; apply findings; `tests/check_settings.R`; commit | done (76a900b): no blocking defect, 4 small findings applied, 94 of 94 checks |
| 18 | The 32 q05 runs, one at a time with the batch driver (plan check, paper snapshot before and after, commit per analysis): 16 at 10x, 16 at 5x; all fits cached; driver scratchpad/run_batch.sh (switches from the id) | batch A (8 at 10x) done 07:12 UTC (c04375d to fc8a0ca); B (8 at 10x) done 09:56 UTC (e86d7c4 to 3d30591); C (8 at 5x) in two halves of four (each 5x run about 14 minutes): C1 done 11:04 UTC (083cc59 to d64446b); C2 done 12:08 UTC (a35e5c0 to dc15d8f); D1 done 13:16 UTC (4572081 to 0192840); D2 done 14:22 UTC (abb63b4 to 2fbc7ed). All 32 done; no paper file changed in any run |
| 19 | Comparison table (64 analyses); issues.md (new issue or I-06/I-52 updates); READMEs; report for the author; independent check | issues.md I-52 paragraph "The q-value cutoff" (b2fd9d6); table with row 20's merge; READMEs, report and the independent check with row 24 |

Restart rule: files survive a container restart, processes do not. If the build agent was lost, its edits are in the working tree (inspect `git diff`), and the session finishes or reverts them. A lost batch is rerun with the same driver command (it skips analyses already done).

Added by the author before going to bed (2026-11-22, about 04:45 UTC), after an evaluation of further variations:
- Counts at 10, 15 and 25 points (25 is the methylKit default) as columns of the comparison table, from `dml_counts_by_threshold_<id>.csv`; no new runs.
- Two new contrasts, each in the full factorial: the ploidy contrast fitted on the 12 control libraries only (I-07; token `ctrl`), and a tested ploidy-by-heat interaction (I-02; token `inter`). Each over the 8 fit settings (filter, coverage, correction): 16 new fits (about 14 h), then on each fit the variants that act after the test (C-to-T exclusion, 50 points, q < 0.05): 128 analyses. The author chose to run everything (about 2 days).

| # | Work | State |
|---|---|---|
| 20 | 25-point (and 10, 15) columns in `compare_analyses.R` | built in a worktree (branch `cmp-thresholds`, cc0eb48; the script stops if a sweep disagrees with the DML counts; all 66 analyses agree); merged 15:38 UTC (2c064da); table of 69 analyses committed (463b40e) |
| 21 | Build the switches `ctrl` and `inter` (not combined with each other) by an agent in a separate git worktree, from the commit of row 17; static tests there; fit tests after the q05 runs; independent review; merge | built on branch `worktree-agent-ac204932127819f20` (8c5d2c9, 6cc644b, b63a90f; 134 settings checks, 25 synthetic checks: ctrl gives A4's p and q exactly on simulated data); independent review: no correctness defect, existing fit fingerprints unchanged; 2 findings applied (6e0c5e4). merged 07:13 UTC (6b22c14; check_settings 134/134); real tests done 08:02 UTC: ctrl (d889352) 2,768 ploidy DML = A4, inter (7c14a71) 16,971 interaction DML = A2, 99 in the draft's shared set; tests/check_contrast_variants_real.R all pass for both; peak memory of the new 10x contrasts under 4 GB at the paper's settings (ctrl must give A4's 2,768 ploidy DML, inter A2's 16,971 interaction DML) with tests/check_contrast_variants_real.R; peak memory of the first 5x run to be measured |
| 22 | The 16 new fits (one at a time; 5x ones by sequence), each with its base analysis; commit per analysis | ctrl and inter at the paper's settings done (row 21); the other 14 with ALLOW_FIT=1, each a single new contrast (plan estimate 12 minutes at 10x, 45 at 5x): F1 (corrected-ctrl, corrected-inter, mnchisq-ctrl) done 15:38 UTC (e5171cd, abba8e6, 9fd7de7; step 05 6.5 to 13.6 minutes, 20 to 30 minutes a run; peak memory 10.2 GB, in step 04); F2 (mnchisq-inter, corrected-mnchisq-ctrl, corrected-mnchisq-inter) done 16:59 UTC (57569d6, 6449536, 39ed2e4; peak 10.1 GB); the eight 5x fits one per task (G1 to G8: cov5-ctrl, cov5-inter, corrected-cov5-ctrl, corrected-cov5-inter, cov5-mnchisq-ctrl, cov5-mnchisq-inter, corrected-cov5-mnchisq-ctrl, corrected-cov5-mnchisq-inter; scratchpad/fits_G<n>.log, memory in scratchpad/mem_G<n>.log): G1 done 17:51 UTC (c9218e5; step 05 31.6 minutes; peak memory of the machine 8.9 GB, of the fit parts 6.4 GB: the 5x memory measurement asked by the review); G2 running since 17:51 UTC. Disk: 22 GB of stale scratch copies deleted 14:30 UTC (30 GB free); if free space falls under 8 GB, the git-ignored per-CpG tables (diffmeth_*, united_cpg_*) of finished analyses are deleted (regenerable) |
| 23 | The 112 remaining variants on the new fits (ctexcl, d50, q05 and their combinations); batch driver; commit per analysis | after 22 |
| 24 | Comparison table of all analyses; issues.md (I-02, I-07, I-06, I-52 updates or new issues); READMEs; report; independent check | after 23 |


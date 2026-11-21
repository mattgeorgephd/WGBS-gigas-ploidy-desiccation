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
| 5 | d50 version of every analysis: legacy, corrected, mnchisq, corrected-mnchisq, ctexcl, cov5, cov5-ctexcl, corrected-cov5-ctexcl, corrected-cov5-mnchisq-ctexcl | batch 1 (the five at 10x) running since 10:33 UTC via scratchpad/run_d50_batch.sh (one at a time; plan check, paper snapshot before and after, commit after each); batch 2 (the four at 5x) after it |
| 6 | Comparison table of all analyses at 20 and 50 points (script, CSV, README table) | script committed (658b730); run after 5 |
| 7 | Results into issues.md, then an independent review | after 6 |
| 8 | Final default run with tests/check_paper_outputs.R | after 7 |
| 9 | `renv.lock`, run from a clean clone, final README pass (task 17) | after 8 |
| 10 | Final error checks: every test in tests/, the outputs of every analysis, documentation consistency; a report for the author | last |

Each run has `tests/snapshot_paper.R` before and after it. Progress is logged in the "State" column and in commits.

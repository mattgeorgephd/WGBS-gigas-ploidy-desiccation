export const meta = {
  name: 'build-dml-sensitivity',
  description: 'Build, review and fix the sensitivity analyses of the differential methylation calls (script only; the full model fits are run afterwards)',
  phases: [
    { title: 'Build', detail: 'write 06_diagnostics/dml_sensitivity.R and test it on a subset of CpGs' },
    { title: 'Review', detail: 'independent statistical and code review on the subset' },
    { title: 'Fix', detail: 'confirmed findings fixed, subset run repeated' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'

const CONTEXT = `
CONTEXT (read all of it before starting)
Repository: ${REPO} (git branch reproducibility-audit). It holds an R-only pipeline that re-derives the results of a manuscript draft on DNA methylation (whole-genome bisulfite sequencing, WGBS) in diploid and triploid Pacific oysters after a heat exposure, together with an audit of that draft (issues.md).
Design of the WGBS data: 22 libraries. 10 are heat-treated animals (5 diploid, 5 triploid; Zymo project zr3534; "ronit" cohort) and 12 are controls (6 diploid, 6 triploid; Zymo project zr3644; "yaamini" cohort). The controls come from a different experiment and sequencing batch than the heat-treated animals, so "heat" is confounded with cohort (issue I-01). Four of the ten heat libraries also received a 45 C acute heat shock (column acute_heat_shock of 03_input_files/metadata/wgbs_libraries.csv; issue I-04).
What the draft did (step 05, 01_pipeline/05_differential_methylation.R; read it and its helpers in 02_R_functions/methylkit_helpers.R): CpGs with at least 10x coverage in all 22 libraries; methylKit::calculateDiffMeth with default settings (logistic regression, no overdispersion correction, SLIM q-values), once with treatment = ploidy and heat as covariate, once with treatment = heat and ploidy as covariate; a DML is a CpG with q < 0.01 and an absolute difference above 20 percentage points; "interaction-DML" are the CpGs in both lists (issue I-02: that is an overlap, not a tested interaction).
Outputs of step 05 you start from (04_output/05_differential_methylation/tables/): united_cpg_<mode>.csv.gz (the tested CpGs: chr, start, end, strand, then coverage<j>, numCs<j>, numTs<j> for j = 1..22, already normalized; the order of j is the row order of read_libraries()); diffmeth_<contrast>_<mode>.csv.gz (every tested CpG with pvalue, qvalue, meth.diff, united_row); dml_{ploidy,heat,shared}_<mode>.csv. Modes: "legacy" (480,579 CpGs; as the draft ran it) and "corrected" (454,555 CpGs; with the 99.9th percentile coverage filter of the Methods).
Read first: issues.md entries I-01, I-02, I-04, I-05, I-06, I-07, I-08, I-09; 01_pipeline/README.md ("How a step is written"); run_config.R; the helpers in 02_R_functions/ (load_packages, input_path, output_path, read_libraries, log_msg, write_table, theme_wgbs, save_figure, mean_percent_methylation, dml_eligible); and 06_diagnostics/README.md and 06_diagnostics/go_random_cpg_control.R as a model of a diagnostics script.
How to run: Rscript is installed (R 4.3.3; methylKit 1.28.0, data.table, ggplot2, patchwork). Export WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 before running.
MACHINE LIMITS: 4 cores and 15 GB of memory shared with other agents working at the same time. A full model fit on all CpGs takes 5 to 10 minutes; the complete set of analyses takes about two hours. YOU MUST NOT RUN THE FULL ANALYSIS. Develop and test on a subset of CpGs only (see the task). Never start a fit on more than 30,000 CpGs. Keep R under 3 GB. Do not run 01_pipeline/04_cpg_coverage.R or 05_differential_methylation.R.
RULES FOR WHAT YOU WRITE: R only. Header block in the style of the existing scripts (Purpose, Run, Reads, Writes, Notes). Paths through here::here and the helpers; no setwd, no absolute paths. Settings from the config list of run_config.R (config$dml_qvalue, config$dml_difference, config$cores, config$modes); do not edit run_config.R. No number typed in as a result. Deterministic: set a seed before anything random. Fail loudly on unexpected input. Plain declarative English in comments and tables; no em dash character (U+2014); the treatment is called "heat", never "desiccation". Cite sources for method choices (authors, year, journal, volume, pages or DOI). Known hazard: inside a data.table call, a function argument or loop variable with the same name as a column refers to the column; name things so that this cannot happen.
Do not run git commands that change state. Do not edit README files or issues.md: return proposed text instead.
`

const SPEC = `Write 06_diagnostics/dml_sensitivity.R: sensitivity analyses of the differential methylation calls. Its purpose is to show the authors how the DML lists change under each defensible alternative that the audit raises, so that they can decide what the paper should report. It is a diagnostics script: nothing in 01_pipeline reads its outputs.

Building the methylKit object. Read united_cpg_<mode>.csv.gz and construct a methylBase object directly (methods::new("methylBase", <data.frame>, sample.ids = ..., assembly = ..., context = "CpG", treatment = ..., coverage.index = ..., numCs.index = ..., numTs.index = ..., destranded = FALSE, resolution = "base"); look at how methylKit::unite builds one and at getSlots("methylBase")). Use methylKit::reorganize() to change the treatment vector or to drop libraries. The library order of the table is the row order of read_libraries().

Analyses. Each is one calculateDiffMeth call; put them in a registry table (id, question it answers, contrast, libraries used, group sizes, covariates, overdispersion, test) that the script loops over, so that adding one is a one-line change.
 A0 primary: exactly step 05 (ploidy with heat as covariate; heat with ploidy as covariate; defaults). This is the equivalence check: on the full data it must reproduce the step 05 DML lists; on the test subset its p-values and differences must equal those of diffmeth_<contrast>_<mode>.csv.gz for the same CpGs (q-values depend on the set of CpGs through SLIM and will differ on a subset).
 A1 overdispersion: the same two contrasts with overdispersion = "MN" and test = "F", and with overdispersion = "MN" and test = "Chisq". The groups are biological replicates, and logistic regression without a dispersion term treats reads as independent draws (issue I-06). Read methylKit's documentation and code to state exactly what each option does, and cite Akalin et al. 2012 (Genome Biology 13:R87) and Wreczycka et al. 2017 (Journal of Biotechnology 261:105-115) for the recommendation only if you have checked what they say.
 A2 interaction, tested: treatment = 1 for triploid heat-treated libraries and 0 otherwise, with ploidy and heat both as covariates. With both main effects in the null model, the test of this indicator is the test of the ploidy by heat interaction (the full model is saturated over the four groups). Run it without overdispersion correction and with "MN" and test = "F". methylKit's meth.diff for this coding is NOT the interaction contrast, so compute the effect size yourself from the pooled group methylation (coverage-weighted, as methylKit's default): (triploid heat - triploid control) - (diploid heat - diploid control), in percentage points. An interaction DML is q < config$dml_qvalue and an absolute interaction contrast above config$dml_difference. Then cross-tabulate with the legacy "shared" loci (dml_shared_<mode>.csv and its class column): how many of them have a significant interaction, with and without the effect-size criterion, and how many interaction DML are not among them. Remember, and say in the header, that under issue I-01 this is a ploidy by cohort interaction.
 A3 without the four libraries that also had the acute heat shock: both contrasts, settings of A0.
 A4 ploidy within one cohort: diploid against triploid among the 12 controls, and among the 10 heat-treated libraries, no covariate; settings of A0 and also "MN" with test = "F". The draft describes its ploidy DML as differences "under control conditions" although they were computed across all 22 libraries (issue I-07); the control-only contrast is what that sentence describes, and both within-cohort contrasts are free of the cohort confound.
 A5 heat within one ploidy: heat against control among diploids, and among triploids; settings of A0. Still a cohort contrast.
 Run every analysis in mode "legacy"; run A0, A1 (MN with F) and A2 in mode "corrected" as well. All on the CpG set of the united table (the 22-library set), so that lists are comparable.
 Dispersion: for a seeded random sample of 5,000 eligible CpGs (dml_eligible(mean_percent_methylation(united), config)), fit a quasi-binomial GLM of the counts on ploidy + heat and record the dispersion estimate (Pearson chi-squared over residual degrees of freedom). Report median, quartiles, and the share above 1.5 and above 3. This says how far the data are from the binomial assumption.

Outputs, all under 06_diagnostics/dml_sensitivity/ (create it): analyses.csv (the registry); dml_counts.csv (per analysis and mode: CpGs tested, DML, hyper, hypo, how many are in the primary list of the same contrast and mode, share of the primary list recovered, share of this list that is in the primary list); interaction_vs_shared_<mode>.csv (the cross-tabulation above); dml_<analysis id>_<mode>.csv (chr, start, pvalue, qvalue, effect size, direction; only CpGs that pass); dispersion_summary.csv; figure_dml_sensitivity (PNG and PDF via save_figure: number of DML per analysis for each contrast, with the part shared with the primary list shown). Model fits are cached as RDS under 06_diagnostics/dml_sensitivity/cache/ (git-ignored), one file per analysis and mode, so that an interrupted run resumes; the script says when it uses a cached fit.

Test mode. When the environment variable WGBS_SENSITIVITY_SUBSET is set to a number N, the script uses a seeded random sample of N CpGs (rows) of each united table, sorted by row. In test mode all outputs and caches go to 06_diagnostics/dml_sensitivity/test_subset/ (git-ignored) instead, and the header of every CSV-producing section must not confuse the two. Use WGBS_SENSITIVITY_SUBSET=20000 for all your testing. The full run will be started by someone else after review.

Also return: a paragraph for 06_diagnostics/README.md describing the script, and proposed text for issues.md (which issues its outputs will inform and how to read each table), without numbers from the subset run, since those are not results.`

const BUILD_SCHEMA = {
  type: 'object',
  properties: {
    files_written: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string' },
    design_decisions: { type: 'array', items: { type: 'string' }, description: 'Each statistical or coding decision you made beyond the specification, with the reason and the source.' },
    validation: { type: 'array', items: { type: 'string' }, description: 'Each check you ran on the subset and its outcome, including the A0 equivalence check with numbers.' },
    caveats: { type: 'array', items: { type: 'string' } },
    estimated_full_run_minutes: { type: 'number' },
    readme_text: { type: 'string' },
    issues_text: { type: 'string' },
  },
  required: ['files_written', 'summary', 'design_decisions', 'validation', 'caveats', 'estimated_full_run_minutes', 'readme_text', 'issues_text'],
}

const REVIEW_SCHEMA = {
  type: 'object',
  properties: {
    runs_clean_on_subset: { type: 'boolean' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          line: { type: 'integer' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          summary: { type: 'string' },
          failure_scenario: { type: 'string' },
          evidence: { type: 'string' },
          suggested_fix: { type: 'string' },
        },
        required: ['severity', 'summary', 'failure_scenario', 'evidence', 'suggested_fix'],
      },
    },
    independent_checks: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
  required: ['runs_clean_on_subset', 'findings', 'independent_checks', 'notes'],
}

phase('Build')
const build = await agent(`${CONTEXT}
Your scratch directory: ${SCRATCH}/build-sensitivity
FILES YOU MAY CREATE OR CHANGE: 06_diagnostics/dml_sensitivity.R and files under 06_diagnostics/dml_sensitivity/ (written by running the script in test mode, which puts them under test_subset/).

YOUR TASK
${SPEC}

When you are done the script must run cleanly in test mode from a fresh R process (WGBS_SENSITIVITY_SUBSET=20000), twice, with identical CSV outputs.`, { label: 'build:dml_sensitivity', phase: 'Build', schema: BUILD_SCHEMA })

phase('Review')
const review = await agent(`${CONTEXT}
Your scratch directory: ${SCRATCH}/review-sensitivity
You may run the script in test mode only (WGBS_SENSITIVITY_SUBSET=20000 Rscript 06_diagnostics/dml_sensitivity.R); do not edit any file.

YOUR TASK: independent statistical and code review of 06_diagnostics/dml_sensitivity.R, written by another agent to this specification:
${SPEC}

The builder's report (claims to check, not facts):
${JSON.stringify(build, null, 1)}

Assume it has defects. Check, with your own code in your scratch directory where a computation decides the point:
1. The methylBase construction: column indices, sample order, that getData() of the object equals the input table, that reorganize() keeps counts aligned with the right libraries.
2. Treatment and covariate coding of every analysis against the library sheet: which libraries are 1 and 0, which are dropped, what the covariates are. For A4 and A5 confirm the subsets and group sizes. For A3 confirm exactly the four acute heat-shock libraries are removed.
3. A0 equivalence on the subset: p-values and meth.diff equal to the step 05 table for the same CpGs (use united_row or chr and start to align).
4. A2: verify by an independent fit that the test is the interaction test: for 200 CpGs of the subset fit glm(cbind(numCs, numTs) ~ ploidy * heat, binomial) yourself and compare the likelihood-ratio p-value of the interaction term with the script's p-value (for the no-overdispersion run and the test methylKit uses there; work out from methylKit's code which reference distribution that is). Verify the interaction contrast formula and its sign convention on a few CpGs by hand.
5. Overdispersion options: does the script's description of "MN" with "F" and with "Chisq" match methylKit's code? Are the citations accurate as far as you can tell?
6. The DML criterion (q below config$dml_qvalue, absolute effect above config$dml_difference: strict or non-strict inequality, and is it the same as step 05?), the overlap statistics, the cross-tabulation with the shared loci.
7. The dispersion estimate: formula, the sample of CpGs, handling of CpGs with zero variance or perfect separation.
8. Test mode: is the sample seeded and identical between runs; do test outputs go only to test_subset/; can a test-mode cache ever be read by a full run or the reverse?
9. Caching: what happens if the united table or the settings change after a fit was cached (stale cache)? Is there a guard?
10. Rules: header complete and accurate, helpers used, no absolute paths, no em dash character, no "desiccation", deterministic outputs (run twice and compare md5 of the CSVs).
Report only real problems with a concrete failure scenario and evidence.`, { label: 'review:dml_sensitivity', phase: 'Review', schema: REVIEW_SCHEMA })

phase('Fix')
const final = await agent(`${CONTEXT}
Your scratch directory: ${SCRATCH}/fix-sensitivity
FILES YOU MAY CREATE OR CHANGE: 06_diagnostics/dml_sensitivity.R and files under 06_diagnostics/dml_sensitivity/test_subset/ (written by running the script in test mode).

YOUR TASK: bring 06_diagnostics/dml_sensitivity.R to a finished state. Specification:
${SPEC}

Builder's report:
${JSON.stringify(build, null, 1)}

Reviewer's report:
${JSON.stringify(review, null, 1)}

For each reviewer finding: confirm it yourself first; fix it if real; if it is not a defect, leave the code and say why. Then run the script in test mode twice from a fresh R process (WGBS_SENSITIVITY_SUBSET=20000), confirm identical CSV outputs, open the figure with the Read tool, and delete nothing outside test_subset/. Do not run the full analysis. In your report give, under validation, each finding with what you did, and each check with its outcome; give final readme_text and issues_text; and give an honest estimate of the full run time on 4 cores based on your subset timings.`, { label: 'fix:dml_sensitivity', phase: 'Fix', schema: BUILD_SCHEMA })

return { build, review, final }

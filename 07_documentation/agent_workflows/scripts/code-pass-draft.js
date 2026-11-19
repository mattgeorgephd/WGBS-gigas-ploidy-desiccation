export const meta = {
  name: 'code-pass',
  description: 'Apply the verified code-review findings to the pipeline, rerun everything, explain every changed output, review the changes adversarially, fix',
  phases: [
    { title: 'Fix', detail: 'four agents with disjoint files: steps 04-05; steps 06-07; steps 08-09; step 10' },
    { title: 'Shared code', detail: 'run_pipeline.R, run_config.R, logging, packages, steps 01, 11, 12; config requests of the other agents' },
    { title: 'Rerun', detail: 'the whole pipeline and all diagnostics in place; every changed file explained' },
    { title: 'Review', detail: 'two independent reviewers: the code diff; the outputs against the documentation' },
    { title: 'Repair', detail: 'fix what the reviewers found and rerun what it touches' },
  ],
}

// args: { base: '<commit hash before the code pass>', findings_dir: '<scratch dir with one text file of findings per group>' }
const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'
const BASE = args.base
const FIND = args.findings_dir

const CHANGE_SCHEMA = {
  type: 'object',
  properties: {
    group: { type: 'string' },
    findings_handled: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          finding: { type: 'string', description: 'the finding as named in your findings file (for example R1.1), or a short name for an extra task' },
          action: { type: 'string', enum: ['fixed', 'already_fixed', 'not_a_defect', 'deferred', 'needs_decision'] },
          what_changed: { type: 'string', description: 'what you changed and where; or why no change' },
          files: { type: 'array', items: { type: 'string' } },
        },
        required: ['finding', 'action', 'what_changed', 'files'],
      },
    },
    outputs_changed: {
      type: 'array',
      description: 'every committed output file (04_output, 05_manuscript, 06_diagnostics) whose content will change or changed because of your edits',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          change: { type: 'string', description: 'what differs: a new column, a corrected label, a number' },
          changes_a_number: { type: 'boolean' },
          old_value: { type: 'string' },
          new_value: { type: 'string' },
          reason: { type: 'string' },
        },
        required: ['file', 'change', 'changes_a_number', 'old_value', 'new_value', 'reason'],
      },
    },
    tests_run: { type: 'array', items: { type: 'object', properties: { what: { type: 'string' }, result: { type: 'string' } }, required: ['what', 'result'] } },
    config_requests: {
      type: 'array',
      description: 'settings you need in run_config.R, which you must not edit',
      items: { type: 'object', properties: { name: { type: 'string' }, value: { type: 'string' }, comment: { type: 'string' }, used_in: { type: 'string' } }, required: ['name', 'value', 'comment', 'used_in'] },
    },
    shared_file_requests: {
      type: 'array',
      description: 'edits to files you do not own (exact old text and new text)',
      items: { type: 'object', properties: { file: { type: 'string' }, old_text: { type: 'string' }, new_text: { type: 'string' }, reason: { type: 'string' } }, required: ['file', 'old_text', 'new_text', 'reason'] },
    },
    issues_text: {
      type: 'array',
      description: 'text for issues.md (which you must not edit): statements that your change makes wrong, and what they should say',
      items: { type: 'object', properties: { issue: { type: 'string' }, current_text: { type: 'string' }, proposed_text: { type: 'string' }, reason: { type: 'string' } }, required: ['issue', 'current_text', 'proposed_text', 'reason'] },
    },
    open_points: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string' },
  },
  required: ['group', 'findings_handled', 'outputs_changed', 'tests_run', 'config_requests', 'shared_file_requests', 'issues_text', 'open_points', 'summary'],
}

const RERUN_SCHEMA = {
  type: 'object',
  properties: {
    steps_run: { type: 'array', items: { type: 'object', properties: { what: { type: 'string' }, status: { type: 'string' }, minutes: { type: 'number' }, peak_memory_gb: { type: 'number' }, note: { type: 'string' } }, required: ['what', 'status', 'minutes', 'peak_memory_gb', 'note'] } },
    changed_files: {
      type: 'array',
      description: 'every tracked file that differs from the base commit after the rerun, and every new untracked file',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          kind: { type: 'string', enum: ['code', 'output_table', 'output_figure', 'documentation', 'other'] },
          what_differs: { type: 'string' },
          changes_a_number: { type: 'boolean' },
          explained_by: { type: 'string', description: 'the finding or task that explains it, or UNEXPLAINED' },
        },
        required: ['file', 'kind', 'what_differs', 'changes_a_number', 'explained_by'],
      },
    },
    legacy_check: { type: 'string', description: 'result of 06_diagnostics/check_legacy_outputs.R: comparisons, matches, each difference' },
    manuscript_values: { type: 'string', description: 'result of step 13: counts by status before and after, and every row whose status or regenerated value changed' },
    second_run_identical: { type: 'string', description: 'which outputs were rerun a second time and whether they were byte-identical' },
    problems: { type: 'array', items: { type: 'string' } },
    issues_text: { type: 'array', items: { type: 'object', properties: { issue: { type: 'string' }, current_text: { type: 'string' }, proposed_text: { type: 'string' }, reason: { type: 'string' } }, required: ['issue', 'current_text', 'proposed_text', 'reason'] } },
    summary: { type: 'string' },
  },
  required: ['steps_run', 'changed_files', 'legacy_check', 'manuscript_values', 'second_run_identical', 'problems', 'issues_text', 'summary'],
}

const REVIEW_SCHEMA = {
  type: 'object',
  properties: {
    review: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          problem: { type: 'string' },
          failure_scenario: { type: 'string' },
          evidence: { type: 'string' },
          proposed_fix: { type: 'string' },
          changes_a_committed_number: { type: 'boolean' },
        },
        required: ['file', 'line', 'severity', 'problem', 'failure_scenario', 'evidence', 'proposed_fix', 'changes_a_committed_number'],
      },
    },
    issues_text: { type: 'array', items: { type: 'object', properties: { issue: { type: 'string' }, current_text: { type: 'string' }, proposed_text: { type: 'string' }, reason: { type: 'string' } }, required: ['issue', 'current_text', 'proposed_text', 'reason'] } },
    verdict: { type: 'string' },
  },
  required: ['review', 'findings', 'issues_text', 'verdict'],
}

const COMMON = `
Project. ${REPO} is an R pipeline that reproduces and audits the analysis of a manuscript draft on DNA methylation (WGBS) in diploid and triploid Pacific oysters. Read ${REPO}/README.md, ${REPO}/01_pipeline/README.md and ${REPO}/02_R_functions/README.md first. Numbered steps live in 01_pipeline/ (01 to 13), shared functions in 02_R_functions/, settings in run_config.R, inputs in 03_input_files/, outputs in 04_output/<step>/{tables,figures}, diagnostics in 06_diagnostics/, the audit log in issues.md. Two analysis modes run side by side: "legacy" (the analysis as the draft ran it; its outputs must keep reproducing the legacy results, see 06_diagnostics/check_legacy_outputs.R and 06_diagnostics/legacy_check_results.csv, 329 of 331 comparisons match) and "corrected".

What this workflow does. Six independent code reviews found defects; a second agent re-examined each finding and confirmed or refuted it. You apply the confirmed findings of your group. Your findings file holds, for each finding, the reviewer's text and the second agent's verdict with its tested fix. Findings marked refuted need no change unless you find the refutation wrong. A confirmed finding may already be fixed in the tree (the reviews are some hours old): check the current file before you edit.

Rules for every agent
- Base commit of this code pass: ${BASE}. "git diff ${BASE} -- <file>" shows what has changed since.
- Edit only the files your task lists as yours. For run_config.R and for any other file you do not own, return the edit you need under config_requests or shared_file_requests. Never edit issues.md; return text for it under issues_text. In README files that several agents touch (01_pipeline/README.md, 04_output/README.md, 02_R_functions/README.md, 06_diagnostics/README.md) change only the rows and sentences about your own steps, and re-read the file immediately before each edit because another agent may have changed it.
- Run no git command that changes state (no add, commit, checkout, stash, reset, clean). The parent session commits.
- Legacy mode must keep reproducing the legacy numbers. If a fix would change a number in legacy mode, stop and think: a fault of the legacy analysis is to be REPRODUCED in legacy mode and documented, not repaired there; the corrected mode is where repairs go. A change to a legacy-mode number is acceptable only when the pipeline, not the legacy analysis, was wrong.
- Every number that changes in a committed output must be listed under outputs_changed with old and new value and the reason. No silent changes.
- Style of code: base R and data.table as in the existing steps; each script starts with a header (Purpose, Run, Reads, Writes, Notes) that must stay true after your edit; comments say why, not what; no em dashes anywhere (use commas or semicolons); treated oysters are "heat" or "heat-treated", the two experiments are the "heat cohort" and the "control cohort". Keep output column names and file names stable unless a finding requires a change; if you rename one, grep for every reader (01_pipeline, 02_R_functions, 06_diagnostics, 03_input_files/manuscript/manuscript_values_map.csv, README files) and report them.
- The machine has 4 cores and 15 GB, and other agents work on it. Use WGBS_CACHE_DIR=/home/user/gannet_cache (the downloaded inputs are there). Check "free -g" before anything that needs more than 2 GB. Only steps 04 and 05 are heavy (up to about 10 GB; step 05 fits models for tens of minutes): the task of each agent says whether it may run them. Never run two heavy jobs at once; use at most WGBS_CORES=2.
- Scratch space: ${SCRATCH}/code-pass/<your group>/ (create it). Keep test scripts and logs there.
- Report honestly: a test you did not run is "not run", a finding you could not fix is "deferred" with the reason.
`

const GROUPS = {
  A: {
    label: 'fix:A-steps-04-05',
    prompt: `${COMMON}
Your group: A, steps 04 and 05, the sample sheet, and the two diagnostics that share their counting rule.
Your files: 01_pipeline/04_cpg_coverage.R, 01_pipeline/05_differential_methylation.R, 02_R_functions/methylkit_helpers.R, 02_R_functions/sample_sheet.R, 06_diagnostics/dml_sensitivity.R, 06_diagnostics/library_methylation_level.R, and the rows about steps 04 and 05 and these two diagnostics in the README files.
Findings file: ${FIND}/group_A.txt (code reviews R1 and parts of R5, R6; the check of 06_diagnostics/library_methylation_level.R; loose ends recorded in issues.md under I-16).

Tasks
1. Apply the confirmed findings. The main ones: (a) the model fit cached by step 05 must carry a fingerprint of everything it depends on (settings used by the fit, sample sheet, the united table or its inputs, methylKit version, the fitting code) and must be refitted when the fingerprint differs or the cache file cannot be read; (b) one rule, in one helper of 02_R_functions/methylkit_helpers.R, for "difference above the threshold", used by steps 05, 06 and 07 (ask groups B for 06 and 07 through shared_file_requests if they must call it) and by both diagnostics, so that a CpG exactly at the threshold is treated the same everywhere and a count retaken from a written table equals the count made in memory; legacy-mode lists must stay exactly the legacy lists; (c) step 04's rounding order against methylKit's; (d) the count of CpGs covered in at least one library (12,918,252 when every line is read) written to a table; (e) read_libraries() must stop on unexpected labels; (f) headers and README rows.
2. Issue I-15 (read it in issues.md): Results 3.1 of the draft gives a genomic distribution of CpGs "at 1x". Extend step 04 to write, for its CpG sets (each coverage threshold, each sample set, both modes), the overlap with the feature tracks of 03_input_files/genome_features using the functions of 02_R_functions/genomic_overlaps.R (owned by group B; do not edit it), counting each CpG once per feature (not overlap events), so that a reader can see which CpG set, if any, gives the draft's percentages. Say in the header what the draft's numbers match.
3. 06_diagnostics/dml_sensitivity.R: write two more tables at the end of a run, from the fits it already has: for every analysis that runs in both modes, the CpGs with q below the cutoff in legacy mode split by whether the coverage filter keeps them (and how many are on the mitochondrial genome NC_001276.1); and the DML counts under the Benjamini-Hochberg adjustment beside those under SLIM. issues.md (I-06, item 6) quotes these from a scratch computation: 157 of 287, 138 mitochondrial, 15 and 4. Then run the script (the fits are cached under 06_diagnostics/dml_sensitivity/cache/, so a run takes a few minutes) and confirm that every table it wrote before is unchanged.
4. 06_diagnostics/library_methylation_level.R: apply the code findings of the independent check (in your findings file), switch it to the shared threshold helper, rerun it, and report every number that changes.

Running. You MAY run step 04 in place (about 8 minutes, 8 GB): nobody else reads its outputs during this phase. You may NOT run step 05 in place: other agents are reading its outputs. To test step 05, copy the working tree without .git, 04_output and 08_archive to your scratch directory (the .Rproj file makes it a project root for here::here()), run step 05 there for both modes with WGBS_CORES=2, and compare every table it writes with the committed table of the same name in ${REPO}/04_output/05_differential_methylation/tables/ (cmp or md5sum). Test the fingerprint: an unchanged rerun must reuse the cache; a changed setting, a changed sample sheet and a truncated cache file must each cause a refit.

Return group "A".`,
  },
  B: {
    label: 'fix:B-steps-06-07',
    prompt: `${COMMON}
Your group: B, steps 06 and 07 and the genome feature tracks.
Your files: 01_pipeline/06_dml_figures.R, 01_pipeline/07_genome_features.R, 02_R_functions/genomic_overlaps.R, 02_R_functions/plot_helpers.R if it exists, 03_input_files/genome_features/README.md, 06_diagnostics/feature_odds_matched.R, and the rows about steps 06 and 07 in the README files.
Findings file: ${FIND}/group_B.txt (code review R2 and its verification; related header points of R6).

Tasks
1. Apply the confirmed findings. The main ones: (a) steps 06 and 07 must not fail or write a meaningless test when a DML set is empty or has fewer than the loci a figure needs: skip that panel or row with a logged message and still write complete, well-formed outputs; (b) the track called TE is the RepeatMasker output of the assembly with 1-based closed coordinates and all repeat classes, read as if it were BED. Find out how the LEGACY analysis read this file (08_archive/code, the notebooks and bedtools calls) and keep that reading in legacy mode if the legacy numbers depend on it; in corrected mode read the coordinates correctly. Describe the track for what it is in the README and in the outputs' labels. Report every count that changes (the verification found 101 CpGs affected) and every statement that names "transposable elements"; (c) the three wrong descriptions in the track README; (d) the header of step 07 and the README must name every file the step writes; (e) "hyper" and "hypo" of the loci in both lists are the sign of the heat difference: say so in the column documentation, the figure captions and the README; (f) the label and legend defects of Figures 2 and 3.
2. 06_diagnostics/feature_odds_matched.R was written by another agent today and reads the same tracks: after your change to the track reading, rerun it and report what changes in its four tables.
3. If group A asks (it owns a new helper for "difference above the threshold" in 02_R_functions/methylkit_helpers.R), steps 06 and 07 must use it. Check ${REPO}/02_R_functions/methylkit_helpers.R when you start and again before you finish; if the helper is there, use it; if it is not there yet, say so under open_points.

Running. You may run steps 06 and 07 in place (minutes, light). Do not run steps 04 or 05. After your edits, legacy-mode tables that the findings do not touch must be byte-identical to the committed ones (compare with "git diff --stat"); run 06_diagnostics/check_legacy_outputs.R and report its count of matches (329 of 331 before the code pass).

Return group "B".`,
  },
  C: {
    label: 'fix:C-steps-08-09',
    prompt: `${COMMON}
Your group: C, steps 08 and 09 (gene annotation and the legacy GO enrichment).
Your files: 01_pipeline/08_gene_annotation.R, 01_pipeline/09_go_enrichment_legacy.R, 02_R_functions/annotation.R, 03_input_files/annotation/README.md, 06_diagnostics/legacy_reader_mechanism.R, 06_diagnostics/universe_file_cut.R, and the rows about steps 08 and 09 in the README files.
Findings file: ${FIND}/group_C.txt (code review R3 and its verification; points of R5 and R6 about steps 08 and 09; points of the integration checks).

Tasks
1. Apply the confirmed findings. The main ones: (a) the what_it_counts definitions of "GO terms" and "Unique genes" that step 09 writes to figure_5_values_legacy.csv are the wording that issues.md withdrew (see issue I-20, the table in its Update paragraph, for the corrected definitions); (b) the header of step 08 on the genes whose name holds a hyphen (5 genes, 24 mRNAs); (c) "Used by step 09 only" in the annotation README; (d) with modes = c("corrected") in run_config.R steps 08 and 09 exit with status 0 and write wrong cross-mode tables: step 09 is a legacy-only step and must say so and stop cleanly (or do nothing) when legacy mode is not configured; step 08 must write correct tables for whichever modes are configured; (e) legacy_universe_by_chromosome.csv and the log line that the reviews mention; (f) headers and README rows that lag behind the code.
2. After your edits rerun steps 08 and 09 in place (minutes) and show that every table is byte-identical to the committed one except where a finding changes it (list those). Run 06_diagnostics/check_legacy_outputs.R and report its count of matches.
3. Test the modes setting without editing run_config.R: copy the working tree without .git and 08_archive to your scratch directory, set modes to c("corrected") in the copy, and run steps 08 and 09 there.

Do not run steps 04 or 05.

Return group "C".`,
  },
  D: {
    label: 'fix:D-step-10',
    prompt: `${COMMON}
Your group: D, step 10 (the gene-level GO enrichment that the audit proposes in place of the legacy one) and its negative control.
Your files: 01_pipeline/10_go_enrichment_gene_level.R, 02_R_functions/go_tests.R, 06_diagnostics/go_random_cpg_control.R, the section on the random-CpG control in 06_diagnostics/README.md, and the rows about step 10 in the README files.
Findings file: ${FIND}/group_D.txt (code review R4 and its verification; the reports C4a and C4b of the first verification, which re-implemented step 10 and proposed changes that issues.md lists under "Status of the pipeline" in issue I-41: read that paragraph).

Tasks
1. Apply the confirmed findings. The main one changes a headline result, so be exact: go_gene_level_object() hands 2017 GO IDs to topGO unchanged and topGO drops every ID that is not a primary ID of the installed GO.db; 395 of the dropped IDs are secondary IDs of current terms. Remap secondary IDs to their primary ID (GO.db) before the tests, record how many IDs and gene-term pairs were remapped and how many remain unknown, and write that to a table. The reviewer found that with the remapping the step's own test calls "organelle" for heat, CC, in both modes (q about 0.018), and that universes and several documented counts shift. Reproduce that and report every changed count.
2. The other findings: pin or check the GO.db release that step 10 uses and record it in the outputs; remove a stale figure_go_enriched from an earlier run when no term is enriched; the retest of the draft's terms must not be hard-wired to legacy mode (stop cleanly or adapt when legacy is not configured); Wallenius p-values above 1; the expected and enrichment columns (unadjusted, rounded Fisher values written next to a Wallenius q-value); the pmax() that does nothing; header and log slips.
3. From issue I-41, "Status of the pipeline": add an exact conditional test beside the Wallenius approximation, and use the logarithm of the number of eligible CpGs as the bias variable if that is what the verification recommended; keep the old columns so that readers can compare. Read the two verification scripts (07_documentation/agent_workflows/verification_scripts/C4a-eligible-theory and C4b-eligible-reimplementation) before you write this.
4. Rerun step 10 in place (both modes) and 06_diagnostics/go_random_cpg_control.R (about 10 minutes). Correct the sentence under the control's result table in 06_diagnostics/README.md that the review names. Step 11 reads no output of step 10; do not run it.
5. Return under issues_text the replacements for every statement in issues.md, README.md and 07_documentation/reproducibility_plan.md that says no term passes at the configured settings (search for "organelle" and "no term" and "No specific").

Do not run steps 04 or 05. You start after group C has finished, so the outputs of step 08 are final.

Return group "D".`,
  },
}

const E_PROMPT = (requests) => `${COMMON}
Your group: E, the code that every step shares, and steps 01, 11, 12 and 13.
Your files: run_pipeline.R, run_config.R, 02_R_functions/logging.R, 02_R_functions/packages.R, 02_R_functions/paths.R, 02_R_functions/legacy_go_db.R, 01_pipeline/01_fetch_inputs.R, 01_pipeline/11_semantic_similarity.R, 01_pipeline/12_qpcr.R, 01_pipeline/13_manuscript_values.R, a new .gitattributes at the repository root, and the README statements about these files (README.md at the root, 01_pipeline/README.md, 02_R_functions/README.md, 03_input_files/README.md, 04_output/README.md).
Findings file: ${FIND}/group_E.txt (code review R5 and its verification; the write_table finding of R3; documentation findings of R6 that concern these files).

Tasks
1. Apply the confirmed findings. The main ones: (a) run_config.R must reject values that cannot work, in one validation block that every step passes through because every step sources this file: modes (non-empty subset of legacy and corrected), cores (WGBS_CORES empty or not a positive whole number must stop at once with a plain message), thresholds; (b) write_table() must write the same bytes whatever the platform and the user's options (line ending fixed to LF, scipen fixed), and .gitattributes must keep git from converting line endings of the tracked tables; confirm on Linux that every committed table is reproduced byte for byte; (c) run_pipeline.R: started from outside the repository it must find the repository or stop with an error, never report success after running nothing; --steps must reject malformed selections; (d) load_packages() must name R.utils when a step reads .gz tables with fread; (e) an interrupted install of the GO.db side library must not block later runs without an explanation; (f) step 01: input_inventory.csv must not depend on the collation locale, and the committed copy must be the current one; (g) the statements that the numbered folder names live only in paths.R.
2. Steps 11 and 12 each open with a block of settings: move the settings into run_config.R (one list per step, with the comments), keep the scripts reading them from config, and show that both steps write byte-identical tables afterwards (step 11 takes about 5 minutes, step 12 under a minute). Both steps read legacy files from 08_archive/ through archive_path(): copy those files into the fitting folder under 03_input_files/ (with a row in that folder's README saying where each came from), point the steps there, and leave 08_archive untouched.
3. Apply the requests of the other groups for run_config.R and for files that none of them owns. Their requests:
${requests}
   Judge each request: apply it, adapt it, or decline it with a reason (report which).
4. Do not run steps 04 or 05. You may run steps 01 (it only checks the cache when the downloads are present), 02, 03, 11, 12 and 13 in place.

Return group "E".`

const RERUN_PROMPT = (summaries) => `${COMMON}
Your job: rerun the whole pipeline and all diagnostics in place after the code pass, and explain every file that changed. You fix nothing; you report. Five agents have edited the code since commit ${BASE}; their summaries:
${summaries}

Steps
1. "git status --short" and "git diff --stat ${BASE}" to see what the agents left. Note every output file that already differs.
2. Run the pipeline in order from the repository root: WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2 Rscript run_pipeline.R (all 13 steps; step 05 will refit its models because its cache now carries a fingerprint, which takes about an hour; run it in the background with nohup and poll its log, and do not start anything else heavy meanwhile). Record minutes per step, and the peak memory of steps 04 and 05 (for example by sampling "ps -o rss" of the R process every 10 seconds into a log).
3. Run every diagnostic of 06_diagnostics/README.md, section "Files", in the order given there, including dml_sensitivity.R (its fits are cached under 06_diagnostics/dml_sensitivity/cache/; if step 05's united tables are unchanged the cache is valid, and the script says what it reuses), library_methylation_level.R, feature_odds_matched.R, locus_nc047563.R, legacy_reader_mechanism.R, coverage_filter_variants.R and check_methylkit_equivalence.R if they exist, go_random_cpg_control.R, table1_group_means.R, cohort_separation.R, universe_file_cut.R, and last check_legacy_outputs.R.
4. "git status --short" and "git diff --stat ${BASE}" again. For EVERY tracked file that differs from ${BASE} and every untracked file: say what differs (for tables: which rows, which columns, old and new values; use R or Python to diff tables, do not eyeball) and which finding or task of the summaries above explains it. Anything you cannot tie to a finding is UNEXPLAINED and goes under problems. session_info.txt files and PDF figures change on every run; list them once as a class.
5. Step 13: compare 05_manuscript/manuscript_values.csv before (git show ${BASE}:05_manuscript/manuscript_values.csv) and after: counts by status, and every row whose status or regenerated value changed, with the reason. A row that turned to an error because a table or column was renamed is a problem to report with the fix for 03_input_files/manuscript/manuscript_values_map.csv.
6. Reproducibility: run steps 06 to 13 a second time and confirm that no tracked table changes between the first and second run.
7. Under issues_text give replacements for statements of issues.md that the changed numbers make wrong (search issues.md for the old values).

Return the result in the schema. Be exact about what you did not run and why.`

const REVIEW_CODE_PROMPT = `${COMMON}
Your job: an adversarial review of the code changes of this code pass. You fix nothing. "git diff ${BASE} -- 01_pipeline 02_R_functions 06_diagnostics run_pipeline.R run_config.R .gitattributes" is the change set (the working tree is not committed; compare the working tree with ${BASE}). For each changed file read the whole file as it now stands, not only the diff.
Look for: a fix that does not do what its finding asked (the findings are in ${FIND}/group_A.txt to group_E.txt); regressions; a legacy-mode number changed where the legacy analysis, not the pipeline, was at fault; a guard that hides an error instead of reporting it; settings read from config that run_config.R does not define, or defined and never read; headers and README rows that no longer match the code; the cache fingerprint of step 05 (try to construct a change that it misses); the shared threshold helper (is it used everywhere a difference is compared with the threshold? grep for comparisons with dml_difference and with the literal 20); run_config.R validation (try bad values in a scratch copy of the tree); run_pipeline.R argument handling (try it); write_table() (bytes); the secondary GO ID remapping of step 10 (is a term counted twice for a gene after remapping? are obsolete IDs handled?).
You may run light steps in a scratch copy of the working tree (copy without .git and 08_archive). Do not run steps 04 or 05. Report each defect with a concrete failure scenario. Return review "code".`

const REVIEW_DOCS_PROMPT = (rerun) => `${COMMON}
Your job: an adversarial check that the documentation and the audit log agree with the outputs as they now stand after the code pass. You fix nothing. The rerun agent reported:
${rerun}

Check
1. Every README under ${REPO} that describes pipeline steps, outputs, inputs or diagnostics (root README.md, 01_pipeline/README.md, 02_R_functions/README.md, 03_input_files/README.md and its subfolders, 04_output/README.md, 05_manuscript/README.md, 06_diagnostics/README.md): each file named must exist, each column described must exist, each number quoted must equal the table it comes from, each run time and memory figure must be the one the rerun measured, each status must be current.
2. The headers (Purpose, Run, Reads, Writes, Notes) of every script in 01_pipeline and 06_diagnostics against what the script reads and writes now.
3. issues.md: for every number that the code pass changed (the rerun report lists them), find every statement of issues.md, README.md and 07_documentation/reproducibility_plan.md that quotes the old value or depends on it, and give the replacement under issues_text. Then take 60 numbers of issues.md at random from different issues, each with a named source file under 04_output or 06_diagnostics, and check them against the files as they now stand.
4. 03_input_files/manuscript/manuscript_values_map.csv against the outputs: every expression must still evaluate (step 13 writes a status per row in 05_manuscript/manuscript_values.csv), and notes that quote numbers must still be right.
Report each mismatch as a finding (file, line, what is wrong, the correct text). Return review "docs".`

const REPAIR_PROMPT = (reviews) => `${COMMON}
Your job: repair what two independent reviewers found after the code pass, then rerun what your repairs touch. You may edit any file under 01_pipeline, 02_R_functions, 06_diagnostics, run_pipeline.R, run_config.R, and the README files; not issues.md (return text for it). The reviewers' findings:
${reviews}

For each finding: check it against the current file, fix it if it holds, and say what you did. Then rerun the affected steps in place (steps 04 and 05 only if a repair changes what they compute; if so run them alone with WGBS_CORES=2 and say so), rerun 06_diagnostics/check_legacy_outputs.R and step 13, and list every output that changed with old and new values. Return group "repair".`

function brief(r) {
  if (!r) return '(no result)'
  const f = (r.findings_handled || []).map(x => `  - ${x.finding}: ${x.action}. ${x.what_changed}`).join('\n')
  const o = (r.outputs_changed || []).map(x => `  - ${x.file}: ${x.change} (${x.old_value} -> ${x.new_value}); ${x.reason}`).join('\n')
  return `Group ${r.group}. ${r.summary}\n Findings:\n${f}\n Outputs changed:\n${o}\n Open points: ${(r.open_points || []).join(' | ')}`
}
function requestsOf(rs) {
  const lines = []
  for (const r of rs.filter(Boolean)) {
    for (const c of r.config_requests || []) lines.push(`   - [group ${r.group}] run_config.R: ${c.name} = ${c.value}  # ${c.comment} (used in ${c.used_in})`)
    for (const s of r.shared_file_requests || []) lines.push(`   - [group ${r.group}] ${s.file}: replace <<${s.old_text}>> with <<${s.new_text}>>; reason: ${s.reason}`)
  }
  return lines.length ? lines.join('\n') : '   (none)'
}

phase('Fix')
// Lane 1: A alone. Lane 2: B, then C, then D (D needs C's outputs of step 08; B and C share README files only).
const laneA = agent(GROUPS.A.prompt, { label: GROUPS.A.label, phase: 'Fix', schema: CHANGE_SCHEMA })
const lane2 = (async () => {
  const b = await agent(GROUPS.B.prompt, { label: GROUPS.B.label, phase: 'Fix', schema: CHANGE_SCHEMA })
  const c = await agent(GROUPS.C.prompt, { label: GROUPS.C.label, phase: 'Fix', schema: CHANGE_SCHEMA })
  const d = await agent(GROUPS.D.prompt, { label: GROUPS.D.label, phase: 'Fix', schema: CHANGE_SCHEMA })
  return [b, c, d]
})()
const a = await laneA
const [b, c, d] = await lane2
log(`Fix phase done: ${[a, b, c, d].filter(Boolean).length} of 4 groups returned`)

phase('Shared code')
const e = await agent(E_PROMPT(requestsOf([a, b, c, d])), { label: 'fix:E-shared-code', phase: 'Shared code', schema: CHANGE_SCHEMA })

phase('Rerun')
const summaries = [a, b, c, d, e].map(brief).join('\n\n')
const rerun = await agent(RERUN_PROMPT(summaries), { label: 'rerun:whole-pipeline', phase: 'Rerun', schema: RERUN_SCHEMA })

phase('Review')
const rerunText = rerun ? JSON.stringify({ changed_files: rerun.changed_files, legacy_check: rerun.legacy_check, manuscript_values: rerun.manuscript_values, problems: rerun.problems, summary: rerun.summary }, null, 1) : '(the rerun agent returned nothing; do the comparison with git diff yourself)'
const reviews = await parallel([
  () => agent(REVIEW_CODE_PROMPT, { label: 'review:code-diff', phase: 'Review', schema: REVIEW_SCHEMA }),
  () => agent(REVIEW_DOCS_PROMPT(rerunText), { label: 'review:docs-and-outputs', phase: 'Review', schema: REVIEW_SCHEMA }),
])

phase('Repair')
const allFindings = reviews.filter(Boolean).flatMap(r => r.findings.map(f => ({ review: r.review, ...f })))
const unexplained = rerun ? rerun.problems : []
let repair = null
if (allFindings.length || unexplained.length) {
  const text = JSON.stringify({ reviewer_findings: allFindings, problems_of_the_rerun: unexplained }, null, 1)
  repair = await agent(REPAIR_PROMPT(text), { label: 'repair', phase: 'Repair', schema: CHANGE_SCHEMA })
} else {
  log('No reviewer finding and no unexplained change: repair stage skipped')
}

return { groups: { a, b, c, d, e }, rerun, reviews, repair }

export const meta = {
  name: 'check-issues-integration',
  description: 'Independent check that issues.md carries every error and caveat of the nine verification reports faithfully, plus a cross-file consistency check',
  phases: [
    { title: 'Check integration', detail: 'one agent per verification report compares it with issues.md' },
    { title: 'Cross-file consistency', detail: 'summary table, corrections table, README, plan, code comments' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const RESULTS = REPO + '/07_documentation/agent_workflows/results/audit-verification.json'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'

const REPORTS = [
  { id: 'C1-universe', issues: 'I-39 (also I-03, I-28)' },
  { id: 'C2-reader', issues: 'I-40, I-43 (also I-39)' },
  { id: 'C3-ontology-release', issues: 'I-37, I-28, I-42, I-43 (the REVIGO point belongs to I-22, which is deliberately not yet updated; report it under not_yet_due, not as an omission)' },
  { id: 'C4a-eligible-theory', issues: 'I-41, I-03 (the point about the number of ploidy DML under label permutation belongs to I-06, which is deliberately not yet updated and is mentioned in I-41 only as a pointer; the code changes it proposes for go_tests.R and go_random_cpg_control.R are deliberately deferred and listed under "Status of the pipeline" in I-41)' },
  { id: 'C4b-eligible-reimplementation', issues: 'I-41 (code changes to step 10 are deliberately deferred and listed under "Status of the pipeline" in I-41)' },
  { id: 'C5-feature-odds', issues: 'I-41, I-14 (the header fix for 07_genome_features.R outputs is deferred to the code-review pass)' },
  { id: 'C6-coverage-filter-locus', issues: 'I-05, I-09, I-38 (the missing script check_methylkit_equivalence.R is deferred; the comment in methylkit_helpers.R was corrected)' },
  { id: 'C7-legacy-go-semantics', issues: 'I-02, I-07, I-20, I-21, I-42, I-45' },
  { id: 'C8-counts-and-correlations', issues: 'I-16, I-44, I-26 (code changes to step 04 are deferred and recorded as loose ends at the end of the I-16 Verification paragraph)' },
]

const SCHEMA = {
  type: 'object',
  properties: {
    report_id: { type: 'string' },
    items_checked: { type: 'integer', description: 'number of errors_found plus missing_caveats items of the report that you compared with issues.md' },
    carried_correctly: { type: 'integer' },
    discrepancies: {
      type: 'array',
      description: 'statements now in issues.md (or README.md, the plan, code comments) that misstate the report or the data: a number changed, a meaning changed, a qualifier lost, a claim stronger than the report supports, or an internal contradiction',
      items: {
        type: 'object',
        properties: {
          where: { type: 'string', description: 'file, issue ID and enough of the sentence to find it' },
          quote: { type: 'string', description: 'the exact text as it stands now' },
          problem: { type: 'string' },
          evidence: { type: 'string', description: 'what you checked and what you found: the report item, and where possible a recomputation from files in the repository' },
          proposed_text: { type: 'string', description: 'exact replacement text, without em dashes' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['where', 'quote', 'problem', 'evidence', 'proposed_text', 'severity'],
      },
    },
    omitted_items: {
      type: 'array',
      description: 'items of the report that issues.md does not carry and that a reader of issues.md needs',
      items: {
        type: 'object',
        properties: {
          item: { type: 'string' },
          why_it_matters: { type: 'string' },
          where_to_add: { type: 'string' },
          proposed_text: { type: 'string' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['item', 'why_it_matters', 'where_to_add', 'proposed_text', 'severity'],
      },
    },
    not_yet_due: { type: 'array', items: { type: 'string' }, description: 'report items that belong to parts flagged as deliberately deferred' },
    notes: { type: 'string' },
  },
  required: ['report_id', 'items_checked', 'carried_correctly', 'discrepancies', 'omitted_items', 'not_yet_due', 'notes'],
}

const COMMON = `
You are checking someone else's editing work. The repository ${REPO} holds an audit of a manuscript on DNA methylation in diploid and triploid Pacific oysters. Its issue log, ${REPO}/issues.md, was first written by one analyst. Nine independent agents then tried to refute its main statements and returned reports. The analyst has just edited issues.md to carry those reports in: paragraphs headed "Verification", statements corrected where they stood, new issues I-43 to I-45, and a table "Corrections after verification" at the end. Your job is to find where that editing went wrong.

Rules
- Read-only. Do not modify any file under ${REPO}. Write scratch files only under ${SCRATCH}/<your own new subdirectory>.
- Do not start long computations. The machine is busy with model fits (3 of its 4 cores). Reading files, grep, and short R or Python snippets on tables that already exist under ${REPO}/04_output and ${REPO}/06_diagnostics are fine (keep any single command under about two minutes and under 2 GB of memory). Do not rerun pipeline steps. Do not read the 22 coverage files.
- Be skeptical of both sides. A report can be wrong too: where issues.md deliberately departs from a report, judge which is right from the data if you can check it cheaply, and say so.
- Style of any text you propose: plain declarative sentences, no em dashes (use commas or semicolons), numbers with their units, no vague wording.
- Your final answer is the structured object. Put everything the analyst needs in it; nobody reads your scratch files unless you name them.
`

function reportPrompt(r) {
  return `${COMMON}
Your report: ${r.id}. The issues it bears on: ${r.issues}.

How to read your report: run
  python3 -c "import json; j=json.load(open('${RESULTS}')); r=[x['result'] for x in j['results'] if x['result'].get('claim_id')=='${r.id}'][0]; print(r['verdict']); print(r['summary']); print('ERRORS'); [print('-', json.dumps(e)) for e in r['errors_found']]; print('CAVEATS'); [print('-', c if isinstance(c,str) else json.dumps(c)) for c in r['missing_caveats']]; print('CHECKS'); [print('-', c if isinstance(c,str) else json.dumps(c)) for c in r['checks_run']]"
The verifying agent's own scripts and small result tables are in ${REPO}/07_documentation/agent_workflows/verification_scripts/${r.id}/ (and its full scratch directory, if it still exists, is ${SCRATCH}/${r.id}/).

What to do
1. Read the whole report (summary, every errors_found item, every missing_caveats item).
2. Read the issues named above in ${REPO}/issues.md in full, including their original text, Update paragraphs, Verification paragraphs, and the rows of the table "Corrections after verification" that concern them. Also read the summary table rows at the top for those issues.
3. For every errors_found item and every missing_caveats item decide: carried correctly; carried with an error (report it under discrepancies, with the exact current text and an exact replacement); omitted although a reader needs it (omitted_items); or belonging to a part flagged above as deliberately deferred (not_yet_due). An item need not be carried word for word; it must not be lost or distorted. Items that concern only the manuscript wording and are low value may be left out without being an omission; use judgement and say what you did in notes.
4. Check every number in the Verification paragraphs and corrected sentences of those issues against the report. Where the number can be recomputed in a minute or two from tables under ${REPO}/04_output or ${REPO}/06_diagnostics or from the saved verification tables, recompute it and say so in evidence. Flag any number that appears in issues.md with a different value in two places.
5. Look for overstatement introduced by the editing: a sentence in issues.md that claims more than the report established, drops a limit the report stated (for example "not run under R 4.1", "one agent only", "inferred"), or presents an agent's untested suggestion as a finding.
6. Check that every statement of 2026-11-19 that was changed in place because of this report has a row in "Corrections after verification", and that each row describes the change accurately.

Return only real problems. If the editing is faithful, say so with empty lists. Severity: high = a reader would be misled about a conclusion; medium = a number or qualifier is wrong; low = wording.`
}

const crossPrompt = `${COMMON}
Your task is consistency across files, not fidelity to one report.

Read ${REPO}/issues.md in full (about 850 lines), then ${REPO}/README.md, ${REPO}/07_documentation/reproducibility_plan.md, ${REPO}/06_diagnostics/README.md, ${REPO}/01_pipeline/README.md, ${REPO}/04_output/README.md, ${REPO}/07_documentation/session_handoff.md, and the header comments of ${REPO}/01_pipeline/07_genome_features.R, 08_gene_annotation.R, 09_go_enrichment_legacy.R, 10_go_enrichment_gene_level.R and of ${REPO}/02_R_functions/methylkit_helpers.R and go_tests.R, and ${REPO}/06_diagnostics/go_random_cpg_control.R.

Check and report (use report_id "cross-file"; put problems under discrepancies, missing things under omitted_items):
1. Summary table of issues.md against the sections: every issue I-01 to I-45 has one row and one section; ID, severity and title agree in substance; the status in the table is not contradicted by the section. List issue sections whose severity or status looks out of date given their own Verification paragraphs (for example an issue whose verification shows it is worse or milder than its rating).
2. Internal contradictions inside issues.md after the edits: the same quantity with two values; a sentence in an older paragraph that a later Verification paragraph contradicts without saying so (the log's convention: text dated 2026-11-18 is left as written and superseded by a Verification paragraph that says so; text of 2026-11-19 is corrected in place). Pay special attention to I-03, I-09, I-16, I-20, I-21, I-37, I-39, I-40, I-41, I-42 and to cross-references between issues (does the issue referred to say what the reference claims?).
3. The table "Corrections after verification": is each row's "as first written" text plausible as an earlier version (you can see earlier versions with: git -C ${REPO} show 16de2c4:issues.md), is the corrected statement what the issue now says, and is any in-place change between 16de2c4 and the working tree missing from the table? Use: git -C ${REPO} diff 16de2c4 -- issues.md
4. Statements that the edits corrected in issues.md but that survive elsewhere in their old form: "no term is enriched" or "no GO term is enriched" stated flatly; "twice as likely" or a twofold gene-body enrichment presented as a result; CpGs that "cannot be DML"; "the third analysis keeps its error rate"; "the enrichment ran under Bioconductor 3.13 (topGO 2.44.0)" stated as established; libraries "about as similar across groups as within"; "collapsed repeat". Search all .md and .R files outside 08_archive and outside 07_documentation/agent_workflows. Note: 01_pipeline/10_go_enrichment_gene_level.R and 01_pipeline/11_semantic_similarity.R contain log or output strings of this kind; report them (they will be changed in a later code pass) with the exact line.
5. Counts quoted about the log itself: number of issues, number bearing on main conclusions, and similar, in README.md, the plan and the handoff file.
6. Formatting: any em dash character in the files listed; broken Markdown tables (rows with a different number of cells from the header) in issues.md; links in issues.md and README.md that point to files that do not exist in the repository.

Return only real problems, each with the exact current text and an exact replacement.`

phase('Check integration')
const reportResults = await parallel(REPORTS.map(r => () =>
  agent(reportPrompt(r), { label: `check:${r.id}`, phase: 'Check integration', schema: SCHEMA })))

phase('Cross-file consistency')
const cross = await agent(crossPrompt, { label: 'check:cross-file', phase: 'Cross-file consistency', schema: SCHEMA })

const all = reportResults.filter(Boolean).concat(cross ? [cross] : [])
log(`${all.length} of ${REPORTS.length + 1} checks returned; ${all.reduce((n, r) => n + r.discrepancies.length, 0)} discrepancies, ${all.reduce((n, r) => n + r.omitted_items.length, 0)} omissions`)
return { checks: all, missing: REPORTS.map(r => r.id).filter(id => !all.some(a => a.report_id === id || a.report_id.startsWith(id.split('-')[0]))) }

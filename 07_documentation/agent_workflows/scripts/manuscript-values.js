export const meta = {
  name: 'manuscript-values',
  description: 'Inventory every quantitative statement of the manuscript draft and map each to the pipeline output that regenerates it',
  phases: [
    { title: 'Extract', detail: 'one reader per manuscript section lists every quantitative statement' },
    { title: 'Map', detail: 'each statement is tied to an R expression over the pipeline outputs and tested' },
    { title: 'Verify', detail: 'an independent agent recomputes each value, checks the mapping, and looks for missed statements' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'
const MANUSCRIPT = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/manuscript/manuscript_gdoc.txt'

const ALL_SECTIONS = {
  methods: { key: 'methods', label: 'Materials and Methods (sections 2 to 2.6)', lines: '17 to 32', idrule: 'M2.0-001 for the opening paragraph of section 2, then M2.1-001, M2.2-001, ... by subsection' },
  results31: { key: 'results31', label: 'Results 3.1', lines: '34 to 37', idrule: 'R3.1-001, R3.1-002, ...' },
  results32_33: { key: 'results32_33', label: 'Results 3.2 and 3.3', lines: '38 to 46', idrule: 'R3.2-001, ... for section 3.2 and R3.3-001, ... for section 3.3' },
  results34_35: { key: 'results34_35', label: 'Results 3.4 and 3.5, and Table 1', lines: '47 to 58 (text) and 160 to 183 (Table 1)', idrule: 'R3.4-001, ..., R3.5-001, ..., and T1-001, ... for the cells of Table 1 (one row per table row and numeric column: start, end, methylation difference; plus one row per table row for the gene ID)' },
  discussion: { key: 'discussion', label: 'Discussion and Conclusions (section 4 and the orphan section 3.5 at the end)', lines: '59 to 115', idrule: 'D4.0-001 for the opening, D4.1-001, D4.2-001, ... by subsection, and D-orphan-001 for lines 114 to 115' },
  legends: { key: 'legends', label: 'Figure legends and supplemental material', lines: '185 to 204', idrule: 'F1-001 for Figure 1, F2-001, ..., FS1-001 for Figure S1' },
}
const wanted = (args && args.sections && args.sections.length) ? args.sections : Object.keys(ALL_SECTIONS)
const SECTIONS = wanted.map(k => ALL_SECTIONS[k]).filter(Boolean)
const EXTRA = (args && args.note) ? args.note : ''

const CONTEXT = `
CONTEXT (read all of it before starting)
Repository: ${REPO} (git branch reproducibility-audit). It holds an R-only pipeline that re-derives the results of a manuscript draft on DNA methylation (whole-genome bisulfite sequencing, WGBS, and a global 5-mC ELISA) in diploid and triploid Pacific oysters after a heat exposure, together with an audit of that draft (issues.md).
The manuscript draft is a text file: ${MANUSCRIPT}. Each line is one paragraph, heading, table row, or figure legend. Use line numbers of that file (sed -n 'Np' or awk) and copy text from it exactly.
Pipeline outputs are CSV tables under 04_output/<step>/tables/. A catalog of every table with its columns is at ${SCRATCH}/output_catalog.txt. The scripts that write them are 01_pipeline/NN_*.R; each has a header block that says what its outputs mean. Read the header of a step before using its tables.
Vocabulary of the pipeline: mode "legacy" = the analysis as it was run for the draft (its numbers should equal the draft's); mode "corrected" = the same with the 99.9th percentile coverage filter that the Methods describe but that was never applied. DML = differentially methylated locus. Sets: "ploidy", "heat", and "shared" (loci in both lists; the draft calls these interaction-DML). "eligible" CpGs = mean methylation between 10% and 90%, the only CpGs that can be DML. Step 09 repeats the draft's GO enrichment exactly; step 10 is the gene-level enrichment the audit recommends (it finds no enriched term).
THE GOAL of this work: a table, 03_input_files/manuscript/manuscript_values_map.csv, with one row per quantitative statement of the draft and an R expression that recomputes that quantity from the pipeline outputs. Step 13 (01_pipeline/13_manuscript_values.R; read its header) evaluates the expressions and reports, for every number in the draft, whether the pipeline reproduces it.
How an expression works: it is one line of R (data.table syntax) that returns exactly one value, with these helpers in scope: tab("<step folder>", "<file>") returns a table of 04_output/<step folder>/tables/ as a data.table; diagnostic("<file>") returns a table of 06_diagnostics/; config is the settings list of run_config.R. Example: tab("05_differential_methylation", "dml_counts_by_mode.csv")[mode == "legacy", ploidy_dml]. A number is compared with the manuscript value within the row's tolerance. For a statement that is an inequality or a ranking (p < 0.05; "the largest share was introns"), write an expression that returns TRUE or FALSE; TRUE means the statement holds. To test expressions, put them one per line in a file in your scratch directory and run: Rscript ${SCRATCH}/eval_expressions.R <file>. Every expression you return must have been run this way.
Machine limits: 4 cores and 15 GB shared with other agents. Do not run pipeline steps. Do not read the large tables whole (united_cpg_*.csv.gz, diffmeth_*.csv.gz) unless a statement cannot be regenerated otherwise.
You are READ-ONLY on the repository: do not create, edit, or delete anything under ${REPO}; no git commands that change state. Other agents are adding steps 11 and 12 while you work; ignore files named 01_pipeline/11_*, 12_* and 04_output/11_*, 12_* unless told otherwise below.
Writing rules for every text field you return: plain declarative English; no em dash character (U+2014); the treatment is "heat" (quote legacy words only inside the verbatim quote field).
${EXTRA}
`

const ROW_PROPS = {
  id: { type: 'string' },
  section: { type: 'string', description: 'For example "Methods 2.5", "Results 3.2", "Table 1", "Figure 2 legend", "Discussion 4.3".' },
  line: { type: 'integer', description: 'Line number in the manuscript text file.' },
  quote: { type: 'string', description: 'The clause that contains the value, copied exactly from the manuscript, at most 250 characters.' },
  item: { type: 'string', description: 'What the value is, in a few plain words, specific enough to recompute it.' },
  manuscript_value: { type: 'string', description: 'The value exactly as written, for example "1,038", "57%", "F = 13.85", "p < 0.01".' },
  manuscript_numeric: { type: ['number', 'null'], description: 'The value as a number (57 for "57%"); null when the statement is not a single number.' },
  kind: { type: 'string', enum: ['count', 'percent', 'statistic', 'p_value', 'parameter', 'design', 'coordinate', 'identifier', 'inequality', 'ranking', 'software_version', 'other'] },
}

const EXTRACT_SCHEMA = {
  type: 'object',
  properties: {
    section_key: { type: 'string' },
    statements: { type: 'array', items: { type: 'object', properties: ROW_PROPS, required: ['id', 'section', 'line', 'quote', 'item', 'manuscript_value', 'manuscript_numeric', 'kind'] } },
    notes: { type: 'string' },
  },
  required: ['section_key', 'statements', 'notes'],
}

const MAPPED_PROPS = Object.assign({}, ROW_PROPS, {
  tolerance: { type: 'number', description: 'Absolute tolerance for a match: 0 for counts and coordinates; half a unit of the last printed digit otherwise.' },
  expression_legacy: { type: 'string', description: 'R expression that recomputes the value as the draft computed it; empty string if it cannot be regenerated.' },
  expression_corrected: { type: 'string', description: 'R expression for the same quantity under the recommended analysis, where that differs; empty string otherwise.' },
  regenerated_legacy: { type: 'string', description: 'What expression_legacy returned when you ran it; empty if no expression.' },
  regenerated_corrected: { type: 'string', description: 'What expression_corrected returned when you ran it; empty if no expression.' },
  status: { type: 'string', enum: ['matches', 'differs', 'not regenerated'] },
  issue: { type: 'string', description: 'IDs of the entries of issues.md that bear on this statement, separated by "; "; empty if none.' },
  note: { type: 'string', description: 'One or two sentences: why it differs, why it cannot be regenerated, or what the reader must know to interpret a match. Empty if nothing to add.' },
})
const MAPPED_REQUIRED = ['id', 'section', 'line', 'quote', 'item', 'manuscript_value', 'manuscript_numeric', 'kind', 'tolerance', 'expression_legacy', 'expression_corrected', 'regenerated_legacy', 'regenerated_corrected', 'status', 'issue', 'note']

const MAP_SCHEMA = {
  type: 'object',
  properties: {
    section_key: { type: 'string' },
    rows: { type: 'array', items: { type: 'object', properties: MAPPED_PROPS, required: MAPPED_REQUIRED } },
    outputs_missing: { type: 'array', items: { type: 'string' }, description: 'Quantities stated in the manuscript that no pipeline output holds, each with what the pipeline would need to write.' },
    notes: { type: 'string' },
  },
  required: ['section_key', 'rows', 'outputs_missing', 'notes'],
}

const VERIFY_SCHEMA = {
  type: 'object',
  properties: {
    section_key: { type: 'string' },
    rows: { type: 'array', items: { type: 'object', properties: MAPPED_PROPS, required: MAPPED_REQUIRED }, description: 'The complete, corrected set of rows for the section: every row you were given (fixed where needed) plus rows for statements that had been missed.' },
    changes: { type: 'array', items: { type: 'object', properties: { id: { type: 'string' }, what: { type: 'string' }, why: { type: 'string' } }, required: ['id', 'what', 'why'] }, description: 'Every row you changed or added, with the reason.' },
    outputs_missing: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
  required: ['section_key', 'rows', 'changes', 'outputs_missing', 'notes'],
}

function extractPrompt(s) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/mv-extract-${s.key}
YOUR TASK: inventory the quantitative statements of one part of the manuscript: ${s.label}, lines ${s.lines} of the text file. Read those lines completely (they are long; print each with sed -n and read all of it).
List EVERY quantitative statement that this study makes about itself: counts of animals, libraries, reads, CpGs, loci, genes, terms; percentages; means; test statistics with their degrees of freedom; p-values and q-values; thresholds, cutoffs and parameters of the analysis; coverage levels; software and database versions; temperatures and durations of the treatment; coordinates and methylation differences in tables; statements of rank or inequality that carry a quantitative claim ("the largest share", "more than", "no significant difference", "p < 0.01"). One row per value: a sentence with three numbers gives three rows, each with the same or an overlapping quote.
Leave out: numbers that describe other studies or general facts, citation years, section and figure numbers, catalogue numbers of reagents, and addresses.
Give each row an ID by this rule, in order of appearance: ${s.idrule}.
The quote must be copied character for character from the text (including symbols such as the multiplication sign in "10×"); keep it short but long enough to identify the value without doubt. When the same number appears twice in the section, make two rows.
Before returning, reread the lines once more and check that no number was skipped. Return section_key "${s.key}".`
}

function mapPrompt(extract, s) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/mv-map-${s.key}
YOUR TASK: for each statement below, taken from ${s.label} of the manuscript (lines ${s.lines}), find the pipeline output that regenerates it and write the expression.
First read: the summary table of issues.md and every issue relevant to this part of the manuscript; the header blocks of the pipeline steps whose outputs you will use; the catalog of tables. Then, for each statement:
1. Read the sentence in context in the manuscript file (line number given) and decide exactly what quantity it states (which samples, which coverage level, which DML set, which denominator).
2. Write expression_legacy: the quantity as the draft's own analysis computed it, from legacy-mode outputs. When the draft's number comes from a known faulty computation that the pipeline replicates (the legacy tables exist for that purpose: for example column legacy_share_pct of step 07 for the location percentages; step 09 for everything about GO terms; mode legacy everywhere), use the replica, so that a match means "the pipeline reproduces how this number came about".
3. Write expression_corrected where the audit recommends a different analysis that gives this quantity a different value: the corrected mode (coverage filter applied) for anything computed from CpGs or DML; shares of loci (pct_of_loci) or the eligible background in place of legacy shares for location statements; step 10 (gene level, universe "eligible", Wallenius, q < config$go_fdr) for enrichment statements. Leave it empty when there is no corrected counterpart, and say why in the note when that is not obvious.
4. Set tolerance: 0 for counts, coordinates and identifiers; half a unit of the last printed digit for everything else (0.5 for "57%", 0.05 for "27.7%", 0.005 for "13.85").
5. Run every expression with the tester and record what it returned (regenerated_legacy, regenerated_corrected). Set status: "matches" if the legacy expression returns the manuscript value within tolerance (or TRUE), "differs" if it returns something else, "not regenerated" if no pipeline output holds the quantity. Do not force a match: if the closest pipeline quantity differs from the manuscript, report "differs" and explain in the note what the pipeline value is and, if you can tell, where the manuscript value came from.
6. A statement is "not regenerated" when it concerns something the pipeline does not compute: wet-laboratory details, sequencing facility numbers, the experiment design as described by the authors, software versions that were not recorded. Say in the note what kind of thing it is and, where the repository holds evidence about it (for example 03_input_files/metadata/wgbs_libraries.csv for numbers of libraries, run_config.R for analysis parameters, 07_documentation/ for upstream processing), put that evidence in the note. Analysis parameters that are settings of run_config.R should be mapped to config (for example config$dml_qvalue), not left unmapped.
7. issue: the IDs in issues.md that bear on the statement.
Rules: keep the IDs, quotes and manuscript values as given unless one is wrong, in which case correct it and say so in the note. Expressions must be self-contained single lines, use only tab(), diagnostic(), config and base R / data.table, and return exactly one value. Prefer the smallest table that holds the value. Do not type a result into an expression.
If a quantity is stated in the manuscript and no output holds it but an output easily could, list it in outputs_missing with the step that should write it.

STATEMENTS (JSON):
${JSON.stringify(extract && extract.statements ? extract.statements : [], null, 1)}

Extractor notes: ${extract && extract.notes ? extract.notes : '(none)'}

Return section_key "${s.key}" and one row per statement, in the same order.`
}

function verifyPrompt(mapped, s) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/mv-verify-${s.key}
YOUR TASK: independent verification of a set of rows for manuscript_values_map.csv covering ${s.label} of the manuscript (lines ${s.lines}). Another agent extracted the statements and wrote the expressions. Assume there are mistakes. Check all of the following and return the complete corrected set.
A. Completeness. Read lines ${s.lines} of the manuscript yourself, in full, and list every quantitative statement (same criteria: every number, percentage, statistic, p-value, threshold, parameter, version, coordinate, and every ranking or inequality that carries a quantitative claim, about this study). Add a fully worked row for each one that is missing, with a new ID that continues the numbering.
B. Quotes. For each row, confirm that the quote occurs verbatim on the stated line (grep -F or an exact string test in R or awk; mind special characters). Fix quotes and line numbers that are off.
C. Semantics. For each row with an expression, decide independently what quantity the sentence states, then recompute it yourself by reading the output files directly with your own code (not by pasting the expression). Does the expression compute that quantity, for the right samples, set, coverage level, mode and denominator? A row that "matches" through the wrong quantity is the worst error: look for it. Check that expression_corrected is the same quantity under the recommended analysis, not a different quantity.
D. Run every expression with the tester and compare with regenerated_legacy and regenerated_corrected; fix the recorded values and the status. Check the tolerance rule (0 for counts, coordinates, identifiers; half a unit of the last printed digit otherwise).
E. Rows marked "not regenerated": could a pipeline output regenerate them after all? Look in the catalog and in 03_input_files (for example the library sheet, the Bismark reports table of step 02, run_config.R). If so, write the expression.
F. Issue references: each ID must exist in issues.md and bear on the statement; add any that are missing.
G. Notes: accurate, specific, and free of the em dash character.
Return every row of the section in its final form (corrected rows, untouched rows, and added rows), and list each change with its reason.

ROWS TO VERIFY (JSON):
${JSON.stringify(mapped && mapped.rows ? mapped.rows : [], null, 1)}

Mapper notes: ${mapped && mapped.notes ? mapped.notes : '(none)'}
Outputs the mapper found missing: ${JSON.stringify(mapped && mapped.outputs_missing ? mapped.outputs_missing : [])}

Return section_key "${s.key}".`
}

const results = await pipeline(
  SECTIONS,
  s => agent(extractPrompt(s), { label: `extract:${s.key}`, phase: 'Extract', schema: EXTRACT_SCHEMA }),
  (extract, s) => agent(mapPrompt(extract, s), { label: `map:${s.key}`, phase: 'Map', schema: MAP_SCHEMA })
    .then(mapped => ({ extract, mapped })),
  (em, s) => agent(verifyPrompt(em.mapped, s), { label: `verify:${s.key}`, phase: 'Verify', schema: VERIFY_SCHEMA })
    .then(verified => ({ n_extracted: em.extract && em.extract.statements ? em.extract.statements.length : 0, n_mapped: em.mapped && em.mapped.rows ? em.mapped.rows.length : 0, verified }))
)

const out = SECTIONS.map((s, i) => ({ section_key: s.key, result: results[i] }))
log(out.map(o => `${o.section_key}: ${o.result && o.result.verified && o.result.verified.rows ? o.result.verified.rows.length : 0} rows`).join('; '))
return out

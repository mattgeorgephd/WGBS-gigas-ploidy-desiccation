export const meta = {
  name: 'build-steps-11-12',
  description: 'Build pipeline steps 11 (semantic similarity of GO terms) and 12 (qPCR), each built, independently reviewed, then fixed',
  phases: [
    { title: 'Build', detail: 'one builder per step writes the script and runs it' },
    { title: 'Review', detail: 'an independent reviewer runs the step and looks for defects' },
    { title: 'Fix', detail: 'confirmed findings are fixed and the step is rerun' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'
const MANUSCRIPT = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/manuscript/manuscript_gdoc.txt'

const CONTEXT = `
CONTEXT (read all of it before starting)
Repository: ${REPO} (git branch reproducibility-audit). It holds an R-only pipeline that re-derives the results of a manuscript draft on DNA methylation (whole-genome bisulfite sequencing) in diploid and triploid Pacific oysters after a heat exposure, together with an audit of that draft. A local text copy of the manuscript is at ${MANUSCRIPT}.
Layout: 01_pipeline/ (numbered R scripts; steps 01 to 10 exist), 02_R_functions/ (helpers sourced by every step), 03_input_files/, 04_output/<step>/{tables,figures}, 06_diagnostics/, 07_documentation/, 08_archive/ (legacy code and outputs of the original analysis; never modify), issues.md (the issue log), run_config.R (all settings), run_pipeline.R.
Read first: 01_pipeline/README.md (section "How a step is written"), run_config.R, the files in 02_R_functions/ (so that you use the existing helpers: load_packages, input_path, output_path, archive_path, start_step, finish_step, log_msg, write_table, theme_wgbs, save_figure), and two existing steps as models of style: 01_pipeline/09_go_enrichment_legacy.R and 01_pipeline/10_go_enrichment_gene_level.R.
Vocabulary: "legacy" = the analysis as run for the draft; "DML" = differentially methylated locus; "shared" = loci in both the ploidy and the heat DML list (the draft calls them interaction-DML); the treatment is called "heat" (legacy files say "desiccation"; never use that word in new text except when quoting a legacy file name or column).
How to run: Rscript is installed (R 4.3.3; packages include data.table, ggplot2, patchwork, ggrepel, topGO, GO.db 3.18.0, rrvgo, GOSemSim, org.Hs.eg.db, simplifyEnrichment, survival, emmeans, multcomp, readxl). Export WGBS_CACHE_DIR=/home/user/gannet_cache before running pipeline code. Run a step as: cd ${REPO} && Rscript 01_pipeline/<script>.R
MACHINE LIMITS: 4 cores, 15 GB of memory, shared with other agents working at the same time. Do not run 01_pipeline/04_cpg_coverage.R or 05_differential_methylation.R. Keep R processes under about 3 GB.
RULES FOR EVERYTHING YOU WRITE INTO THE REPOSITORY
- R only. No Python, shell, or other languages in the pipeline.
- A step script begins with a comment header in the exact style of the existing steps: Purpose, Inputs, Outputs, Manuscript items, Notes. Then the same four preamble lines the other steps use (library(here); source run_config.R; source every file in 02_R_functions; load_packages(...)), then step <- "<NN_name>"; t0 <- start_step(step, "..."), and it ends with finish_step(step, t0).
- All paths through the helpers (here::here underneath). No setwd, no absolute paths, no rstudioapi.
- Settings come from the config list in run_config.R (the blocks for steps 11 and 12 are already there). Do not edit run_config.R; if you need a setting that is not there, define it in a clearly marked block at the top of your script and report it.
- No number is typed into the script as a result. Values quoted from the manuscript or from legacy files, used for comparison, are allowed and must be labelled as such.
- Outputs go only to 04_output/<step>/tables (CSV via write_table; gzip only if larger than a few MB) and 04_output/<step>/figures (via save_figure, which writes PNG and PDF). Output must be deterministic: set a seed before anything random and avoid algorithms whose result changes between runs.
- The script must run from a fresh R process and must fail loudly (stop with a clear message) rather than produce a wrong or partial result. It must also cope with an empty input where that can legitimately happen.
- Writing style in comments, messages, and tables: plain declarative English, short sentences, no em dash character (U+2014) anywhere, no marketing words. State what a number is. Cite sources for methods (authors, year, journal, volume, page or DOI) where a method choice is made.
- data.table hazard known in this code base: inside DT[...], a function argument or loop variable with the same name as a column refers to the column. Name variables so that this cannot happen, and assert row counts after joins.
FILES YOU MAY CREATE OR CHANGE: only the ones listed in your task. Do not create files in 02_R_functions/ (keep step-specific functions inside your step script). Do not edit README files, issues.md, or run_config.R: return the text you propose for them instead. Do not run git commands that change state (no add, commit, checkout, stash, push). Scratch work goes under your scratch directory (given below).
`

const STEPS = [
  {
    id: '11_semantic_similarity',
    script: '01_pipeline/11_semantic_similarity.R',
    files: '01_pipeline/11_semantic_similarity.R and everything under 04_output/11_semantic_similarity/ (written by running the script)',
    task: `Build pipeline step 11: semantic similarity of GO terms (manuscript Figures 6, 7 and 8 and the text of Results 3.4 that describes them).

Background you need.
- The draft's Figures 6 to 8 are "semantic similarity plots of significantly enriched GO terms", biological process only, for ploidy-DML (Figure 6), heat-DML (Figure 7) and interaction-DML (Figure 8). The Methods say they were made with the REVIGO web server using the SimRel measure. A web server cannot be part of an R-only pipeline, so this step redoes the computation in R with the Bioconductor package rrvgo (Sayols 2023, MicroPublication Biology, doi:10.17912/micropub.biology.000811), which implements the same approach: pairwise semantic similarity between terms, grouping of terms more similar than a threshold, and a two-dimensional layout of the similarity matrix.
- The terms that went into those figures are the BP terms with Fisher p < 0.01 in the legacy enrichment. Step 09 regenerates them exactly: 04_output/09_go_enrichment_legacy/tables/legacy_fisher_<set>_<ontology>.csv for set in ploidy, heat, shared and ontology in BP, CC, MF; use the rows with below_p_cutoff == TRUE (383, 856 and 341 BP terms). Those tables were computed under the 2021 ontology (GO.db 3.13.0); this step uses the GO.db of the R installation (3.18.0), so some term IDs no longer exist: drop them, count them, and report the counts in a table.
- Read issues I-03, I-21, I-22, I-39, I-40 and I-41 in issues.md first. The audit found that the legacy enrichment is unsound and that, at gene level against genes that could hold a DML, no GO term is enriched (step 10; 04_output/10_go_enrichment_gene_level/tables/go_enriched_terms.csv is empty apart from its header). So this step has two jobs: (1) regenerate the draft's three figures in R from the legacy term lists, so that they can be traced; (2) run the same procedure on the terms the gene-level analysis calls enriched, which at present is nothing: the script must handle that case by writing a one-row table that says so and drawing no figure for it, without failing.
- What survives of the REVIGO runs: 08_archive/code/Revigo.R, 08_archive/code/OLD/Revigo.R and 08_archive/figures/revigo_figure/Revigo-PLOIDY-BP.R are plotting scripts exported by the server with the data embedded (term_ID, description, frequency, plot_X, plot_Y, log_size, value, uniqueness, dispensability), and 08_archive/plots/ and 08_archive/figures/ hold images. Work out which term list each export corresponds to (compare its term IDs with the p < 0.01 lists of step 09 for the nine set and ontology combinations) and report it.

What the step must produce.
1. For each DML set and ontology of the legacy lists: the similarity matrix among the listed terms (rrvgo::calculateSimMatrix with the organism database, measure, and threshold given in config$semantic_similarity), the reduction into clusters (rrvgo::reduceSimMatrix with scores = -log10 p), and a table semantic_clusters_legacy_<set>_<ontology>.csv with one row per term: GO ID, term, p, score, cluster, cluster representative (parent), cluster size, and the two layout coordinates. Use a deterministic layout (classical multidimensional scaling of 1 - similarity, stats::cmdscale), not UMAP. Record how many terms were dropped because the ontology or the organism database does not hold them.
2. Figures 6, 7 and 8 regenerated: BP for ploidy, heat, shared; your own ggplot built from the coordinates (points sized by cluster or term size, coloured by cluster, the config$semantic_similarity$n_labels highest-scoring cluster representatives labelled with ggrepel using a fixed seed). File names figure_6_legacy, figure_7_legacy, figure_8_legacy. Make the label rule explicit in the header: labels are chosen by rule, not by hand (issue I-22 says the draft's labels were hand-picked).
3. A redraw of each surviving REVIGO export from its embedded data, for comparison with the R version (figure_revigo_export_<name>), and a table revigo_exports.csv saying which export matches which term list, how many of its terms are in that list, and how many of the list's terms are in the export.
4. A table named_terms_results_3_4.csv that takes every GO term named in Results 3.4 of the manuscript (the three paragraphs that begin "Semantic clustering of enriched GO terms for ploidy-DML", "Enriched GO terms for heat-DML", and "Semantic clustering of interaction-DML") and, for the DML set that paragraph is about, gives: the phrase as written in the manuscript, the GO ID and official term name you matched it to (match on term names in GO.db, exact where possible; record the match type as exact, close, or not found; do not guess silently), whether the term is in the legacy p < 0.01 list of that set, its rank and p-value there, the number of transcripts annotated and of interest, the number of distinct genes behind the transcripts of interest (column genes_significant of the step 09 table), and its gene-level result from step 10 (04_output/10_go_enrichment_gene_level/tables/go_gene_level_<set>_legacy.csv, universe "eligible": genes annotated, genes with a DML, p_wallenius, q_wallenius, or a note when it was not tested).
5. The gene-level counterpart described above (empty case handled).
6. A short table semantic_summary.csv: per set and ontology, terms in, terms dropped and why, clusters out.
Checks to make before you finish: the script runs twice in a row with identical CSV outputs; every figure opened and inspected (use the Read tool on the PNG) so that labels are legible and nothing overlaps badly; the counts in semantic_summary.csv add up.`,
  },
  {
    id: '12_qpcr',
    script: '01_pipeline/12_qpcr.R',
    files: '01_pipeline/12_qpcr.R and everything under 04_output/12_qpcr/ (written by running the script)',
    task: `Build pipeline step 12: the qPCR analysis.

Background you need.
- Read 03_input_files/qpcr/README.md, section 3 and Appendix C of 07_documentation/audit/elisa_qpcr_audit.md (an earlier audit of these data, with numbers you must be able to reproduce), and issue I-29 in issues.md. Then find every statement in the manuscript that rests on qPCR (search the text copy for qPCR, HSP90, COX1, mitochond, expression): they are in the Discussion only; the Methods and Results do not describe qPCR at all.
- Inputs: 03_input_files/qpcr/raw/qpcr_data_<gene>.csv for nine genes (actin, dnmt1, hat, hathap2, hif1a, hsc70, hsp90, mbd2, mecp2): columns Sample and "Cq Value", 80 rows each = 40 animals x 2 technical replicates; non-detects are the text N/A; sample IDs like D_N_N01 = ploidy (D diploid, T triploid), treatment (N control, D heat), then a letter for the acute heat shock (N none, H 1 h at 45 C) joined to the animal number. And 03_input_files/qpcr/raw/qpcr_data_consolidated.csv: one row per animal for the 32 animals without acute heat shock, IDs like D_C_N_01 (second token C control, D heat), one column per gene, built by hand in a spreadsheet; it also holds SOD, ATPsynthetase and COX1, for which no replicate-level file exists.
- Legacy analysis (08_archive/gene_expression/scripts/dct_analysis.R, dct_analysis_NA-to-45.R, 08_archive/ELISA/dct_analysis_NA-to-45_COX1convert.Rmd and its knitted .md with the printed ANOVA tables): actin is the only reference; expression = 2^-(Cq target - Cq actin) per animal, natural log, then aov(log expression ~ Ploidy + Desiccation + Ploidy:Desiccation) and TukeyHSD per gene on the 32 animals; non-detects either dropped or replaced by Cq 45; no correction across genes; COX1 analysed as raw Cq or 2^-Cq without a reference.

What the step must produce.
1. A per-animal table built in R from the replicate files (qpcr_per_animal.csv): animal, ploidy, treatment, acute heat shock flag, and per gene the mean Cq of the detected replicates, the number of detected replicates, and the replicate difference. Parse the sample IDs with explicit checks (stop if an ID does not parse or a group has an unexpected size). Then a check table comparing it with the hand-built consolidated file for the 32 animals x 9 genes (qpcr_consolidated_check.csv: largest absolute difference per gene, number of cells that differ, and which). SOD, ATPsynthetase and COX1 can only come from the consolidated file; carry them with a flag saying no replicate data exist.
2. The legacy analysis reproduced for every gene, in both legacy variants (non-detects dropped; non-detects set to config$qpcr$legacy_nondetect), with a table of F and p for ploidy, treatment and interaction (qpcr_anova_legacy.csv), compared with the numbers printed in the knitted legacy .md and in the audit note (for example HSP90 with non-detects at 45: ploidy F(1,28) = 4.90, p = 0.035; treatment F = 6.39, p = 0.017). Report agreement gene by gene in a column.
3. An analysis you can defend, next to the legacy one (qpcr_models.csv), on delta Cq (target minus reference), per gene, with the model delta Cq ~ ploidy * treatment on the 32 animals: (a) complete cases; (b) a censored regression that treats a non-detect as "Cq above config$qpcr$max_cycles" (survival::survreg with a Gaussian distribution on delta Cq, right-censored; say in the header why substitution by 45 is not a sound treatment of non-detects and cite McCall et al. 2014, Bioinformatics 30:2310, doi:10.1093/bioinformatics/btu239, if you use that argument); Benjamini-Hochberg across genes for each term; effect sizes as fold changes with 95% confidence intervals; number of non-detects per gene and group. Decide and state what to do with genes that are mostly non-detects (MeCP2 is a non-detect in 75% of wells): reporting that a gene cannot be analysed is an acceptable outcome.
4. Reference gene diagnostics (qpcr_reference_gene.csv): actin Cq by group, a test of whether actin Cq differs by ploidy or treatment (if the only reference differs between groups, normalization moves the result), and a list of animals whose actin Cq marks a failed sample (the audit note names D_C_N_05 with actin Cq 35.26). Show how each conclusion changes when failed samples are excluded.
5. COX1: what can and cannot be said. It was measured without a nuclear reference on the same template, so "mitochondrial abundance" is not estimable; report the raw Cq by group and the tests the audit note lists, labelled for what they are.
6. A table qpcr_manuscript_statements.csv: each qPCR-based statement of the manuscript quoted verbatim with its location, the numbers from this step that bear on it, and a verdict in plain words (supported, partly supported with the condition, not supported, cannot be evaluated), with one sentence of reasoning.
7. Figures: figure_qpcr_expression (all genes with replicate-level data: relative expression per animal on a log2 axis by ploidy and treatment, every animal shown, non-detects marked, group means), and figure_qpcr_cox1_actin (raw Cq of COX1 and of actin by group). Colours from config$colours. Inspect each PNG with the Read tool.
8. The 8 animals with the acute heat shock (numbers 09, 10, 19, 20 of each ploidy) are in the replicate files but were left out of the legacy analysis. Keep them in the per-animal table with their flag, leave them out of the models by default, and say so in the header.
Checks to make before you finish: every number in section 3.3 of the audit note that your tables should contain is reproduced (list them with the regenerated value); the script runs twice in a row with identical CSV outputs.`,
  },
]

const BUILD_SCHEMA = {
  type: 'object',
  properties: {
    step: { type: 'string' },
    files_written: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string', description: 'What the step does and the main results, with numbers.' },
    key_results: { type: 'array', items: { type: 'string' } },
    validation: { type: 'array', items: { type: 'string' }, description: 'Each check you ran and its outcome.' },
    caveats: { type: 'array', items: { type: 'string' } },
    settings_needed: { type: 'string', description: 'Any setting you defined in the script that should move to run_config.R; empty if none.' },
    readme_row: { type: 'string', description: 'Proposed one-row description of the step for the table in 01_pipeline/README.md.' },
    issues_text: { type: 'string', description: 'Proposed text for issues.md: updates to existing issues or new issues that this step supports, with numbers and file references.' },
    run_time_minutes: { type: 'number' },
  },
  required: ['step', 'files_written', 'summary', 'key_results', 'validation', 'caveats', 'settings_needed', 'readme_row', 'issues_text'],
}

const REVIEW_SCHEMA = {
  type: 'object',
  properties: {
    step: { type: 'string' },
    runs_clean: { type: 'boolean' },
    deterministic: { type: 'boolean' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          summary: { type: 'string' },
          failure_scenario: { type: 'string' },
          evidence: { type: 'string' },
          suggested_fix: { type: 'string' },
        },
        required: ['file', 'severity', 'summary', 'failure_scenario', 'evidence', 'suggested_fix'],
      },
    },
    independent_checks: { type: 'array', items: { type: 'string' }, description: 'Key numbers you recomputed by other means, with both values.' },
    notes: { type: 'string' },
  },
  required: ['step', 'runs_clean', 'deterministic', 'findings', 'independent_checks', 'notes'],
}

function buildPrompt(s) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/build-${s.id}
FILES YOU MAY CREATE OR CHANGE: ${s.files}

YOUR TASK
${s.task}

When you are done, the script must exist, run cleanly from a fresh R process, and its outputs must be in place. Return step "${s.id}" with an honest account: what works, what you verified and how, and what you could not do.`
}

function reviewPrompt(build, s) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/review-${s.id}
You may run the step (cd ${REPO} && Rscript ${s.script}); running it rewrites its own outputs under 04_output/${s.id}/, which is expected. Do not edit any file.

YOUR TASK: independent review of a pipeline step that another agent has just written. You did not write it and you should assume it has defects. The specification the builder was given:
${s.task}

The builder's own report (treat as claims to check, not as facts):
${JSON.stringify(build, null, 1)}

Do all of the following. (1) Run the script from a fresh R process; record whether it exits cleanly and how long it takes. Run it a second time and compare the CSV outputs byte for byte (md5sum) to test determinism. (2) Read the whole script for correctness defects: wrong joins, index or factor-level misalignment, data.table scoping errors, mislabelled columns or axes, statistics that do not do what the comment says, silent dropping of rows, a legitimate empty input that would crash it. (3) Recompute at least five key numbers of its output tables by independent means (your own short R code in your scratch directory, not a copy of the builder's logic) and report both values. (4) Open every PNG it wrote with the Read tool and check that each figure shows what its file name and the header claim, with legible labels, correct legends, and no overplotted or cut-off text. (5) Check the rules: header block complete and accurate (every input and output named exists), helpers used, no absolute paths, no setwd, no typed results, no em dash character (search the script and its CSV outputs for U+2014), the word "desiccation" only where a legacy name is quoted, sources cited for method choices. (6) Check that every item of the specification was delivered, and list any that was not.
Report only real problems, each with a concrete failure scenario and evidence. Return step "${s.id}".`
}

function fixPrompt(review, s, build) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/fix-${s.id}
FILES YOU MAY CREATE OR CHANGE: ${s.files}

YOUR TASK: a pipeline step was written by one agent and reviewed by another. Bring it to a finished state. The specification:
${s.task}

The builder's report:
${JSON.stringify(build, null, 1)}

The reviewer's report:
${JSON.stringify(review, null, 1)}

For each reviewer finding: reproduce or confirm it yourself first; fix it if it is real; if you conclude it is not a defect, leave the code and say why. Also deliver any item of the specification the reviewer found missing. Then run the script from a fresh R process twice, confirm identical CSV outputs, and open each PNG with the Read tool to confirm the figures are right. Do not rewrite parts that work.
Return step "${s.id}" with the final state: in summary and key_results give the results as they now stand (rerun numbers, not the builder's), in validation list each reviewer finding with what you did about it and each check you ran, and give final proposed readme_row and issues_text.`
}

const results = await pipeline(
  STEPS,
  s => agent(buildPrompt(s), { label: `build:${s.id}`, phase: 'Build', schema: BUILD_SCHEMA }),
  (build, s) => agent(reviewPrompt(build, s), { label: `review:${s.id}`, phase: 'Review', schema: REVIEW_SCHEMA })
    .then(review => ({ build, review })),
  (br, s) => agent(fixPrompt(br.review, s, br.build), { label: `fix:${s.id}`, phase: 'Fix', schema: BUILD_SCHEMA })
    .then(final => ({ build: br.build, review: br.review, final }))
)
return STEPS.map((s, i) => ({ step: s.id, result: results[i] }))

export const meta = {
  name: 'check-sensitivity-text',
  description: 'Independent, adversarial check of the sensitivity results written into issues.md (I-06, I-49 and five updates), of the new diagnostic script, and of the citations used',
  phases: [
    { title: 'Check', detail: 'five agents: I-06 numbers and reading; I-49 numbers and alternative explanations; code review of library_methylation_level.R; cross-file consistency; literature' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'
const MANUSCRIPT = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/manuscript/manuscript_gdoc.txt'

const SCHEMA = {
  type: 'object',
  properties: {
    check_id: { type: 'string' },
    statements_checked: { type: 'integer', description: 'number of separate numbers or statements you compared with the data' },
    confirmed: { type: 'integer' },
    discrepancies: {
      type: 'array',
      description: 'statements that are wrong as written: a number that the tables or your recomputation do not give, a wrong file or column reference, an internal contradiction, a stale statement',
      items: {
        type: 'object',
        properties: {
          where: { type: 'string', description: 'file, issue ID or section, and enough of the sentence to find it' },
          quote: { type: 'string', description: 'the exact text as it stands now' },
          problem: { type: 'string' },
          evidence: { type: 'string', description: 'what you computed or read, with file names and values' },
          proposed_text: { type: 'string', description: 'exact replacement text, without em dashes' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['where', 'quote', 'problem', 'evidence', 'proposed_text', 'severity'],
      },
    },
    challenges: {
      type: 'array',
      description: 'for each interpretive claim you tried to refute: the objection or alternative explanation, the test you ran, and the verdict',
      items: {
        type: 'object',
        properties: {
          claim: { type: 'string', description: 'the claim, quoted or closely paraphrased, with its location' },
          objection: { type: 'string', description: 'the alternative explanation or objection you tested' },
          test_done: { type: 'string', description: 'what you computed, on which files, with which method' },
          result: { type: 'string', description: 'the numbers you obtained' },
          verdict: { type: 'string', enum: ['claim_stands', 'claim_needs_qualification', 'claim_refuted', 'could_not_test'] },
          proposed_text: { type: 'string', description: 'replacement or added text if the verdict is not claim_stands; empty otherwise; no em dashes' },
        },
        required: ['claim', 'objection', 'test_done', 'result', 'verdict', 'proposed_text'],
      },
    },
    additions: {
      type: 'array',
      description: 'things a reader needs that the text does not say',
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
    code_findings: {
      type: 'array',
      description: 'defects of code (only for the code review)',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          problem: { type: 'string' },
          failure_scenario: { type: 'string', description: 'concrete input or state that gives a wrong output, a crash, or a misleading table; or "none at current settings" with the setting that would trigger it' },
          changes_a_committed_number: { type: 'boolean' },
          proposed_fix: { type: 'string' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['file', 'line', 'problem', 'failure_scenario', 'changes_a_committed_number', 'proposed_fix', 'severity'],
      },
    },
    citations: {
      type: 'array',
      description: 'only for the literature check',
      items: {
        type: 'object',
        properties: {
          citation_as_written: { type: 'string' },
          where_used: { type: 'string' },
          metadata_verified: { type: 'boolean' },
          source_checked: { type: 'string', description: 'URL or API endpoint that answered, or what refused' },
          correct_metadata: { type: 'string', description: 'authors, year, title, journal, volume, pages or article number, DOI, as the source gives them' },
          content_supports_use: { type: 'string', enum: ['yes', 'partly', 'no', 'could_not_check'] },
          excerpt: { type: 'string', description: 'verbatim sentence(s) from the abstract or text that bear on the use, at most 600 characters; empty if not accessible' },
          note: { type: 'string' },
        },
        required: ['citation_as_written', 'where_used', 'metadata_verified', 'source_checked', 'correct_metadata', 'content_supports_use', 'excerpt', 'note'],
      },
    },
    notes: { type: 'string', description: 'anything else, including what you could not check and why' },
  },
  required: ['check_id', 'statements_checked', 'confirmed', 'discrepancies', 'challenges', 'additions', 'notes'],
}

const COMMON = `
You are an independent checker whose job is to REFUTE. The repository ${REPO} holds an audit of a manuscript draft on DNA methylation (whole-genome bisulfite sequencing, WGBS) in diploid and triploid Pacific oysters, heat-treated or not. The audit's issue log is ${REPO}/issues.md. An analyst has just written the results of a set of sensitivity analyses into it. Your default assumption is that the analyst made mistakes: numbers copied wrongly, readings stronger than the evidence, alternative explanations not considered. Find them. If after honest effort a statement holds, say so.

Background you need
- 22 WGBS libraries. "Heat cohort": 10 libraries of Zymo project zr3534 (zr3534_1 to 5 diploid, 6 to 10 triploid; all heat-treated; libraries 4, 5, 9, 10 also had an acute heat shock). "Control cohort": 12 libraries of project zr3644 (zr3644_1 to 6 diploid, 13 to 18 triploid; untreated animals of a different experiment). Heat treatment is therefore completely confounded with experiment and sequencing batch (issue I-01). Sample sheet: ${REPO}/03_input_files/metadata/wgbs_libraries.csv.
- The draft's analysis: methylKit calculateDiffMeth without overdispersion correction, ploidy tested with heat as covariate and heat tested with ploidy as covariate, a DML = q < 0.01 and difference above 20 points. "legacy" mode = the CpG set as the draft ran it (480,579 CpGs with at least 10 reads in all 22 libraries); "corrected" mode = with the 99.9th percentile coverage filter the Methods describe (454,555 CpGs). "Eligible" CpGs = mean methylation of the 22 libraries between 10% and 90%.
- Step outputs are under ${REPO}/04_output/<step>/tables/. The united count tables are 04_output/05_differential_methylation/tables/united_cpg_<mode>.csv.gz (columns chr, start, end, strand, then coverage<j>, numCs<j>, numTs<j> for library j = 1..22 in the row order of the sample sheet; normalized counts as methylKit tests them). Per-CpG results of the draft's two contrasts: diffmeth_{ploidy,heat}_<mode>.csv.gz in the same folder.
- Sensitivity analyses: script ${REPO}/06_diagnostics/dml_sensitivity.R, outputs in ${REPO}/06_diagnostics/dml_sensitivity/ (analyses.csv describes every fit; dml_counts.csv; permutation_reference.csv; dispersion_summary.csv; one dml_<analysis>_<mode>.csv per fit; cached fits in cache/fit_<analysis>_<mode>.rds, each a list whose element $fit is a data.table with pvalue, qvalue, meth.diff, rows in the order of the united table of that mode). Documentation: ${REPO}/06_diagnostics/README.md.
- Library-level diagnostic: script ${REPO}/06_diagnostics/library_methylation_level.R, outputs in ${REPO}/06_diagnostics/library_methylation_level/ (11 CSV files, described in the script header and in 06_diagnostics/README.md).
- R 4.3.3 with methylKit 1.28.0 and data.table is installed. Run R from the repository root.

Rules
- Read-only: do not modify, create or delete any file under ${REPO}, and run no git command that changes state. Write your own scripts and scratch files only under ${SCRATCH}/<a new subdirectory named after your check>. Keep your scripts there: they are a record.
- The machine has 4 cores and 15 GB and other agents are working on it. Use one thread (data.table::setDTthreads(1); no parallel workers). Before a computation that needs more than 1 GB, run "free -g" and wait if less than 4 GB is available. Never run two R processes at once. Keep each computation under about 10 minutes and 4 GB. Reading a united table takes about 20 seconds and 1 GB.
- Other agents may be writing unrelated files under ${REPO}/06_diagnostics (legacy_reader_*, coverage_filter_*, check_methylkit_*, feature_odds_*, locus_nc047563*). Ignore those.
- The manuscript draft is unpublished and is not in the repository. A text copy for reference is at ${MANUSCRIPT}; do not copy it anywhere and quote at most a clause from it.
- In anything you send to an outside service (web requests), include only what is needed to look up a published paper (a DOI, a title, author names). No file contents, no project or user details.
- Proposed text must be plain, exact, and ready to paste: full sentences, no em dashes (use commas or semicolons), treated oysters are "heat" or "heat-treated", numbers with the same rounding style as the surrounding text, and the file that holds each number named.
- Report every discrepancy, however small, with severity: high = a reader would draw a wrong conclusion about the manuscript; medium = a number or reference is wrong or a needed qualifier is missing; low = wording, rounding, style.
`

const CHECKS = [
  {
    id: 'S1-I06-numbers-and-reading',
    prompt: `${COMMON}
Your check: the section of ${REPO}/issues.md headed "### I-06 (S1) No overdispersion correction in the differential methylation model", from the paragraph "**Rating.**" to the end of that issue (the paragraphs before it are older text; read them for context).

Tasks
(a) Numbers. Recompute every number and every table cell of that text from the files in 06_diagnostics/dml_sensitivity/ and report each mismatch. Item 6 of the text takes some counts "from a scratch computation on the cached fits" (CpGs with q < 0.01 by whether the coverage filter keeps them, how many are mitochondrial (chromosome NC_001276.1), the q-value ranges of the ten CpGs that leave the list, the lists under Benjamini-Hochberg): recompute those yourself from cache/fit_A1_heat_mn_f_legacy.rds and cache/fit_A1_heat_mn_f_corrected.rds. Check the per-library values quoted for NC_047565.1:35,030,701 against the united table.
(b) Quotations. Verify the two quotations from the methylKit documentation word for word against the installed package (help page of calculateDiffMeth; vignette under system.file("doc", package = "methylKit")).
(c) The reading. Try to refute each of the following, with a computation where one is possible:
  1. Item 1, "The test statistic of the uncorrected model is thus inflated about threefold at a typical CpG." Is the dispersion estimator of the script the quantity by which methylKit's correction divides the statistic (read methylKit's source: methylKit:::logReg or the function calculateDiffMeth dispatches to)? Is the comparison value 1.18 right?
  2. Item 2, "With the correction, the ploidy contrast of the draft and the ploidy contrast among the controls have nearly the p-value distribution of groups that do not differ." This assumes the corrected F test holds its level on these data. Test that: take a random sample of about 3,000 eligible CpGs from the legacy united table, reassign the ploidy labels within cohort about 30 times, and for each relabelling compute the corrected p-value exactly as methylKit does for overdispersion = "MN", test = "F" with heat as covariate (reproduce its computation in base R after reading its source, and confirm on the observed labels that your p-values equal those of cache/fit_A1_ploidy_mn_f_legacy.rds at the sampled CpGs). Report the shares of p below 0.05, 0.01 and 0.001 under relabelling. If the test is anticonservative or conservative under relabelling, say what that does to the reading of the table in item 2 and of the row for the observed ploidy labels.
  3. Item 5, the relabelling reference: is the comparison fair (same statistic, same libraries, cohort kept)? Is "more than one such grouping in ten yields more than the real one" what the table says? Does anything in permutation_reference.csv contradict the text?
  4. Item 7: recompute the F and chi-squared statistics needed for a single call from the residual degrees of freedom in analyses.csv. Is the statement about what a short corrected list means fair?
  5. "What follows for the manuscript" and "Options". Argue the other side: what would a defender of the draft say (for example that a 20-point threshold with q < 0.01 is common practice in the bivalve methylation literature, that the corrected F test has little power with 11 animals per group, that pooling across cohorts adjusted by a covariate is legitimate)? Does the text answer those points fairly, or does it claim more than the evidence gives (look hard at "cannot be reported as differences between diploid and triploid oysters" and "This removes the basis of Results 3.2 to 3.5")? Propose exact rewording where it overreaches or where it is too weak.
(d) Anything in the sensitivity tables that the text should report and does not (for example a row of dml_counts.csv that contradicts the summary, the share_of_list_in_primary values, the number of interaction DML by sign).

Return check_id "S1-I06-numbers-and-reading".`,
  },
  {
    id: 'S2-I49-evidence-and-alternatives',
    prompt: `${COMMON}
Your check: the issue "### I-49 (S1) The group differences in the WGBS data are carried by five libraries with a higher overall methylation level and by five mutually similar triploid libraries" in ${REPO}/issues.md (near the end of the file, before "## Corrections after verification"). It was written from the tables in ${REPO}/06_diagnostics/library_methylation_level/. This issue makes the strongest interpretive claims of the audit, so it gets the hardest check.

Tasks
(a) Numbers. Compare every number in I-49 with the tables. Then recompute the key ones with your own code, not by calling the script's functions: the five libraries and the gap; the counts by sign for the ploidy contrast inside each cohort (legacy mode; pooled reads per group; a CpG counts when the absolute difference exceeds 20 points); the correlation of the two within-cohort differences across CpGs; the 21 counts of the two-against-five enumeration in the heat cohort without its three high libraries; the mean pairwise correlation of the five heat-cohort triploids and its rank among all 252 subsets of five (correlation matrix: 04_output/04_cpg_coverage/tables/library_correlation_10x_legacy.csv).
(b) Alternative explanations. For each, state the test you ran and the numbers:
  1. Is the "level" an artifact of sequencing depth, duplication, mapping efficiency or coverage at the tested CpGs? Use 04_output/02_alignment_qc/tables/bismark_by_library.csv, 04_output/04_cpg_coverage/tables/library_coverage.csv and the per-library coverage columns of the united table. Does the level correlate with any of them, within cohort?
  2. Is the mutual similarity of the heat-cohort triploids (item 6) explained by sampling noise, that is, by those libraries having more reads at the tested CpGs? Pairwise Pearson correlation of percent methylation falls with binomial sampling noise. Test it: thin the counts of every library of the heat cohort binomially to a common coverage per CpG (for example the minimum over the ten libraries at that CpG), recompute the pairwise correlations and the rank of the triploid subset; repeat with a few seeds. Also correlate pairwise r with the coverage of the two libraries.
  3. Is the similarity instead a consequence of the level (the three high diploids lower the diploid correlations)? Check the correlations among the seven libraries that are not high: are the triploid pairs still the most similar? Is the pair D11, D12 (zr3534_1, zr3534_2) unusual?
  4. The gap rule. Would other label-free rules pick the same five libraries (for example a two-component mixture, k-means with k = 2, the top quartile)? In the 24 libraries of the control project the genome-wide values form a continuum with a tail (bismark summary: 03_input_files/bismark_reports/061021-big_bismark_summary_report.txt, rows zr3644_1 to zr3644_24): does "five libraries stand apart" survive that, or should it be worded as "the five highest"?
  5. Item 4, the slope analysis: can the association between a CpG's slope on the level (estimated in one cohort) and the ploidy difference of the other cohort arise without any level effect, for example because CpGs with high variance among animals have both large slopes and large differences? The claim is about sign. Test it by recomputing the decile table after permuting the level values among the libraries of the cohort in which the slope is estimated (several permutations), and report whether the asymmetry (47% against 1% in the top tenth) disappears.
  6. The draft's own hypothesis is that heat changed methylation differently in diploids and triploids. Could that produce exactly what I-49 describes (three of five heat-treated diploids with a raised level; uniform triploids)? List what in the data speaks for and against it, and check whether any sentence of I-49, of the update of I-01 or of the "What follows" paragraph of I-06 excludes that hypothesis more firmly than the data allow. Quote such sentences and propose replacements.
  7. Item 2: is the argument about non-CpG methylation and bisulfite conversion sound? Check the numbers in bismark_by_library.csv.
  8. The ELISA paragraph under "Why it matters": check the group means and the six paired animals against 04_output/03_global_methylation/tables/ (elisa_percent_5mC.csv, group_summary.csv, pct_mCpG_by_library.csv). Is "The two assays agree in direction" fair, given that three of the eight ELISA diploid heat animals are flagged above the top standard?
(c) Wording. Is the title of I-49 supported ("are carried by")? Is the severity S1 justified by the definitions at the top of issues.md? Is anything in "What would settle it" or "Needs from authors" wrong or missing?

Return check_id "S2-I49-evidence-and-alternatives".`,
  },
  {
    id: 'S3-code-review-library-level',
    prompt: `${COMMON}
Your check: a code review of ${REPO}/06_diagnostics/library_methylation_level.R (about 520 lines, eight parts), written today by the analyst and not yet reviewed. Its outputs are the 13 tables in ${REPO}/06_diagnostics/library_methylation_level/. Helper functions it uses are in ${REPO}/02_R_functions/ and settings in ${REPO}/run_config.R.

Tasks
1. Read every line. For each of the eight parts decide whether the code computes what the header says. Report defects under code_findings with a concrete failure scenario. Look in particular at:
   - relabel_reference(): the construction of the groupings (the numbers of groupings must be 252, 924, 21, 210, 4,410 and 100 for the six comparisons), the handling of the observed grouping and its mirror image, the chunked matrix arithmetic, and the claim in the comment that dropping CpGs whose library values span no more than the threshold changes no count (prove it or give a counterexample; think about the tie tolerance and about libraries outside the comparison);
   - above(), at_threshold() and the tie tolerance: can a CpG be counted or missed wrongly because of floating-point error, given that pooled percentages are ratios of sums of non-negative numbers with at most a few thousand reads? Are the counts consistent between the vectorised matrix path and the rowSums path?
   - above_largest_gap(): ties, and what happens if the largest gap is at the bottom of the distribution;
   - the slope computation and the decile table (cut() with quantile breaks: duplicated breaks, boundary handling);
   - part 7: is the enumeration right, and is "largest_mean_r_of_another_subset" what its name says;
   - part 8: the definition of "apart", the single-row case, the join on chr and start;
   - round_cols(): does rounding to six significant digits alter any integer-valued column or any p-value in a way that matters;
   - the p-value definitions (share of groupings at or above the observed, the observed included): are they valid exact permutation p-values, and are they described correctly in the header and in ${REPO}/06_diagnostics/README.md (section "What the group differences consist of")?
2. Reproduce. Copy the script to your scratch directory, change only the line that sets out_dir so that it points into your scratch directory, and run the copy from the repository root (about 2 minutes, under 3 GB; check free memory first). Compare every table with the committed one (md5 or all.equal). Run it a second time and confirm the two runs are identical.
3. Recompute independently, with different code (for example data.table by-group operations instead of matrices, and fisher/wilcox functions of base R): (i) contrast_direction.csv rows for the two within-cohort ploidy contrasts in legacy mode, including n_exactly_at_threshold; (ii) relabelling_reference.csv row "acute heat shock inside the heat cohort" in legacy mode (100 groupings); (iii) level_group_tests.csv for pct_mCpG_genome; (iv) pair_similarity.csv for the heat cohort.
4. Compare with the other diagnostic: 06_diagnostics/dml_sensitivity/permutation_reference.csv reports 3,350 CpGs for the control cohort and medians of 4,019.5 and 2,912.5 where this script reports 3,349, 4,018.5 and 2,912.5. The header attributes the differences to CpGs exactly at 20 points. Verify that explanation exactly: find the CpGs involved and show how each script treats them (read the counting code of 06_diagnostics/dml_sensitivity.R).
5. Check the header (Purpose, Reads, Writes, Notes) against the code and the files, line by line, and the run time and memory it states against what you measured (use /usr/bin/time if installed, else R's proc.time and gc).
6. Say whether any committed number in the 13 tables is wrong, and list every number of the tables that would change under your proposed fixes.

Report code defects under code_findings; wrong or misleading statements in the header or README under discrepancies; and missing checks or outputs under additions. Return check_id "S3-code-review-library-level".`,
  },
  {
    id: 'S4-updates-and-cross-file',
    prompt: `${COMMON}
Your check: the smaller pieces of text written in the same pass, and consistency across files. State of the files: since those pieces were written, the two long issues of the pass, I-06 and I-49, have each been checked by another agent and rewritten (their corrections are listed at the end of issues.md), and the paragraphs named below received small corrections from those checks. Read everything as it stands now; where a paragraph named below still disagrees with the rewritten I-06 or I-49, that is a finding.

Part 1. In ${REPO}/issues.md, the paragraphs that begin "**Update 2026-11-19: what the within-cohort analyses show" (issue I-01), "**Update 2026-11-19: the interaction, tested" (I-02), "**Update 2026-11-19: sensitivity (computed).**" (I-04), "**Update 2026-11-19: the control-only analysis" (I-07) and the paragraph "**Update 2026-11-19 (computed; \`06_diagnostics/library_methylation_level/corrected_dml_library_values.csv\`" (I-08). Recompute every number from the tables named (06_diagnostics/dml_sensitivity/*.csv and 06_diagnostics/library_methylation_level/*.csv) and report mismatches. Check that each file and column reference exists and holds what the text says. Check each inference: does the conclusion follow from the numbers given (for example in I-04, "The data give no reason to drop them"; in I-07, "The draft's 'baseline' list is thus mostly a list of differences between the diploids and triploids of the heat experiment"; in I-02, "half of the 'interaction-DML' are not interaction loci under the draft's own model and thresholds")?

Part 2. Consistency. Compare with each other and with the tables:
  - the Summary table at the top of issues.md (rows I-06 and I-49; severities; statuses) against the issue headings and texts;
  - the header paragraph "**Updates.**" of issues.md (fourth pass);
  - ${REPO}/README.md, the paragraph "Where the work stands" (49 problems, nine bearing on main conclusions; the nine numbered statements): count the issues and the S1 issues in issues.md yourself;
  - ${REPO}/06_diagnostics/README.md: the row for library_methylation_level.R, the section "Result of the sensitivity analyses" (every number of its table and bullets against dml_counts.csv and permutation_reference.csv) and the section "What the group differences consist of" (every file description against the actual columns of the 11 tables);
  - ${REPO}/07_documentation/reproducibility_plan.md, section 7, decisions 1, 2, 3, 15 and 16, and the row of phase 5 in section 5: numbers and whether the question put to the authors is the one the evidence raises;
  - ${REPO}/07_documentation/session_handoff.md: State table and "Next actions" against the above.
Part 3. Stale or contradicted statements. The new results contradict or outdate statements written earlier. Search issues.md (all issues), the READMEs under ${REPO}/01_pipeline, ${REPO}/04_output, ${REPO}/05_manuscript and ${REPO}/06_diagnostics, the plan, and the headers of the scripts under ${REPO}/01_pipeline and ${REPO}/06_diagnostics for statements such as: the sensitivity analyses being in progress or planned; "The ploidy contrast is balanced across the two experiments ... so batch does not align with ploidy" (I-01, now qualified by its update: is the qualification sufficient and placed where a reader of the original paragraph will see it?); promises that "the pipeline will produce" something that now exists or that the results make moot; statements that a ploidy difference exists or that ploidy DML are a valid list; statements in I-41, I-05, I-09, I-44 or I-45 that refer to I-06 for something I-06 does not say, or that the new results change. Quote each and propose the replacement or the added sentence.
Part 4. The manuscript side. Using the text copy of the draft named in the rules, check the "Where" lines of I-06 and I-49 and the statements in those issues about what the draft says (for example "F(1,18) = 14.49", the list of sections and figures that depend on the ploidy DML). Report wrong section or figure references.

Return check_id "S4-updates-and-cross-file".`,
  },
  {
    id: 'S5-literature',
    prompt: `${COMMON}
Your check: the published sources that the new text relies on, and what the literature says about the alternative explanations in issue I-49. You need outbound web access; try the Crossref REST API (https://api.crossref.org/works/<DOI>), the Europe PMC REST API (https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:<DOI>&resultType=core&format=json, and its fullTextXML endpoint for open-access articles), the bioRxiv API (https://api.biorxiv.org/details/biorxiv/<DOI>), DataCite (https://api.datacite.org/dois/<DOI>), doi.org, and publisher pages. Earlier attempts in this project to open the article of Wreczycka et al. at the publisher and at bioRxiv were refused (403 and 429); try the APIs above before giving up, and wait between requests. If a host is blocked by the network policy, say so and move on; do not try to get around it.

Part 1. Verify each citation below: exact metadata, and whether the paper says what the text uses it for. Give a verbatim excerpt from the abstract or text when you can reach it.
  1. Wreczycka K, et al. 2017. Strategies for analyzing bisulfite sequencing data. J Biotechnol 261:105, doi:10.1016/j.jbiotec.2017.08.007. Issue I-06 of ${REPO}/issues.md (paragraph "Why it matters") says that the methylKit authors' own benchmark "shows a markedly higher false positive rate without the correction and recommends overdispersion = "MN" with a chi-squared test". What does the paper actually report about methylKit with and without overdispersion correction, with the F test and the chi-squared test, and what does it recommend? This is the most important item of your check.
  2. Feng H, Conneely KN, Wu H. 2014. Nucleic Acids Res 42:e69, doi:10.1093/nar/gku154 (DSS; used in I-06 "Options" and in 06_diagnostics/README.md).
  3. Park Y, Wu H. 2016. Bioinformatics 32:1446, doi:10.1093/bioinformatics/btw026 (DSS for general designs; I-02, I-06).
  4. Korthauer K, et al. 2019. Biostatistics 20:367-383, doi:10.1093/biostatistics/kxy007 (dmrseq; 06_diagnostics/README.md).
  5. Phipson B, Smyth GK. 2010. Stat Appl Genet Mol Biol 9:Article 39, doi:10.2202/1544-6115.1585 (permutation p-values; 06_diagnostics/README.md).
  6. McCullagh P, Nelder JA. 1989. Generalized Linear Models, 2nd edition, Chapman and Hall (the overdispersion correction; I-06). Metadata only.
  7. Leek JT, et al. 2010. Nat Rev Genet 11:733, doi:10.1038/nrg2825 (batch effects; I-01).
  8. Gelman A, Stern H. 2006. Am Stat 60:328, doi:10.1198/000313006X152649, and Nieuwenhuis S, et al. 2011. Nat Neurosci 14:1105, doi:10.1038/nn.2886 (I-02).
  9. Liu Y, et al. 2012. Genome Biol 13:R61, doi:10.1186/gb-2012-13-7-r61 (C to T variants in bisulfite data; I-08).
Part 2. Find and verify sources for statements that I-49 currently makes without a citation, or that would let it be more precise. For each, give full metadata with DOI and a verbatim excerpt, and say plainly if you found nothing:
  a. Library preparation can shift the absolute, genome-wide methylation level measured by WGBS between libraries (for example bisulfite-induced fragment loss or amplification bias). A candidate is Olova N, et al. 2018, Genome Biol 19:33, doi:10.1186/s13059-018-1408-2; verify what it shows.
  b. DNA methylation patterns in oysters (Crassostrea gigas or other oysters) differ between families or follow genotype. Candidates: Olson CE, Roberts SB, "Indication of family-specific DNA methylation patterns in developing oysters" (bioRxiv; find the DOI and year), and Silliman K, Spencer LH, White SJ, Roberts SB, 2023, "Epigenetic and genetic population structure is coupled in a marine invertebrate", Genome Biol Evol (find the DOI). Search for others.
  c. Published comparisons of DNA methylation between diploid and triploid oysters or other bivalves, and what they report about global level and about variation among individuals.
  d. Guidance that methylation analyses with biological replicates need a model with overdispersion or a beta-binomial model, and that related or co-reared animals are not independent replicates (pseudoreplication). Candidates: Hurlbert SH 1984, Ecol Monogr 54:187, doi:10.2307/1942661; Lazic SE 2010, BMC Neurosci 11:5, doi:10.1186/1471-2202-11-5. Verify.
  f. Apparent methylation of mitochondrial DNA in bisulfite sequencing as an artifact of incomplete conversion of circular or structured mitochondrial DNA, and how it depends on library preparation (issue I-01 now reports 0.31% to 0.58% pooled methylation at mitochondrial CpGs in the ten heat libraries against 1.98% to 3.64% in the twelve control libraries). A candidate is Mechta M, et al. 2017, Front Genet 8:166, doi:10.3389/fgene.2017.00166; verify what it shows and search for others.
  e. Whether commercial triploid Pacific oysters are typically produced as one or a few crosses (tetraploid male by diploid female), which bears on whether the triploids of an experiment are likely to be related. Give a primary source if one exists (for example Guo X, DeBrosse GA, Allen SK 1996, Aquaculture 142:149, doi:10.1016/0044-8486(95)01243-5; verify).
Part 3. Read the paragraph "**Correction to "Why it matters" (2026-11-19).**" of I-06. It now quotes an account of the benchmark that one agent wrote into the header of ${REPO}/06_diagnostics/dml_sensitivity.R (note "Which of the two tests", about lines 107 to 115): lowest specificity and best overall accuracy for the default test without the correction, results close to DSS for the F test with the correction, the chi-squared test with the correction named as a good compromise, three samples per group in the simulation. Confirm or refute each of these four points from the paper or the preprint (doi:10.1101/109512), with verbatim excerpts. Then write the exact replacement for the original sentence on the benchmark (plain statement of what the paper reports and recommends, with the citation), or state that it remains unverified and why.

Put each source under citations. Put proposed text for issues.md under additions (with where_to_add), and anything in the current text that a source contradicts under discrepancies. Return check_id "S5-literature".`,
  },
]

phase('Check')
const results = await parallel(CHECKS.map(c => () =>
  agent(c.prompt, { label: `check:${c.id}`, phase: 'Check', schema: SCHEMA })
))
const ok = results.filter(Boolean)
log(`${ok.length} of ${CHECKS.length} checks returned`)
return {
  checks: ok.map(r => ({
    check_id: r.check_id,
    statements_checked: r.statements_checked,
    confirmed: r.confirmed,
    discrepancies: r.discrepancies.length,
    high: r.discrepancies.filter(d => d.severity === 'high').length,
    challenges_not_standing: r.challenges.filter(c => c.verdict !== 'claim_stands').length,
    additions: r.additions.length,
    code_findings: (r.code_findings || []).length,
    citations: (r.citations || []).length,
  })),
  full: ok,
}

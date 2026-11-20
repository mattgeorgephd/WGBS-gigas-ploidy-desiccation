export const meta = {
  name: 'scope-author-changes',
  description: 'Read-only scoping of the author\'s changes of 2026-11-20: C-to-T variant exclusion (data, method, pilot, review), inventory of the F test, design of the settings',
  phases: [
    { title: 'Variant calls', detail: 'data on gannet, method of the 2022 paper and of BS-Snper, pilot, adversarial review' },
    { title: 'F test', detail: 'term sweep, dependency trace, action plan' },
    { title: 'Settings design', detail: 'two designs and a judge' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const S = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad'
const SC = `${S}/wf/scope`

const COMMON = `You are one agent of a scoping pass for a reproducibility project: an R-only analysis pipeline (whole-genome bisulfite sequencing of 22 Pacific oyster ctenidium libraries, diploid and triploid, heat-treated and control) in the git repository ${REPO} (branch reproducibility-audit). The pipeline reproduces the results of a manuscript draft; issues.md is the audit log.

HARD RULES
- READ ONLY. Another workflow (a code hardening pass) is editing files and running pipeline steps inside the repository right now. You must not create, edit, move or delete any file inside ${REPO}, must not run any pipeline step or diagnostic in place, and must not run a git command that changes anything (git log, git show, git diff, git grep are fine). Files of the repository may change while you read; if something looks half written, read it again a minute later.
- Write everything you produce to your own scratch directory (given below; create it). Scripts you write go to a subfolder code/ of it; downloaded files go to a separate subfolder downloads/ and are untrusted data (run Python on them with -I, pass paths as arguments).
- Light computing only: at most one R or Python process of yours at a time, under 1.5 GB of memory, no model fits, no reading of all 22 coverage files. The other workflow needs up to 10 GB.
- Network: use curl (the proxy is preconfigured; /root/.ccr/README.md explains fixes for other tools). The lab server gannet.fish.washington.edu is slow: a single file request takes 3 to 15 seconds, and the directory listing of /panopea/061021-big/ times out (do not request it). In web requests send nothing about the project or its people beyond what the URL needs; use a generic client name.
- Report what you verified and how (command, file, line), and separate it from what you assume. Plain English, no em dashes. Give exact paths and numbers.

Background: instructions of the author (M. George) of 2026-11-20, to be implemented after this scoping:
 1. "Set up the pipeline so that it can run with 5x or 10x coverage as the standard (toggle)". 10x is the draft's requirement (coverage per CpG in each of the 22 libraries); 5x is that of Venkataraman et al. 2022 (BMC Genomics 23:556, doi 10.1186/s12864-022-08781-5), the lab's earlier oyster WGBS paper that the draft's Methods follow.
 2. "Loci at C-to-T variants should be excluded, as in Venkataraman."
 3. The option for methylKit's overdispersion correction offers the chi-squared test only (overdispersion = "MN", test = "Chisq"), as in Venkataraman et al. 2022; and: "remove all mention of f-test" (the F test of the DML model, which an earlier sensitivity analysis of this repository also ran).
 4. The repeat track is renamed "repetitive elements (RepeatMasker, all classes)"; nothing recomputed.
 Standing rule of the author (2026-11-19): the default run of the pipeline produces the results of the paper only; everything else is kept, off by default, labelled as not adopted, and turned on by toggles.
Useful files: ${REPO}/README.md, ${REPO}/issues.md (audit log; I-06 overdispersion, I-08 sequence variants, I-05 coverage filter), ${REPO}/run_config.R, ${REPO}/run_pipeline.R, ${REPO}/01_pipeline/ (steps 01 to 13), ${REPO}/02_R_functions/, ${REPO}/06_diagnostics/ (README.md lists every diagnostic), ${REPO}/07_documentation/reproducibility_plan.md and session_handoff.md, ${REPO}/07_documentation/upstream_processing/ (job scripts of the alignment run). The download cache of the pipeline is /home/user/gannet_cache (coverage files in cov_files/, the reference genome in genome/cgigas_uk_roslin_v1_genomic-mito.fa).`

const REPORT_SCHEMA = {
  type: 'object',
  properties: {
    summary: { type: 'string', description: 'Plain summary of what you found, with the numbers that matter' },
    facts: { type: 'array', items: { type: 'object', properties: {
      claim: { type: 'string' }, evidence: { type: 'string', description: 'command, file and line, URL, or quote' }, verified: { type: 'boolean' } }, required: ['claim', 'evidence', 'verified'] } },
    recommendation: { type: 'string' },
    risks: { type: 'array', items: { type: 'string' } },
    open_points: { type: 'array', items: { type: 'string' }, description: 'what could not be settled, and decisions that belong to the author' },
    files_written: { type: 'array', items: { type: 'string' }, description: 'absolute paths in your scratch directory' },
  },
  required: ['summary', 'facts', 'recommendation', 'open_points', 'files_written'],
}

// ---------------------------------------------------------------- Variant calls
const V1 = `${COMMON}
Your scratch directory: ${SC}/V1-data
Your task: establish what data exist for calling sequence variants in the 22 libraries, and how an R script can read them.

Known so far (verified by the parent session with HEAD requests on 2026-11-20): https://gannet.fish.washington.edu/panopea/061021-big/zr3534_1_R1_bismark_bt2_pe.deduplicated.sorted.bam exists (5,669,455,726 bytes, accept-ranges: bytes) with its index .bam.bai (1,610,472 bytes); zr3644_1_R1_bismark_bt2_pe.deduplicated.sorted.bam exists (4,304,427,733 bytes). These are the sorted, deduplicated Bismark alignments of the June to July 2021 run that produced the coverage files of this project (${REPO}/07_documentation/upstream_processing/README.md and 061021-big_job.sh; Bismark 0.22.3, --non_directional, paired end).

Do this
1. The 22 libraries are listed in the sample sheet of the repository (find it under ${REPO}/03_input_files/metadata/ and through ${REPO}/02_R_functions/sample_sheet.R) and in /home/user/gannet_cache/cov_files/samples.txt. For each, HEAD the sorted BAM and its .bai and record URL, size, last-modified. Write the table to your scratch directory as bam_files.csv. Sum the sizes.
2. Look for variant calls that already exist for these libraries: HEAD requests for the usual BS-Snper output names in that directory and in the project directory (the manifest of the project directory is ${REPO}/07_documentation/audit/gannet_manifest_2026-11-18.tsv; search it for snp, vcf, bam), and in the OSF archive of Venkataraman et al. 2022 only as a pattern for file names. Say plainly if none exist.
3. Test whether R can read regions of a BAM on the server without downloading the file: Rsamtools is installed (check Rhtslib's build for libcurl support). Try scanBamHeader and scanBam on the URL with a region (ScanBamParam(which = ...)), with the index fetched to your scratch directory or by URL. htslib uses libcurl, which needs the proxy and its CA bundle (environment variables such as https_proxy and CURL_CA_BUNDLE=/root/.ccr/ca-bundle.crt; see /root/.ccr/README.md). If it works, time 200 region queries of 2 bases each, spread over the genome, on one library, and report queries per minute. If it does not work, say exactly why (error text) and try the alternative of fetching byte ranges yourself: the .bai gives the BGZF virtual offsets of each region, so a script can fetch only the blocks it needs; say whether that is practical in R alone.
4. Measure download speed: fetch one 300 MB byte range of one BAM (curl -r) to your downloads folder, time it, delete it. From it estimate the time to download each BAM in turn. The session has about 30 GB of free disk; say what plan fits (one BAM at a time, deleted after use).
5. From a few thousand alignments of one heat library and one control library (region reads or the partial download), tabulate: FLAG values, MAPQ, read lengths, the tags present (XR, XG, XM), the four combinations of XR and XG (the libraries are non-directional), the share of read pairs whose mates overlap each other, and whether mate pairs are properly flagged. This is needed to assign each read to the DNA strand it came from.
6. Say what else an R script would need to turn these alignments into per-position base counts by strand of origin at a given list of CpG positions (Rsamtools::pileup with ScanBamParam and flags, or scanBam plus own counting; the tag XG is not usable in pileup directly, so say how strand of origin can be had: by FLAG as BS-Snper does, by splitting on the tag with a filter function, or otherwise). Do not run a pilot over many loci; another agent does that after you.
Return the schema. In facts include one fact per item above.`

const V2 = `${COMMON}
Your scratch directory: ${SC}/V2-method
Your task: establish exactly how Venkataraman et al. 2022 identified and excluded C-to-T variants, and exactly what rule BS-Snper applies, so that the exclusion can be reproduced for this project "as in Venkataraman".

Sources to read in full (primary sources only; quote the lines that matter)
- The article: Venkataraman et al. 2022, BMC Genomics 23:556 (PMC9351233). Full text XML: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9351233/fullTextXML . Its Methods, section "SNP identification", say: sorted deduplicated BAM files were merged (samtools merge) to identify SNPs; variants were also identified in individual BAM files; default BS-Snper settings except a minimum coverage of 5x; "The C->T SNPs were used to exclude loci of interest". Results: 13,234,183 unique SNPs "in individual or merged BAM files", of which 300,278 C->T SNPs; 11,238,223 CpG loci at 5x became 10,939,918 after removal; 315 of 1,599 candidate DML overlapped C->T SNPs and were removed.
- Its code. The article's archive is doi.org/10.17605/OSF.IO/YGCTB (Open Science Framework, "Oyster Gonad Methylation"), which mirrors a GitHub repository of the Roberts lab (look for the repository of that project, for example RobertsLab/project-gigas-oa-meth, and its code/ folder: a notebook named 07-BS-SNPer.ipynb or similar, and the notebooks or R scripts that use its output for the methylation characterisation and for the DML). Use the OSF API (https://api.osf.io/v2/) or raw.githubusercontent.com with curl. Find and quote: the BS-Snper command lines and parameters; how "C->T SNPs" were taken from the output (which columns and values; reference C with alternative T only, or also reference G with alternative A, which is the same change on the other strand; homozygous and heterozygous; which FILTER values); whether calls from the merged file, from individual files, or their union were used; at which stage loci were excluded (from the coverage files before methylKit, from the united table, or from the DML list) and with which command (for example bedtools intersect -v) and coordinate convention; whether a CpG was excluded when the variant lies at its C, at its G, or either.
- BS-Snper itself: source and README at github.com/hellbelly/BS-Snper (BS-Snper.pl and the C source it calls), and Gao et al. 2015, Bioinformatics 31:4006 (doi 10.1093/bioinformatics/btv507) if an open copy can be read. Establish: the default values of every parameter (--minhetfreq, --minhomfreq, --minquali, --mincover, --maxcover, --minread2, --errorate, --mapvalue and any other); the decision rule in words and in code terms (which reads count as Watson and Crick, by which FLAG values; which strand's reads are used for which reference base; how a C-to-T variant is told from an unmethylated cytosine; how heterozygous and homozygous calls are made; what is written to the VCF and to the candidates file); and whether the rule assumes directional libraries (the libraries here and in the 2022 study are non-directional, Bismark --non_directional).
Then
- Write the rule as a specification that an R script could implement for a given list of CpG positions and a set of BAM files: inputs, counts needed per position and strand, thresholds, output columns. State every place where an R implementation would differ from BS-Snper (for example base quality handling, overlapping mates, indels, pooled calling over a merged BAM against per-library calling) and what the difference could do to the calls.
- State the alternative of running BS-Snper itself (Perl and C; is perl installed here, does it need compiling, what would it need as input) and what "as in Venkataraman" strictly requires.
- The animals here are diploid and triploid. Say what the thresholds mean for a triploid (one variant allele of three gives an expected frequency of one third; two of three give two thirds) and whether BS-Snper's defaults would call such sites.
Write your specification to spec.md in your scratch directory and return the schema; put the quotes from the code in facts with file names and line numbers or cell numbers.`

const V3 = (v1, v2) => `${COMMON}
Your scratch directory: ${SC}/V3-pilot
Your task: a pilot of the C-to-T variant exclusion on real data, and from it a design for a pipeline step, with measured costs. Two agents prepared the ground; their reports follow. Their files are in ${SC}/V1-data and ${SC}/V2-method (read spec.md there).

REPORT OF THE DATA AGENT
${v1}

REPORT OF THE METHOD AGENT
${v2}

Do this
1. Take as test loci: every draft DML (the committed tables ${REPO}/04_output/05_differential_methylation/tables/dml_ploidy_legacy.csv and dml_heat_legacy.csv: 1,038 and 2,351 loci; columns chr, start, end where start is the 0-based position of the C of the CpG and the CpG occupies the two 1-based positions start + 1 and start + 2) and 3,000 CpGs drawn at random (fixed seed) from the analysed set (united_cpg_legacy.csv.gz in the same folder, 480,579 rows; read only the first three columns).
2. For four libraries (two heat-treated and two control, one diploid and one triploid of each; take them from the sample sheet), obtain for each test locus the counts of A, C, G, T at the C position and at the G position, separately for reads that come from each DNA strand, by the access method the data agent found to work (remote region reads if they work at a useful speed; otherwise download one BAM at a time into your scratch directory and delete it after use; the session has about 30 GB of free disk). Write the R code so that it would run unchanged over all 22 libraries and any list of loci. Record the time per library.
3. Apply the rule of the method agent's specification and list the loci that it calls C-to-T variants (at the C) or G-to-A variants (at the G), per library and pooled over the four libraries. Report: how many of the 1,038 and 2,351 DML and of the 3,000 random CpGs are called, the per-library allele frequencies at the called loci, and for the called DML the methylation values of the libraries (from united_cpg_legacy.csv.gz: numCs and numTs columns in library order; the sample sheet gives the order). Print the raw strand-specific counts of 15 called and 15 uncalled loci so that a reader can check the logic by eye.
4. Check the pilot against something independent: at a locus that is a C-to-T variant in a library, the reads of the strand that carries the cytosine must show T (Bismark counts it unmethylated) and the reads of the other strand must show A at the paired position where the reference has G. Verify on the raw reads of five called loci that this is what the alignments show. Also compare, at the test loci, the total C and T counts of your pileup with numCs and numTs of the united table for the same library (they will differ where the pipeline normalised coverage or where base quality filters differ; say how closely they agree and why not exactly).
5. From the measured times, estimate the cost of the full job for 22 libraries: (a) at the candidate DML only, (b) at all analysed CpGs at 10x (480,579) and at 5x (2,297,236), (c) genome-wide as BS-Snper does. Say which of these "as in Venkataraman" needs (they removed variant loci from all characterised CpGs and from the DML).
6. Propose the design of a pipeline step in R (name and number among the steps of ${REPO}/01_pipeline/, inputs fetched from gannet, the table it writes, what it caches, what could be committed to the repository so that a user without the alignments can still run the exclusion, the settings it needs in run_config.R, and where steps 04 and 05 would apply the exclusion). Say what it would cost to run BS-Snper itself instead and what would then be stored.
Return the schema; put the measured numbers in facts.`

const V4 = (v3) => `${COMMON}
Your scratch directory: ${SC}/V4-review
Your task: an adversarial review of a pilot and a design for excluding CpG loci at C-to-T sequence variants. You fix nothing; you try to break it. The pilot's files are in ${SC}/V3-pilot (code, tables, logs); the method specification is ${SC}/V2-method/spec.md; the data notes are in ${SC}/V1-data. The pilot agent reported:

${v3}

Try to refute it
1. Strand logic. The libraries are non-directional and paired end. For a C-to-T variant at a reference C, which reads are informative? Work it out from first principles for all four Bismark strand types (OT, CTOT, OB, CTOB; tags XR and XG; FLAG values of read 1 and read 2) and check that the pilot's code assigns every read to the right class, including read 2 of each pair. Take ten loci and print the raw reads (position, FLAG, XR, XG, base at the locus) to check the counts by hand.
2. Double counting where the two mates of a pair overlap at the locus; soft clipping, insertions and deletions near the locus; MAPQ and base quality filters; duplicate flags; positions at contig ends.
3. Thresholds: does the pilot apply the BS-Snper defaults as the specification states them, with minimum coverage 5? What do the thresholds do for triploid animals? What happens at loci with fewer than the minimum number of informative reads in a library: are they treated as "no variant" (risky) or as "not assessable"? How many of the test loci are not assessable in at least one library?
4. Per-library calls against pooled calls: the 2022 study called variants on a merged BAM and on individual BAMs. Does the pilot's pooled rule match that? Could a variant present in a few animals be missed, or a sequencing error be called?
5. Biological sanity: are the called loci plausible? Among the DML called as variants, do the methylation values of the libraries sit at genotype-like levels (near 0, one third, one half, two thirds, 1)? Among random CpGs, is the share called in line with the 2.7% of CpG loci that Venkataraman et al. 2022 removed (298,305 of 11,238,223) and the 19.7% of their candidate DML?
6. The design of the pipeline step: can a user reproduce it; what breaks if the server is slow or a download is cut; is anything too large to commit; does the default run of the pipeline stay the paper's analysis.
Run small checks of your own on the pilot's outputs and, where needed, on a few region reads. Return the schema: in summary give your verdict (sound, sound with changes, or not sound) and under risks every defect with a concrete failure scenario and the change that would fix it.`

// ---------------------------------------------------------------- F test
const F1 = `${COMMON}
Your scratch directory: ${SC}/F1-terms
Your task: find every place in the repository that mentions or uses the F test of methylKit's differential methylation model, by searching for terms. The author's instruction is "remove all mention of f-test"; the model with the overdispersion correction is to be offered with the chi-squared test only.

Background: 06_diagnostics/dml_sensitivity.R fitted variants of the DML model, among them methylKit's calculateDiffMeth with overdispersion = "MN" and test = "F" (analysis ids ending in _mn_f) and with test = "Chisq" (ids ending in _mn_chisq). Results and text built on them are spread over 06_diagnostics/, issues.md, README files, 07_documentation/ and the table of manuscript values. The draft's own analyses of variance (ELISA of global methylation, qPCR) also report F statistics; those are results of the paper and are NOT meant, but must be listed separately so that nobody removes them by mistake.

Do this
1. Search every tracked file except .git and 08_archive (git ls-files gives the list; include CSV, TSV, JSON, logs, R, md, txt; for PDF and PNG figures say whether their generating code draws an F-test variant) for these and any further patterns you think of: "F test", "F-test", "F tests", test = "F", test="F", "F statistic", "F(" followed by numbers, "_mn_f", "mn_f", "MN + F", "F and Chisq", "quasi", "default test", "Chisq" (to see the neighbours), "the corrected test", "corrected model", "with the correction", "overdispersion".
2. Write the complete list to hits.csv in your scratch directory: file, line number, the line (trimmed to 300 characters), class (one of: DML model F test, named explicitly; DML model with correction, test not named but could be the F test; chi-squared variant; analysis of variance of the draft or another F statistic; unrelated), and for the first two classes which output file or analysis id the statement or code rests on if the line shows it.
3. Summarise by file and class (counts), and list separately: code that would have to change (script, lines), output files that exist only for the F variant (to delete), output tables with rows or columns for the F variant (to regenerate or filter), figures that show it, and text passages (file, heading or issue number).
4. Records of earlier agent work are kept under 07_documentation/agent_workflows/ (results JSON, scripts, verification scripts and logs). List which of them mention the F test and how often; they are historical records, so do not propose edits to them, only report.
Return the schema with the counts in summary and facts; files_written must include hits.csv.`

const F2 = `${COMMON}
Your scratch directory: ${SC}/F2-dependencies
Your task: trace what depends on the F test of methylKit's differential methylation model, starting from the code, so that statements that rest on it without naming it are found. The author's instruction is "remove all mention of f-test"; the model with the overdispersion correction is to be kept with the chi-squared test only (overdispersion = "MN", test = "Chisq").

Do this
1. Read 06_diagnostics/dml_sensitivity.R in full, and every other script of 06_diagnostics/ and 01_pipeline/ that fits or reads a model with overdispersion or reads the outputs of dml_sensitivity (06_diagnostics/README.md lists the diagnostics; library_methylation_level.R reads lists of "corrected" DML). For each analysis id (A0 to A5 and their variants) record: the call (covariates, overdispersion, test), the libraries used, the analysis modes (legacy, corrected), the output files it writes, and the numbers it gives (06_diagnostics/dml_sensitivity/dml_counts.csv, analyses.csv, permutation_reference.csv, dispersion_summary.csv, interaction_vs_shared_*.csv). Note how methylKit behaves: with test = "F" the F test is applied only where the dispersion estimate is above 1, and a chi-squared test otherwise (methylKit:::logReg; verify in the installed package, version 1.28.0).
2. List every output file, table row, table column and figure panel that comes from an F-test fit, directly or through a later script (for example 06_diagnostics/library_methylation_level/corrected_dml_library_values.csv and corrected_dml_group_apart.csv).
3. Then read issues.md in full (it is long; the issues I-01, I-02, I-06, I-07, I-08, I-49 and the corrections table near the end matter most), README.md, 06_diagnostics/README.md, 07_documentation/reproducibility_plan.md, 07_documentation/session_handoff.md, 03_input_files/manuscript/manuscript_values_map.csv (columns note and expression_*), and find EVERY statement that rests on a number or a list from an F-test fit, whether or not it names the test. Typical forms: "with the correction", "the corrected model", "the corrected test", counts such as 2 and 16 DML, 24 and 287 significant, 125 CpGs, an F statistic above 84. For each: file, line or issue and paragraph, the statement (trimmed), the output it rests on, and whether a chi-squared counterpart exists already (the legacy-mode fits A1_ploidy_mn_chisq and A1_heat_mn_chisq exist; say what else exists) or would need a new fit (say which: analysis, libraries, mode).
4. For each statement propose one of: replace the numbers with those of the chi-squared fit (give the new numbers where the fit exists); needs a new chi-squared fit, then rewrite; withdraw (the statement only makes sense for the F test). Say which conclusions of the audit log would change in substance if the chi-squared fit is used instead (for example 14 and 185 DML instead of 2 and 16).
5. Estimate the computing needed for the new chi-squared fits from the timings in the logs or session_info of the diagnostics (fits are cached under 06_diagnostics/dml_sensitivity/cache/; say whether a cached fit could be re-tested with another test without refitting, which would make the change cheap: look at what the cache stores).
Write statements.csv and analyses.csv to your scratch directory and return the schema.`

const F3 = (f1, f2) => `${COMMON}
Your scratch directory: ${SC}/F3-plan
Your task: merge two inventories into one action plan for removing the F test of methylKit's DML model from the repository and keeping the model with the overdispersion correction with the chi-squared test only. Their files are in ${SC}/F1-terms (hits.csv) and ${SC}/F2-dependencies (statements.csv, analyses.csv). Their reports:

TERM SWEEP
${f1}

DEPENDENCY TRACE
${f2}

Do this
1. Cross-check the two: every explicit hit of the sweep must be covered by the trace or explained, and every statement of the trace must be locatable (open the files and confirm line by line for at least 40 entries chosen across files; fix wrong line numbers).
2. Write plan.md in your scratch directory with, in this order: (a) code changes by file and line (dml_sensitivity.R and every script that reads its F-test outputs), keeping the chi-squared variants and adding those that are needed as replacements; (b) the fits to run, with the cheapest way to get them and the estimated minutes; (c) output files to delete, and tables and figures to regenerate; (d) text changes by file and paragraph, each with the action (replace with chi-squared numbers, rewrite after a new fit, withdraw) and, where the numbers exist already, the proposed new wording; (e) rows of 03_input_files/manuscript/manuscript_values_map.csv to change; (f) what must NOT be touched (F statistics of the draft's analyses of variance and anything else that only looks similar), listed by file and line; (g) the historical records under 07_documentation/agent_workflows/ that mention the F test, with the options (leave, delete the files concerned, delete the folder) and what each option loses, without choosing; (h) a check that can be run afterwards to prove that no mention is left (the exact search) and that every number in the rewritten text equals its source table.
3. Say plainly which conclusions of issues.md change in substance when the chi-squared variant replaces the F variant, and which cannot be kept at all.
Return the schema; put the plan's path in files_written and its main points in summary.`

// ---------------------------------------------------------------- Settings design
const DESIGN_BRIEF = `Design the settings of the pipeline so that:
- The default run produces the results of the paper and nothing else ("legacy" analysis: at least 10 reads per CpG in each of the 22 libraries, no upper coverage filter, the first line of each coverage file dropped as the legacy reader did, methylKit's logistic regression without overdispersion correction, DML at q < 0.01 and a difference of more than 20 percentage points). Its output files keep their present names and bytes (file names carry the mode: _legacy).
- Each of the following is off by default and is turned on by a setting in run_config.R and by a switch of run_pipeline.R; its outputs are labelled as optional and can never overwrite or be mistaken for the paper's outputs:
  (a) the coverage filter that the Methods describe (now the mode "corrected": 99.9th percentile filter and the first line read);
  (b) minimum coverage 5x instead of 10x as the standard of the run (2,297,236 CpGs instead of 480,579 in the legacy reading);
  (c) the DML model with methylKit's overdispersion correction and the chi-squared test (overdispersion = "MN", test = "Chisq"); no other test is offered;
  (d) exclusion of CpG loci at C-to-T sequence variants (a new step will supply the list of loci; the author wrote that such loci "should be excluded", so it may become the standard later: make its default a one-line choice);
  (e) the gene-level GO enrichment (step 10);
  (f) the diagnostics of 06_diagnostics/.
  Combinations must work (for example 5x with the correction and the exclusion).
- Steps 06 to 13 work for whichever analyses were run; step 13 (table of manuscript values, expressions per row in 03_input_files/manuscript/manuscript_values_map.csv) reports optional values as not run when their outputs are absent; tables that span analyses (for example dml_counts_by_mode.csv, annotation_summary.csv, cpg_counts_by_coverage.csv) do not lose or mix rows when only some analyses are run.
- The cached fit of step 05 (it now carries a fingerprint of settings and inputs) stays valid per analysis.
- A reader of the repository can tell from a file's name or folder which analysis it belongs to, and the READMEs can say it in one table.
Read the present code first: run_config.R, run_pipeline.R, 02_R_functions/ (paths.R, logging.R, methylkit_helpers.R), every step of 01_pipeline/, 01_pipeline/README.md and 04_output/README.md (the hardening pass is editing some of these right now; note the time you read them). Count how many places would change. Give: the settings and switches with exact names and defaults; the naming of outputs with examples; how each step loops over analyses; what happens to existing committed outputs of the "corrected" mode and of step 10; the changes by file; the run time and memory to expect for a 5x run (step 05 at 10x takes about 18 minutes per mode at two cores and 6.5 GB; say how you scale it); the risks to the byte-identical reproduction of the paper's outputs and how a test would catch them.`

const D1 = `${COMMON}
Your scratch directory: ${SC}/D1-design
Your task: a design proposal. Your angle: the smallest change to the existing code and file names; extend what is there (the modes vector and the mode tag in file names) rather than replace it.
${DESIGN_BRIEF}
Write design.md to your scratch directory and return the schema (recommendation = the design in brief; facts = what you verified in the code).`

const D2 = `${COMMON}
Your scratch directory: ${SC}/D2-design
Your task: a design proposal. Your angle: the clearest separation for a user: the paper's results in one place, each optional analysis in its own clearly labelled place, and one obvious command per choice; accept more code change if it buys clarity and safety.
${DESIGN_BRIEF}
Write design.md to your scratch directory and return the schema (recommendation = the design in brief; facts = what you verified in the code).`

const D3 = (d1, d2) => `${COMMON}
Your scratch directory: ${SC}/D3-judge
Your task: judge two design proposals for the settings of the pipeline and write the specification to build from. The proposals are in ${SC}/D1-design/design.md and ${SC}/D2-design/design.md; their authors reported:

PROPOSAL 1 (smallest change)
${d1}

PROPOSAL 2 (clearest separation)
${d2}

The brief both worked from:
${DESIGN_BRIEF}

Do this
1. Check each proposal against the present code (read the files it names) and against every point of the brief. List what each gets wrong or leaves out, with file and line.
2. Score them on: safety of the paper's outputs (names and bytes unchanged, impossible to overwrite), clarity for a user, amount and risk of code change, how combinations are handled, step 13 and the cross-analysis tables, cache validity.
3. Write build_spec.md in your scratch directory: one design (the better proposal with the best ideas of the other), precise enough to build from: settings and switches (names, values, defaults), naming rules with examples for every step's outputs, the change list by file, the order of work, and the tests that prove (i) a default run reproduces every committed legacy output byte for byte, (ii) no optional run changes a legacy output, (iii) each toggle alone and three combinations run through, (iv) step 13 reports optional values as not run when absent.
4. List the decisions that belong to the author (for example names of settings, whether the C-to-T exclusion is on by default), each with your recommendation and the reason.
Return the schema with the specification's path in files_written.`

function txt(r) {
  if (!r) return '(no report: the agent returned nothing)'
  const facts = (r.facts || []).map(f => `  - ${f.verified ? '[verified]' : '[not verified]'} ${f.claim} (${f.evidence})`).join('\n')
  return `SUMMARY: ${r.summary}\nFACTS:\n${facts}\nRECOMMENDATION: ${r.recommendation}\nRISKS: ${(r.risks || []).join(' | ')}\nOPEN POINTS: ${(r.open_points || []).join(' | ')}\nFILES: ${(r.files_written || []).join(', ')}`
}

const variantLane = (async () => {
  const [v1, v2] = await parallel([
    () => agent(V1, { label: 'variants:data-on-gannet', phase: 'Variant calls', schema: REPORT_SCHEMA }),
    () => agent(V2, { label: 'variants:method-2022-and-bs-snper', phase: 'Variant calls', schema: REPORT_SCHEMA }),
  ])
  const v3 = await agent(V3(txt(v1), txt(v2)), { label: 'variants:pilot-and-design', phase: 'Variant calls', schema: REPORT_SCHEMA })
  const v4 = await agent(V4(txt(v3)), { label: 'variants:adversarial-review', phase: 'Variant calls', schema: REPORT_SCHEMA })
  return { v1, v2, v3, v4 }
})()

const fLane = (async () => {
  const [f1, f2] = await parallel([
    () => agent(F1, { label: 'ftest:term-sweep', phase: 'F test', schema: REPORT_SCHEMA }),
    () => agent(F2, { label: 'ftest:dependency-trace', phase: 'F test', schema: REPORT_SCHEMA }),
  ])
  const f3 = await agent(F3(txt(f1), txt(f2)), { label: 'ftest:action-plan', phase: 'F test', schema: REPORT_SCHEMA })
  return { f1, f2, f3 }
})()

const dLane = (async () => {
  const [d1, d2] = await parallel([
    () => agent(D1, { label: 'design:smallest-change', phase: 'Settings design', schema: REPORT_SCHEMA }),
    () => agent(D2, { label: 'design:clearest-separation', phase: 'Settings design', schema: REPORT_SCHEMA }),
  ])
  const d3 = await agent(D3(txt(d1), txt(d2)), { label: 'design:judge', phase: 'Settings design', schema: REPORT_SCHEMA })
  return { d1, d2, d3 }
})()

const [variants, ftest, design] = await Promise.all([variantLane, fLane, dLane])
return { variants, ftest, design }

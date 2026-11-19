export const meta = {
  name: 'audit-verification',
  description: 'Adversarially verify the audit claims in issues.md and review the pipeline code built so far (read-only on the repository)',
  phases: [
    { title: 'Verify claims', detail: 'one skeptic per major claim, each trying to refute it from the data' },
    { title: 'Review code', detail: 'reviewers by pipeline area, plus a run of steps 06 to 10 from a clean checkout' },
    { title: 'Check findings', detail: 'each reviewer finding list re-examined by an independent verifier' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'
const MANUSCRIPT = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/manuscript/manuscript_gdoc.txt'

const CONTEXT = `
CONTEXT (read all of it before starting)
Repository: ${REPO} (git branch reproducibility-audit). It holds an R-only pipeline that re-derives the results of a manuscript draft on DNA methylation (whole-genome bisulfite sequencing, WGBS) in diploid and triploid Pacific oysters after a heat exposure, together with an audit of that draft.
Layout: 01_pipeline/ (numbered steps 01 to 10, each an R script with a header block), 02_R_functions/ (helpers sourced by every step), 03_input_files/, 04_output/<step>/{tables,figures}, 06_diagnostics/, 07_documentation/, 08_archive/ (the legacy code and outputs of the original analysis; never modified), issues.md (the issue log, the main statement of the audit's claims), run_config.R (all settings), run_pipeline.R (runs steps in order).
Read first: the summary table of issues.md and the issue(s) you are assigned, 01_pipeline/README.md, and the header block of every script you examine.
Vocabulary: "legacy" mode = the analysis as it was run for the draft; "corrected" mode = with the 99.9th percentile coverage filter applied. "DML" = differentially methylated locus (CpG with q < 0.01 and at least 20 percentage points difference). "shared" = loci in both the ploidy and the heat DML list (the draft calls them interaction-DML). "eligible" CpG = mean methylation across the 22 libraries between 10% and 90%. 22 libraries: 10 heat-treated (Zymo project zr3534) and 12 controls (zr3644), half diploid and half triploid in each.
A local text copy of the manuscript draft is at ${MANUSCRIPT}.
How to run R: Rscript is installed (R 4.3.3 with data.table, GenomicRanges, methylKit 1.28.0, topGO 2.54.0, goseq, readxl). Before running any pipeline code: export WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=2. The 22 Bismark coverage files are in /home/user/gannet_cache/cov_files (about 440 MB each; tab separated, no header: sequence, 0-based start, end, percent methylation, methylated count, unmethylated count).
HARD LIMITS: this machine has 4 cores and 15 GB of memory, and several other agents are working on it at the same time. Do NOT run 01_pipeline/04_cpg_coverage.R or 01_pipeline/05_differential_methylation.R (30 minutes, 10 GB). Never hold more than one coverage file in memory; prefer awk or grep for targeted extraction from them. Keep every R process under about 3 GB.
You are READ-ONLY on the repository: do not create, edit, or delete anything under ${REPO}, and do not run git commands that change its state (no add, commit, checkout, stash, reset, push). Other agents are adding files named 01_pipeline/11_*, 01_pipeline/12_*, 04_output/11_*, 04_output/12_* while you work: ignore those.
Write your scratch scripts and outputs only under your own scratch directory (given below); create it first.
Standards: report only what you checked yourself. Give numbers. Say plainly when something could not be determined. Do not soften a refutation and do not invent a problem to have something to report.
`

const CLAIM_SCHEMA = {
  type: 'object',
  properties: {
    claim_id: { type: 'string' },
    verdict: { type: 'string', enum: ['confirmed', 'partly_confirmed', 'refuted', 'undetermined'] },
    summary: { type: 'string', description: 'What you found, in plain sentences with the numbers.' },
    checks_run: { type: 'array', items: { type: 'string' }, description: 'Each independent check you performed and its result.' },
    errors_found: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          where: { type: 'string', description: 'File and line or section, e.g. issues.md I-39, or 01_pipeline/09_go_enrichment_legacy.R:120' },
          problem: { type: 'string' },
          correction: { type: 'string', description: 'The corrected statement, number, or code.' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['where', 'problem', 'correction', 'severity'],
      },
    },
    missing_caveats: { type: 'array', items: { type: 'string' }, description: 'Limits or alternative explanations the issue text should state and does not.' },
  },
  required: ['claim_id', 'verdict', 'summary', 'checks_run', 'errors_found', 'missing_caveats'],
}

const FINDINGS_SCHEMA = {
  type: 'object',
  properties: {
    area: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          summary: { type: 'string' },
          failure_scenario: { type: 'string', description: 'Concrete inputs or state that give a wrong output, a crash, or a misleading statement.' },
          evidence: { type: 'string', description: 'What you ran or read that shows it.' },
          suggested_fix: { type: 'string' },
        },
        required: ['file', 'severity', 'summary', 'failure_scenario', 'evidence', 'suggested_fix'],
      },
    },
    notes: { type: 'string', description: 'Anything that checked out and is worth recording, and anything you could not check.' },
  },
  required: ['area', 'findings', 'notes'],
}

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    area: { type: 'string' },
    verdicts: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          summary: { type: 'string' },
          verdict: { type: 'string', enum: ['confirmed', 'plausible', 'refuted'] },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          reason: { type: 'string', description: 'What you ran or read that decides it.' },
          fix: { type: 'string', description: 'The fix you recommend, as concretely as you can; empty if refuted.' },
        },
        required: ['file', 'summary', 'verdict', 'severity', 'reason', 'fix'],
      },
    },
  },
  required: ['area', 'verdicts'],
}

const CLAIMS = [
  {
    id: 'C1-universe',
    text: `Issue I-39. The gene universe of the legacy GO enrichment, 03_input_files/annotation/geneid2go-union1x-BOTH.tab.gz (27,469 transcripts), lacks 15,715 of the 43,184 annotated mRNAs. 15,666 of the missing ones form one unbroken run of accession numbers (XM_011455141.3 to XM_034461018.1). The share of annotated mRNAs present per linkage group runs from 18% to 96%. The loss is an artifact (a failed write in the legacy notebook 08_archive/code/4_WGBS_GOterm-Annotation-BOTH.ipynb), not a consequence of the intended filter (transcripts overlapping a CpG with a read in at least one library). Of the genes that hold a DML and have GO terms, the universe contains 73% (ploidy), 75% (heat), 78% (shared).
Try to refute it. In particular: (1) test the alternative explanation that the missing transcripts simply had no covered CpG: use 04_output/08_gene_annotation/tables/transcript_table.csv.gz and gene_table_legacy.csv.gz (n_cpg_tested counts CpGs with 10x in all libraries per gene), and for a direct 1x check take 40 genes from the missing block and look for covered CpGs inside them in one coverage file with awk. (2) Is "one unbroken run" true, and under which ordering of the IDs (byte order, locale order)? Would a different but equally natural ordering break it up? (3) Recompute the per-linkage-group shares (04_output/09_go_enrichment_legacy/tables/legacy_universe_by_chromosome.csv) and the 73/75/78% figures yourself. (4) Read the notebook (it is JSON; the saved cell outputs matter) and judge whether "failed write" is supported or only one of several possible causes; say what the evidence does and does not show.`,
  },
  {
    id: 'C2-reader',
    text: `Issue I-40. The legacy script 08_archive/code/5_WGBS_Functional-enrichement-analysis-BOTH.R read the three DML annotation tables with read.table() defaults; apostrophes in GO term names made it return 10,979 of 46,626 lines (ploidy), 81,109 of 106,621 (heat), 6,941 of 7,029 (shared). Transcripts kept: 623 of 2,348; 3,444 of 4,606; 379 of 379. Genes kept: 171 of 622; 949 of 1,237; 117 of 117. Combined with the incomplete universe, genes that entered the test: 168 of 622 (27%), 610 of 1,237 (49%), 91 of 117 (78%) (04_output/09_go_enrichment_legacy/tables/legacy_funnel.csv).
Try to refute it, using the legacy tables in 03_input_files/annotation/legacy/ (gz, ten columns, no header). (1) Is the apostrophe really the mechanism? Read with quote = "", with quote = "\\"" only, and with defaults, and compare. Does comment.char or fill play a part? (2) Would R 4.1 (2021) have behaved the same as R 4.3.3 for this call? Check the read.table and scan documentation and NEWS for relevant changes. (3) Recompute every number above independently of the pipeline code (do not source 02_R_functions/annotation.R). (4) The issue says which rows were lost is arbitrary. Examine which rows and genes survive: is there a pattern that matters biologically or statistically (for example whole genes lost versus single rows)? (5) Is the claim "the script's own comments record 10979, 81109, 6941" accurate? Quote the lines.`,
  },
  {
    id: 'C3-ontology-release',
    text: `Issues I-37 (update), I-28, and the header of 01_pipeline/09_go_enrichment_legacy.R. The nine archived topGO tables (08_archive/GOterms/BOTH-{PLOIDY,HEAT,COMBO}/{PLOIDY,HEAT,COMBO}-{BP,CC,MF}-FisherTestResults.csv) are reproduced exactly (every term, with identical annotated and significant counts) by step 09 with GO.db 3.13.0, and not with GO.db 3.14.0 or 3.18.0; hence the legacy enrichment ran under Bioconductor 3.13 (topGO 2.44.0), not with topGO 2.34.0 as the Methods say.
Try to refute it. (1) Compare 04_output/09_go_enrichment_legacy/tables/legacy_fisher_<set>_<ontology>.csv (sets ploidy, heat, shared) with the archive yourself. The archive CSVs were written unquoted and term names contain commas, so parse from both ends of each line. Compare p-values too wherever the archive value is parseable (it is printed to two or three significant digits, sometimes as "< 1e-30"). (2) Re-run the Fisher test for at least the BP ontology of one DML set under GO.db 3.14.0 (side library /home/user/gannet_cache/godb/lib314; put it first on .libPaths() before loading topGO, in a fresh R process) and under the default GO.db 3.18.0, to confirm that they do not reproduce the archive. Follow the logic of 01_pipeline/09_go_enrichment_legacy.R; you may source 02_R_functions for helpers. (3) Challenge the inference from "GO.db 3.13.0" to "Bioconductor 3.13, topGO 2.44.0": could topGO 2.34.0 have been used with that GO.db? What do the legacy scripts and any saved session files or comments in 08_archive say about R, Bioconductor, and package versions? (4) Check that the 2.34.0 statement is really in the manuscript Methods.`,
  },
  {
    id: 'C4a-eligible-theory',
    text: `Issue I-41 (read it in full) and 01_pipeline/10_go_enrichment_gene_level.R, 02_R_functions/go_tests.R, 02_R_functions/methylkit_helpers.R (dml_eligible, mean_percent_methylation). The argument: (a) a CpG can be a DML at the 20-point threshold only if its mean methylation across the 22 libraries lies between 10% and 90%; (b) so the fair background for "where are DML" and "which genes hold DML" is the eligible CpGs, and genes with an eligible CpG; (c) defining the background by the pooled mean does not bias the comparison (independent filtering, Bourgon et al. 2010); (d) with the eligible gene universe and the goseq Wallenius test with the number of eligible CpGs per gene as the bias variable, no GO term is enriched at Benjamini-Hochberg q < 0.05, and this is the right primary analysis.
You are a statistician refereeing this argument. Attack each link. (a) methylKit's meth.diff is a difference of coverage-weighted pooled group proportions (read methylKit's calculateDiffMeth and the pooling code), the groups are 11 and 11 libraries for ploidy and 10 and 12 for heat (check 03_input_files/metadata/wgbs_libraries.csv), and mean_percent_methylation() is an unweighted mean of per-library percentages. Is the 10% to 90% bound exact, approximate, or wrong? Quantify with the data: how many DML fall outside, and how far (04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz, dml_*_legacy.csv). (c) Is the overall mean independent of the test statistic under the null hypothesis for this logistic regression with a covariate, or can conditioning on it bias the comparison, for example through unequal group sizes, or through the systematic coverage and methylation differences between the two cohorts? (b) Is restricting the universe to genes with eligible CpGs over-conditioning that could hide real biology? What is the strongest counter-argument an author could make, and is it right? (d) Is the Wallenius approximation of goseq sound when the bias is this strong (the probability that a gene holds a DML rises from 5% to 69% with its number of eligible CpGs)? Is Benjamini-Hochberg across GO terms appropriate given the dependence between terms? Is "no term enriched" robust to go_min_genes (5), to the eligible range, and to using the smallest q rather than a count? Run numerical checks where they decide a point, using 04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz, gene_go_map.csv.gz and 04_output/10_go_enrichment_gene_level/tables/*.csv. Also assess the negative control 06_diagnostics/go_random_cpg_control.R: does it show what issues.md says it shows, and what would make it misleading? Finally, check each literature citation in I-41 and I-03 for accuracy as far as you can from your own knowledge, and flag any you are not sure of.`,
  },
  {
    id: 'C4b-eligible-reimplementation',
    text: `Issue I-41, final claim: against the genes that could hold a DML, and adjusted for the number of eligible CpGs per gene, no GO term is enriched at FDR 0.05 for the ploidy or the heat DML (step 10: 01_pipeline/10_go_enrichment_gene_level.R, results in 04_output/10_go_enrichment_gene_level/tables/).
Test it by INDEPENDENT RE-IMPLEMENTATION. Do not call go_gene_level_test() or goseq. Implement two different tests of GO enrichment among genes holding DML that account for the number of eligible CpGs per gene: (i) logistic regression per term: gene holds a DML ~ smooth function of log(number of eligible CpGs) + term membership, one-sided for enrichment, Benjamini-Hochberg across terms with at least 5 annotated genes (fit the null model once and use score tests, or fit per term for the terms with the smallest unadjusted p-values and bound the rest); (ii) a permutation test: permute the DML labels among the eligible CpGs keeping the number of DML fixed (at least 1,000 permutations, vectorized), recompute which genes hold at least one DML, and compare the observed number of DML genes per term with its permutation distribution; adjust across terms. Legacy mode only; ploidy and heat; ontologies BP, CC, MF.
Inputs: 04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz (the tested CpGs with per-library counts, columns numCs<j>, numTs<j>), dml_ploidy_legacy.csv, dml_heat_legacy.csv (column united_row is the row of the united table), 04_output/08_gene_annotation/tables/gene_table_legacy.csv.gz and gene_go_map.csv.gz (direct annotations only; you need ancestors too: build a topGOdata object as go_gene_level_object() in 02_R_functions/go_tests.R does and use genesInTerm(), or use GO.db). Helper functions you may source from 02_R_functions: mean_percent_methylation, dml_eligible, cpg_ranges, overlap_pairs, read_genes. Protein-coding genes only; a DML belongs to every gene whose body contains it.
Report, per DML set and ontology: terms tested, terms with q < 0.05 under each of your two methods, the smallest q, and the five strongest terms with their counts. State whether your result agrees with step 10 (zero terms). If either of your methods finds enriched terms, list them all and investigate why the Wallenius test did not.`,
  },
  {
    id: 'C5-feature-odds',
    text: `Issue I-41 (location table) and the I-14 update. Odds ratios of a DML lying in a genome feature, against all tested CpGs and against eligible CpGs, legacy mode: gene body 7.1 and 2.5 (ploidy), 6.1 and 2.2 (heat); intron 1.9 and 0.80, 2.0 and 0.87; CDS 1.5 and 1.7, 1.4 and 1.6; transposable element 0.77 and 0.54, 0.83 and 0.57; intergenic 0.11 and 0.30, 0.14 and 0.37. Also: step 07 reproduces every DML location percentage of Results 3.3 of the manuscript through the column legacy_share_pct of 04_output/07_genome_features/tables/feature_overlap_by_track_legacy.csv.
Recompute independently: write your own overlap code with GenomicRanges, reading the gz tracks in 03_input_files/genome_features/ directly (GFF is 1-based closed; BED and the CpG positions are 0-based starts; a CpG covers two bases). Do not call 01_pipeline/07_genome_features.R or 02_R_functions/genomic_overlaps.R, though read them afterwards to compare. Inputs: 04_output/05_differential_methylation/tables/united_cpg_legacy.csv.gz, dml_ploidy_legacy.csv, dml_heat_legacy.csv. Then: (1) how is the 2x2 table of step 07's Fisher test set up (is a DML compared with the other CpGs of the background, and are DML that are not themselves eligible handled sensibly)? Is that the right comparison? (2) Do the conclusions (gene-body enrichment a third as large; introns not enriched; CDS still enriched) hold if the eligible range is 15% to 85% or 5% to 95%? (3) Compare the percentages quoted in Results 3.3 of the manuscript (search the text copy for "Relative to this background" and the sentences around it) with the regenerated legacy_share_pct values, one by one, and list every mismatch. (4) Are neighbouring CpGs treated as independent in these tests, and does the issue text say so?`,
  },
  {
    id: 'C6-coverage-filter-locus',
    text: `Issues I-05, I-09 and I-38. (I-05) methylKit::filterByCoverage takes hi.perc; the legacy call passed high.perc, so no upper filter ran; with the filter there are 454,555 tested CpGs and 952, 2,273 and 159 DML (ploidy, heat, shared). (I-09) All 54 tested CpGs in NC_047563.1:39,633,323-39,635,520 are removed by the 99.9th percentile filter, because each is at or above the cutoff in at least one library; the excess coverage suggests a collapsed repeat or a copy-number difference. (I-38) methRead() defaults to header = TRUE, so line 1 of every coverage file was dropped by the legacy import.
Verify. (1) Read methylKit's own code and argument lists (args(), getMethod(), showMethods) for filterByCoverage and methRead; confirm what happens to an unknown argument such as high.perc. (2) For the locus: extract the region from each of the 22 coverage files with awk, one file at a time (/home/user/gannet_cache/cov_files/<library>_R1.CpG_report.merged_CpG_evidence.cov), and compare coverage there with each library's cutoff in 04_output/04_cpg_coverage/tables/library_coverage.csv (column upper_cutoff). Is every one of the 54 CpGs at or above the cutoff in at least one library? In how many libraries each? Is the excess present in both cohorts and both ploidies, or only in some libraries (this bears on the copy-number reading)? What is the ratio of coverage in the region to the library median? What does the flanking 5 kb look like? Is the interpretation in I-09 justified, overstated, or one of several? (3) Does the pipeline's filter (02_R_functions/methylkit_helpers.R, read_filtered_counts and upper_coverage_cutoff in coverage_io.R: remove a CpG from a library when its coverage is at or above the 99.9th percentile computed among CpGs with at least 2 reads, then keep CpGs present in all 22 libraries) match methylKit's filterByCoverage(lo.count = 10, hi.perc = 99.9) followed by unite() exactly? Read methylKit's code for the comparison operators and the quantile type. (4) Check the DML counts by mode (04_output/05_differential_methylation/tables/dml_counts_by_mode.csv) against the DML files, and the table added to I-05 in issues.md against 04_output/07_genome_features/tables/dml_region_summary_<mode>.csv and 04_output/08_gene_annotation/tables/annotation_summary.csv.`,
  },
  {
    id: 'C7-legacy-go-semantics',
    text: `Issues I-20 (update), I-21 (update), I-02 (update), and step 09 (01_pipeline/09_go_enrichment_legacy.R with outputs in 04_output/09_go_enrichment_legacy/tables/). Claims: the table in I-20 states what each of the four panels of Figure 5 counts; I-21 states how "DML within genes with enriched GO terms" (Figure 9, Table 1) was computed; the methylation differences of the interaction rows are those of the heat contrast.
Verify against the legacy code that produced each number: 08_archive/code/5_WGBS_Functional-enrichement-analysis-BOTH.R and 08_archive/code/gene-plots.Rmd (the values of Figures 5 and 9 are typed there), and 08_archive/code/2_WGBS_Methylkit-all-gigas.R for how the interaction table was built. Is each definition in the I-20 table exactly right? Is step 09 a faithful replica of the legacy logic (compare line by line where it matters: the left joins, the unique() calls, >= 0 and <= 0 for hyper and hypo)? Then compare 04_output/09_go_enrichment_legacy/tables/table_1_legacy.csv with Table 1 of the manuscript (text copy; search "Table 1"), row by row, and check every gene-level statement in Results 3.5 against it: "kinesin-like protein unc-104 (+42%)", "ataxin-2-like protein (-21%)", "CAD protein (-31%)", "nesprin-1 (+21%)", "dihydropyrimidinase (+31%)", "F-box DNA helicase 1 (-21%)", "neurexin-1 (-32%)", "MAP kinase-activating death domain protein (+22%)", the percentages "CC 49%, MF 31%, BP 20%" and "MF 68%, CC 28%, BP 4%", and the interaction counts "2 hypermethylated and 1 hypomethylated ... 8 and 8 ... 5 and 3". List every mismatch between the manuscript and the regenerated values, and every place where the manuscript's wording describes something other than what the number is.`,
  },
  {
    id: 'C8-counts-and-correlations',
    text: `Issue I-16 (update) and step 04 (01_pipeline/04_cpg_coverage.R, outputs in 04_output/04_cpg_coverage/tables/). Claims: the reference holds 13,245,836 CpG sites; CpGs with a read in at least one library number 12,918,252; the manuscript's 12,794,074 cannot be traced; pairwise Pearson correlations are 0.79 to 0.90 at 1x and 0.84 to 0.94 at 10x, against "0.87 to 0.93" in Results 3.1; within-cohort and between-cohort median correlations are 0.862 and 0.856 at 1x.
(1) Try to find where 12,794,074 and "0.87 to 0.93" came from. Search 08_archive thoroughly: R scripts, Rmd, notebooks (JSON, including saved outputs), md and txt files, csv, and every sheet of every xlsx workbook (readxl). Search for the digit strings with and without separators and for nearby values. Test cheap hypotheses for the correlation range with the saved matrices 04_output/04_cpg_coverage/tables/library_correlation_{1x,10x}_{legacy,corrected}.csv: within one cohort only, within one ploidy only, a subset of libraries, Spearman cannot be tested from these so say so. (2) Check the arithmetic in the I-16 update (96.6% and 3.6%). (3) Review how step 04 computes the correlation and the PCA (percent methylation from counts scaled like methylKit's normalizeCoverage; Pearson; PCA on CpGs with SD above the median, centred and scaled) against what the legacy script 08_archive/code/1_coverage_analysis.R did (getCorrelation, PCASamples and their arguments); are they the same computation? (4) Confirm 13,245,836 by counting CG dinucleotides in the reference with a method of your own (the reference is /home/user/gannet_cache/genome/cgigas_uk_roslin_v1_genomic-mito.fa, 656 MB; use awk or a streaming count, not R, and mind CG pairs that span a line break) and say whether the count should include lower-case (soft-masked) sequence and both strands.`,
  },
]

const REVIEWS = [
  {
    id: 'R1-coverage-dml',
    text: `Area: 01_pipeline/04_cpg_coverage.R, 01_pipeline/05_differential_methylation.R, 02_R_functions/coverage_io.R, 02_R_functions/methylkit_helpers.R, 02_R_functions/sample_sheet.R, 03_input_files/metadata/wgbs_libraries.csv.
Look for correctness bugs: wrong joins or index alignment (the code matches CpGs by a numeric key and by row position), library order assumptions, off-by-one coordinate errors, the methylKit-style scaling and rounding, the legacy first-line handling, the treatment vectors and covariates passed to calculateDiffMeth, the definition of hyper and hypo, the shared set and its classes, the threshold sweep, and whether "corrected" mode really differs from "legacy" only in what the header says. Compare with the legacy scripts 08_archive/code/1_coverage_analysis.R and 08_archive/code/2_WGBS_Methylkit-all-gigas.R: does the new code compute what the legacy code computed? You cannot run steps 04 and 05, but you can test their functions on small synthetic inputs or on the first 200,000 lines of two coverage files (head), and you can check their committed outputs in 04_output/04_cpg_coverage and 04_output/05_differential_methylation for internal consistency (for example that dml files agree with diffmeth files at the stated thresholds, that united_row indexes the united table, that counts by threshold are monotone).`,
  },
  {
    id: 'R2-figures-features',
    text: `Area: 01_pipeline/06_dml_figures.R, 01_pipeline/07_genome_features.R, 02_R_functions/genomic_overlaps.R, 02_R_functions/plot_theme.R, the feature tracks in 03_input_files/genome_features/ and their README.
Look for correctness bugs: coordinate conventions (GFF 1-based closed, BED 0-based half-open, CpG two bases wide), strand handling, the exclusive class assignment and its priority, the legacy_share_pct formula, the Fisher and chi-squared tests (table construction, which CpGs are in the comparison row), the region grouping (neighbours within 1 kb, across sequence boundaries?), counts per linkage group, the PCA (same filter as methylKit::PCASamples?), heat map inputs, and the labels and legends of Figures 2, 3 and 4 (open the PNGs in 04_output/06_dml_figures/figures and 04_output/07_genome_features/figures and check that what is drawn matches the tables). You may run both steps (each takes under a minute and about 2 GB) but only into a copy: copy the repository working tree to your scratch directory with rsync, excluding .git and 04_output/05_differential_methylation/cache, and run there. Check that rerunning reproduces the committed CSV outputs byte for byte.`,
  },
  {
    id: 'R3-annotation-legacy-go',
    text: `Area: 01_pipeline/08_gene_annotation.R, 01_pipeline/09_go_enrichment_legacy.R, 02_R_functions/annotation.R, 02_R_functions/legacy_go_db.R, 03_input_files/annotation/ (README, build_legacy_go_terms.R, the maps).
Look for correctness bugs: GFF attribute parsing (genes or transcripts whose attributes break the regular expressions; URL decoding of products), genes sharing an ID, the transcript name repair (true_transcript_ids), data.table scoping errors (a function argument or loop variable with the same name as a column is a known hazard in this code base), merges that silently multiply or drop rows (allow.cartesian), the per-gene summaries, the n_cpg_tested and n_cpg_eligible counts, whether the gene-level GO map really is the union over a gene's mRNAs, and in step 09 whether updateGenes() reuse can leak state between DML sets, whether genes_significant is computed from the right transcript set, and whether the side-library mechanism can load the wrong GO.db (for example when GO.db is already loaded, or when the cache is missing). You may run both steps in a copy of the working tree made with rsync into your scratch directory (exclude .git and 04_output/05_differential_methylation/cache); step 08 takes 20 seconds and step 09 three minutes. Check that rerunning reproduces the committed CSV and gz outputs byte for byte.`,
  },
  {
    id: 'R4-gene-level-go',
    text: `Area: 01_pipeline/10_go_enrichment_gene_level.R, 02_R_functions/go_tests.R, 06_diagnostics/go_random_cpg_control.R, and their outputs (04_output/10_go_enrichment_gene_level/tables, 06_diagnostics/go_random_cpg_control*.csv).
Look for correctness bugs in the statistics as implemented: is the universe handed to topGO the one the header describes; does updateGenes() on a shared object give the same result as a fresh object (test it for one case); are the genes passed to goseq::nullp() and goseq::goseq() aligned with the bias vector; does use_genes_without_cat = FALSE do what the code assumes; is the Wallenius p-value one-sided for over-representation; are Benjamini-Hochberg q-values computed over the right set of terms; is "tested" (at least 5 annotated genes) applied before adjustment; is the handling of DML genes without an eligible CpG (pmax with any_dml) sound or does it leak outcome information into the universe; in the random control, are the random CpGs mapped to genes the same way real DML are, and is the comparison with the real DML like for like (same universe, same genes filter). Read goseq's source for nullp and goseq to confirm. You may run step 10 in a copy of the working tree made with rsync into your scratch directory (exclude .git and 04_output/05_differential_methylation/cache); it takes three minutes. Do not rerun the random control in full (10 minutes); test its functions on one draw.`,
  },
  {
    id: 'R5-clean-checkout',
    text: `Area: whether the repository runs as its documentation says, from a clean checkout. Make a clone into your scratch directory: git clone --quiet --branch reproducibility-audit ${REPO} <scratch>/clone (cloning is a read of the source repository and is allowed). In the clone, with WGBS_CACHE_DIR=/home/user/gannet_cache and WGBS_CORES=2 exported, run: Rscript run_pipeline.R --list ; then Rscript run_pipeline.R --steps 2:3 ; then Rscript run_pipeline.R --steps 6:10 (about 8 minutes; do NOT run steps 4 and 5; step 1 may be run, it only verifies checksums for 2 minutes). Then, inside the clone, use git status --short and git diff --stat to see which tracked outputs changed. PDF files and session_info.txt carry time stamps and are expected to change; every CSV, TSV and gz table should be byte-identical. Report each tracked table that changed, with a diagnosis (non-determinism, dependence on run order, locale, package behaviour). Then run Rscript 06_diagnostics/check_legacy_outputs.R in the clone and confirm that 329 of 331 comparisons match and that the two differences are the documented ones.
Also review run_pipeline.R, run_config.R and 02_R_functions/{paths,packages,logging,fetch_gannet,legacy_go_db,plot_theme}.R for: failure modes of argument parsing; what happens when the cache is absent or a download is interrupted; portability to Windows and macOS (forked workers, gz reading, file.rename across devices, sort order and locale, fonts in figures); hidden dependence on this machine (absolute paths, environment variables, packages that 07_documentation/cloud_environment_setup.sh does not install: compare library() and loadNamespace() and pkg:: uses across 01_pipeline and 02_R_functions and 06_diagnostics with that script); unseeded randomness. Check every claim about running the pipeline in README.md and 01_pipeline/README.md against what you observed (run times, memory, which steps need which downloads).`,
  },
  {
    id: 'R6-documentation',
    text: `Area: the accuracy of the written record: README.md, issues.md (the dated updates of 2026-11-19 and issues I-38 to I-42 above all), 01_pipeline/README.md, 02_R_functions/README.md, 03_input_files/README.md and 03_input_files/annotation/README.md, 04_output/README.md, 06_diagnostics/README.md, 07_documentation/reproducibility_plan.md, and the header blocks of the scripts in 01_pipeline.
Check every number in those texts against the output tables it comes from (04_output/*/tables, 06_diagnostics/*.csv); every file name and path against the file system; every function listed in 02_R_functions/README.md against the code; every relative link; every statement of the form "step N does X" against the script. Check internal consistency between documents (for example the same quantity given as two different numbers, issue statuses in the summary table against the issue text, step numbers that changed when the pipeline went from twelve to thirteen steps, references to "step 09" or "phase 4" that now point at the wrong thing). Check the writing rules of this project: no em dash character anywhere in new text (the dash U+2014; search for it); the treatment is called "heat" or "heat-treated", never "desiccation", except inside 08_archive, in quoted legacy names, and in the repository name; claims carry a source. Use grep and small R or awk checks; cite file and line for every finding. Severity high = a wrong number or statement that would mislead an author acting on it; medium = inconsistency or stale text; low = wording.`,
  },
]

function claimPrompt(c) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/${c.id}
YOUR TASK: adversarial verification of one claim of the audit. Default to skepticism: your job is to find out whether the claim is wrong, overstated, or missing a caveat, by checking it against the data and code yourself. A claim you could not break after real attempts is "confirmed"; one that holds in substance but has wrong details is "partly_confirmed" with the details in errors_found; one whose central statement is wrong is "refuted".

THE CLAIM
${c.text}

Return claim_id "${c.id}". In checks_run list each check with its numeric result. In errors_found give every wrong number, wrong statement, or code error you found, with the correction. In missing_caveats list limits or alternative explanations that the issue text should state.`
}

function reviewPrompt(r) {
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/${r.id}
YOUR TASK: code and documentation review of one area. Find real defects: things that give a wrong number, a wrong figure, a crash on a plausible input, a result that changes from run to run, or a statement that is false. Do not report style preferences. For every finding give the concrete failure scenario and the evidence (what you ran or read). If the area is sound, say so and return few or no findings; an empty list is an acceptable answer.

${r.text}

Return area "${r.id}".`
}

function verifyPrompt(review, r) {
  const list = (review && review.findings) ? review.findings : []
  return `${CONTEXT}
Your scratch directory: ${SCRATCH}/${r.id}-verify
YOUR TASK: an independent reviewer examined one area of this repository and reported the findings below. Re-examine each one yourself and decide: "confirmed" (you reproduced the defect or read the code and it is certainly there), "plausible" (likely real but you could not demonstrate it), or "refuted" (the code or data show it is not a defect; say why). Default to refuted when you cannot find the defect after a real attempt. For confirmed and plausible findings give the most concrete fix you can, as code or exact replacement text. You may run code under the same limits as above, in your scratch directory.

Area reviewed: ${r.id}
Scope given to the reviewer: ${r.text}

FINDINGS TO CHECK (JSON):
${JSON.stringify(list, null, 1)}

Reviewer notes: ${review && review.notes ? review.notes : '(none)'}

Return area "${r.id}" and one verdict per finding, in the same order.`
}

phase('Verify claims')
const claimsPromise = parallel(CLAIMS.map(c => () =>
  agent(claimPrompt(c), { label: `claim:${c.id}`, phase: 'Verify claims', schema: CLAIM_SCHEMA })))

const reviewsPromise = pipeline(
  REVIEWS,
  r => agent(reviewPrompt(r), { label: `review:${r.id}`, phase: 'Review code', schema: FINDINGS_SCHEMA }),
  (review, r) => {
    if (!review || !review.findings || review.findings.length === 0) return { review, verdicts: { area: r.id, verdicts: [] } }
    return agent(verifyPrompt(review, r), { label: `verify:${r.id}`, phase: 'Check findings', schema: VERDICT_SCHEMA })
      .then(v => ({ review, verdicts: v }))
  }
)

const claims = (await claimsPromise)
const reviews = (await reviewsPromise)

const claimResults = CLAIMS.map((c, i) => claims[i] || { claim_id: c.id, verdict: 'undetermined', summary: 'agent returned nothing', checks_run: [], errors_found: [], missing_caveats: [] })
const reviewResults = REVIEWS.map((r, i) => reviews[i] || { review: { area: r.id, findings: [], notes: 'agent returned nothing' }, verdicts: { area: r.id, verdicts: [] } })
log(`claims: ${claimResults.map(c => c.claim_id + '=' + c.verdict).join(', ')}`)
log(`review findings: ${reviewResults.map(r => (r.review && r.review.area) + ':' + ((r.review && r.review.findings) ? r.review.findings.length : 0)).join(', ')}`)
return { claims: claimResults, reviews: reviewResults }

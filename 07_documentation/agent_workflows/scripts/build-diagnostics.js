export const meta = {
  name: 'build-diagnostics',
  description: 'Turn five one-off verification computations into diagnostics scripts that run from the repository, each built, independently reviewed, and fixed',
  phases: [
    { title: 'Build', detail: 'one agent writes each script under 06_diagnostics and runs it' },
    { title: 'Review', detail: 'a second agent reruns it, recomputes key numbers independently, reads the code' },
    { title: 'Fix', detail: 'confirmed findings fixed; final run' },
  ],
}

const REPO = '/home/user/WGBS-gigas-ploidy-desiccation'
const SCRATCH = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad/wf'
const VER = REPO + '/07_documentation/agent_workflows/verification_scripts'

const CONTEXT = `
CONTEXT (read all of it before starting)
Repository: ${REPO} (git branch reproducibility-audit). It holds an R-only pipeline that re-derives the results of a manuscript draft on DNA methylation (whole-genome bisulfite sequencing) in diploid and triploid Pacific oysters after a heat exposure, and an audit of that draft (issues.md, 45 issues).
Design: 22 libraries. 10 are heat-treated animals (5 diploid, 5 triploid; "heat" cohort) and 12 are controls (6 diploid, 6 triploid; "control" cohort). The two cohorts come from different experiments and sequencing batches, so heat is confounded with cohort (issue I-01). The sample sheet is 03_input_files/metadata/wgbs_libraries.csv, read with read_libraries(); its row order is the library order of every matrix in the pipeline.
What the draft did (step 05): CpGs with at least 10 reads in all 22 libraries (480,579; mode "legacy"); methylKit::calculateDiffMeth with defaults, once for ploidy with heat as covariate and once for heat with ploidy as covariate; a DML is a CpG with q < 0.01 and an absolute difference above 20 percentage points (1,038 ploidy DML, 2,351 heat DML, 200 loci in both lists). Mode "corrected" applies the 99.9th percentile coverage filter that the draft's Methods describe (454,555 CpGs; 952, 2,273, 159).
Outputs you start from: 04_output/05_differential_methylation/tables/ holds united_cpg_<mode>.csv.gz (the tested CpGs: chr, start, end, strand, then coverage<j>, numCs<j>, numTs<j> for j = 1..22, NORMALIZED counts, in read_libraries() order), diffmeth_<contrast>_<mode>.csv.gz (every tested CpG with pvalue, qvalue, meth.diff, united_row), dml_{ploidy,heat,shared}_<mode>.csv, library_filtering_<mode>.csv. 04_output/04_cpg_coverage/tables/library_coverage.csv has per-library medians and cutoffs. 04_output/07_genome_features/tables/ has the feature overlap and association tables; 04_output/08_gene_annotation/tables/ the gene tables; 04_output/09_go_enrichment_legacy/tables/ the legacy enrichment.
The 22 raw coverage files (Bismark merged CpG evidence, no header, 6 columns: chr, start, end, percent, count methylated, count unmethylated; about 450 MB and 10 million lines each) are in /home/user/gannet_cache/cov_files/ and are reached in code through fetch_gannet_file() or gannet_cache_dir(config, "cov_files") with the file names of the cov_file column of the sample sheet (see how 01_pipeline/04_cpg_coverage.R does it, and read_cov() in 02_R_functions/coverage_io.R).
Why this task exists: independent agents re-derived the audit's statements and found more. Their numbers now stand in issues.md in paragraphs headed "Verification", but the code behind them ran in a scratch directory and cannot be rerun from the repository. Your script makes one group of those numbers reproducible. The verifying agents' own scripts and small result tables are in ${VER}/<ID>/ (their full scratch directories, if still present, are ${SCRATCH}/<ID>/). Use them as a starting point and as a cross-check, not as something to copy blindly: they contain absolute paths and exploratory dead ends, and an agent can be wrong.
Read first: the issue named in your task in ${REPO}/issues.md (whole entry); 01_pipeline/README.md section "How a step is written"; run_config.R; the helper files in 02_R_functions/ (load_packages, input_path, output_path, read_libraries, log_msg, write_table, gannet_cache_dir, read_cov, upper_coverage_cutoff, read_filtered_counts, scale_counts_like_methylkit, mean_percent_methylation, dml_eligible, read_feature_track, cpg_ranges, overlaps_track, overlap_pairs, read_dml_annotation, use_legacy_go_db); and, as models of a diagnostics script, 06_diagnostics/table1_group_means.R and 06_diagnostics/cohort_separation.R, and 06_diagnostics/README.md.
How to run: Rscript is installed (R 4.3.3; data.table, methylKit 1.28.0, GenomicRanges, Biostrings, topGO, goseq, GO.db, ggplot2, patchwork, curl). Export WGBS_CACHE_DIR=/home/user/gannet_cache WGBS_CORES=1 before running anything.
MACHINE LIMITS (strict; the machine has 4 cores and 15 GB shared with other agents and a long model-fitting job): one R process at a time, single-threaded (data.table::setDTthreads(1) while you develop; no parallel, no mclapply, WGBS_CORES=1), under 3 GB of memory. Read the raw coverage files one at a time, keep only the rows or the summary you need, and free the rest before the next file. Do not run any script of 01_pipeline. Do not run 06_diagnostics/dml_sensitivity.R (it is running; do not touch its folder 06_diagnostics/dml_sensitivity/). If a run would take more than about 25 minutes, redesign it or reduce the number of permutations and say so.
RULES FOR WHAT YOU WRITE: R only. One script, 06_diagnostics/<name>.R, with a header block in the style of the two model scripts (Purpose, Run with its approximate run time, Reads, Writes, Notes). Outputs are CSV files next to the script, named <name>*.csv (small: no table above about 2 MB; no per-CpG dumps of the whole genome). Paths through here::here and the helpers; no setwd, no absolute paths, nothing that points at a scratch directory. Settings from the config list of run_config.R; do not edit run_config.R, the helpers, any 01_pipeline script, any README, or issues.md (return proposed text instead). No number typed in as a result: every number the script reports is computed by it; expected values may appear only in stopifnot-style checks of inputs or in comments. Deterministic: set.seed before anything random and say how many draws. Fail loudly on unexpected input. Plain declarative English in comments, column names and messages; no em dash character (U+2014); the treatment is "heat", never "desiccation". Cite sources for method choices (authors, year, journal, volume, pages or DOI) only if you have checked what they say. Known hazard: inside a data.table call, a function argument or loop variable with the same name as a column refers to the column; name things so that this cannot happen.
Do not run git commands that change state. Scratch work goes under ${SCRATCH}/<your own new subdirectory>.
A departure from a target number is not a failure. If your script, checked carefully, gives a different number from the one in issues.md, report both and explain which is right and why. Do not tune the method to hit a target.
`

const ITEMS = [
  {
    key: 'feature_odds_matched',
    issue: 'I-41 (location of DML), and its corrections',
    basis: 'C5-feature-odds and C4a-eligible-theory',
    spec: `Write 06_diagnostics/feature_odds_matched.R. Question: are DML more likely than comparable CpGs to lie in gene bodies, coding sequence, introns, transposable elements, or intergenic sequence? "Comparable" is the point. Step 07 (01_pipeline/07_genome_features.R; read it) compares DML with all tested CpGs and with "eligible" CpGs (mean methylation across libraries between 10% and 90%). Verification showed the eligible background is still unfair: inside it the chance of exceeding 20 points by chance depends on mean methylation and on how much animals of the same group differ.
The script computes, for mode legacy and mode corrected, for the ploidy and heat DML lists, and for the tracks gene body, CDS, intron, exon UTR, lncRNA, transposable element, intergenic (same track definitions and overlap convention as step 07; reuse its helpers):
 (a) the odds ratio against all tested CpGs and against eligible CpGs, with Fisher 95% intervals; these must equal step 07's feature_association_tests_<mode>.csv, and the script checks that they do;
 (b) the Mantel-Haenszel odds ratio with 95% interval against eligible CpGs stratified by 5-point bins of mean methylation;
 (c) the same, stratified by methylation bin and decile of within-group standard deviation, where the within-group SD of a CpG is the pooled standard deviation of per-library percent methylation around the means of the four groups (ploidy by cohort);
 (d) sensitivity of (a, eligible) to the range: 5% to 95% and 15% to 85%;
 (e) the transposable-element odds ratio against eligible CpGs separately in the top quartile of mean coverage and in the rest;
 (f) within genes, the odds ratio of CDS against intron;
 (g) a block version of the interval for (a, eligible): the same odds ratios with intervals from a bootstrap over 100 kb blocks (or a jackknife; say which and why), because neighbouring CpGs are not independent;
 (h) for ploidy only, a label-permutation reference: shuffle the ploidy labels within cohort (5 and 5 in the heat cohort, 6 and 6 in the control cohort), recompute the read-pooled difference triploid minus diploid at every tested CpG from the normalized counts (methylKit's definition of meth.diff: pooled methylated reads over pooled reads per group), call the CpGs whose absolute difference exceeds config$dml_difference, and compute the odds ratios of (a, eligible) for that list. 1,000 permutations with a fixed seed. Report per track the median, the 2.5th and 97.5th percentiles, and the share of permutations whose odds ratio is at least the real one; also the number of CpGs called per permutation (median, 5th and 95th percentiles, share of permutations reaching the real number of ploidy DML), and the share of eligible CpGs called by mean-methylation bin, and inside against outside gene bodies. State plainly in the header that the permutation applies the difference rule only, not the q-value rule, and give the evidence that the q-value rule is not binding in the real data (count it from diffmeth_<contrast>_<mode>.csv.gz). For heat no such permutation exists, because heat coincides with cohort; say so in the output, do not invent one.
Also write the median within-group SD of DML and of other eligible CpGs.
Targets from the verification (legacy mode; ploidy then heat): eligible gene body 2.52 and 2.20, CDS 1.68 and 1.57, intron 0.80 and 0.87, transposable element 0.54 and 0.57, intergenic 0.30 and 0.37; within methylation bins gene body 2.0 and 1.7; within bin and SD decile gene body 1.13 (0.76 to 1.70) and 1.01 (0.76 to 1.34), CDS 1.25 and 1.17, intron 0.81 and 0.87, transposable element 1.05 and 1.02, intergenic 0.80 and 0.92; range 15-85 gene body 2.10 and 1.84, range 5-95 3.59 and 3.14; transposable element 0.20 and 0.29 in the top coverage quartile, 0.81 and 0.83 in the rest; CDS against intron within genes 1.60 and 1.47; permutation: gene body median 2.68 (1.44 to 4.53), the real 2.52 reached in 58% of permutations, CDS 1.58, transposable element 0.44, intergenic 0.30; CpGs called per permutation median 751 to 762 (5th to 95th percentile about 510 to 1,320), 12% to 13% of permutations reaching 1,038; share called 0.30% at mean 10-20% and 2.70% at 50-60%; 2.55 times higher inside gene bodies; median within-group SD 22 points for DML and 15 for other eligible CpGs. Stratification details of the verifying agents may differ from yours; where a number differs, say why.
Outputs: feature_odds_matched.csv (long format: mode, dml_set, track, comparison, odds_ratio, ci_low, ci_high, n strata used, notes), feature_odds_label_permutation.csv (per track summary), feature_odds_label_permutation_counts.csv (CpGs called: summary and by methylation bin), all under 06_diagnostics/.`,
  },
  {
    key: 'locus_nc047563',
    issue: 'I-09 (the locus on NC_047563.1) and I-05',
    basis: 'C6-coverage-filter-locus (locus.R, locus_meth.R, loo.R, flank.R, kmers.R, locus_per_library.csv, locus_per_cpg.csv)',
    spec: `Write 06_diagnostics/locus_nc047563.R. One region, NC_047563.1:39,633,323-39,635,520, holds 53 of the 1,038 ploidy DML, 39 of the 2,351 heat DML and 39 of the 200 loci in both lists (legacy mode), and all of the draft's "lncRNA" signal. Every one of its 54 tested CpGs is removed by the 99.9th percentile coverage filter. The script documents what the locus is, from the raw coverage files and the step 05 tables:
 (a) the extent of the block of excess coverage: per-library raw coverage of every CpG between about 39,625,000 and 39,645,000, expressed relative to the library median coverage (library_coverage.csv) and to the library's upper cutoff; determine the block boundaries by a stated rule and report them, the number of CpG sites and of tested CpGs inside, and the coverage of the flanks;
 (b) per library, for the tested CpGs of the DML region: median raw coverage, its ratio to the library median and to the cutoff, how many of the CpGs are at or above the cutoff, and the read-pooled percent methylation; which libraries carry the excess (state the rule), and their distribution over ploidy and cohort with Fisher exact tests;
 (c) leave-one-library-out: for each library in turn, the number of CpGs in the region whose read-pooled difference (normalized counts of united_cpg_legacy.csv.gz) still exceeds config$dml_difference for the ploidy contrast and for the heat contrast when that library is left out; state that this applies the difference rule only and count from diffmeth how many CpGs in the region fail the q-value rule in the full data;
 (d) animal-level methylation of the region by ploidy: medians, Wilcoxon rank-sum test overall, within each cohort, and overall without the most extreme library;
 (e) the genes overlapping the block (gene track of 03_input_files/genome_features) and whether the block overlaps the repeat track;
 (f) optional, only if cheap and within the memory limit: GC content of the block and flanks and the share of 30-mers of the block that occur once in the reference (the reference FASTA is in the cache under genome/; reading the whole genome into memory may exceed the limit, in which case leave this out and say so).
Targets from the verification: block NC_047563.1:39,633,065-39,636,662 (3,599 bp, 82 CpG sites, 56 tested CpGs); flanks at 0.56 to 1.72 times the library median; 13 of 22 libraries carry the excess (median coverage 10.6 to 35.5 times the library median; 4 of 10 heat, 9 of 12 control; 5 of 11 diploid, 8 of 11 triploid; Fisher p 0.39 for ploidy and 0.19 for cohort), the other 9 at 2.0 to 5.9 times with no CpG at the cutoff; each of the 54 CpGs at or above the cutoff in 9 to 13 libraries; without zr3534_4 the region has 1 ploidy DML instead of 53 and 0 heat DML instead of 39; zr3534_4 at 34 times its library median, 84% methylated, 39% of the heat group's reads; methylation median 34.1% in diploids against 3.2% in triploids (Wilcoxon p 0.007; 0.012 without zr3534_4; heat cohort 0.008; control cohort 0.09); gene LOC105341636 (lncRNA), not in the repeat track; 42 of 60 sampled 30-mers unique; GC 42.5% against about 32%.
Outputs under 06_diagnostics/: locus_nc047563_per_library.csv, locus_nc047563_leave_one_out.csv, locus_nc047563_summary.csv (one row per fact: quantity, value, how defined), and locus_nc047563_coverage_profile.csv (per CpG of the window: position, tested or not, median relative coverage among the libraries that carry the excess and among the others).`,
  },
  {
    key: 'legacy_reader_mechanism',
    issue: 'I-40 (the file-reading fault) with I-39',
    basis: 'C2-reader',
    spec: `Write 06_diagnostics/legacy_reader_mechanism.R. The legacy enrichment script read three DML annotation tables with read.table() at default quoting, and GO term names contain apostrophes; most of the ploidy table and a quarter of the heat table were skipped. Step 09 (01_pipeline/09_go_enrichment_legacy.R; read it, and read_dml_annotation() in 02_R_functions/annotation.R) replicates the faulty read and writes legacy_reader_loss.csv and legacy_funnel.csv. This script documents the mechanism and its consequence, for the three tables in 03_input_files/annotation/legacy/ (ploidy, heat, shared):
 (a) lines, transcripts and genes in the file, returned by the default read, and returned with quote = "" and with quote = "\\"" only; any warnings R gives for each read (capture and report them);
 (b) the number of lines with an odd number of apostrophes outside a double-quoted field, the number of merged rows they produce, and the proof by intervention: after replacing the apostrophes on those lines only, the default read returns every line;
 (c) what is lost: transcripts and genes that vanish entirely, and the share of skipped lines that belong to them;
 (d) where the survivors are: for each table, the longest run of surviving transcripts in accession order, and the distribution of surviving and lost transcripts and genes over the ten linkage groups (join to 04_output/08_gene_annotation/tables/transcript_table.csv.gz); for ploidy, the share of tested transcripts of interest (survivors that are also in the legacy universe) on NC_047559.1 to NC_047563.1 against the share of the universe there;
 (e) overlap with the universe gap of I-39: genes admitted by the universe alone, by the reader alone, and by both;
 (f) survival by number of transcripts per gene in the table (one, two to five, more than five);
 (g) the position-only control for ploidy BP: take as "transcripts of interest" every transcript of the legacy universe whose accession lies in the surviving accession run, with no DML information, run the legacy BP test exactly as step 09 does (transcripts as genes, legacy universe with its repeated names, classic Fisher, GO.db 3.13.0 through use_legacy_go_db(config), which must be called before topGO is loaded), and count terms with p < 0.01; compare with random lists of the same size from the universe (at least 20 seeded draws; median and range) and with the real legacy result (383 terms), and report how many of the ten legacy top BP terms for ploidy are among the terms with p < 0.01 in the position-only run.
Targets from the verification: rows 10,979 of 46,626, 81,109 of 106,621, 6,941 of 7,029; transcripts 623 of 2,348, 3,444 of 4,606, 379 of 379; genes 171 of 622, 949 of 1,237, 117 of 117; odd-apostrophe lines 27, 88, 10; merged rows 14, 44, 5; with the apostrophe as the only quote character the ploidy read returns 10,981 rows; R warns "EOF within quoted string" for the ploidy table only; 1,725 ploidy transcripts and 451 genes vanish (heat 1,162 and 288); 603 of the 623 surviving ploidy transcripts are one accession run (XM_034462509.1 to XM_034479694.1); 601 of the 612 tested transcripts of interest (98.2%) on NC_047559.1 to NC_047563.1, which hold 62.1% of the universe; genes admitted by universe alone 456, 930, 91, by reader alone 171, 949, 117, by both 168, 610, 91; heat genes with one transcript survive in 64% and with more than five in 96% (ploidy 20% and 37%); position-only control 1,041 of 12,720 BP terms with p < 0.01 against a median of 52.5 for random lists, nine of the ten legacy top terms among them.
Outputs under 06_diagnostics/: legacy_reader_mechanism.csv (long format: table, quantity, value, definition), legacy_reader_by_linkage_group.csv, legacy_reader_position_control.csv.`,
  },
  {
    key: 'coverage_filter_variants',
    issue: 'I-05 (the coverage filter that was never applied)',
    basis: 'C6-coverage-filter-locus (cutoffs.R, pctile.R, sens.R, pooled.R, filter_rule_sensitivity.csv, cutoffs_mine.csv)',
    spec: `Write 06_diagnostics/coverage_filter_variants.R. The draft's Methods say CpGs "in the 99.9th percentile of coverage" were removed; the call misspelled the argument and no upper filter ran. The pipeline's corrected mode applies methylKit's filterByCoverage(hi.perc = 99.9) per library, with the percentile taken over CpGs with at least config$import_min_coverage (2) reads. The Methods do not fix that choice. This script shows how much the choice matters:
 (a) per library, the 99.9th percentile cutoff (R quantile type 7, as methylKit) with the percentile taken over CpGs with at least 1 read, at least 2 reads, and at least 10 reads, computed from a coverage histogram of the raw file; and one pooled cutoff over all 22 libraries (CpGs with at least 2 reads);
 (b) for the 480,579 tested CpGs of legacy mode, the raw coverage in each library, and from it, for each rule below, the number of CpGs kept, removed, and the number of legacy ploidy DML, heat DML and loci in both lists that are at kept CpGs: no upper filter; per-library cutoff with base 2 (remove if at or above the cutoff in any library; this must reproduce the pipeline's corrected mode exactly: 454,555 CpGs and the positions of united_cpg_corrected.csv.gz, and the corrected DML lists); the same with removal only above the cutoff; base 1; base 10; the pooled cutoff; and three variants of the "any library" rule with base 2: at least 2 libraries, more than half of the libraries, all 22 libraries;
 (c) how the removed CpGs are distributed: number of libraries at or above the cutoff per removed CpG (exactly 1, 2 to 21, all 22), share in the repeat track, number on the mitochondrial sequence, share lying in runs of 10 or more removed CpGs;
 (d) the evidence that DML counts at kept CpGs are exact for the alternative rules: the q-value rule is not binding (from diffmeth_<contrast>_<mode>.csv.gz count CpGs above the difference threshold that fail q < config$dml_qvalue, in both modes), so the lists change only by deletion.
Raw coverage of the tested CpGs: read each of the 22 files once (one at a time), build the histogram for (a) and keep the coverage at the tested positions for (b). Legacy mode loses the first line of each file; none of those CpGs is tested, but handle the files the way the helpers do.
Targets from the verification: kept CpGs and DML (ploidy, heat, both): no filter 480,579 and 1,038, 2,351, 200; base 2 454,555 and 952, 2,273, 159; above the cutoff only 455,253 and 952, 2,273, 159; base 1 451,328 and 949, 2,267, 159; base 10 469,453 and 959, 2,292, 159; pooled cutoff 88 reads, 450,065 and 949, 2,268, 159; at least 2 libraries 462,296 and 957, 2,285, 159; more than half 475,175 and 1,002, 2,322, 176; all 22 478,611 and 1,038, 2,351, 200. Removed CpGs (base 2): 26,024, of which 7,741 exceed in exactly one library and 1,968 in all 22; 43.4% in the repeat track against 21.7% of retained CpGs; 363 mitochondrial; 62.3% in runs of 10 or more. q-value rule: fails for 0 of 1,038 ploidy and 1 of 2,352 heat candidates in legacy mode, 0 and 0 in corrected mode.
Outputs under 06_diagnostics/: coverage_filter_variants.csv (one row per rule), coverage_filter_cutoffs.csv (per library: the three cutoffs and the pooled one), coverage_filter_removed_cpgs.csv (the summary of c), coverage_filter_qvalue_rule.csv.`,
  },
  {
    key: 'check_methylkit_equivalence',
    issue: 'I-05 and I-38, and the comment at the top of 02_R_functions/methylkit_helpers.R',
    basis: 'R1-coverage-dml (t05_equiv.R, t04_native.R) and C6-coverage-filter-locus (toy_equiv.R, mk_real.R)',
    spec: `Write 06_diagnostics/check_methylkit_equivalence.R. Step 05 does not call methylKit's methRead, filterByCoverage, normalizeCoverage and unite on the 22 full libraries, because that route needs about 9 GB; it applies the same rules one library at a time with helpers in 02_R_functions/methylkit_helpers.R (read_filtered_counts, scale_counts_like_methylkit, make_methylraw) and then calls calculateDiffMeth. The claim to be made repeatable: the helper route gives exactly what methylKit's own route gives, in both modes. The script:
 (a) builds a regional extract of all 22 raw coverage files: a fixed set of regions chosen to exercise every rule (for example the first 6 Mb of NC_047559.1, 300 kb around NC_047563.1:39,633,323 where coverage is extreme, the mitochondrial sequence NC_001276.1, and some unplaced scaffolds), INCLUDING the first line of each file so that the first-line rule of legacy mode is exercised; the extract is written to a temporary directory (tempdir()), one file per library, in the original format;
 (b) native route, legacy mode: methRead with the legacy arguments (header left at its default TRUE, which consumes the first line; mincov = config$import_min_coverage), filterByCoverage(lo.count = config$min_coverage) with no upper filter, normalizeCoverage, unite, then calculateDiffMeth for ploidy with heat as covariate and for heat with ploidy as covariate, as step 05 does; native route, corrected mode: methRead(header = FALSE), filterByCoverage(lo.count = config$min_coverage, hi.perc = config$upper_coverage_percentile), then the same;
 (c) helper route on the same extract files, through the same functions step 05 uses (read 01_pipeline/05_differential_methylation.R to see exactly how it chains them; call the helpers, do not reimplement them);
 (d) compares, per mode: the set of CpGs in the united table; every count; the per-library cutoffs, kept rows and scaling factors; meth.diff and p-values of both contrasts (largest absolute and relative difference); the DML called at the configured thresholds on the extract. Writes one row per comparison with the number of items compared, the number that differ and the largest difference, and stops with an error if anything differs beyond floating-point tolerance (state the tolerance);
 (e) also runs methylKit::normalizeCoverage against scale_counts_like_methylkit on synthetic counts built to contain many exact .5 ties, and reports the number of values that differ.
On the extract the percentile cutoffs differ from the full-data cutoffs; that is fine, both routes see the same files. Keep the extract small enough that the native route stays under 3 GB and the whole script under about 15 minutes; say what share of the tested CpGs the extract covers.
Targets from the verification (on one agent's extract; yours will differ in size): united tables identical in both modes (6,744 and 6,419 CpGs on that extract, zero count differences), same rows kept per library, same scaling factors, identical meth.diff, same DML at thresholds 5, 10 and 20; scale_counts_like_methylkit bit-identical to normalizeCoverage over 80 synthetic libraries with 20,129 ties.
Also report, but do not fix (other files are out of your scope): 01_pipeline/04_cpg_coverage.R scales counts for its correlation and PCA in a different order of operations from methylKit (it rounds (new_cov * nc) / cv where methylKit rounds new_cov * (nc / cv)); include in the script a small check that counts how many values on the extract differ between that order and methylKit's, so the size of the effect is on record.
Output under 06_diagnostics/: check_methylkit_equivalence.csv.`,
  },
]

const BUILD_SCHEMA = {
  type: 'object',
  properties: {
    script: { type: 'string' },
    outputs: { type: 'array', items: { type: 'string' } },
    run_time_minutes: { type: 'number' },
    peak_memory_gb: { type: 'number' },
    summary: { type: 'string' },
    targets: {
      type: 'array',
      description: 'every target number of the task: what the script obtained and whether it agrees',
      items: {
        type: 'object',
        properties: { quantity: { type: 'string' }, target: { type: 'string' }, obtained: { type: 'string' }, agrees: { type: 'boolean' }, comment: { type: 'string' } },
        required: ['quantity', 'target', 'obtained', 'agrees'],
      },
    },
    design_decisions: { type: 'array', items: { type: 'string' } },
    left_out: { type: 'array', items: { type: 'string' }, description: 'parts of the specification not implemented, each with the reason' },
    caveats: { type: 'array', items: { type: 'string' } },
    readme_row: { type: 'string', description: 'one Markdown table row for 06_diagnostics/README.md: | `file` | purpose, run time, issues | outputs |' },
  },
  required: ['script', 'outputs', 'run_time_minutes', 'summary', 'targets', 'design_decisions', 'left_out', 'caveats', 'readme_row'],
}

const REVIEW_SCHEMA = {
  type: 'object',
  properties: {
    runs_clean: { type: 'boolean' },
    outputs_reproduced_on_rerun: { type: 'boolean', description: 'a second run of the script gave byte-identical CSV files' },
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
    numbers_recomputed: { type: 'array', items: { type: 'string' }, description: 'each number you recomputed with your own code, with both values' },
    statements_in_issues_md_affected: { type: 'array', items: { type: 'string' }, description: 'sentences of issues.md that the script output contradicts or refines, with the exact current text and the value the script gives' },
    verdict: { type: 'string' },
  },
  required: ['runs_clean', 'outputs_reproduced_on_rerun', 'findings', 'numbers_recomputed', 'statements_in_issues_md_affected', 'verdict'],
}

const FIX_SCHEMA = {
  type: 'object',
  properties: {
    script: { type: 'string' },
    outputs: { type: 'array', items: { type: 'string' } },
    changes_made: { type: 'array', items: { type: 'string' } },
    findings_not_applied: { type: 'array', items: { type: 'string' }, description: 'each review finding you did not apply, with the reason' },
    final_run_clean: { type: 'boolean' },
    run_time_minutes: { type: 'number' },
    targets: {
      type: 'array',
      items: {
        type: 'object',
        properties: { quantity: { type: 'string' }, target: { type: 'string' }, obtained: { type: 'string' }, agrees: { type: 'boolean' }, comment: { type: 'string' } },
        required: ['quantity', 'target', 'obtained', 'agrees'],
      },
    },
    issues_md_changes: {
      type: 'array',
      description: 'for the analyst who maintains issues.md: each sentence or number there that should now cite this script, or that the script shows to be wrong',
      items: {
        type: 'object',
        properties: { issue: { type: 'string' }, current_text: { type: 'string' }, proposed_text: { type: 'string' }, reason: { type: 'string' } },
        required: ['issue', 'current_text', 'proposed_text', 'reason'],
      },
    },
    readme_row: { type: 'string' },
    notes: { type: 'string' },
  },
  required: ['script', 'outputs', 'changes_made', 'findings_not_applied', 'final_run_clean', 'run_time_minutes', 'targets', 'issues_md_changes', 'readme_row', 'notes'],
}

const buildPrompt = it => `${CONTEXT}
YOUR TASK (build). Issue: ${it.issue}. Verifying agents whose work is your starting point: ${it.basis} (folders under ${VER}/).
${it.spec}
Work in this order: read; design; write the script; run it from the repository root with Rscript; check every target; read your own output tables as a skeptical reader would. Run the finished script twice and confirm the CSV files are byte-identical between runs. Leave the script and its outputs in 06_diagnostics/. Return the structured report; list every target with what you obtained.`

const reviewPrompt = (it, build) => `${CONTEXT}
YOUR TASK (independent review). Another agent has just written 06_diagnostics/${it.key}.R for this purpose:
${it.spec}
Its own report: ${JSON.stringify(build)}
You did not write it and you should assume it contains mistakes. Do not modify the script or its outputs (work on copies in your scratch directory if you need to change anything to test it).
1. Read the script line by line against the specification and against the helpers it calls. Look for: a statistic that is not what its column name says; a comparison that silently passes (for example a check written so that it cannot fail); wrong strata, wrong group assignment, wrong library order; off-by-one in coordinates (the coverage files and united tables use the 0-based start of the C; BED tracks are 0-based half-open; GFF is 1-based closed; see 02_R_functions/genomic_overlaps.R); numbers typed in instead of computed; results that depend on thread count or on an unset seed; data.table scoping mistakes; paths that would not work from a fresh clone; violation of the memory or style rules.
2. Run it (copy its current outputs aside first), then run it again, and compare the CSV files of the two runs and with the copies.
3. Recompute with your own code, written without looking at how the script does it, at least five of the key numbers it reports, chosen to cover different parts of the script. Use a different route where you can (base R instead of data.table, direct counting instead of the helper).
4. Compare the script's outputs with the statements in issues.md that it is meant to support (the Verification paragraphs of the issue named above) and list every sentence that the output contradicts or sharpens.
Report only real problems, each with evidence. Severity: high = a reported number is wrong or a check cannot fail; medium = a number is right but its label, definition or interval is misleading, or the script breaks on a legitimate rerun; low = style, wording, robustness.`

const fixPrompt = (it, build, review) => `${CONTEXT}
YOUR TASK (fix and finish). The script 06_diagnostics/${it.key}.R was written for this purpose:
${it.spec}
Builder's report: ${JSON.stringify(build)}
Independent review: ${JSON.stringify(review)}
For each finding of the review: check it yourself first. Apply the fix if the finding is right; if it is wrong or not worth applying, leave the code and say why under findings_not_applied. Do not make changes the review did not ask for, except to repair something you break. Then run the script from the repository root, twice, confirm the outputs are byte-identical between the two runs, and re-check every target. Remove output files in 06_diagnostics/ that the final script no longer writes (only files named ${it.key}*.csv or the names listed in the specification). Finally, read the Verification paragraphs of the issue in ${REPO}/issues.md and list under issues_md_changes every sentence or number that should now cite this script's output, or that the output shows to be wrong, with exact current text and exact proposed text (plain sentences, no em dashes). Do not edit issues.md or any README yourself.`

const results = await pipeline(
  ITEMS,
  it => agent(buildPrompt(it), { label: `build:${it.key}`, phase: 'Build', schema: BUILD_SCHEMA }),
  (build, it) => build ? agent(reviewPrompt(it, build), { label: `review:${it.key}`, phase: 'Review', schema: REVIEW_SCHEMA }).then(review => ({ build, review })) : null,
  (prev, it) => (prev && prev.review) ? agent(fixPrompt(it, prev.build, prev.review), { label: `fix:${it.key}`, phase: 'Fix', schema: FIX_SCHEMA }).then(fix => ({ key: it.key, build: prev.build, review: prev.review, fix })) : (prev ? { key: it.key, build: prev.build, review: null, fix: null } : null)
)
const done = results.filter(Boolean)
log(`${done.filter(r => r.fix && r.fix.final_run_clean).length} of ${ITEMS.length} diagnostics finished clean`)
return { diagnostics: done, not_finished: ITEMS.map(i => i.key).filter(k => !done.some(r => r.key === k && r.fix)) }

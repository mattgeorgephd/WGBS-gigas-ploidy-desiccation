# 06_diagnostics

Checks that the pipeline reproduces what it should, and controls for the analyses that replace legacy ones. These scripts are not part of the analysis: no analysis step reads their outputs. Step 13, which ties the manuscript's numbers to outputs, may quote them.

## Files

| File | Purpose | Output |
|---|---|---|
| `check_legacy_outputs.R` | Compares the legacy-mode outputs of the pipeline with the outputs of the legacy analysis (`08_archive/`, and the legacy tables in `03_input_files/annotation/legacy/` that exist only on gannet). One row per comparison. | `legacy_check_results.csv` |
| `go_random_cpg_control.R` | Negative control for the GO enrichment: CpGs drawn at random from the eligible CpGs are passed through the gene-level tests as if they were DML. About 10 minutes. | `go_random_cpg_control.csv` (every draw), `go_random_cpg_control_summary.csv` |
| `dml_sensitivity.R` | Sensitivity of the DML calls of step 05. The per-CpG tests are repeated under each alternative that the audit raises (overdispersion correction, a tested interaction, without the acute heat shock libraries, within one cohort, within one ploidy), and each DML list is compared with the list of the draft. It also writes what a list size must be read with: a permutation reference for the ploidy contrasts, the share of small p-values, the test statistic that a call needs, and the dispersion of the counts next to that of binomial reads. One methylKit fit per analysis and mode (22 fits): about 2 hours on 4 free cores; with `WGBS_CORES=3` on a shared machine a fit took 8 minutes on average, which makes about 3 hours. Fits are cached, an interrupted run resumes, and a run that finds every fit in the cache takes about 10 minutes. Issues I-01, I-02, I-04, I-05, I-06, I-07. | Folder `dml_sensitivity/` (tables, DML lists, figure; see "Sensitivity of the DML calls") |
| `cohort_separation.R` | Do the libraries group by cohort or by ploidy before any locus is tested? Median correlation within against between groups, and the number of libraries a single threshold on each principal component assigns to their group, each with a permutation p-value (10,000 label permutations). Reads the step 04 tables. Half a minute. Issues I-16 and I-44. | `cohort_separation_correlation.csv`, `cohort_separation_pca.csv` |
| `universe_file_cut.R` | Shows where the gene universe of the legacy GO enrichment lost its run of 15,666 transcripts (15,653 of them to the cut): fetches 4 KB of the 20 GB table on gannet around the place where a line is cut and a line from far later in the sort order follows, and checks that the two transcripts at the cut bound the gap in the universe. A few seconds; needs gannet. Issue I-39. | `universe_file_cut.csv`, `universe_file_cut_context.txt` |
| `table1_group_means.R` | What the "methylation difference" values of Table 1 are made of: read-pooled percent methylation of the four groups at each of its 20 loci, the ploidy difference within each cohort, the cohort difference within each ploidy, and the per-library values. Descriptive; no test. Under a minute. Issues I-02 and I-07. | `table1_group_means.csv`, `table1_library_values.csv` |
| `library_methylation_level.R` | What the differences between groups of libraries consist of. Overall methylation level of each library; the ploidy contrast inside each cohort by direction; whether a CpG's dependence on the level predicts the contrast of the other cohort; counts with one library left out; complete enumeration of regroupings, with and without the five libraries of highest level, and for the acute heat shock; similarity of the libraries of each group; library values at the CpGs called with the overdispersion correction. Tests nothing at single CpGs. About 3 minutes, 3 GB; part 8 needs the outputs of `dml_sensitivity.R`. Issues I-49, I-01, I-02, I-04, I-07, I-08. Written by the audit session; independent review pending. | `library_methylation_level/` (13 tables, described below) |
| `fastq_read_names.R` | Sequencer run, flowcell and lane of each library, from the read names at the start of each trimmed FASTQ file on gannet (first 10,000 reads of read 1, streamed; about a minute; needs the network). Methods 2.2 says one lane of a single flowcell; the ten heat libraries are run 163, flowcell HLTJMDSXY, lane 1, and the twelve control libraries run 173, flowcell HMGVGDSXY, lane 2. Issues I-01 and I-28 | `fastq_read_names.csv` |
| `feature_odds_matched.R` | Location of DML against matched and label-permuted backgrounds (I-41). Written by an agent; review pending | `feature_odds_matched.csv`, `feature_odds_matched_within_group_sd.csv`, `feature_odds_label_permutation.csv`, `feature_odds_label_permutation_counts.csv` |
| `locus_nc047563.R` | What the locus of I-09 is: extent of the excess coverage, per-library coverage, leave-one-library-out. Needs the coverage files and the reference genome. Written by an agent; review pending | `locus_nc047563_summary.csv`, `locus_nc047563_per_library.csv`, `locus_nc047563_leave_one_out.csv`, `locus_nc047563_coverage_profile.csv` |
| `legacy_reader_mechanism.R` | How the legacy reader lost annotation lines and what the loss did to the test (I-40, I-39, I-43). Written by an agent; review pending | `legacy_reader_mechanism.csv`, `legacy_reader_by_linkage_group.csv`, `legacy_reader_position_control.csv` |
| `coverage_filter_variants.R` | How far the CpG and DML counts depend on the reading of the upper coverage filter (I-05). Needs the coverage files. Written by an agent; review pending | `coverage_filter_variants.csv`, `coverage_filter_cutoffs.csv`, `coverage_filter_removed_cpgs.csv`, `coverage_filter_qvalue_rule.csv` |

Run from the repository root, after the pipeline steps they check:

```sh
Rscript 06_diagnostics/check_legacy_outputs.R
Rscript 06_diagnostics/go_random_cpg_control.R
Rscript 06_diagnostics/table1_group_means.R
Rscript 06_diagnostics/cohort_separation.R
Rscript 06_diagnostics/universe_file_cut.R
Rscript 06_diagnostics/dml_sensitivity.R                                 # full run, after step 05; hours, see below
Rscript 06_diagnostics/library_methylation_level.R                       # after steps 02 to 05 and, for its part 8, after dml_sensitivity.R
Rscript 06_diagnostics/fastq_read_names.R                                # needs the network; reads the start of 22 FASTQ files on gannet
Rscript 06_diagnostics/feature_odds_matched.R                            # review pending
Rscript 06_diagnostics/locus_nc047563.R                                  # review pending; needs the coverage files and the genome
Rscript 06_diagnostics/legacy_reader_mechanism.R                         # review pending
Rscript 06_diagnostics/coverage_filter_variants.R                        # review pending; needs the coverage files
WGBS_SENSITIVITY_SUBSET=20000 Rscript 06_diagnostics/dml_sensitivity.R   # test of the code on a random sample of CpGs, 10 to 15 minutes
```

## Result of the legacy comparison

`legacy_check_results.csv` holds 331 comparisons; 329 match. Checks whose pipeline outputs are absent are reported as "not run" instead of failing.

| What is compared | Comparisons | Result |
|---|---|---|
| ELISA %5-mC per animal and the two-way ANOVA | 5 | match |
| %mCpG two-way ANOVA (interaction F) | 1 | 14.43 against 14.49; the legacy sheet mixed two alignment runs (issue I-13), so exact agreement is not expected |
| CpGs covered in all libraries of a set, 5 sets by 5 coverage levels | 25 | match |
| Ploidy DML, heat DML, loci in both lists: counts, loci, differences, q-values, row positions | 19 | match (q-values to five significant digits) |
| DML counts at ten difference thresholds | 60 | 59 match; the other is a swap of two typed values in the legacy workbook (issue I-19) |
| DML overlapping each of 13 feature tracks, 6 DML sets | 78 | match |
| Transcript-to-gene map; the three DML annotation tables, line for line | 17 | match |
| Nine topGO tables: terms, counts per term, top-ten lists, terms with p < 0.01 | 45 | match |
| DML in genes with a top-ten term (the gene lists behind Figure 9 and Table 1) | 18 | match |
| Figure 5 values | 36 | match |
| Figure 9 values | 27 | match |

Three conditions are needed for the GO tables to match, and all three are part of step 09: the ontology release GO.db 3.13.0, the legacy reader of the annotation tables (issue I-40), and the universe names taken with their repeats, one per line of the map (issue I-43).

## Result of the random-CpG control

Number of GO terms called enriched at FDR 0.05 (BP, CC, MF), for random CpGs (median of 10 draws) and for the real DML:

| Analysis | Random, ploidy-sized list | Real ploidy DML | Random, heat-sized list | Real heat DML |
|---|---|---|---|---|
| Genes with a tested CpG as universe, Fisher | 37, 40, 39 | 102, 48, 52 | 129.5, 67.5, 63.5 | 197, 82, 78 |
| Genes with an eligible CpG as universe, Fisher | 49.5, 12, 1.5 | 0, 2, 8 | 87.5, 22, 11.5 | 22, 0, 31 |
| Genes with an eligible CpG as universe, Wallenius (primary analysis of step 10) | 0, 0, 0 | 0, 0, 0 | 0, 0, 0 | 0, 0, 0 |

Random CpGs are "enriched" for dozens to hundreds of terms against the tested universe (15 to 200 terms per draw and ontology) and for none under the third analysis. Under the second analysis the number varies widely (medians from 1.5 to 87.5; 0 to 194 terms per draw and ontology). Enrichment under the first two therefore reflects which genes are methylated and how many methylated CpGs they hold, not ploidy or heat. See issue I-41.

**Limits of this control.** CpGs drawn uniformly from the eligible CpGs are the one null under which the number of eligible CpGs per gene is exactly the right adjustment, so this control is easy on the third analysis and cannot tell a calibrated test from a conservative one. Independent verification found the Wallenius analysis conservative under this null (a term in under 1% of runs), at the nominal level when ploidy labels are shuffled within cohort (a term in 4% to 6% of 250 permutations), and above it (6% to 17%) for heat-sized lists drawn in proportion to each CpG's chance of passing the 20-point rule. A label-permutation control and a spike-in positive control are to be added here (plan, phase 5); until then those figures come from `07_documentation/agent_workflows/verification_scripts/C4a-eligible-theory`.

## Sensitivity of the DML calls

`dml_sensitivity.R` shows how the DML lists of step 05 change under each defensible alternative to the analysis of the draft, so that the authors can decide what the paper should report. Each analysis is one entry of a registry in the script and one `methylKit::calculateDiffMeth` call on the CpGs of `united_cpg_<mode>.csv.gz`. All analyses therefore test the same CpGs, and a DML is defined as in step 05 (`config$dml_qvalue`, `config$dml_difference`).

| Analysis | What it changes | Issues |
|---|---|---|
| A0 | Nothing: ploidy with heat as covariate and heat with ploidy as covariate, as in step 05. It is fitted first, and the script stops unless it reproduces the p-values, q-values, both DML lists, and the loci in both lists of step 05 | |
| A1 | Overdispersion correction (`overdispersion = "MN"`), once with the F test and once with the chi-squared test | I-06 |
| A2 | A tested ploidy by heat interaction, without and with the correction, set against the loci in both primary lists. Because heat is confounded with cohort, it is a ploidy by cohort interaction | I-02, I-01 |
| A3 | Both contrasts without the heat libraries that also had the acute heat shock | I-04 |
| A4 | Ploidy among the controls alone and among the heat-treated libraries alone, without and with the correction | I-07, I-01 |
| A5 | Heat against control among diploids alone and among triploids alone. Still a contrast between cohorts | I-01 |

All analyses run in `legacy` mode. A0, A1 with the F test, and A2 also run in `corrected` mode (I-05).

Files in `dml_sensitivity/`. Every CSV has a column `run_scope` that says whether it comes from a full run or a test run.

| File | Content |
|---|---|
| `analyses.csv` | The registry: question, tested indicator, libraries, group sizes, covariates, overdispersion, test, residual degrees of freedom, reference distribution of the test statistic, definition of the effect size |
| `equivalence_with_step05.csv` | A0 against step 05: largest differences in p-values, q-values, and effect sizes, and whether the lists are identical. It is written as soon as a mode has been checked, also when the check fails and the run stops |
| `dml_counts.csv` | One row per analysis and mode. The list: `n_cpg_tested`, `n_significant` (q below the cutoff), `n_dml`, `n_hyper`, `n_hypo`. Its overlap with the primary list of the same contrast: `n_primary`, `n_in_primary`, `n_in_primary_same_direction`, `share_of_primary_recovered`, `share_of_list_in_primary`. What the size of the list must be read with: `n_effect_above_cutoff`, `n_cpg_eligible`, `median_effect_size_eligible`, `share_eligible_p_below_0_05`, `_0_01`, `_0_001`, `p_needed_for_a_single_call`, `statistic_needed_for_a_single_call` |
| `dml_<analysis id>_<mode>.csv` | The CpGs that pass: position, p-value, q-value, effect size, direction. The interaction lists also hold the four group percentages |
| `interaction_vs_shared_<mode>.csv` | The loci in both primary lists, by class, and all other CpGs: how many have a significant tested interaction, with and without the effect-size criterion |
| `permutation_reference.csv` | One row per contrast, set of libraries, and mode: the number of CpGs above the difference cutoff with the observed labels (`n_effect_above_cutoff`), and the distribution of that number when the ploidy labels are permuted among the libraries of each cohort (`n_above_cutoff_min`, `_q05`, `_q25`, `_median`, `_q75`, `_q95`, `_max`), with the number of labellings that exist and that were used and the share of them that give at least the observed number. No model is fitted. The rows of the heat contrasts say that no reference exists |
| `dispersion_summary.csv` | Dispersion of the read counts (Pearson chi-squared over residual degrees of freedom, model ploidy + heat) at 5,000 random eligible CpGs: quartiles and the shares above 1, 1.5, and 3. Three rows per mode: the observed counts; binomial reads simulated at the same CpGs and scaled and rounded as in step 05; the same reads not scaled |
| `figure_dml_sensitivity.png`, `.pdf` | DML per analysis for each contrast and mode, with the part that is also in the primary list |
| `session_info.txt` | R and package versions of the run |

### How to read the size of a list

A long list is not by itself evidence of differences between the groups, and a short list is not evidence of their absence. Read each row of `dml_counts.csv` with the following.

1. **`n_effect_above_cutoff` against `n_dml`.** Where the two are equal, the q-value criterion removes nothing: the list is the set of CpGs whose pooled group difference exceeds the cutoff. The difference between two small groups is noisier than the difference between two large ones, so that set is larger by chance in the analyses on fewer libraries (A3, A4, A5) and for the interaction (A2, four groups of five or six libraries).
2. **`permutation_reference.csv`.** It shows how many CpGs exceed the cutoff when the ploidy labels are permuted among the libraries of each cohort. Where at most 1,000 labellings exist, all are counted (924 for the controls, 252 for the heat-treated libraries); otherwise the observed labelling and 999 random ones. `share_of_labellings_at_least_observed` is a permutation p-value in which the observed labelling counts as one labelling (Phipson and Smyth 2010, Stat Appl Genet Mol Biol 9:Article 39, doi:10.2202/1544-6115.1585). A share that is not small means that the observed labels give no more CpGs above the cutoff than random labels do. A small share means that the diploid and the triploid libraries of the analysis differ (for the interaction: that their difference is not the same in the two cohorts). It does not say which CpGs differ, and it cannot separate ploidy from anything else that separates the diploid from the triploid libraries of a cohort, such as a preparation batch. The smallest share that can occur is 1 over the number of labellings used, and 2 over that number where all labellings are used and the two groups of a cohort have the same size, because exchanging the two labels gives the same count. The heat contrasts have no reference: every heat-treated library is of one cohort and every control of the other, so the heat label cannot be permuted (I-01).
3. **The p-values.** `share_eligible_p_below_0_05`, `_0_01`, `_0_001` are the shares of the eligible CpGs (mean methylation across libraries between half the difference cutoff and 100 minus half of it) with a p-value below 0.05, 0.01, and 0.001. Where the groups do not differ, a test that holds its level gives about these shares. In the analyses with the F test of the correction (`_mn_f`) a share well above these values indicates differences; a share slightly above them does not, because under relabelled ploidy labels that test gives 3.5% to 9.3% of p-values below 0.05 in nine relabellings of ten. The chi-squared form of the correction (`_mn_chisq`) does not hold its level on these libraries (7.2%, 2.2% and 0.47% under relabelling), so its shares cannot be read this way (`issues.md`, I-06, items 2 and 4). In the analyses without it a larger share also arises from the dispersion alone (see 6).
4. **`statistic_needed_for_a_single_call`.** A CpG whose p-value is below `p_needed_for_a_single_call` (the q-value cutoff divided by the number of CpGs tested) has a q-value below the cutoff whatever the other CpGs do; where k CpGs are called together, a p-value k times as large is enough. The column gives the test statistic at that p-value for the test of the analysis (`reference_distribution` in `analyses.csv`). For the F test it rises steeply as the residual degrees of freedom fall. A corrected list that is short while the shares of 3 are well above their nominal values therefore means that the groups differ but that single CpGs cannot be called from this number of animals. A corrected list that is short while the shares are at their nominal values gives no evidence of differences at single CpGs. The F threshold applies only where the dispersion estimate of a CpG is above 1. Where it is not, methylKit sets the estimate to 1 and uses the chi-squared test, so the "corrected" p-value is the uncorrected one and a statistic of 31.4 is enough (`analyses.csv`, column `reference_distribution`). Most CpGs of the corrected lists are of this kind, and relabelled libraries get such calls too (`issues.md`, I-06, item 7), so a corrected list is not a list with 1% false discoveries.
5. **`median_effect_size_eligible` with `n_hyper` and `n_hypo`.** A median away from 0 means that one group of libraries is shifted as a whole at CpGs with intermediate methylation. Many CpGs then pass the difference cutoff in one direction, and a list that is nearly all hyper or nearly all hypo reflects that shift and not loci that stand out from the others.
6. **`dispersion_summary.csv`.** The observed row is the dispersion of the counts that are tested. The model of the draft refers its statistic to a chi-squared distribution, which assumes a dispersion of 1. The second row shows what binomial reads give after step 05 has multiplied the counts of a library by its scaling factor and rounded (about 1.18): the scaling alone makes the uncorrected test somewhat anticonservative. The third row shows that the estimator gives about 1 for binomial counts as sequenced. The test statistics of the analyses without the correction are too large by about the observed dispersion (3.7-fold at the median eligible CpG: a factor of about 1.2 from the scaling and about 3.1 from variation among libraries).

The analyses of this script are tests of single CpGs with the variation among animals either ignored (as in the draft) or estimated per CpG from 10 to 22 libraries. Methods made for few animals share information between CpGs when they estimate that variation, or test regions (DSS: Feng, Conneely and Wu 2014, Nucleic Acids Res 42:e69, doi:10.1093/nar/gku154; dmrseq: Korthauer et al. 2019, Biostatistics 20:367-383, doi:10.1093/biostatistics/kxy007). They are not run here; see I-06.

### Test mode

`WGBS_SENSITIVITY_SUBSET=<N>` runs the same code on a seeded random sample of N CpGs of each united table and writes everything to `dml_sensitivity/test_subset/`. That folder and the `cache/` folders are git-ignored. The files of a test run test the code and are not results: q-values depend on the set of CpGs tested. The variable selects test mode as soon as it is set; a value that is empty, not a whole number, or larger than a united table stops the run before any fit.

## Result of the sensitivity analyses

Full run of 2026-11-19 (22 fits; 16 in legacy mode, 6 in corrected mode; A0 reproduces step 05 exactly in both). The reading is in `issues.md`, I-06, with the consequences for the interaction (I-02), the acute heat shock libraries (I-04), the control-only contrast (I-07) and the confound (I-01). In short, legacy mode, with corrected mode in parentheses where it was run:

| Analysis | Libraries | DML, uncorrected | DML with the correction, F test | Share of eligible CpGs with p < 0.05, uncorrected and corrected |
|---|---|---|---|---|
| Ploidy, heat as covariate (the draft) | 22 | 1,038 (952) | 2 (2); 14 with the chi-squared test | 33% and 5.7% |
| Heat, ploidy as covariate (the draft) | 22 | 2,351 (2,273) | 16 (6); 185 with the chi-squared test | 43% and 11% |
| Ploidy by heat interaction | 22 | 16,971 (16,257) | 2 (2) | 50% and 20% |
| Ploidy without the acute heat shock libraries | 18 | 1,202 | not run | 30% |
| Heat without the acute heat shock libraries | 18 | 3,545 | not run | 39% |
| Ploidy among the controls | 12 | 2,768 | 2 | 37% and 6.7% |
| Ploidy among the heat-treated | 10 | 12,405 | 125 | 50% and 18% |
| Heat among the diploids | 11 | 12,487 | not run | 49% |
| Heat among the triploids | 11 | 6,632 | not run | 46% |

- Dispersion at eligible CpGs: median 3.67 for the observed counts, 1.18 for binomial reads scaled as in step 05 (corrected mode 3.59 and 1.19).
- Relabelling reference, CpGs above 20 points with the observed ploidy labels against labels reassigned within cohort: all libraries 1,038 against a median of 749, reached by 113 of 1,000 relabellings (corrected mode 952 against 706.5, 137 of 1,000); without the acute heat shock libraries 1,214 against 1,220.5 (505 of 1,000); controls only 3,350 against 2,912.5 (242 of all 924); heat-treated only 13,028 against 4,019.5 (2 of all 252, the smallest share possible); interaction 23,503 against 9,379 (1 of 1,000; corrected mode 22,846 against 9,017, 1 of 1,000).
- Of the draft's 200 loci in both lists, 99 pass the tested interaction without the correction and none with it (58 and none of the 159 in corrected mode).
- What an independent check added (I-06; scripts in `07_documentation/agent_workflows/verification_scripts/S1-I06-numbers-and-reading/`): the corrected F test holds its level on average under relabelling (5.6%, 1.2% and 0.14% of p-values below 0.05, 0.01 and 0.001), the chi-squared form of the correction does not (7.2%, 2.2%, 0.47%), and in the far tail neither is calibrated: relabelled libraries get a median of 0, 2 and 10 corrected calls (22 libraries, controls, heat-treated), nearly all at CpGs whose dispersion estimate is not above 1, where methylKit applies no correction. To be written by the script in the code pass.

What the script does not yet write, and what a statement in I-06 therefore takes from a scratch computation on the cached fits: the number of CpGs with q < 0.01 by whether the coverage filter keeps them (157 of the 287 corrected heat CpGs of legacy mode are removed by the filter, 138 of them mitochondrial), and the lists under the Benjamini-Hochberg adjustment (15 and 4 heat CpGs). To be added to the script in the code pass.

## What the group differences consist of

`library_methylation_level.R` describes the libraries one by one and the contrasts between groups of them (I-49). Files in `library_methylation_level/`:

| File | Content |
|---|---|
| `library_level.csv` | One row per library: genome-wide percent CpG, CHG and CHH methylation (Bismark, step 02), mean percent methylation over the tested and over the eligible CpGs in each mode, principal components 1 to 4 of step 04, the ELISA value of the animal where one exists, and `high_level` (above the largest gap in the sorted genome-wide values; `gap_rule` states the gap) |
| `level_group_tests.csv` | Exact Wilcoxon rank-sum tests of each level measure: ploidy inside each cohort, cohort, ploidy over all libraries, acute heat shock inside the heat cohort; ranges, medians, and whether the two groups are completely separated |
| `level_correlations.csv` | Genome-wide level against the principal components, the means over tested and eligible CpGs, CHH methylation (all libraries and by cohort) and the ELISA value (six animals) |
| `contrast_direction.csv` | For eight contrasts and both modes: CpGs whose read-pooled difference exceeds the DML threshold, by sign; CpGs exactly at the threshold; mean, median and quartiles of the difference over the eligible CpGs |
| `contrast_agreement.csv` | The ploidy difference of the heat cohort against that of the control cohort across CpGs: Pearson and Spearman correlation, CpGs above the threshold in both, how many of these have the same sign, and the overlap expected if the two sets were unrelated; with all libraries and without the high ones |
| `level_slope_prediction.csv`, `level_slope_deciles.csv` | Slope of each CpG's methylation on the genome-wide level, estimated among the libraries of one cohort, against the ploidy difference of the other cohort: correlations, and by tenth of slope the share of eligible CpGs above the threshold in each direction |
| `leave_one_library_out.csv` | CpGs above the threshold for the ploidy contrast inside each cohort with each library left out in turn |
| `relabelling_reference.csv` | Complete enumeration: CpGs above the threshold for the observed grouping against every regrouping of the same libraries into groups of the same sizes (for the comparison of both cohorts the labels are moved within cohort, and for the acute heat shock within ploidy; column `labels_moved`). Six comparisons per mode; `largest_other_grouping` excludes the observed grouping and its mirror image |
| `relabelling_small_sets.csv` | Every grouping, with its count, for the comparisons that have at most 252 groupings |
| `pair_similarity.csv` | Mean pairwise correlation among the diploids and among the triploids of each cohort, against all subsets of the same size from that cohort |
| `corrected_dml_library_values.csv`, `corrected_dml_group_apart.csv` | Percent methylation of every library at the CpGs that `dml_sensitivity.R` calls with the overdispersion correction, the range of each group, and whether one group lies outside the range of all other libraries |

Four points for reading these tables.

1. **Thresholds and ties.** A CpG counts as above the threshold when the absolute difference exceeds 20 points by more than 1e-9. A difference of exactly 20 points is not counted; `contrast_direction.csv` reports such CpGs separately (20 for the ploidy contrast inside the heat cohort). `dml_sensitivity.R` compares floating-point values without that margin, so its counts can differ from these by a few CpGs (3,350 against 3,349 for the control cohort; medians 4,019.5 against 4,018.5 for the heat cohort). The code pass will give both scripts one rule.
2. **Complete enumeration.** Where all regroupings are counted, the share at or above the observed count is an exact p-value for the hypothesis that the grouping is unrelated to methylation. It cannot fall below (groupings as extreme as the observed one) divided by (all groupings): 2 of 252 for five against five, 1 of 21 for two against five.
3. **Exploratory rows.** The five high libraries are found by a rule that ignores the group labels. Setting them aside is nevertheless a choice made after seeing the data, and the rows "without the high libraries" must be read as exploration.
4. **What an independent check found about reading these tables** (issue I-49; scripts in `07_documentation/agent_workflows/verification_scripts/S2-I49-evidence-and-alternatives/`). `pair_similarity.csv` cannot be compared between cohorts as it stands: the control libraries have fewer reads, and a correlation falls with sampling noise; at equal read depth the heat triploids are as similar to each other as ordinary control libraries are. `level_slope_deciles.csv` shows an association of sign that is in effect the contrast between the two high control libraries and the other ten, not a graded dependence on the level. `level_correlations.csv`: the correlations of the level with the first principal component and with the mean over eligible CpGs are those of the five high libraries against the other 17. `relabelling_reference.csv`: "2 of 252" for the heat cohort mostly restates that its three high libraries are all diploid; the next regrouping is 206 CpGs behind. The computations behind these statements are to be added to the script in the code pass.

## Still to be written

| File | Purpose |
|---|---|
| `check_methylkit_equivalence.R` | Shows that the low-memory import of step 05 gives the same united table as methylKit's own `methRead`, `filterByCoverage`, `normalizeCoverage`, `unite`. Checked once outside the repository: identical united tables on a small end-to-end example in both modes, and identical cutoffs, kept rows and removed positions on four real libraries; in legacy mode the DML lists also equal the archived lists of the legacy analysis (`check_legacy_outputs.R`). No comparison on the complete 22 libraries in corrected mode is on record; the script will make the check repeatable. |
| `check_modes.R` | Side-by-side table of legacy and corrected results for every reported quantity. |

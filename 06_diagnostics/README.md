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

Run from the repository root, after the pipeline steps they check:

```sh
Rscript 06_diagnostics/check_legacy_outputs.R
Rscript 06_diagnostics/go_random_cpg_control.R
Rscript 06_diagnostics/table1_group_means.R
Rscript 06_diagnostics/cohort_separation.R
Rscript 06_diagnostics/universe_file_cut.R
Rscript 06_diagnostics/dml_sensitivity.R                                 # full run, after step 05; hours, see below
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

Random CpGs are "enriched" for dozens to hundreds of terms under the first two analyses and for none under the third. Enrichment under the first two therefore reflects which genes are methylated and how many methylated CpGs they hold, not ploidy or heat. See issue I-41.

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
3. **The p-values.** `share_eligible_p_below_0_05`, `_0_01`, `_0_001` are the shares of the eligible CpGs (mean methylation across libraries between half the difference cutoff and 100 minus half of it) with a p-value below 0.05, 0.01, and 0.001. Where the groups do not differ, a test that holds its level gives about these shares. In the analyses with the correction (`_mn_f`, `_mn_chisq`) a larger share indicates differences. In the analyses without it a larger share also arises from the dispersion alone (see 6).
4. **`statistic_needed_for_a_single_call`.** A CpG whose p-value is below `p_needed_for_a_single_call` (the q-value cutoff divided by the number of CpGs tested) has a q-value below the cutoff whatever the other CpGs do; where k CpGs are called together, a p-value k times as large is enough. The column gives the test statistic at that p-value for the test of the analysis (`reference_distribution` in `analyses.csv`). For the F test it rises steeply as the residual degrees of freedom fall. A corrected list that is short while the shares of 3 are well above their nominal values therefore means that the groups differ but that single CpGs cannot be called from this number of animals. A corrected list that is short while the shares are at their nominal values gives no evidence of differences at single CpGs.
5. **`median_effect_size_eligible` with `n_hyper` and `n_hypo`.** A median away from 0 means that one group of libraries is shifted as a whole at CpGs with intermediate methylation. Many CpGs then pass the difference cutoff in one direction, and a list that is nearly all hyper or nearly all hypo reflects that shift and not loci that stand out from the others.
6. **`dispersion_summary.csv`.** The observed row is the dispersion of the counts that are tested. The model of the draft assumes the dispersion of the second row, not exactly 1: step 05 multiplies the counts of a library by its scaling factor and rounds, which raises the dispersion of binomial reads somewhat above 1. The third row shows that the estimator gives about 1 for binomial counts as sequenced. The test statistics of the analyses without the correction are too large by about the observed dispersion.

The analyses of this script are tests of single CpGs with the variation among animals either ignored (as in the draft) or estimated per CpG from 10 to 22 libraries. Methods made for few animals share information between CpGs when they estimate that variation, or test regions (DSS: Feng, Conneely and Wu 2014, Nucleic Acids Res 42:e69, doi:10.1093/nar/gku154; dmrseq: Korthauer et al. 2019, Biostatistics 20:367-383, doi:10.1093/biostatistics/kxy007). They are not run here; see I-06.

### Test mode

`WGBS_SENSITIVITY_SUBSET=<N>` runs the same code on a seeded random sample of N CpGs of each united table and writes everything to `dml_sensitivity/test_subset/`. That folder and the `cache/` folders are git-ignored. The files of a test run test the code and are not results: q-values depend on the set of CpGs tested. The variable selects test mode as soon as it is set; a value that is empty, not a whole number, or larger than a united table stops the run before any fit.

## Result of the sensitivity analyses

The full run was in progress when this section was written (all 16 analyses of legacy mode fitted; corrected mode fitting). The result will be summarized here and in `issues.md` (I-02, I-04, I-06, I-07) when the run has written its tables.

## Still to be written

| File | Purpose |
|---|---|
| `check_methylkit_equivalence.R` | Shows that the low-memory import of step 05 gives the same united table as methylKit's own `methRead`, `filterByCoverage`, `normalizeCoverage`, `unite`. Confirmed once outside the repository (identical DML lists in both modes); the script will make the check repeatable. |
| `check_modes.R` | Side-by-side table of legacy and corrected results for every reported quantity. |

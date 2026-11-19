# 06_diagnostics

Checks that the pipeline reproduces what it should, and controls for the analyses that replace legacy ones. These scripts are not part of the analysis: no analysis step reads their outputs. Step 13, which ties the manuscript's numbers to outputs, may quote them.

## Files

| File | Purpose | Output |
|---|---|---|
| `check_legacy_outputs.R` | Compares the legacy-mode outputs of the pipeline with the outputs of the legacy analysis (`08_archive/`, and the legacy tables in `03_input_files/annotation/legacy/` that exist only on gannet). One row per comparison. | `legacy_check_results.csv` |
| `go_random_cpg_control.R` | Negative control for the GO enrichment: CpGs drawn at random from the eligible CpGs are passed through the gene-level tests as if they were DML. About 10 minutes. | `go_random_cpg_control.csv` (every draw), `go_random_cpg_control_summary.csv` |
| `cohort_separation.R` | Do the libraries group by cohort or by ploidy before any locus is tested? Median correlation within against between groups, and the number of libraries a single threshold on each principal component assigns to their group, each with a permutation p-value (10,000 label permutations). Reads the step 04 tables. Half a minute. Issues I-16 and I-44. | `cohort_separation_correlation.csv`, `cohort_separation_pca.csv` |
| `table1_group_means.R` | What the "methylation difference" values of Table 1 are made of: read-pooled percent methylation of the four groups at each of its 20 loci, the ploidy difference within each cohort, the cohort difference within each ploidy, and the per-library values. Descriptive; no test. Under a minute. Issues I-02 and I-07. | `table1_group_means.csv`, `table1_library_values.csv` |

Run from the repository root, after the pipeline steps they check:

```sh
Rscript 06_diagnostics/check_legacy_outputs.R
Rscript 06_diagnostics/go_random_cpg_control.R
Rscript 06_diagnostics/table1_group_means.R
Rscript 06_diagnostics/cohort_separation.R
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

Two conditions are needed for the GO tables to match, and both are part of step 09: the ontology release GO.db 3.13.0, and the legacy reader of the annotation tables (issue I-40).

## Result of the random-CpG control

Number of GO terms called enriched at FDR 0.05 (BP, CC, MF), for random CpGs (median of 10 draws) and for the real DML:

| Analysis | Random, ploidy-sized list | Real ploidy DML | Random, heat-sized list | Real heat DML |
|---|---|---|---|---|
| Genes with a tested CpG as universe, Fisher | 37, 40, 39 | 102, 48, 52 | 129.5, 67.5, 63.5 | 197, 82, 78 |
| Genes with an eligible CpG as universe, Fisher | 49.5, 12, 1.5 | 0, 2, 8 | 87.5, 22, 11.5 | 22, 0, 31 |
| Genes with an eligible CpG as universe, Wallenius (primary analysis of step 10) | 0, 0, 0 | 0, 0, 0 | 0, 0, 0 | 0, 0, 0 |

Random CpGs are "enriched" for dozens to hundreds of terms under the first two analyses and for none under the third. Enrichment under the first two therefore reflects which genes are methylated and how many methylated CpGs they hold, not ploidy or heat. See issue I-41.

**Limits of this control.** CpGs drawn uniformly from the eligible CpGs are the one null under which the number of eligible CpGs per gene is exactly the right adjustment, so this control is easy on the third analysis and cannot tell a calibrated test from a conservative one. Independent verification found the Wallenius analysis conservative under this null (a term in under 1% of runs), at the nominal level when ploidy labels are shuffled within cohort (a term in 4% to 6% of 250 permutations), and above it (6% to 17%) for heat-sized lists drawn in proportion to each CpG's chance of passing the 20-point rule. A label-permutation control and a spike-in positive control are to be added here (plan, phase 5); until then those figures come from `07_documentation/agent_workflows/verification_scripts/C4a-eligible-theory`.

## Still to be written

| File | Purpose |
|---|---|
| `check_methylkit_equivalence.R` | Shows that the low-memory import of step 05 gives the same united table as methylKit's own `methRead`, `filterByCoverage`, `normalizeCoverage`, `unite`. Confirmed once outside the repository (identical DML lists in both modes); the script will make the check repeatable. |
| `check_modes.R` | Side-by-side table of legacy and corrected results for every reported quantity. |

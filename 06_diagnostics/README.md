# 06_diagnostics

Checks that the pipeline reproduces what it should, and controls for the analyses that replace legacy ones. These scripts are not part of the analysis: nothing in `01_pipeline/` reads their outputs.

## Files

| File | Purpose | Output |
|---|---|---|
| `check_legacy_outputs.R` | Compares the legacy-mode outputs of the pipeline with the outputs of the legacy analysis (`08_archive/`, and the legacy tables in `03_input_files/annotation/legacy/` that exist only on gannet). One row per comparison. | `legacy_check_results.csv` |
| `go_random_cpg_control.R` | Negative control for the GO enrichment: CpGs drawn at random from the eligible CpGs are passed through the gene-level tests as if they were DML. About 10 minutes. | `go_random_cpg_control.csv` (every draw), `go_random_cpg_control_summary.csv` |

Run from the repository root, after the pipeline steps they check:

```sh
Rscript 06_diagnostics/check_legacy_outputs.R
Rscript 06_diagnostics/go_random_cpg_control.R
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

## Still to be written

| File | Purpose |
|---|---|
| `check_methylkit_equivalence.R` | Shows that the low-memory import of step 05 gives the same united table as methylKit's own `methRead`, `filterByCoverage`, `normalizeCoverage`, `unite`. Confirmed once outside the repository (identical DML lists in both modes); the script will make the check repeatable. |
| `check_modes.R` | Side-by-side table of legacy and corrected results for every reported quantity. |

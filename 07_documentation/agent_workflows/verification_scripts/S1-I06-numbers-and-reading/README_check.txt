Check S1-I06-numbers-and-reading: scripts and outputs (read-only on the repository; one R process at a time, one thread).
01_methylkit_source.R/.out          source of methylKit 1.28.0 logReg, calculateDiffMeth, SLIM, QValuesfun
02_quotations.R/.out                help page of calculateDiffMeth as installed; vignette file list
03_numbers_from_cache.R/.log/.csv   every count and share of dml_counts.csv recomputed from the cached fits and united tables
04_item6.R/.log                     item 6: 157 of 287, 138 mitochondrial, q ranges, Benjamini-Hochberg lists
05_locus_rows.txt                   united-table rows of NC_047565.1:35030701 and NC_047562.1:40685928
06_sample_and_parity.R/.log         3,000 eligible CpGs; my base R version of logReg against the cache and against methylKit:::logReg
07_relabel_level_pooled.R, _a/_b    corrected F, uncorrected and MN chi-squared tests under 150 + 150 relabellings (two CpG samples)
08_relabel_level_within_cohort.R    same inside each cohort (heat: all 125 other groupings; control: 150)
09_permutation_reference.R          item 5 recomputed with the labellings of the audit script; decomposition; hybrid labellings
10_combine_relabel.R/.log           runs a and b combined
11_dml_corrected_p_under_relabelling.R  corrected p-values at the CpGs above 20 points, real labels against 150 relabellings
12_dispersion_reference.R, 12b_binomial_null.R  the 1.18 reference; what scaling alone does to the draft's test
13_what_the_text_omits.R/.log       within-cohort breakdown of the 1,038; corrected lists; q < 0.05 and 0.1; Benjamini-Hochberg
14_correction_not_applied.R/.log    corrected calls whose p-value equals the uncorrected one
15_*.R/.log/.csv                    length of the corrected lists under relabelling (heat cohort, pooled, control cohort)
17_crosstab_permutation.R/.log      which relabellings reach the observed count
18_, 19_                            relabelled shares by kind of relabelling; within-cohort regroupings by placement of the high-level libraries
20_phi_at_corrected_calls.R/.log    dispersion estimates at the corrected calls
21_stratified_enumeration.R/.log    all 4,410 labellings with the five high-level libraries fixed
16_extract_preprint_text.py         parser prepared for the bioRxiv page of Wreczycka et al.; the page answered once (HTTP 200) and then HTTP 429, so nothing was extracted
issues_I06_snapshot_2.txt/_3.txt    the text of I-06 as checked (md5 of issues.md 4419afee681fa01cbcd5faaa73c7be94)

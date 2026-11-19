Check S4-updates-and-cross-file: scripts and logs (read-only on the repository; one R process of mine at a time, one thread).
issues.md as checked: md5 ce80a4f948edd1cdbe346c976bc5f6e2 (modified 2026-11-19 15:58 UTC).
00_inspect_cache.R                 structure of a cached fit
01_recompute_united.R / .log       numbers of the updates of I-01, I-02, I-04, I-07, I-08 from united_cpg_legacy.csv.gz
                                   (N = numCs + numTs, as the repository scripts do) and from the cached fits; writes
                                   united_legacy_matrices.RData for the later scripts (267 MB, deleted after the check; rerun 01 to get it back)
02_alternatives.R / .log           interaction contrast without the five high libraries against all 4,410 regroupings;
                                   acute heat shock: main contrast, acute by ploidy contrast, within-ploidy pairs, 3 against 3;
                                   share of a pooled list above 20 points inside each cohort, real labels against 300 relabellings
02_pooled_list_within_cohort_relabelled.csv   the 300 relabellings of the last point
03_group_apart.R / .log            "one group apart": table counts recomputed; all sets of five heat libraries over all tested CpGs
04_mito_and_spike.R / .log         mitochondrial CpGs (ratio, coverage); spike-in of a common ploidy difference against the sign count
05_mito_coverage_tests.R / .log    relative mitochondrial coverage by cohort and ploidy
06_figure1b_dedup.R / .log         interaction contrast of Figure 1B on deduplicated counts and on counts before deduplication
07_gaps_17.R / .log                margins of the 17 "apart" CpGs of the heat-cohort corrected list

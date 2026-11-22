# Comparison of all analyses

Written by `01_pipeline/optional/compare_analyses.R` from `analyses_included.csv` and `05_differential_methylation/tables/dml_counts_by_mode.csv` in this folder; do not edit by hand. Only `legacy` is the paper's analysis; every other row is an optional analysis, not adopted. A DML has a SLIM q-value below the cutoff of its row (column Q: 0.01 in the paper, 0.05 with option `q05`) and an absolute difference of more than the threshold of its row (percentage points). "Upper filter" is the 99.9th percentile coverage filter that the draft describes but did not apply (I-05). Shared: loci in both the ploidy and the heat list; opposite: different directions in the two lists. Contrasts: the paper's (ploidy with heat as covariate and heat with ploidy as covariate, all 22 libraries), `ctrl` (the ploidy contrast on the 12 control libraries alone; I-07), or `inter` (the paper's, and a tested ploidy by heat interaction whose DML list is in the column Interaction: its DML (hyper/hypo) and how many of them are among the shared loci; I-02). The row marked V2022 has the settings that the Methods of Venkataraman et al. 2022 state for coverage, the test, variant removal and the threshold, at q < 0.01 (the variant rule counts reads with their full CIGAR, I-08).

| Analysis | Coverage | Upper filter | Overdispersion | Contrasts | C-to-T loci | Threshold | Q | CpGs tested | Ploidy DML (hyper/hypo) | Heat DML (hyper/hypo) | Shared (both hyper/both hypo/opposite) | Opposite % | Interaction DML (hyper/hypo; in shared) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy (paper) | 10x | no | none | paper's | kept | 20 | 0.01 | 480,579 | 1,038 (442/596) | 2,351 (1,150/1,201) | 200 (28/49/123) | 61.5 |  |
| ctexcl | 10x | no | none | paper's | excluded | 20 | 0.01 | 452,159 | 717 (302/415) | 1,667 (826/841) | 139 (16/26/97) | 69.8 |  |
| mnchisq | 10x | no | MN, chi-squared | paper's | kept | 20 | 0.01 | 480,579 | 14 (9/5) | 185 (95/90) | 4 (1/2/1) | 25.0 |  |
| mnchisq-ctexcl | 10x | no | MN, chi-squared | paper's | excluded | 20 | 0.01 | 452,159 | 12 (7/5) | 142 (75/67) | 4 (1/2/1) | 25.0 |  |
| corrected | 10x | yes | none | paper's | kept | 20 | 0.01 | 454,555 | 952 (430/522) | 2,273 (1,089/1,184) | 159 (27/49/83) | 52.2 |  |
| corrected-ctexcl | 10x | yes | none | paper's | excluded | 20 | 0.01 | 427,349 | 635 (292/343) | 1,596 (772/824) | 100 (15/26/59) | 59.0 |  |
| corrected-mnchisq | 10x | yes | MN, chi-squared | paper's | kept | 20 | 0.01 | 454,555 | 14 (9/5) | 157 (84/73) | 4 (1/2/1) | 25.0 |  |
| corrected-mnchisq-ctexcl | 10x | yes | MN, chi-squared | paper's | excluded | 20 | 0.01 | 427,349 | 12 (7/5) | 118 (65/53) | 4 (1/2/1) | 25.0 |  |
| cov5 | 5x | no | none | paper's | kept | 20 | 0.01 | 2,297,236 | 9,942 (4,408/5,534) | 21,361 (10,090/11,271) | 1,784 (318/485/981) | 55.0 |  |
| cov5-ctexcl | 5x | no | none | paper's | excluded | 20 | 0.01 | 2,153,017 | 7,046 (3,157/3,889) | 15,530 (7,265/8,265) | 1,219 (228/311/680) | 55.8 |  |
| cov5-mnchisq | 5x | no | MN, chi-squared | paper's | kept | 20 | 0.01 | 2,297,236 | 126 (58/68) | 1,191 (511/680) | 27 (10/12/5) | 18.5 |  |
| cov5-mnchisq-ctexcl | 5x | no | MN, chi-squared | paper's | excluded | 20 | 0.01 | 2,153,017 | 96 (48/48) | 923 (382/541) | 19 (8/8/3) | 15.8 |  |
| corrected-cov5 | 5x | yes | none | paper's | kept | 20 | 0.01 | 2,265,759 | 9,791 (4,380/5,411) | 21,200 (9,987/11,213) | 1,711 (312/482/917) | 53.6 |  |
| corrected-cov5-ctexcl | 5x | yes | none | paper's | excluded | 20 | 0.01 | 2,122,969 | 6,902 (3,133/3,769) | 15,376 (7,168/8,208) | 1,148 (222/308/618) | 53.8 |  |
| corrected-cov5-mnchisq | 5x | yes | MN, chi-squared | paper's | kept | 20 | 0.01 | 2,265,759 | 126 (58/68) | 1,133 (494/639) | 27 (10/12/5) | 18.5 |  |
| corrected-cov5-mnchisq-ctexcl | 5x | yes | MN, chi-squared | paper's | excluded | 20 | 0.01 | 2,122,969 | 96 (48/48) | 875 (368/507) | 19 (8/8/3) | 15.8 |  |
| d50 | 10x | no | none | paper's | kept | 50 | 0.01 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| ctexcl-d50 | 10x | no | none | paper's | excluded | 50 | 0.01 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| mnchisq-d50 | 10x | no | MN, chi-squared | paper's | kept | 50 | 0.01 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| mnchisq-ctexcl-d50 | 10x | no | MN, chi-squared | paper's | excluded | 50 | 0.01 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| corrected-d50 | 10x | yes | none | paper's | kept | 50 | 0.01 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| corrected-ctexcl-d50 | 10x | yes | none | paper's | excluded | 50 | 0.01 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| corrected-mnchisq-d50 | 10x | yes | MN, chi-squared | paper's | kept | 50 | 0.01 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| corrected-mnchisq-ctexcl-d50 | 10x | yes | MN, chi-squared | paper's | excluded | 50 | 0.01 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| cov5-d50 | 5x | no | none | paper's | kept | 50 | 0.01 | 2,297,236 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |  |
| cov5-ctexcl-d50 | 5x | no | none | paper's | excluded | 50 | 0.01 | 2,153,017 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |  |
| cov5-mnchisq-d50 | 5x | no | MN, chi-squared | paper's | kept | 50 | 0.01 | 2,297,236 | 1 (0/1) | 63 (30/33) | 0 (0/0/0) |  |  |
| cov5-mnchisq-ctexcl-d50 | 5x | no | MN, chi-squared | paper's | excluded | 50 | 0.01 | 2,153,017 | 1 (0/1) | 48 (21/27) | 0 (0/0/0) |  |  |
| corrected-cov5-d50 | 5x | yes | none | paper's | kept | 50 | 0.01 | 2,265,759 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |  |
| corrected-cov5-ctexcl-d50 | 5x | yes | none | paper's | excluded | 50 | 0.01 | 2,122,969 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |  |
| corrected-cov5-mnchisq-d50 | 5x | yes | MN, chi-squared | paper's | kept | 50 | 0.01 | 2,265,759 | 1 (0/1) | 62 (30/32) | 0 (0/0/0) |  |  |
| corrected-cov5-mnchisq-ctexcl-d50 (V2022) | 5x | yes | MN, chi-squared | paper's | excluded | 50 | 0.01 | 2,122,969 | 1 (0/1) | 47 (21/26) | 0 (0/0/0) |  |  |
| q05 | 10x | no | none | paper's | kept | 20 | 0.05 | 480,579 | 1,038 (442/596) | 2,352 (1,150/1,202) | 201 (28/50/123) | 61.2 |  |
| ctexcl-q05 | 10x | no | none | paper's | excluded | 20 | 0.05 | 452,159 | 717 (302/415) | 1,668 (826/842) | 140 (16/27/97) | 69.3 |  |
| mnchisq-q05 | 10x | no | MN, chi-squared | paper's | kept | 20 | 0.05 | 480,579 | 37 (19/18) | 410 (191/219) | 6 (1/3/2) | 33.3 |  |
| mnchisq-ctexcl-q05 | 10x | no | MN, chi-squared | paper's | excluded | 20 | 0.05 | 452,159 | 28 (14/14) | 314 (153/161) | 6 (1/3/2) | 33.3 |  |
| corrected-q05 | 10x | yes | none | paper's | kept | 20 | 0.05 | 454,555 | 952 (430/522) | 2,273 (1,089/1,184) | 159 (27/49/83) | 52.2 |  |
| corrected-ctexcl-q05 | 10x | yes | none | paper's | excluded | 20 | 0.05 | 427,349 | 635 (292/343) | 1,596 (772/824) | 100 (15/26/59) | 59.0 |  |
| corrected-mnchisq-q05 | 10x | yes | MN, chi-squared | paper's | kept | 20 | 0.05 | 454,555 | 34 (17/17) | 355 (166/189) | 6 (1/3/2) | 33.3 |  |
| corrected-mnchisq-ctexcl-q05 | 10x | yes | MN, chi-squared | paper's | excluded | 20 | 0.05 | 427,349 | 25 (12/13) | 267 (131/136) | 6 (1/3/2) | 33.3 |  |
| cov5-q05 | 5x | no | none | paper's | kept | 20 | 0.05 | 2,297,236 | 10,335 (4,571/5,764) | 21,604 (10,216/11,388) | 1,891 (334/515/1,042) | 55.1 |  |
| cov5-ctexcl-q05 | 5x | no | none | paper's | excluded | 20 | 0.05 | 2,153,017 | 7,357 (3,287/4,070) | 15,717 (7,359/8,358) | 1,297 (238/335/724) | 55.8 |  |
| cov5-mnchisq-q05 | 5x | no | MN, chi-squared | paper's | kept | 20 | 0.05 | 2,297,236 | 321 (141/180) | 2,834 (1,203/1,631) | 45 (12/19/14) | 31.1 |  |
| cov5-mnchisq-ctexcl-q05 | 5x | no | MN, chi-squared | paper's | excluded | 20 | 0.05 | 2,153,017 | 240 (109/131) | 2,195 (915/1,280) | 32 (10/12/10) | 31.2 |  |
| corrected-cov5-q05 | 5x | yes | none | paper's | kept | 20 | 0.05 | 2,265,759 | 10,192 (4,546/5,646) | 21,453 (10,117/11,336) | 1,820 (330/512/978) | 53.7 |  |
| corrected-cov5-ctexcl-q05 | 5x | yes | none | paper's | excluded | 20 | 0.05 | 2,122,969 | 7,220 (3,266/3,954) | 15,572 (7,266/8,306) | 1,228 (234/332/662) | 53.9 |  |
| corrected-cov5-mnchisq-q05 | 5x | yes | MN, chi-squared | paper's | kept | 20 | 0.05 | 2,265,759 | 313 (137/176) | 2,737 (1,163/1,574) | 42 (11/19/12) | 28.6 |  |
| corrected-cov5-mnchisq-ctexcl-q05 | 5x | yes | MN, chi-squared | paper's | excluded | 20 | 0.05 | 2,122,969 | 232 (105/127) | 2,116 (883/1,233) | 29 (9/12/8) | 27.6 |  |
| d50-q05 | 10x | no | none | paper's | kept | 50 | 0.05 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| ctexcl-d50-q05 | 10x | no | none | paper's | excluded | 50 | 0.05 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| mnchisq-d50-q05 | 10x | no | MN, chi-squared | paper's | kept | 50 | 0.05 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| mnchisq-ctexcl-d50-q05 | 10x | no | MN, chi-squared | paper's | excluded | 50 | 0.05 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| corrected-d50-q05 | 10x | yes | none | paper's | kept | 50 | 0.05 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| corrected-ctexcl-d50-q05 | 10x | yes | none | paper's | excluded | 50 | 0.05 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| corrected-mnchisq-d50-q05 | 10x | yes | MN, chi-squared | paper's | kept | 50 | 0.05 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |  |
| corrected-mnchisq-ctexcl-d50-q05 | 10x | yes | MN, chi-squared | paper's | excluded | 50 | 0.05 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |  |
| cov5-d50-q05 | 5x | no | none | paper's | kept | 50 | 0.05 | 2,297,236 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |  |
| cov5-ctexcl-d50-q05 | 5x | no | none | paper's | excluded | 50 | 0.05 | 2,153,017 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |  |
| cov5-mnchisq-d50-q05 | 5x | no | MN, chi-squared | paper's | kept | 50 | 0.05 | 2,297,236 | 3 (1/2) | 83 (37/46) | 0 (0/0/0) |  |  |
| cov5-mnchisq-ctexcl-d50-q05 | 5x | no | MN, chi-squared | paper's | excluded | 50 | 0.05 | 2,153,017 | 2 (0/2) | 60 (26/34) | 0 (0/0/0) |  |  |
| corrected-cov5-d50-q05 | 5x | yes | none | paper's | kept | 50 | 0.05 | 2,265,759 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |  |
| corrected-cov5-ctexcl-d50-q05 | 5x | yes | none | paper's | excluded | 50 | 0.05 | 2,122,969 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |  |
| corrected-cov5-mnchisq-d50-q05 | 5x | yes | MN, chi-squared | paper's | kept | 50 | 0.05 | 2,265,759 | 3 (1/2) | 82 (36/46) | 0 (0/0/0) |  |  |
| corrected-cov5-mnchisq-ctexcl-d50-q05 | 5x | yes | MN, chi-squared | paper's | excluded | 50 | 0.05 | 2,122,969 | 2 (0/2) | 59 (25/34) | 0 (0/0/0) |  |  |
| ctrl | 10x | no | none | ploidy in controls only (ctrl) | kept | 20 | 0.01 | 480,579 | 2,768 (2,211/557) | 2,351 (1,150/1,201) | 194 (69/30/95) | 49.0 |  |
| mnchisq-ctrl | 10x | no | MN, chi-squared | ploidy in controls only (ctrl) | kept | 20 | 0.01 | 480,579 | 10 (6/4) | 185 (95/90) | 0 (0/0/0) |  |  |
| corrected-ctrl | 10x | yes | none | ploidy in controls only (ctrl) | kept | 20 | 0.01 | 454,555 | 2,628 (2,150/478) | 2,273 (1,089/1,184) | 187 (67/28/92) | 49.2 |  |
| inter | 10x | no | none | paper's and interaction (inter) | kept | 20 | 0.01 | 480,579 | 1,038 (442/596) | 2,351 (1,150/1,201) | 200 (28/49/123) | 61.5 | 16,971 (1,746/15,225; 99) |
| corrected-inter | 10x | yes | none | paper's and interaction (inter) | kept | 20 | 0.01 | 454,555 | 952 (430/522) | 2,273 (1,089/1,184) | 159 (27/49/83) | 52.2 | 16,257 (1,550/14,707; 58) |

## Counts at other thresholds

Ploidy, heat and shared DML at differences of more than 10, 15, 20 and 25 points (25 is the default of methylKit `getMethylDiff`), from `dml_counts_by_threshold_<id>.csv` of each analysis, at the q-value cutoff of its row. Only the 20-point rows are listed: the sweep does not depend on the threshold an analysis reports, so each `d50` analysis has the sweep of its 20-point twin. The 20-point counts equal the table above (checked by the script).

| Analysis | Q | Ploidy at 10 / 15 / 20 / 25 | Heat at 10 / 15 / 20 / 25 | Shared at 10 / 15 / 20 / 25 |
|---|---|---|---|---|
| legacy | 0.01 | 5,352 / 2,872 / 1,038 / 372 | 10,747 / 5,910 / 2,351 / 928 | 1,623 / 738 / 200 / 50 |
| ctexcl | 0.01 | 4,143 / 2,063 / 717 / 250 | 8,607 / 4,468 / 1,667 / 628 | 1,213 / 506 / 139 / 34 |
| mnchisq | 0.01 | 68 / 26 / 14 / 7 | 398 / 259 / 185 / 123 | 29 / 6 / 4 / 0 |
| mnchisq-ctexcl | 0.01 | 59 / 23 / 12 / 5 | 338 / 208 / 142 / 91 | 25 / 5 / 4 / 0 |
| corrected | 0.01 | 4,837 / 2,665 / 952 / 318 | 9,950 / 5,671 / 2,273 / 903 | 1,474 / 665 / 159 / 41 |
| corrected-ctexcl | 0.01 | 3,694 / 1,879 / 635 / 198 | 7,890 / 4,249 / 1,596 / 605 | 1,085 / 438 / 100 / 25 |
| corrected-mnchisq | 0.01 | 62 / 25 / 14 / 7 | 312 / 214 / 157 / 108 | 26 / 5 / 4 / 0 |
| corrected-mnchisq-ctexcl | 0.01 | 54 / 22 / 12 / 5 | 262 / 169 / 118 / 78 | 22 / 4 / 4 / 0 |
| cov5 | 0.01 | 28,979 / 19,822 / 9,942 / 3,654 | 59,176 / 41,224 / 21,361 / 9,251 | 8,422 / 4,858 / 1,784 / 418 |
| cov5-ctexcl | 0.01 | 22,084 / 14,463 / 7,046 / 2,525 | 46,294 / 31,068 / 15,530 / 6,485 | 6,161 / 3,371 / 1,219 / 269 |
| cov5-mnchisq | 0.01 | 431 / 223 / 126 / 70 | 2,187 / 1,582 / 1,191 / 879 | 203 / 72 / 27 / 9 |
| cov5-mnchisq-ctexcl | 0.01 | 357 / 176 / 96 / 51 | 1,821 / 1,271 / 923 / 668 | 161 / 54 / 19 / 6 |
| corrected-cov5 | 0.01 | 28,189 / 19,444 / 9,791 / 3,568 | 57,942 / 40,722 / 21,200 / 9,194 | 8,161 / 4,709 / 1,711 / 399 |
| corrected-cov5-ctexcl | 0.01 | 21,373 / 14,121 / 6,902 / 2,442 | 45,181 / 30,621 / 15,376 / 6,430 | 5,930 / 3,236 / 1,148 / 250 |
| corrected-cov5-mnchisq | 0.01 | 425 / 220 / 126 / 70 | 2,080 / 1,509 / 1,133 / 838 | 201 / 71 / 27 / 9 |
| corrected-cov5-mnchisq-ctexcl | 0.01 | 351 / 173 / 96 / 51 | 1,727 / 1,210 / 875 / 634 | 159 / 53 / 19 / 6 |
| q05 | 0.05 | 7,678 / 3,164 / 1,038 / 372 | 13,805 / 6,102 / 2,352 / 928 | 2,572 / 833 / 201 / 50 |
| ctexcl-q05 | 0.05 | 5,990 / 2,274 / 717 / 250 | 11,100 / 4,615 / 1,668 / 628 | 1,918 / 565 / 140 / 34 |
| mnchisq-q05 | 0.05 | 136 / 63 / 37 / 16 | 936 / 610 / 410 / 259 | 43 / 11 / 6 / 0 |
| mnchisq-ctexcl-q05 | 0.05 | 114 / 48 / 28 / 10 | 799 / 492 / 314 / 189 | 37 / 8 / 6 / 0 |
| corrected-q05 | 0.05 | 7,124 / 2,991 / 952 / 318 | 13,003 / 5,890 / 2,273 / 903 | 2,390 / 764 / 159 / 41 |
| corrected-ctexcl-q05 | 0.05 | 5,496 / 2,115 / 635 / 198 | 10,373 / 4,419 / 1,596 / 605 | 1,756 / 498 / 100 / 25 |
| corrected-mnchisq-q05 | 0.05 | 125 / 58 / 34 / 15 | 766 / 518 / 355 / 227 | 42 / 11 / 6 / 0 |
| corrected-mnchisq-ctexcl-q05 | 0.05 | 107 / 45 / 25 / 9 | 644 / 412 / 267 / 163 | 36 / 8 / 6 / 0 |
| cov5-q05 | 0.05 | 45,803 / 26,008 / 10,335 / 3,656 | 82,221 / 48,288 / 21,604 / 9,251 | 14,978 / 6,789 / 1,891 / 420 |
| cov5-ctexcl-q05 | 0.05 | 35,448 / 19,163 / 7,357 / 2,527 | 65,094 / 36,650 / 15,717 / 6,485 | 11,069 / 4,739 / 1,297 / 271 |
| cov5-mnchisq-q05 | 0.05 | 904 / 516 / 321 / 185 | 5,180 / 3,806 / 2,834 / 1,967 | 326 / 114 / 45 / 15 |
| cov5-mnchisq-ctexcl-q05 | 0.05 | 748 / 405 / 240 / 136 | 4,297 / 3,041 / 2,195 / 1,482 | 267 / 86 / 32 / 9 |
| corrected-cov5-q05 | 0.05 | 44,898 / 25,667 / 10,192 / 3,570 | 80,959 / 47,879 / 21,453 / 9,194 | 14,642 / 6,655 / 1,820 / 401 |
| corrected-cov5-ctexcl-q05 | 0.05 | 34,647 / 18,853 / 7,220 / 2,444 | 63,954 / 36,279 / 15,572 / 6,430 | 10,775 / 4,616 / 1,228 / 252 |
| corrected-cov5-mnchisq-q05 | 0.05 | 887 / 505 / 313 / 179 | 4,956 / 3,667 / 2,737 / 1,905 | 321 / 111 / 42 / 14 |
| corrected-cov5-mnchisq-ctexcl-q05 | 0.05 | 735 / 395 / 232 / 130 | 4,105 / 2,929 / 2,116 / 1,432 | 262 / 83 / 29 / 8 |
| ctrl | 0.01 | 7,577 / 5,008 / 2,768 / 1,156 | 10,747 / 5,910 / 2,351 / 928 | 1,087 / 559 / 194 / 50 |
| mnchisq-ctrl | 0.01 | 25 / 16 / 10 / 8 | 398 / 259 / 185 / 123 | 0 / 0 / 0 / 0 |
| corrected-ctrl | 0.01 | 6,794 / 4,649 / 2,628 / 1,098 | 9,950 / 5,671 / 2,273 / 903 | 936 / 527 / 187 / 50 |
| inter | 0.01 | 5,352 / 2,872 / 1,038 / 372 | 10,747 / 5,910 / 2,351 / 928 | 1,623 / 738 / 200 / 50 |
| corrected-inter | 0.01 | 4,837 / 2,665 / 952 / 318 | 9,950 / 5,671 / 2,273 / 903 | 1,474 / 665 / 159 / 41 |


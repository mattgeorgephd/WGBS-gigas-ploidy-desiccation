# Comparison of all analyses

Written by `01_pipeline/optional/compare_analyses.R` from `analyses_included.csv` and `05_differential_methylation/tables/dml_counts_by_mode.csv` in this folder; do not edit by hand. Only `legacy` is the paper's analysis; every other row is an optional analysis, not adopted. A DML has SLIM q < 0.01 and an absolute difference of more than the threshold (percentage points). "Upper filter" is the 99.9th percentile coverage filter that the draft describes but did not apply (I-05). Shared: loci in both the ploidy and the heat list; opposite: different directions in the two lists. The row marked V2022 has the settings that the Methods of Venkataraman et al. 2022 state for coverage, the test, variant removal and the threshold (the variant rule counts reads with their full CIGAR, I-08).

| Analysis | Coverage | Upper filter | Overdispersion | C-to-T loci | Threshold | CpGs tested | Ploidy DML (hyper/hypo) | Heat DML (hyper/hypo) | Shared (both hyper/both hypo/opposite) | Opposite % |
|---|---|---|---|---|---|---|---|---|---|---|
| legacy (paper) | 10x | no | none | kept | 20 | 480,579 | 1,038 (442/596) | 2,351 (1,150/1,201) | 200 (28/49/123) | 61.5 |
| ctexcl | 10x | no | none | excluded | 20 | 452,159 | 717 (302/415) | 1,667 (826/841) | 139 (16/26/97) | 69.8 |
| mnchisq | 10x | no | MN, chi-squared | kept | 20 | 480,579 | 14 (9/5) | 185 (95/90) | 4 (1/2/1) | 25.0 |
| mnchisq-ctexcl | 10x | no | MN, chi-squared | excluded | 20 | 452,159 | 12 (7/5) | 142 (75/67) | 4 (1/2/1) | 25.0 |
| corrected | 10x | yes | none | kept | 20 | 454,555 | 952 (430/522) | 2,273 (1,089/1,184) | 159 (27/49/83) | 52.2 |
| corrected-ctexcl | 10x | yes | none | excluded | 20 | 427,349 | 635 (292/343) | 1,596 (772/824) | 100 (15/26/59) | 59.0 |
| corrected-mnchisq | 10x | yes | MN, chi-squared | kept | 20 | 454,555 | 14 (9/5) | 157 (84/73) | 4 (1/2/1) | 25.0 |
| corrected-mnchisq-ctexcl | 10x | yes | MN, chi-squared | excluded | 20 | 427,349 | 12 (7/5) | 118 (65/53) | 4 (1/2/1) | 25.0 |
| cov5 | 5x | no | none | kept | 20 | 2,297,236 | 9,942 (4,408/5,534) | 21,361 (10,090/11,271) | 1,784 (318/485/981) | 55.0 |
| cov5-ctexcl | 5x | no | none | excluded | 20 | 2,153,017 | 7,046 (3,157/3,889) | 15,530 (7,265/8,265) | 1,219 (228/311/680) | 55.8 |
| cov5-mnchisq | 5x | no | MN, chi-squared | kept | 20 | 2,297,236 | 126 (58/68) | 1,191 (511/680) | 27 (10/12/5) | 18.5 |
| cov5-mnchisq-ctexcl | 5x | no | MN, chi-squared | excluded | 20 | 2,153,017 | 96 (48/48) | 923 (382/541) | 19 (8/8/3) | 15.8 |
| corrected-cov5 | 5x | yes | none | kept | 20 | 2,265,759 | 9,791 (4,380/5,411) | 21,200 (9,987/11,213) | 1,711 (312/482/917) | 53.6 |
| corrected-cov5-ctexcl | 5x | yes | none | excluded | 20 | 2,122,969 | 6,902 (3,133/3,769) | 15,376 (7,168/8,208) | 1,148 (222/308/618) | 53.8 |
| corrected-cov5-mnchisq | 5x | yes | MN, chi-squared | kept | 20 | 2,265,759 | 126 (58/68) | 1,133 (494/639) | 27 (10/12/5) | 18.5 |
| corrected-cov5-mnchisq-ctexcl | 5x | yes | MN, chi-squared | excluded | 20 | 2,122,969 | 96 (48/48) | 875 (368/507) | 19 (8/8/3) | 15.8 |
| d50 | 10x | no | none | kept | 50 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| ctexcl-d50 | 10x | no | none | excluded | 50 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| mnchisq-d50 | 10x | no | MN, chi-squared | kept | 50 | 480,579 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| mnchisq-ctexcl-d50 | 10x | no | MN, chi-squared | excluded | 50 | 452,159 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| corrected-d50 | 10x | yes | none | kept | 50 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| corrected-ctexcl-d50 | 10x | yes | none | excluded | 50 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| corrected-mnchisq-d50 | 10x | yes | MN, chi-squared | kept | 50 | 454,555 | 1 (0/1) | 4 (3/1) | 0 (0/0/0) |  |
| corrected-mnchisq-ctexcl-d50 | 10x | yes | MN, chi-squared | excluded | 50 | 427,349 | 1 (0/1) | 2 (1/1) | 0 (0/0/0) |  |
| cov5-d50 | 5x | no | none | kept | 50 | 2,297,236 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |
| cov5-ctexcl-d50 | 5x | no | none | excluded | 50 | 2,153,017 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |
| cov5-mnchisq-d50 | 5x | no | MN, chi-squared | kept | 50 | 2,297,236 | 1 (0/1) | 63 (30/33) | 0 (0/0/0) |  |
| cov5-mnchisq-ctexcl-d50 | 5x | no | MN, chi-squared | excluded | 50 | 2,153,017 | 1 (0/1) | 48 (21/27) | 0 (0/0/0) |  |
| corrected-cov5-d50 | 5x | yes | none | kept | 50 | 2,265,759 | 6 (2/4) | 92 (41/51) | 0 (0/0/0) |  |
| corrected-cov5-ctexcl-d50 | 5x | yes | none | excluded | 50 | 2,122,969 | 4 (1/3) | 63 (28/35) | 0 (0/0/0) |  |
| corrected-cov5-mnchisq-d50 | 5x | yes | MN, chi-squared | kept | 50 | 2,265,759 | 1 (0/1) | 62 (30/32) | 0 (0/0/0) |  |
| corrected-cov5-mnchisq-ctexcl-d50 (V2022) | 5x | yes | MN, chi-squared | excluded | 50 | 2,122,969 | 1 (0/1) | 47 (21/26) | 0 (0/0/0) |  |


# qpcr

Raw quantification cycle (Cq) values from the qPCR assays on ctenidia of animals in the heat experiment. Whether and how these data enter the manuscript is undecided; see issues I-29, I-46 and I-47 in [`issues.md`](../../issues.md). The files are analysed by `01_pipeline/12_qpcr.R`; its tables and figures are in `04_output/12_qpcr/`.

## Files

| File | Contents |
|---|---|
| `raw/qpcr_data_<gene>.csv` | One file per gene for actin, dnmt1, hat, hathap2, hif1a, hsc70, hsp90, mbd2, mecp2. Two columns, `Sample` and `Cq Value`; 80 rows: 40 animals by 2 technical replicates. Sample IDs have the form `D_N_N01` (ploidy, treatment, acute heat shock, animal number). Non-detects are the text `N/A`. The header of the actin file has a trailing space after `Sample`. |
| `raw/qpcr_data_consolidated.csv` | One row per animal for the 32 animals without acute heat shock, one column per gene. Each value is the mean of the two technical replicates, ignoring non-detects; `NA` where both failed. It also holds SOD, ATP synthetase, and COX1, for which no per-gene replicate file exists. Sample IDs have the form `D_C_N_01`; the second letter is `C` for control and `D` for the heat treatment (legacy "desiccation"). This file was assembled in a spreadsheet, not by code. |
| `legacy/dct_analysis_NA-to-45_COX1convert.md` | The knitted output of the legacy qPCR analysis (non-detects set to Cq 45; COX1 as 2^-Cq): the R code with the printed ANOVA and TukeyHSD tables of 11 genes. Step 12 repeats that analysis and compares its numbers with the printed ones. A copy, byte for byte, of `08_archive/ELISA/dct_analysis_NA-to-45_COX1convert.md` (md5 `9e9710a6642c43c787bf8e964538892a`), kept here so that no step reads the archive; the figures the file links to are in the archive only. |

## Things to know

- Actin is the only reference gene, and its Cq ranges from 21.9 to 35.3 across animals.
- Several targets sit near the detection limit (median Cq 37 to 38); MeCP2 is a non-detect in 75% of wells.
- COX1 has no nuclear reference measured on the same template, so these data do not give mitochondrial abundance relative to the nuclear genome.
- Details, per-gene statistics, and the verbatim consolidated table are in `07_documentation/audit/elisa_qpcr_audit.md`.

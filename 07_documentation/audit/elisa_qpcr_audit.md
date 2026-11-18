# Audit of the ELISA 5-mC, percent CpG methylation and qPCR analyses

Paths in this document refer to the repository layout BEFORE the reorganization, at commit 523d734 (branch master). After the reorganization, see 07_documentation/reproducibility_plan.md for where each file moved.

## Summary

This audit covers the three small-data analyses behind the manuscript: the global 5-mC ELISA (N=32, Figure 1A), the Bismark percent CpG methylation comparison (N=22, Figure 1B), and the qPCR data cited in the Discussion, together with every script, spreadsheet and report that feeds them. Both reported two-way ANOVA results (ELISA and percent mCpG) are reproduced exactly by `code/5-mC & 5mCpG plots.Rmd` and match `stats/*.csv`. The ELISA statistics, however, are on arcsine-transformed data (untransformed, ploidy and the interaction are significant at 0.05), and Figure 1A silently drops three of the eight diploid heat-treated samples because of an axis limit. In the percent CpG analysis the treatment is completely confounded with sequencing project: the 12 controls come from a different dataset, and 4 of the 10 heat-treated libraries are from animals that also received the acute heat shock and carry wrong names in the spreadsheet. For qPCR, nothing in the repository measures mitochondrial abundance relative to a nuclear gene, the COX1 pattern is not statistically significant, and the HSP90 ploidy effect depends on how one failed sample is handled.

## How this was checked

- The checks reported here did not modify any repository file.
- Every number below was recomputed independently in Python from the committed data files (plate-reader exports, spreadsheets, CSV tables and Bismark reports) and compared with the R output already committed in the repository (`stats/*.csv` and the knitted `.md` files under `ELISA/`).
- The statistics were then repeated in R 4.3.3 (readxl 1.4.3, agricolae 1.3.8, car 3.1.2, broom 1.0.5, ggplot2 3.4.4). The Python results, the R results and the committed outputs agree:
  - both HSD files are reproduced byte for byte in R, and both ANOVA files to at least 12 significant digits in Python and in R;
  - the 2019 ELISA knit output is reproduced in both: max SD 0.01871462, slope 0.6104974, R-squared 0.9978712, model p 0.0007834074, and all six Tukey p-values;
  - all 11 qPCR ANOVA tables in `ELISA/dct_analysis_NA-to-45_COX1convert.md` are reproduced in Python (HSP90 and COX1 also in R), as are its ten "Removed N rows" plot warnings;
  - the boxplot statistics in section 1.5 were confirmed with ggplot2 itself.
- The plate-reader `.xls` files were read two ways (a spreadsheet converter and R's readxl), with identical values.
- No repository script was run exactly as committed, because several cannot run (Appendix E). Their statistical steps were re-implemented line by line.
- One step could not be repeated: the Johnson transformation in `ELISA/GlobalDNAMeth_Polyploids.R:210-229`. The Johnson package was not available and nothing was installed for the audit.
- Section 2.5 also uses four files that the audit retrieved from the lab server: the two alignment job scripts and the March summary report under `07_documentation/upstream_processing/`, and the June summary report at `03_input_files/bismark_reports/061021-big_bismark_summary_report.txt`. They did not exist at commit 523d734, so they are cited at their current locations.

## Conventions used in this document

- **Verified** means read directly or recomputed. **Inference** is an interpretation that the files support but do not prove. **Undetermined** cannot be answered from the repository.
- Citations are `path:line` at commit 523d734.
- **Treatment names.** The treatment is 24 h of aerial emersion at 27 C. The legacy files call it "desiccation" (for example the `Desiccation` factor in the qPCR scripts and the `Dessication` column of `gene_expression/data/zr3534_wgbs_info.csv`). This document calls it "heat" and the animals "heat-treated", and keeps the legacy word only inside verbatim code and labels. The additional 1 h exposure at 45 C is called "acute heat shock" (legacy names include `heat_shock`, `Heat_Stress` and the `H` code in sample IDs).
- **Sample names** follow `README.md:26-40`: 01-08 control; 09-10 control plus acute heat shock; 11-18 heat; 19-20 heat plus acute heat shock. Prefix D is diploid and T is triploid; the suffix -C is ctenidia.
- **Figure 1A** is the ELISA panel and **Figure 1B** the percent mCpG panel of `figures/final_figures/Figure-1-methylation-percentage.png`.
- The project-wide issue log is `issues.md` at the repository root. This document holds the detailed evidence for the ELISA, percent mCpG and qPCR analyses.

## 1. Global 5-mC ELISA (N=32, Figure 1A)

### 1.1 Computational path (verified)

There are two paths. The manuscript statistics come from path B.

**Path A: the scripted analysis of 2019 and 2021** (`ELISA/GlobalDNAMeth_Polyploids.R` and `ELISA/GlobalDNAMeth_Polyploids.Rmd`).

1. Raw absorbance is `ELISA/docs/20190206_FirstReadingSTRIGG.xls`, sheet "plate list", columns 3 (Well) and 6 ("Absorbance @ 450 (1.0s) (A)") (`ELISA/GlobalDNAMeth_Polyploids.R:28`, `ELISA/GlobalDNAMeth_Polyploids.R:33`). Readings 2 and 3 are used only for a stability check (`ELISA/GlobalDNAMeth_Polyploids.R:36-41`); the maximum per-well SD across the three readings is 0.01871462 (`ELISA/GlobalDNAMeth_Polyploids.md:69`).
2. The plate map is `ELISA/docs/20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx`, sheet "PlateMap_list" (`ELISA/GlobalDNAMeth_Polyploids.R:46`; `ELISA/GlobalDNAMeth_Polyploids.Rmd:46` uses `sheet = 3`).
3. Wells labelled `no_sample` (16), `contam` (A03, B03, C03) and `no_ab` (H08) are excluded (`ELISA/GlobalDNAMeth_Polyploids.R:49-51`). That leaves 76 wells: 14 standard wells and 62 sample wells.
4. Duplicate wells are averaged per sample (`ELISA/GlobalDNAMeth_Polyploids.R:58`). D02-C (well D03) and T08-C (well G08) have a single well each.
5. The standard curve uses the mean OD per standard, with the nominal percentage assigned by row position (`ELISA/GlobalDNAMeth_Polyploids.R:71`). The 5% and 10% points are dropped as plateauing (`ELISA/GlobalDNAMeth_Polyploids.R:79-83`). The fit on 0, 0.1, 0.2, 0.5 and 1% gives slope 0.6104974 and R-squared 0.9978712 (`ELISA/GlobalDNAMeth_Polyploids.md:161`, `ELISA/GlobalDNAMeth_Polyploids.md:171`).
6. The negative-control (NC) OD is the mean of two wells, 0.2240971 (`ELISA/GlobalDNAMeth_Polyploids.R:116`).
7. %5-mC = (mean OD - NC OD) / (slope x 25) x 100, where 25 is the assumed ng of input DNA (`ELISA/GlobalDNAMeth_Polyploids.R:117-120`). This is the kit formula.
8. Group labels come from the sample name: the first letter gives ploidy; a number above 10 is `heat_stress` and below 10 is `control` (`ELISA/GlobalDNAMeth_Polyploids.R:126-139`).
9. Statistics: untransformed `aov(perc_meth ~ ploidy + treatment + ploidy:treatment)` and `TukeyHSD` (`ELISA/GlobalDNAMeth_Polyploids.R:194`, `ELISA/GlobalDNAMeth_Polyploids.R:200`; `ELISA/GlobalDNAMeth_Polyploids.Rmd:178`, `ELISA/GlobalDNAMeth_Polyploids.Rmd:187`).
10. The `.R` version then runs a Shapiro test, applies a Johnson transformation if p<0.05 (`ELISA/GlobalDNAMeth_Polyploids.R:210-225`), and repeats the ANOVA (`ELISA/GlobalDNAMeth_Polyploids.R:229`). It stops at `ELISA/GlobalDNAMeth_Polyploids.R:231-233`, which reference objects that do not exist (`MN`, `DMA`, `kPa`); in R the first of these lines fails with "invalid 'envir' argument of type 'closure'".

**Path B: the manuscript statistics, February 2026.**

1. `data/ELISA-percent-mCPG-summary.xlsx`, sheet `raw` (A1:K63, 62 well rows), holds Excel formulas:
   - G = AVERAGE of the duplicate wells (or the single well);
   - I = G - 0.224; J = 25*0.6105; K = I/J*100.
   - The constants 0.224 and 0.6105 are typed in, rounded from 0.2240971 and 0.6104974.
2. Sheet `summary` (A1:D33: ploidy, treatment, sample, percent) holds 32 static numbers equal to `raw` column K (maximum difference 9.8e-15). They differ from the R-script values by at most 0.00063 percentage points.
3. The `raw` sheet is a copy of the `samp_mn_sd` block (columns O-W) of `ELISA/docs/MethylFlashELISA_StdCrv_and_Data_Analysis_Feb2019.xlsx`, with ploidy and treatment columns added. Its absorbances equal the reading-1 values for the 62 sample wells (maximum difference 4.7e-15).
4. `code/5-mC & 5mCpG plots.Rmd` reads sheet `summary` (`code/5-mC & 5mCpG plots.Rmd:47`) and then:
   - runs `shapiro.test` on the pooled percent (`code/5-mC & 5mCpG plots.Rmd:84`);
   - sets response = `asin(sqrt(percent/100))` (`code/5-mC & 5mCpG plots.Rmd:89`);
   - fits `aov(response ~ ploidy * treatment)` (`code/5-mC & 5mCpG plots.Rmd:101`);
   - fits `aov(response ~ interaction(treatment, ploidy))` and runs `agricolae::HSD.test` (`code/5-mC & 5mCpG plots.Rmd:107-111`);
   - writes `stats/AOV_5-mC.csv` and `stats/HSD_test_5-mC.csv` (`code/5-mC & 5mCpG plots.Rmd:116-117`).
5. There is no outlier removal, no log and no Johnson transformation in this path. Duplicate wells are averaged upstream, in the spreadsheet.

### 1.2 Match to the manuscript (verified, recomputed)

| Term | Manuscript | `stats/AOV_5-mC.csv` (identical to the recomputation) |
|---|---|---|
| treatment | F(1,28)=13.85, p<0.001 | F=13.8511287, p=0.000881540 |
| ploidy | F(1,28)=3.93, p=0.057 | F=3.9254685, p=0.0574573 |
| interaction | F(1,28)=3.91, p=0.058 | F=3.9069827, p=0.0580107 |

The residual df is 28. `stats/HSD_test_5-mC.csv` gives heat.diploid "a" (0.19734) and "b" for heat.triploid (0.14828), control.diploid (0.12670) and control.triploid (0.12664).

- **Scale.** The "response" means in the HSD file are on the arcsine scale. Back-transformed they are 3.84, 2.18, 1.60 and 1.60%.
- **Tukey p-values.** heat.diploid vs heat.triploid 0.043, vs control.diploid 0.0021, vs control.triploid 0.0021; all other pairs at least 0.61.
- **Figure letters.** The letters in Figure 1A are the same grouping relabelled. They are typed text boxes in `figures/final_figures/Figure-1-methylation-percentage.pptx`, not code output.

### 1.3 The result depends on the transformation (verified)

| Analysis | ploidy | treatment | interaction |
|---|---|---|---|
| Untransformed (`ELISA/GlobalDNAMeth_Polyploids.R:194`, `ELISA/GlobalDNAMeth_Polyploids.md:267-294`) | F=5.02, p=0.0331 | F=12.63, p=0.00137 | F=4.84, p=0.0363 |
| Arcsine square root (manuscript) | F=3.93, p=0.0575 | F=13.85, p=0.00088 | F=3.91, p=0.0580 |
| Johnson (`ELISA/GlobalDNAMeth_Polyploids.R:210-229`) | not recorded | not recorded | not recorded |

- The older figure `ELISA/figures/5mC_figure.png` (its `.pptx` is dated 2021-03-30) prints the untransformed p-values: 0.033 for ploidy, 0.001 for treatment and 0.036 for the interaction.
- Untransformed Tukey p-values (`ELISA/GlobalDNAMeth_Polyploids.md:281-294`, reproduced): diploid heat vs triploid heat 0.0195; diploid heat vs diploid control 0.00188; diploid heat vs triploid control 0.00174; the other three pairs 0.77 or higher.
- The transformation is triggered by a Shapiro-Wilk test on the pooled raw response (W=0.8805, p=0.0020). On the residuals of the two-way model, on the raw scale, Shapiro-Wilk gives p=0.59. After the arcsine transformation the pooled test gives p=0.196.
- Raw-scale variances are heterogeneous: Brown-Forsythe p=0.043, and the largest group variance is 6.9 times the smallest (4.1 times after arcsine; Brown-Forsythe p=0.33).
- `ELISA/GlobalDNAMeth_Polyploids.R` would apply the Johnson transformation (same Shapiro p-value), but that result is not saved anywhere and could not be recomputed in the audit.

### 1.4 Samples (verified)

There are exactly 32 samples, 8 per group, all ctenidia. No animals given the acute heat shock (09, 10, 19, 20) were assayed. Per-sample values are in Appendix A.

| Group | Samples | Mean %5-mC | SD | Median |
|---|---|---|---|---|
| Diploid control | D01-C to D08-C | 1.698 | 0.941 | 1.452 |
| Diploid heat | D11-C to D18-C | 4.019 | 1.773 | 3.980 |
| Triploid control | T01-C to T08-C | 1.681 | 0.853 | 1.627 |
| Triploid heat | T11-C to T18-C | 2.227 | 0.676 | 2.175 |

### 1.5 Figure 1A drops three of the eight diploid heat samples (verified)

- `ylim(0.5, 5)` (`code/5-mC & 5mCpG plots.Rmd:60`) removes D14-C, D15-C and D16-C (5.21, 5.22 and 6.96%) before the boxplot statistics are computed. ggplot2 reports "Removed 3 rows containing non-finite values".
- The diploid heat box is therefore drawn from n=5:

| Statistic | Drawn in Figure 1A (n=5) | All data (n=8) |
|---|---|---|
| Lower hinge | 2.286 | 2.812 |
| Median | 2.987 | 3.980 |
| Upper hinge | 3.594 | 5.208 |
| Maximum | 4.367 | 6.958 |

- The same applies to `plots/methylation_boxplot.png` and `plots/methylation_combined_facet.png`. Figure 1B is unaffected: no percent mCpG value lies outside its limits of 9 to 12.5 (`code/5-mC & 5mCpG plots.Rmd:130`).
- The ANOVA uses all 8 values, so the statistics and the figure describe different data.
- The older figure `ELISA/figures/5mC_figure.png` uses `limits = c(0, 3)` (`ELISA/GlobalDNAMeth_Polyploids.R:155`) and drops 9 of 32 values, 5 of them diploid heat: D02-C 3.616, D12-C 3.594, D13-C 4.367, D14-C 5.205, D15-C 5.218, D16-C 6.958, T04-C 3.387, T12-C 3.004 and T16-C 3.157. Its diploid heat box is drawn from n=3 (median 2.29).
- The earlier knit figure `ELISA/GlobalDNAMeth_Polyploids_files/figure-markdown_github/unnamed-chunk-14-1.png` shows the full data.

### 1.6 Manual, non-code steps (verified)

- `data/ELISA-percent-mCPG-summary.xlsx` is hand-built: typed constants and group labels, and pasted values in `summary`. No code writes it.
- The plate map and its exclusion flags are typed in `ELISA/docs/20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx`.
- The Figure 1 letters and axis labels are typed in PowerPoint.

### 1.7 Open points

- **Per-sample DNA input (undetermined).** The analysis assumes 25 ng for every sample. The sheet `Brconc.` of `ELISA/docs/MethylFlashELISA_StdCrv_and_Data_Analysis_Feb2019.xlsx` lists a per-sample "concentration" (2.0 to 4.4; exactly 2.0 for 15 of 32, including all 8 diploid controls) and "DNA added" = concentration x 5 (10 to 22 ng). The repository does not say what these values are.
  - If they are the true inputs and the input is set per sample (slope 0.6105), the recomputation gives treatment p=0.049 raw and 0.036 arcsine, ploidy p=0.12 to 0.13, and interaction p=0.27 to 0.39.
  - The correlation between "DNA added" and the published %5-mC is weak (r=0.10).
- **Absolute scale (inference).** The plate map labels the standards "5 ul". The kit manual `ELISA/docs/MethylFlashGlobalDNAELISA_P-1030.pdf` specifies 2 ul (100 ng) of positive control per well and says the control mass must equal the sample mass. The absolute %5-mC scale is therefore uncertain by a constant factor. A constant factor does not change the untransformed F tests, but it does change the arcsine-transformed ones and the reported levels.
- **Extrapolation (verified).** Four diploid heat samples (D13-C to D16-C, mean OD 0.89 to 1.29) lie above the top standard used in the fit (OD 0.826). The sheet `polynom` of the same workbook holds an unused polynomial fit.
- **`contam` wells (undetermined).** Why A03, B03 and C03 were flagged is not documented. The only well of D02-C (D03) is adjacent to them and is the diploid control outlier (3.62%).

## 2. Percent CpG methylation (N=22, Figure 1B)

### 2.1 Source of the values (verified)

- The values are `data/bismark-percent-mCPG-summary.xlsx`, sheet `raw`, column H `%mCpG` = `ROUND(F/SUM(F:G)*100,2)`, where F is `methylated_in_CpG` and G is `unmethylated_C_in_CpG`.
- All 44 typed counts equal lines 24 and 29 of the corresponding report in `data/PE_reports/` ("Total methylated C's in CpG context" and "Total unmethylated C's in CpG context").
- The 2-decimal values agree with the reports' printed 1-decimal "C methylated in CpG context" (line 34) for 22 of 22 libraries (maximum difference 0.05). They are not literally the same numbers.
- The counts come from the alignment report, so they are before deduplication.

### 2.2 Script, model and match to the manuscript (verified, recomputed)

- `code/5-mC & 5mCpG plots.Rmd:146-189` runs `shapiro.test` on the pooled %mCpG (`code/5-mC & 5mCpG plots.Rmd:156`): W=0.868, p=0.0073.
- That triggers `asin(sqrt(x/100))` (`code/5-mC & 5mCpG plots.Rmd:161`), then `aov(response ~ ploidy * treatment)` (`code/5-mC & 5mCpG plots.Rmd:174`), then the HSD test (`code/5-mC & 5mCpG plots.Rmd:179-183`) and `write.table` (`code/5-mC & 5mCpG plots.Rmd:188-189`).

| Term | Manuscript | `stats/AOV_%mCpG.csv` (identical to the recomputation) |
|---|---|---|
| interaction | F(1,18)=14.49, p=0.0013 | F=14.486105, p=0.00129380 |
| ploidy | p=0.86 | F=0.0313, p=0.8615 |
| treatment | p=0.42 | F=0.6767, p=0.4215 |

- HSD letters (`stats/HSD_%mCpG.csv`): heat.diploid a; control.triploid ab; heat.triploid b; control.diploid b. The letters in Figure 1B are the same grouping relabelled.
- Tukey p-values: heat.diploid vs control.diploid 0.020; heat.diploid vs heat.triploid 0.041; control.triploid vs control.diploid 0.106; control.triploid vs heat.triploid 0.188; the other two pairs 0.78 or higher.
- Untransformed, the interaction is F=14.31, p=0.00136. After the transformation the pooled Shapiro p-value is still 0.0094, and the code proceeds regardless. On the residuals of the two-way model Shapiro-Wilk gives p=0.41 on the raw scale and p=0.44 after the transformation.
- `code/1_WGBS_PE_report_analysis.R` is not the source: its input is missing at the referenced path, it has n=10, and its models are `CpG ~ ploidy` and `CpG ~ ploidy + heat_shock` (`code/1_WGBS_PE_report_analysis.R:87`, `code/1_WGBS_PE_report_analysis.R:155`).

### 2.3 Groups and what the Bismark reports reveal (verified)

Group means (SD, n): diploid control 10.085 (0.321, 6); diploid heat 11.006 (0.591, 5); triploid control 10.737 (0.619, 6); triploid heat 10.142 (0.165, 5). The per-library table is Appendix B.

Lines 1 to 3 of each report in `data/PE_reports/` give:

- **Bismark version:** v0.21.0 for the 9 files with "fastp-trim.20201202" in the name; v0.22.3 for the other 13.
- **Bowtie2 options:** `-q --score-min L,0,-0.6 -p 8 --reorder --ignore-quals --no-mixed --no-discordant --dovetail --maxins 500` (v0.21.0); the v0.22.3 runs use `-p 4`.
- **Strand setting:** `--non_directional` in all 22, against the genome folder `/gscratch/srlab/sr320/data/Cgig-genome/roslin_M/`.
- **Inputs:** read from `.../data/cg/` for the fastp-named files and `.../data/cg-big/` for the others.
- **With "fastp-trim.20201202" in the name:** zr3534_1, 2, 3, 4, 5, 6, 8, 9, 10. **Without:** zr3534_7 and all 12 zr3644 libraries. The missing suffix does not mean untrimmed reads (section 2.5).

### 2.4 Problems affecting interpretation

1. **Treatment is completely confounded with sequencing project (verified).** All 10 heat libraries are project zr3534 and all 12 controls are project zr3644.
   - Depth differs between projects: 66 to 91 million read pairs (zr3534) versus 36 to 58 million (zr3644).
   - The duplicate rate differs: 16.5 to 20.6% of aligned pairs (zr3534) versus 10.2 to 12.6% (zr3644), from the June 2021 run (section 2.5).
   - The interaction is a pair of opposite-sign ploidy differences in two projects: in zr3534, diploid minus triploid = +0.864 (t(8)=3.15, p=0.014); in zr3644, -0.652 (t(10)=-2.29, p=0.045).
2. **The controls belong to a different dataset.**
   - Verified: `metadata.xlsx`, sheet `list`, column G marks zr3534 `ronit` and zr3644 `yaamini`. Sheet `all` labels zr3644_7 to 12 and 19 to 24 `OA`. `code/2_WGBS_Methylkit-all-gigas.R:71` downloads from a `project-gigas-oa-meth` directory. `data/bismark-percent-mCPG-summary.xlsx` has no library names for zr3644 (cells B12:B23 are blank).
   - Inference: the 12 controls are the control arm of a separate ocean-acidification project, not animals D01 to D08 and T01 to T08 of this experiment.
   - Undetermined: whether those animals share cohort, date or handling with this experiment.
3. **Four library names are wrong and the heat level is mixed (verified).**
   - `gene_expression/data/zr3534_wgbs_info.csv:5-6` and `gene_expression/data/zr3534_wgbs_info.csv:10-11` give zr3534_4, 5, 9 and 10 as D19-C, D20-C, T19-C and T20-C: heat plus the acute heat shock. The spreadsheet names them D14-C, D15-C, T14-C and T15-C (cells B5, B6, B10, B11).
   - `metadata.xlsx`, sheet `list3` (E5:E6, E10:E11), and the `heat_shock` column of the older workbook `old_WGBS/WGBS/data/PE_reports/percent_methylation_summary.xlsx` (F5:F6, F10:F11) agree with the sample sheet. Only the typed names disagree. The D11, D12, D13, D14, D15 pattern suggests spreadsheet autofill (inference).
   - The heat level in Figure 1B is therefore 6 heat-only libraries and 4 heat plus acute heat shock libraries.
   - By the sample sheet, only 6 animals have both an ELISA and a WGBS value (D11-C, D12-C, D13-C, T11-C, T12-C, T13-C); across them the two measures correlate at r=0.74 (p=0.090).
4. **The spreadsheet mixes two alignment runs (verified, section 2.5).** Nine heat values come from the March 2021 run (Bismark 0.21.0) and the other 13 from the June 2021 run (Bismark 0.22.3). This applies only to the hand-built spreadsheet behind Figure 1B, not to the differential methylation analysis, and its numerical effect is negligible.
5. **Stale `summary` sheet (verified).** The sheet `summary` of `data/bismark-percent-mCPG-summary.xlsx` (10.63, 11.05, 10.07, 9.80) cannot be derived from `raw`. `code/5-mC & 5mCpG plots.Rmd` does not read it.

### 2.5 The alignment runs and trimming

**What the committed job script does (verified).** `old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh` performs no trimming. Despite "fastp" in its file name it contains no fastp command. It aligns reads that are already trimmed (`*_R1.fastp-trim.20201202.fq.gz`; `old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh:40-48`) with Bismark 0.21.0, bowtie2 2.3.4.1 and samtools 1.9 (`old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh:21-23`), then deduplicates, extracts methylation and builds coverage files. A condensed listing is in Appendix H. The fastp version and parameters are recorded nowhere in the repository; the only traces are the "fastp-trim.20201202" file-name suffix and a commented file name with a 20201206 suffix (`old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh:38`).

**How this script relates to the two upstream runs (verified against copies retrieved from the lab server).**

- The committed script is the job script of the March 2021 alignment run (server path `/panopea/030521-ronrosM/job.sh`, copy at `07_documentation/upstream_processing/030521-ronrosM_job.sh`). Every command line is identical. The two files differ only in white space: 31 blank lines carry a single space in the server copy, one line has a trailing space, and the server copy ends with two extra blank lines.
- A second, later run aligned all libraries together (server path `/panopea/061021-big/job.sh`, June to July 2021 per the server listing, copy at `07_documentation/upstream_processing/061021-big_job.sh`). It uses Bismark 0.22.3, bowtie2 2.3.4.1 and samtools 1.9 (`07_documentation/upstream_processing/061021-big_job.sh:21-23`), reads from a `cg-big` directory named `<sample>_R1.fq.gz` (`07_documentation/upstream_processing/061021-big_job.sh:24`), and the same `-score_min L,0,-0.6 --non_directional` settings with `-p 4` (`07_documentation/upstream_processing/061021-big_job.sh:60-68`).
- The two Bismark summary reports (`07_documentation/upstream_processing/030521-ronrosM_bismark_summary_report.txt` for March and `03_input_files/bismark_reports/061021-big_bismark_summary_report.txt` for June) show identical read-pair totals for all ten zr3534 libraries in both runs (10 of 10). The June run therefore used the same fastp-trimmed reads under shorter names.
- The read-pair totals of all 22 reports in `data/PE_reports/` equal the June-run totals (22 of 22).
- The repository itself points to the June run: the download log in `code/4_WGBS_GOterm-Annotation-BOTH.ipynb:1178` and `code/4_WGBS_GOterm-Annotation-BOTH.ipynb:1439` lists a `061021-big` directory holding short-named reports for all ten zr3534 libraries as well as the zr3644 libraries.
- The March run did align zr3534_7 (66,080,050 pairs; `07_documentation/upstream_processing/030521-ronrosM_bismark_summary_report.txt:9`). The repository simply lacks its March report and holds the June one instead.
- The merged CpG coverage files used for all differential methylation analyses are the June (Bismark 0.22.3) outputs for all 22 libraries. This was checked by the wider audit against the server, and it is consistent with `metadata.xlsx`, sheet `list`, whose 22 paths all use the short June-run file names (the older `code/2_WGBS_Methylkit.R:58-67` used the March names).
- **Conclusion.** The mixed-version problem applies only to the hand-built percent mCpG spreadsheet behind Figure 1B, not to the differential methylation analysis.

**How much the mixing matters (verified from the two summary reports).**

- For the same library, the June run gives a deduplicated percent mCpG that is 0.010 to 0.013 percentage points higher than the March run (10 of 10 libraries), with 67,221 to 100,006 fewer aligned pairs (about 0.17%). Group differences in Figure 1B are 0.6 to 0.9 points.
- Deduplicated values are lower than the spreadsheet values by a mean of 0.179 points (zr3534) and 0.079 points (zr3644).
- **Robustness check.** Using the June-run deduplicated counts for all 22 libraries (one pipeline, after deduplication), the two-way ANOVA on the arcsine scale gives interaction F(1,18)=16.42, p=0.00075; ploidy p=0.87; treatment p=0.74. Untransformed, the interaction is F=16.16, p=0.00080. The HSD lettering is unchanged (heat.diploid a; control.triploid ab; control.diploid b; heat.triploid b). Group means are 10.007 (diploid control), 10.820 (diploid heat), 10.657 (triploid control) and 9.970 (triploid heat). Per-library values are in Appendix B.
- The Figure 1B statistic is therefore robust to the pipeline inconsistency. The design problems in section 2.4 remain.

**Still undetermined:** the fastp version and parameters, and whether the zr3644 reads were trimmed with the same settings (they were only in the June run, so no cross-run comparison is possible).

## 3. qPCR

### 3.1 What exists (verified)

- **Genes** in `gene_expression/data/qpcr_ct_values/qpcr_data_consolidated.csv` (32 rows x 13 columns): Actin, HSC70, DNMT1, MBD2, MeCP2, HIF1A, HATHaP2, HAT, HSP90, SOD, ATPsynthetase, COX1.
- **Gene-list discrepancies:**
  - SOD, ATPsynthetase and COX1 are in the consolidated file. They were added as columns in commits fb94b04 and ba4ccdd, and no per-gene replicate file exists for them.
  - `gene_expression/scripts/dct_analysis.R` handles 8 targets and runs no MBD2 ANOVA.
  - `gene_expression/scripts/dct_analysis_NA-to-45.R`, `ELISA/Facet_wrap_expression_vals.R` and `ELISA/dct_analysis_NA-to-45_COX1convert.Rmd` handle all 11.
- **Per-gene files** (`qpcr_data_actin.csv`, `_dnmt1`, `_hat`, `_hathap2`, `_hif1a`, `_hsc70`, `_hsp90`, `_mbd2`, `_mecp2` in `gene_expression/data/qpcr_ct_values/`):
  - 80 rows each: 40 animals x 2 technical replicates.
  - IDs look like `D_N_N01`: ploidy, treatment code (`N` control, `D` heat), then the acute heat shock code and animal number (`N01`, `H09`). They include the 8 animals given the acute heat shock (H09, H10, H19, H20 per ploidy), which the consolidated file omits.
  - Non-detects are the string `N/A`. Counts per 80 wells: MeCP2 60, HIF1A 28, MBD2 23, HAT 17, DNMT1 15, HSP90 3, HATHaP2 2, Actin 1, HSC70 0.
  - The actin header is "Sample " with a trailing space.
- **Consolidated file:**
  - IDs look like `D_C_N_01`: ploidy, treatment code (`C` control, `D` heat), acute heat shock code (always `N`), animal number.
  - Each value is the mean of the available technical replicates, verified for all 9 genes x 32 animals. A single replicate is used when the other is `N/A`.
  - Missing is `NA` (44 cells; MeCP2 18 of 32).
  - No code builds this file. Commit 1b3e175 ("replaced #DIV0 with NA") indicates it was made in a spreadsheet.
- **Method:**
  - Actin is the only normalizer.
  - Expression is 2^-(Cq target - Cq Actin) per animal, with no calibrator, then the natural log.
  - Per gene, `aov(<gene>log ~ Ploidy + Desiccation + Ploidy:Desiccation)` plus `TukeyHSD` (`gene_expression/scripts/dct_analysis.R:44-77`; `gene_expression/scripts/dct_analysis_NA-to-45.R:55-108`; `ELISA/dct_analysis_NA-to-45_COX1convert.Rmd:60-113`).
  - There is no correction across genes.
- **Non-detects:**
  - `gene_expression/scripts/dct_analysis.R` leaves `NA`, so `aov` drops those animals.
  - The "NA to 45" scripts set every `NA` to 45 before normalization (`gene_expression/scripts/dct_analysis_NA-to-45.R:24`; `ELISA/dct_analysis_NA-to-45_COX1convert.Rmd:26`; `ELISA/Facet_wrap_expression_vals.R:25`).
- **Figures:**
  - `gene_expression/DNMT1_ploidy.tiff` and `gene_expression/HSP90_ploidy.tiff` (`gene_expression/scripts/dct_analysis.R:122`, `gene_expression/scripts/dct_analysis.R:225`).
  - Eight PNGs in `gene_expression/analyses/` (DNMT1, HATHaP2, HAT, HIF1A, HSC70, HSP90, MBD2, MeCP2). No code writes them; they use an older naming scheme.
  - Eleven knit PNGs under `ELISA/dct_analysis_NA-to-45_COX1convert_files/figure-markdown_github/`.
  - `ELISA/Facet_wrap_expression_vals.R` plots to screen only.

### 3.2 Mitochondrial abundance (verified)

Nothing in the repository quantifies COX1 relative to a nuclear gene.

- COX1 is used as a raw Cq (`gene_expression/scripts/dct_analysis_NA-to-45.R:51` takes `log()` of the raw Cq; plot at `gene_expression/scripts/dct_analysis_NA-to-45.R:225-234`) or as 2^-Cq (`ELISA/dct_analysis_NA-to-45_COX1convert.Rmd:39`), with the axis label "COX1 DNA Abundance" (`ELISA/dct_analysis_NA-to-45_COX1convert.Rmd:236`).
- Undetermined: the template (DNA or cDNA), primers, input mass and efficiencies.
- COX1 Cq (20.845 to 22.440) is uncorrelated with Actin Cq across animals (r=0.002). This suggests a different template preparation (inference).

### 3.3 Do the data support the Discussion statements?

All values were computed from the consolidated file. dCq is Cq(HSP90) minus Cq(Actin) per animal; relative expression is 2^-dCq; the geometric mean is 2^-(mean dCq).

**HSP90**, with `NA` set to 45, n=8 per group. Complete-case values are in parentheses where they differ.

| Group | Mean dCq | Geometric mean | Median | Arithmetic mean |
|---|---|---|---|---|
| Diploid control | 8.394 (8.201, n=7) | 0.00297 (0.00340) | 0.00317 | 0.00345 (0.00378) |
| Diploid heat | 6.982 | 0.00791 | 0.00563 | 0.01818 |
| Triploid control | 7.133 | 0.00713 | 0.00746 | 0.00845 |
| Triploid heat | 6.118 | 0.01440 | 0.00920 | 0.02752 |

- **`NA` set to 45** (`ELISA/dct_analysis_NA-to-45_COX1convert.md:436-440`, reproduced): ploidy F(1,28)=4.90, p=0.035 (triploid over diploid = 2.09-fold); treatment F=6.39, p=0.017 (heat over control = 2.32-fold); interaction p=0.68. Tukey, triploid control vs diploid control: p=0.27.
- **Complete case** (the analysis in `gene_expression/scripts/dct_analysis.R`): ploidy F(1,27)=3.61, p=0.068; treatment F=5.22, p=0.030; interaction p=0.84. Tukey, triploid control vs diploid control: p=0.44.
- **Control-only unadjusted t-test** (computed for the audit): p=0.016 with `NA` set to 45; p=0.033 complete case.
- **The single non-detect** is D_C_N_05, whose Actin Cq is 35.26 (a failed sample). Imputing 45 for it is what moves the ploidy p-value from 0.068 to 0.035.
- **Assessment.** A roughly 2-fold higher HSP90 in triploids is consistent across variants, but its significance depends on the analysis. The better-supported effect is induction by heat in both ploidies, with no interaction. The data do not show a pattern that is independent of the heat treatment and specific to triploids.

**COX1 raw Cq**, n=8 per group:

| Group | Mean | SD | Median |
|---|---|---|---|
| Diploid control | 21.543 | 0.253 | 21.538 |
| Diploid heat | 21.358 | 0.499 | 21.198 |
| Triploid control | 21.356 | 0.272 | 21.423 |
| Triploid heat | 21.644 | 0.339 | 21.605 |

- Change with heat: diploid -0.185 cycles (2^-Cq x 1.14); triploid +0.288 cycles (x 0.82).
- ANOVA on log(2^-Cq) (`ELISA/dct_analysis_NA-to-45_COX1convert.md:652-656`, reproduced): ploidy p=0.69, treatment p=0.69, interaction F(1,28)=3.56, p=0.070.
- All Tukey pairs have p of 0.38 or higher. Within-ploidy t-tests (computed for the audit): diploid p=0.37, triploid p=0.083.
- **Assessment.** The direction matches the Discussion statement, but nothing is significant, the shifts are under 0.3 cycles, and no nuclear reference is used.

## 4. Bugs and statistical problems

**Directly affecting a manuscript number or figure**

1. **Percent mCpG design.** Treatment is fully confounded with sequencing project and dataset; the controls come from another dataset (section 2.4, items 1 and 2). The interaction F(1,18)=14.49 cannot be separated from a between-project difference.
2. **Percent mCpG labels.** Four libraries from animals given the acute heat shock carry wrong names and are pooled into the heat level; 12 libraries are unnamed (section 2.4, item 3).
3. **Figure 1A drops data.** `ylim(0.5, 5)` removes three diploid heat samples before the boxplot is computed (section 1.5).
4. **Transformation choice decides significance.** ELISA ploidy and interaction go from p=0.033 and 0.036 untransformed to 0.057 and 0.058 after arcsine (section 1.3).
   - Three variants coexist in the repository: untransformed, Johnson and arcsine.
   - The normality test is applied to the pooled response rather than to the model residuals.
   - The HSD "response" means are on the transformed scale and unlabelled.
5. **qPCR imputation.** With `NA` set to 45 and actin normalization, the failed sample D_C_N_05 becomes the highest "expression" value for DNMT1 and HAT although both are non-detects. It also drives the HSP90 ploidy p-value (section 3.3).
6. **qPCR quality.**
   - Actin Cq spans 21.9 to 35.3 and tends to differ by ploidy (p=0.078).
   - Five targets have a median Cq of 37 to 38.
   - MeCP2 is 75% non-detect at well level.
7. **COX1.** It is unnormalized, and `gene_expression/scripts/dct_analysis_NA-to-45.R:51` analyses `log(raw Cq)`.
8. **ELISA input assumptions.** A constant 25 ng input versus the `Brconc.` sheet; the unknown DNA mass in the standard wells; four samples extrapolated above the fitted range (section 1.7).
9. **Percent mCpG spreadsheet mixes two alignment runs.** Negligible numerically, and limited to Figure 1B (section 2.5).

**Code that cannot run as committed**

10. `gene_expression/scripts/dct_analysis.R` does not parse: there are git conflict markers at `gene_expression/scripts/dct_analysis.R:160-164` (R reports "unexpected input" at line 160).
11. `ELISA/GlobalDNAMeth_Polyploids.R:231-233` references undefined objects and stops with an error.
12. Paths and environment:
    - absolute paths to other repositories (`ELISA/GlobalDNAMeth_Polyploids.Rmd`, `ELISA/GlobalDNAMeth_Polyploids_pilot.Rmd`, `ELISA/dct_analysis_NA-to-45_COX1convert.Rmd`);
    - a missing pilot plate map;
    - an rstudioapi requirement (`ELISA/GlobalDNAMeth_Polyploids.R`);
    - dependence on the working directory (`code/5-mC & 5mCpG plots.Rmd`, both `gene_expression/scripts/` files, `ELISA/Facet_wrap_expression_vals.R`);
    - unconditional `install.packages` calls in the qPCR scripts; `code/5-mC & 5mCpG plots.Rmd:20-22` installs any package on its list that is missing, including `bestNormalize`, which the code never uses.
13. **Package-version dependence.** `ELISA/GlobalDNAMeth_Polyploids.Rmd:182` prints `glance(aov)$p.value`. With broom 1.0.5 `glance()` on an `aov` returns no `p.value` column, so the line that printed 0.0007834074 in 2019 (R 3.4.4, `ELISA/GlobalDNAMeth_Polyploids.md:12`, `ELISA/GlobalDNAMeth_Polyploids.md:258`) now returns nothing.

**Other**

14. **ELISA plate layout.** Groups occupy contiguous strips (3-4 diploid control, 5-6 diploid heat, 7-8 triploid control, 9-10 triploid heat; `ELISA/docs/20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx`, sheet `PlateMap`, rows 15-23), on a single plate. Position is confounded with group.
15. **Labels and gaps in `gene_expression/scripts/dct_analysis.R`.**
    - X-axis labels are swapped for the current CSV coding: `labels=c("Desiccated + Elevated Temp.","Control")` at lines 92, 138, 149, 161/163, 175, 186 and 197, while `C` (control) sorts first. They were correct before commit 6bffd5d changed the coding. The two saved TIFFs are unaffected.
    - The MBD2 plot is duplicated (`gene_expression/scripts/dct_analysis.R:130-150`).
    - There is no MBD2 ANOVA (`gene_expression/scripts/dct_analysis.R:44-77` covers seven genes).
16. **qPCR plots drop data.** Every qPCR boxplot with `ylim` drops observations.
    - `ELISA/dct_analysis_NA-to-45_COX1convert.md:662-698` reports the removed rows: HSC70 2, DNMT1 1, MBD2 4, MeCP2 3, HIF1A 3, HATHaP2 3, HAT 4, HSP90 5, SOD 4, ATPsynthetase 2.
    - HSP90 loses 5 values, of which 3 are triploid heat.
    - `gene_expression/HSP90_ploidy.tiff` (`gene_expression/scripts/dct_analysis.R:202`, limits 0 to 0.035) drops the three highest values (D_D_N_11 0.1015, T_D_N_11 0.0998, T_D_N_12 0.0602) and the non-detect.
    - `gene_expression/scripts/dct_analysis_NA-to-45.R:233`, `ylim(c(21,22))`, drops 4 COX1 values.
    - The axis label names the quantity as an expression delta Ct, but the values plotted are 2^-dCq.
17. **Label ambiguity.** `heat` and `heat_stress` mean three different things across files: the 27 C treatment (ELISA scripts and both data workbooks), the 45 C acute heat shock (`ELISA/GlobalDNAMeth_Polyploids_pilot.Rmd`), and a mixture of the two (zr3534 libraries labelled `heat`).
18. **Housekeeping.**
    - The `stats/*.csv` files are space-delimited `write.table` output, not CSV.
    - Percent mCpG is rounded to 2 decimals before analysis; unrounded count ratios give F=14.39 (p=0.00133).
    - Office lock files (`~$...xlsx`) are committed.
19. **Checked and not found.** There is no pseudo-replication from technical replicates. ELISA duplicate wells and qPCR duplicates are averaged per animal before every ANOVA.

## Appendix A. Per-sample ELISA values

Mean OD is from reading 1 (`ELISA/docs/20190206_FirstReadingSTRIGG.xls`). %5-mC is the value in sheet `summary` of `data/ELISA-percent-mCPG-summary.xlsx`, which the manuscript statistics use.

| Sample | Group | Wells | Mean OD | SD of duplicates | %5-mC |
|---|---|---|---|---|---|
| D01-C | diploid control | 2 | 0.3981 | 0.0304 | 1.141 |
| D02-C | diploid control | 1 | 0.7760 | n/a | 3.616 |
| D03-C | diploid control | 2 | 0.3852 | 0.0400 | 1.056 |
| D04-C | diploid control | 2 | 0.5105 | 0.0688 | 1.877 |
| D05-C | diploid control | 2 | 0.4846 | 0.0066 | 1.707 |
| D06-C | diploid control | 2 | 0.4065 | 0.0138 | 1.196 |
| D07-C | diploid control | 2 | 0.3233 | 0.0145 | 0.650 |
| D08-C | diploid control | 2 | 0.5808 | 0.0459 | 2.338 |
| D11-C | diploid heat | 2 | 0.6798 | 0.0444 | 2.987 |
| D12-C | diploid heat | 2 | 0.7725 | 0.0604 | 3.594 |
| D13-C | diploid heat | 2 | 0.8905 | 0.0122 | 4.367 |
| D14-C | diploid heat | 2 | 1.0184 | 0.0993 | 5.205 |
| D15-C | diploid heat | 2 | 1.0205 | 0.0393 | 5.218 |
| D16-C | diploid heat | 2 | 1.2859 | 0.1859 | 6.958 |
| D17-C | diploid heat | 2 | 0.5730 | 0.0754 | 2.286 |
| D18-C | diploid heat | 2 | 0.4589 | 0.0489 | 1.539 |
| T01-C | triploid control | 2 | 0.4926 | 0.0495 | 1.760 |
| T02-C | triploid control | 2 | 0.3752 | 0.0413 | 0.991 |
| T03-C | triploid control | 2 | 0.3619 | 0.1030 | 0.904 |
| T04-C | triploid control | 2 | 0.7410 | 0.1972 | 3.387 |
| T05-C | triploid control | 2 | 0.5099 | 0.0162 | 1.873 |
| T06-C | triploid control | 2 | 0.4519 | 0.0250 | 1.493 |
| T07-C | triploid control | 2 | 0.5620 | 0.0510 | 2.214 |
| T08-C | triploid control | 1 | 0.3498 | n/a | 0.824 |
| T11-C | triploid heat | 2 | 0.6204 | 0.0029 | 2.597 |
| T12-C | triploid heat | 2 | 0.6825 | 0.0800 | 3.004 |
| T13-C | triploid heat | 2 | 0.4995 | 0.0386 | 1.805 |
| T14-C | triploid heat | 2 | 0.4268 | 0.0432 | 1.328 |
| T15-C | triploid heat | 2 | 0.6010 | 0.0365 | 2.470 |
| T16-C | triploid heat | 2 | 0.7059 | 0.0024 | 3.157 |
| T17-C | triploid heat | 2 | 0.5109 | 0.0146 | 1.880 |
| T18-C | triploid heat | 2 | 0.4644 | 0.0765 | 1.575 |

Standards (mean OD of two wells; `ELISA/GlobalDNAMeth_Polyploids.md:110-116`):

| Standard | Nominal %5-mC | Mean OD | Used in the fit |
|---|---|---|---|
| 5 ul NC | 0 | 0.2241 | yes |
| 5 ul 0.1 % PC | 0.1 | 0.2639 | yes |
| 5 ul 0.2 % PC | 0.2 | 0.3224 | yes |
| 5 ul 0.5% PC | 0.5 | 0.5024 | yes |
| 5 ul 1% PC | 1 | 0.8259 | yes |
| 5 ul 5% PC | 5 | 1.7441 | no |
| 5 ul PC | 10 | 2.1136 | no |

## Appendix B. Per-library percent mCpG and Bismark paired-end report table

"Label (xlsx)" is the name typed in `data/bismark-percent-mCPG-summary.xlsx`; "Label (sample sheet)" is the name in `gene_expression/data/zr3534_wgbs_info.csv`. Report percentages are as printed in `data/PE_reports/`. The last column is from the June 2021 summary report (`03_input_files/bismark_reports/061021-big_bismark_summary_report.txt`), after deduplication.

| seqID | Label (xlsx) | Label (sample sheet) | Ploidy | Treatment (xlsx) | fastp in file name | Bismark | Pairs analysed | Unique alignments | Map eff % | CpG / CHG / CHH % | %mCpG (xlsx) | %mCpG, June run, deduplicated |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zr3534_1 | D11-C | D11-C | diploid | heat | yes | v0.21.0 | 76,522,716 | 47,127,724 | 61.6 | 10.5 / 0.5 / 0.5 | 10.46 | 10.235 |
| zr3534_2 | D12-C | D12-C | diploid | heat | yes | v0.21.0 | 75,413,899 | 46,140,925 | 61.2 | 10.4 / 0.5 / 0.5 | 10.35 | 10.209 |
| zr3534_3 | D13-C | D13-C | diploid | heat | yes | v0.21.0 | 78,091,157 | 47,623,705 | 61.0 | 11.8 / 0.5 / 0.5 | 11.76 | 11.551 |
| zr3534_4 | D14-C | D19-C (acute heat shock) | diploid | heat | yes | v0.21.0 | 74,728,031 | 45,908,181 | 61.4 | 11.2 / 0.5 / 0.5 | 11.21 | 11.018 |
| zr3534_5 | D15-C | D20-C (acute heat shock) | diploid | heat | yes | v0.21.0 | 84,046,138 | 51,605,460 | 61.4 | 11.2 / 0.5 / 0.5 | 11.25 | 11.087 |
| zr3534_6 | T11-C | T11-C | triploid | heat | yes | v0.21.0 | 75,321,697 | 45,769,591 | 60.8 | 10.3 / 0.5 / 0.5 | 10.33 | 10.186 |
| zr3534_7 | T12-C | T12-C | triploid | heat | no | v0.22.3 | 66,080,050 | 40,302,825 | 61.0 | 10.1 / 0.5 / 0.5 | 10.14 | 9.930 |
| zr3534_8 | T13-C | T13-C | triploid | heat | yes | v0.21.0 | 90,605,411 | 55,543,581 | 61.3 | 10.3 / 0.5 / 0.5 | 10.28 | 10.081 |
| zr3534_9 | T14-C | T19-C (acute heat shock) | triploid | heat | yes | v0.21.0 | 73,776,588 | 45,591,578 | 61.8 | 10.0 / 0.5 / 0.5 | 10.01 | 9.837 |
| zr3534_10 | T15-C | T20-C (acute heat shock) | triploid | heat | yes | v0.21.0 | 77,576,026 | 48,007,674 | 61.9 | 10.0 / 0.4 / 0.5 | 9.95 | 9.813 |
| zr3644_1 | blank | not listed | diploid | control | no | v0.22.3 | 50,613,363 | 31,400,091 | 62.0 | 9.5 / 0.5 / 0.5 | 9.48 | 9.460 |
| zr3644_2 | blank | not listed | diploid | control | no | v0.22.3 | 54,652,876 | 33,337,869 | 61.0 | 10.0 / 0.5 / 0.5 | 10.04 | 10.030 |
| zr3644_3 | blank | not listed | diploid | control | no | v0.22.3 | 38,671,370 | 23,557,386 | 60.9 | 10.3 / 0.5 / 0.5 | 10.25 | 10.164 |
| zr3644_4 | blank | not listed | diploid | control | no | v0.22.3 | 49,072,984 | 29,921,502 | 61.0 | 10.1 / 0.5 / 0.5 | 10.08 | 9.947 |
| zr3644_5 | blank | not listed | diploid | control | no | v0.22.3 | 43,865,548 | 26,784,959 | 61.1 | 10.4 / 0.5 / 0.5 | 10.35 | 10.234 |
| zr3644_6 | blank | not listed | diploid | control | no | v0.22.3 | 52,883,485 | 32,436,346 | 61.3 | 10.3 / 0.5 / 0.5 | 10.31 | 10.204 |
| zr3644_13 | blank | not listed | triploid | control | no | v0.22.3 | 52,950,721 | 32,442,876 | 61.3 | 10.4 / 0.5 / 0.5 | 10.38 | 10.300 |
| zr3644_14 | blank | not listed | triploid | control | no | v0.22.3 | 36,184,563 | 22,328,832 | 61.7 | 11.8 / 0.5 / 0.5 | 11.78 | 11.615 |
| zr3644_15 | blank | not listed | triploid | control | no | v0.22.3 | 58,102,004 | 35,572,288 | 61.2 | 11.2 / 0.5 / 0.5 | 11.18 | 11.036 |
| zr3644_16 | blank | not listed | triploid | control | no | v0.22.3 | 50,134,077 | 30,903,797 | 61.6 | 10.5 / 0.5 / 0.5 | 10.49 | 10.410 |
| zr3644_17 | blank | not listed | triploid | control | no | v0.22.3 | 57,127,374 | 35,222,770 | 61.7 | 10.1 / 0.5 / 0.5 | 10.13 | 10.181 |
| zr3644_18 | blank | not listed | triploid | control | no | v0.22.3 | 50,675,924 | 31,440,847 | 62.0 | 10.5 / 0.5 / 0.5 | 10.46 | 10.402 |

## Appendix C. qPCR consolidated table

`gene_expression/data/qpcr_ct_values/qpcr_data_consolidated.csv`, verbatim. Sample IDs are ploidy (`D`, `T`), treatment code (`C` control, `D` heat), acute heat shock code (`N` none) and animal number.

```
Sample,Actin,HSC70,DNMT1,MBD2,MeCP2,HIF1A,HATHaP2,HAT,HSP90,SOD,ATPsynthetase,COX1
D_C_N_01,23.445,31.805,37.44,38.47,36.82,36.3,32.46,36.53,31.535,27.465,33.515,21.975
D_C_N_02,23.885,32.515,36.895,37.805,NA,NA,31.11,36.99,31.075,25.715,34.385,21.115
D_C_N_03,26.565,33.02,38.31,39.58,NA,38.53,34.46,NA,34.87,28.29,36.45,21.41
D_C_N_04,24.64,32.175,35.915,39.45,NA,38.355,35.195,38.47,33.905,28.425,34.62,21.42
D_C_N_05,35.255,37.19,NA,NA,39.16,NA,NA,NA,NA,NA,NA,21.715
D_C_N_06,28.15,34.505,38.89,NA,37.93,39.67,36.015,NA,36.45,31.655,38.955,21.595
D_C_N_07,27.24,32.7,37.99,39.94,NA,37.98,35.08,39.8,34.68,30.03,35.795,21.48
D_C_N_08,24.255,31.96,35.935,38.02,NA,36.315,34.565,38.125,33.07,28.345,34.89,21.63
D_D_N_11,30.015,36.215,NA,39.22,NA,NA,31.96,39.88,33.315,27.265,39.22,21.155
D_D_N_12,24.205,31.59,36.71,37.915,37.4,37.59,30.205,36.315,30.55,26.515,34.29,20.845
D_D_N_13,25,32.1,36.08,38.37,37.56,36.81,32.93,38.77,33.05,28.99,35.675,21.705
D_D_N_14,25.475,32.25,37.515,NA,39.29,39.36,33.815,37.44,32.115,28.06,36.06,21.22
D_D_N_15,25.125,31.745,35.83,37.84,38.58,38.45,33.555,38.015,32.935,27.325,35.205,21.25
D_D_N_16,31.465,36.47,NA,NA,NA,39.65,38.865,NA,39.81,31.42,NA,21.07
D_D_N_17,24.88,32.68,35.845,39.17,NA,38.09,33.045,38.28,33.045,28.035,34.975,21.175
D_D_N_18,25.025,31.295,37.155,38.195,39.42,39.15,34.165,37.64,32.225,27.83,35.51,22.44
T_C_N_01,23.87,31.06,37.925,37.105,36.49,37.985,31.26,35.395,29.625,25.27,33.355,21.56
T_C_N_02,21.875,30.315,38.035,36.75,NA,37.345,28.7,34.085,28.335,24.82,32.58,21.14
T_C_N_03,23.625,30.765,35.445,38.075,NA,37.025,32.46,36.525,30.37,26.325,33.555,20.94
T_C_N_04,22.43,30.1,34.87,36.58,NA,34.565,31.86,36.31,31.525,26.035,32.595,21.63
T_C_N_05,28.29,34.58,NA,NA,NA,39.67,34.71,39.2,35.19,29.7,37.765,21.605
T_C_N_06,24.47,30.965,36.215,39.07,35.43,37.275,33.385,37.8,31.8,26.745,35.55,21.105
T_C_N_07,25.005,32.53,37.715,38.605,NA,37.88,33.11,38.765,32.525,27.315,35.35,21.585
T_C_N_08,23.945,29.8,35.37,37.505,NA,36.91,33.13,36.58,31.2,26.91,33.685,21.285
T_D_N_11,25.565,32.07,NA,37.225,NA,39.045,32.945,38.05,28.89,33.66,37.555,21.54
T_D_N_12,25.95,33.4,38.855,38.02,38.03,NA,31.005,39.325,30.005,26.05,37.13,21.27
T_D_N_13,26.455,32.865,37.49,38.815,37.89,38.995,33.285,39.605,33.135,27.805,37.125,21.61
T_D_N_14,26.95,33.275,37.375,39.66,NA,NA,33.71,39.44,32.275,27.66,37.1,21.6
T_D_N_15,25.84,33.1,37.28,39.59,39.69,NA,34.125,38.665,32.695,28.795,36.045,21.215
T_D_N_16,24.91,32.74,37.355,38.26,NA,38.57,34.38,37.925,32.94,27.47,36.245,22.245
T_D_N_17,24.365,31.795,36.76,NA,NA,37.55,33.09,39.34,32.165,26.41,34.79,21.965
T_D_N_18,24.125,30.69,35.66,37.42,36.53,37.215,33.34,36.15,30.995,26.42,33.9,21.705
```

## Appendix D. Catalog of `ELISA/docs/`

| File | Sheets (dimensions) and headers | Contents |
|---|---|---|
| `20190206_FirstReadingSTRIGG.xls`, `20190206_2ndReadingSTRIGG.xls`, `20190206_3rdReadingSTRIGG.xls` | `plate list` (A1:F97): Plate, Repeat, Well, Type, Time, Absorbance @ 450 (1.0s) (A). `Plate` (A1:L15). `Protocol` (A1:A37). `Notes` (A1:A4) | **Raw plate absorbance for the N=32 run**: 96 wells, measured 2/6/2019 at 7:08:55, 7:12:19 and 7:15:48 PM. `Plate` is the 8 x 12 grid, `Protocol` the reader settings. No sample IDs. Only the First reading is used |
| `20181018_FirstReadingSTRIGG.xls`, `20181018_2ndReadingSTRIGG.xls`, `20181018_3rdReadingSTRIGG.xls` | `List ; Plates 1 - 1` (A1:F97, same headers), `Plate`, `Protocol`, `Notes` | Pilot plate of 10/18/2018; raw absorbance; no sample IDs |
| `20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx` | `Concentrations` (declared A1:X1004; data A1:B33): Sample, DNA Concentration (ng/uL). `PlateMap` (A1:L23). `PlateMap_list` (A1:C97): Well, Type, Sample | Stock DNA concentrations for the 32 samples (6.48 to 452 ng/uL). Planned and "ACTUAL" plate layouts. **Well-to-sample map**: IDs D01-C to D08-C, D11-C to D18-C, T01-C to T08-C, T11-C to T18-C, seven standards, and the flags `no_sample`, `contam`, `no_ab` |
| `MethylFlashELISA_StdCrv_and_Data_Analysis_Feb2019.xlsx` | 8 sheets, listed in the next table | **Standard curve and per-sample %5-mC** (exploratory workbook) |
| `MethylFlashELISA_StdCrv_and_Data_Analysis.xlsx` | `DataAndStdCrv` (A1:R41); `MethylationCalcs` (A2:U46) | Pilot. Plate map, three readings and control tables; notes "NC OD = 0.203, slope = 0.285, ng sample DNA = 100"; percentage table for D1, D2, D9, D10, T1, T2, T9, T10, 4G, 4C, 4Ms, 4M, two sea-lice samples and two blanks. Single wells, mantle tissue |
| `MethylFlashGlobalDNAELISA_P-1030.pdf` | not a spreadsheet | Kit manual, with typed margin notes about the pilot run |
| `~$20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx` | not a spreadsheet | 165-byte Office lock file |

Sheets of `MethylFlashELISA_StdCrv_and_Data_Analysis_Feb2019.xlsx`:

| Sheet | Dimensions | Headers and contents |
|---|---|---|
| `Sheet1` | A1:D97 | Well, Abs, Type, Sample: all 96 wells of reading 1 with plate-map labels |
| `stdcrv` | A1:L15 | Per-well standards plus mean and SD formulas; 10 charts. NC mean 0.2241 |
| `polynom` | A1:AC71 | Quadratic-fit percentage (a=-0.0751, b=0.6845); not used downstream |
| `practice` | A1:U140 | Worked polynomial fits and four linear fits with typed slopes 0.6105, 0.1863, 0.4941, 0.7456 |
| `example` | A1:H4 | Three-row polynomial example |
| `samp_mn_sd` | A1:AD80 | Well, Abs, Type, Sample, mean, sd, numerator, denominator, perc. Columns O-W use slope 0.6105 and are **the per-sample %5-mC values used in the manuscript**; columns A-I use slope 0.1863. Group mean and SD at AB5:AD9 |
| `compare` | A1:AB63 | February versus pilot comparisons |
| `Brconc.` | A1:I63 | Abs, Sample, mean, sd, concentration, DNA added, numerator, denominator, perc (see section 1.7) |

The pilot analysis (`ELISA/GlobalDNAMeth_Polyploids_pilot.Rmd`) uses single wells, mantle tissue, 100 ng input and n=2 per group; in it `heat_stress` means the acute heat shock applied to control animals (D9, D10, T9, T10). Its own text notes that an ANOVA is not appropriate with two points per group (`ELISA/GlobalDNAMeth_Polyploids_pilot.Rmd:216`).

## Appendix E. Per-script inventory

"Exists" was checked at commit 523d734. No lock file or `sessionInfo()` output exists for any of these scripts.

| Script | Path issues, `setwd`, rstudioapi | Packages | Inputs (exist at the referenced path?) | Outputs (in the repository?) |
|---|---|---|---|---|
| `ELISA/GlobalDNAMeth_Polyploids.R` | `rm(list=ls())` (line 7); rstudioapi plus `setwd` (lines 17-19; fails outside RStudio); `setwd` at lines 25, 178, 179 | readxl, ggplot2, tidyverse, RColorBrewer, broom, rstudioapi, Johnson, agricolae, nlme, multcomp | Three `20190206_*ReadingSTRIGG.xls`, sheet "plate list" (lines 28-30): yes. Plate map, sheet "PlateMap_list" (line 46): yes | `ELISA/figures/[boxplot]5mC_comparison.tiff` (line 182): no |
| `ELISA/GlobalDNAMeth_Polyploids.Rmd` | Home-relative absolute paths to another repository (lines 18-20, 46: `~/Documents/GitHub/C_gigas/Polyploids/docs/...`) | readxl, dplyr, ggplot2, broom | Do not resolve here; same-named files are in `ELISA/docs/` | `ELISA/GlobalDNAMeth_Polyploids.md` plus 4 PNGs: yes |
| `ELISA/GlobalDNAMeth_Polyploids_pilot.Rmd` | Same absolute-path pattern (lines 17-19, 46) | readxl, dplyr, ggplot2, broom | `20181018_*.xls`: exist in `ELISA/docs/`. `MethylFlashGlobalDNAMethylation_20181017.xlsx` (line 46): nowhere in the repository | `ELISA/GlobalDNAMeth_Polyploids_pilot.md` plus 5 PNGs: yes |
| `ELISA/Facet_wrap_expression_vals.R` | Unconditional `install.packages` x 4 (lines 4-7) | plyr, ggplot2, splitstackshape, tidyr | `data/qpcr_ct_values/qpcr_data_consolidated.csv` (line 16): exists only under `gene_expression/` | none |
| `ELISA/dct_analysis_NA-to-45_COX1convert.Rmd` | Absolute path (line 17: `~/Documents/GitHub/gigas_ploidy_ronits/...`) | plyr, ggplot2, splitstackshape | Does not resolve; the consolidated file in `gene_expression/data/qpcr_ct_values/` reproduces its output | `ELISA/dct_analysis_NA-to-45_COX1convert.md` plus 11 PNGs: yes |
| `code/5-mC & 5mCpG plots.Rmd` | `rm(list=ls())` (line 14); auto-installs missing packages (lines 21-22); no `setwd` or rstudioapi | readxl, tidyverse, scales, bestNormalize, agricolae, nlme, multcomp, patchwork (line 20; scales, bestNormalize, nlme and multcomp are unused) | `data/ELISA-percent-mCPG-summary.xlsx`, sheet "summary" (line 47), and `data/bismark-percent-mCPG-summary.xlsx`, sheet "raw" (line 48): exist from the repository root, not from `code/` | 3 PNGs in `plots/` (lines 66, 136, 209) and 4 files in `stats/` (lines 116-117, 188-189): yes. No knitted output |
| `gene_expression/scripts/dct_analysis.R` | Unconditional `install.packages` x 3 (lines 4-6) | plyr, ggplot2, splitstackshape | Relative CSV path (line 14): resolves only with the working directory `gene_expression/` | `gene_expression/DNMT1_ploidy.tiff` (line 122), `gene_expression/HSP90_ploidy.tiff` (line 225): yes |
| `gene_expression/scripts/dct_analysis_NA-to-45.R` | Unconditional `install.packages` x 3 (lines 5-7) | plyr, ggplot2, splitstackshape | Same relative path (line 15) | none |
| `old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh` | Cluster absolute paths (lines 16, 21-25, 27, 78) | Bismark 0.21.0, bowtie2 2.3.4.1, samtools 1.9, multiqc | Reads and genome: not in the repository | 9 paired-end reports in the repository (zr3534_1 to 6, 8 to 10; identical copies in `data/PE_reports/` and `old_WGBS/WGBS/data/PE_reports/`) |

Additional notes:

- `code/5-mC & 5mCpG plots.Rmd` must have been run with the repository root as working directory. Knitting from `code/` with default settings would fail at `code/5-mC & 5mCpG plots.Rmd:47`.
- `gene_expression/scripts/dct_analysis.R` does not parse as committed: there are git conflict markers at `gene_expression/scripts/dct_analysis.R:160-164`.
- `ELISA/GlobalDNAMeth_Polyploids.R` stops with an error at `ELISA/GlobalDNAMeth_Polyploids.R:231-233`.
- `code/1_WGBS_PE_report_analysis.R` (related, older) needs rstudioapi and `setwd` (`code/1_WGBS_PE_report_analysis.R:15-19`) and reads `data/PE_reports/percent_methylation_summary.xlsx` (`code/1_WGBS_PE_report_analysis.R:24-25`), which exists only under `old_WGBS/WGBS/data/PE_reports/`.

## Appendix F. `stats/*.csv`

All four files are verbatim, and all were reproduced by recomputation.

```
AOV_5-mC.csv
"Df" "Sum Sq" "Mean Sq" "F value" "Pr(>F)"
"ploidy          " 1 0.00482559048460679 0.00482559048460679 3.92546851356249 0.0574572812803601
"treatment       " 1 0.0170272350127791 0.0170272350127791 13.8511287124149 0.000881540398059937
"ploidy:treatment" 1 0.00480286582872205 0.00480286582872205 3.9069827092984 0.0580107135318376
"Residuals       " 28 0.0344204858865031 0.00122930306737511 NA NA

HSD_test_5-mC.csv
"response" "groups"
"heat.diploid" 0.197338207328138 "a"
"heat.triploid" 0.148275892602682 "b"
"control.diploid" 0.126701365523232 "b"
"control.triploid" 0.126643468092149 "b"

AOV_%mCpG.csv
"Df" "Sum Sq" "Mean Sq" "F value" "Pr(>F)"
"ploidy          " 1 1.77993192375073e-06 1.77993192375073e-06 0.0313241828842046 0.861495338917365
"treatment       " 1 3.84519752806729e-05 3.84519752806729e-05 0.676698187092797 0.421490521341217
"ploidy:treatment" 1 0.000823142963045614 0.000823142963045614 14.4861049853825 0.00129379504931095
"Residuals       " 18 0.0010228127816119 5.68229323117721e-05 NA NA

HSD_%mCpG.csv
"response" "groups"
"heat.diploid" 0.338072113786742 "a"
"control.triploid" 0.333731364363768 "ab"
"heat.triploid" 0.324101934697854 "b"
"control.diploid" 0.32313249436849 "b"
```

## Appendix G. The two data workbooks

- **`data/ELISA-percent-mCPG-summary.xlsx`**
  - `raw` (A1:K63): Well, Abs, Type, ploidy, treatment, Sample, mean, sd, numerator, denominator, perc. It has 62 well rows; the formulas are in section 1.1.
  - `summary` (A1:D33): the 32 values in Appendix A.
- **`data/bismark-percent-mCPG-summary.xlsx`**
  - `raw` (A1:H23): seqID, library_name, tissue (all ctenidia), ploidy, treatment, methylated_in_CpG, unmethylated_C_in_CpG, %mCpG. The 22 rows are in Appendix B.
  - `summary` (declared A1:P11; data A1:D5): ploidy, heat_shock, CpG, se, with four stale rows (diploid no 10.633, yes 11.05; triploid no 10.067, yes 9.8).

## Appendix H. Alignment job script

Condensed listing of `old_WGBS/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh`; every command line is quoted with its line number.

```
3-16   #SBATCH job-name=ron-rosM, account/partition=coenv, nodes=1, time=15-00:00:00, mem=100G,
       mail-user=<lab address>, chdir=/gscratch/scrubbed/sr320/030521-ronrosM
21     bismark_dir="/gscratch/srlab/programs/Bismark-0.21.0"
22     bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/"
23     samtools="/gscratch/srlab/programs/samtools-1.9/samtools"
24     reads_dir="/gscratch/srlab/sr320/data/cg/"
25     genome_folder="/gscratch/srlab/sr320/data/Cgig-genome/roslin_M/"
27     source /gscratch/srlab/programs/scripts/paths.sh
31-35  ${bismark_dir}/bismark_genome_preparation --verbose --parallel 28 --path_to_aligner ${bowtie2_dir} ${genome_folder}
38     #/zr3644_11_R2.fastp-trim.20201206.fq.gz
40-48  find ${reads_dir}*_R1.fastp-trim.20201202.fq.gz | xargs basename -s _R1.fastp-trim.20201202.fq.gz |
       xargs -I{} ${bismark_dir}/bismark --path_to_bowtie ${bowtie2_dir} -genome ${genome_folder} -p 8
       -score_min L,0,-0.6 --non_directional -1 ${reads_dir}{}_R1.fastp-trim.20201202.fq.gz -2 ${reads_dir}{}_R2.fastp-trim.20201202.fq.gz
52-57  find *.bam | xargs basename -s .bam | xargs -I{} ${bismark_dir}/deduplicate_bismark --bam --paired {}.bam
61-65  ${bismark_dir}/bismark_methylation_extractor --bedGraph --counts --scaffolds --multicore 28 --buffer_size 75% *deduplicated.bam
71     ${bismark_dir}/bismark2report
75     ${bismark_dir}/bismark2summary
78     /gscratch/srlab/programs/anaconda3/bin/multiqc .
83-87  find *deduplicated.bam | xargs basename -s .bam | xargs -I{} ${samtools} sort --threads 28 {}.bam -o {}.sorted.bam
92-95  find *.sorted.bam | xargs basename -s .sorted.bam | xargs -I{} ${samtools} index -@ 28 {}.sorted.bam
101-108 find *deduplicated.bismark.cov.gz | xargs basename -s _bismark_bt2_pe.deduplicated.bismark.cov.gz |
       xargs -I{} ${bismark_dir}/coverage2cytosine --genome_folder ${genome_folder} -o {} --merge_CpG --zero_based {}_bismark_bt2_pe.deduplicated.bismark.cov.gz
113-145 four awk loops over *merged_CpG_evidence.cov: rows with $5+$6 >= 10 (or >= 5) written to
       _10x.bedgraph / _5x.bedgraph (cols 1-4) and _10x.tab / _5x.tab (cols 1-6)
```

- Tool versions: Bismark 0.21.0, bowtie2 2.3.4.1, samtools 1.9; the multiqc version is unspecified.
- The script only matches the 20201202 (zr3534) files, so it cannot have produced the 13 reports made with Bismark 0.22.3. Those come from the June 2021 run (section 2.5).

## Questions for the authors

Only points that cannot be determined from the repository are listed.

1. **ELISA input mass.** What DNA mass was actually loaded per sample well, versus the 25 ng assumed? What does the `Brconc.` sheet represent (when and how were those concentrations measured, and is 2.0 a lower reporting limit)?
2. **ELISA standards.** What DNA mass did the "5 ul" standard wells contain, and was it matched to the sample mass as the kit manual requires?
3. **ELISA exclusions.** Why were wells A03, B03 and C03 flagged `contam`, and is the neighbouring single well of D02-C (D03) trustworthy?
4. **ELISA analysis choice.** Which analysis is intended for the manuscript: untransformed, Johnson or arcsine? If the Johnson transformation is still wanted, its result needs to be supplied, since it is not saved anywhere.
5. **Control libraries.** Which animals are zr3644_1 to 6 and zr3644_13 to 18, and how were they held and sampled relative to this experiment (cohort, dates, tank conditions, tissue handling)?
6. **Acute heat shock libraries.** Can you confirm that zr3534_4, 5, 9 and 10 are D19, D20, T19 and T20, as the sample sheet says, and whether pooling them with the heat-only animals is intended?
7. **Trimming.** Which fastp version and parameters were used for both projects, and were the zr3644 reads trimmed with the same settings?
8. **COX1.** Was the template DNA or cDNA, with what primers and input mass, and was a nuclear reference gene measured on the same template? Do replicate-level Cq values exist for SOD, ATP synthetase and COX1?
9. **qPCR statements in the Discussion.** Which analysis and figure were the HSP90 and mitochondrial statements based on?
10. **Software versions.** Which R and package versions produced `stats/*.csv` and the Figure 1 panels in February 2026?
11. **Pilot plate map.** Is `MethylFlashGlobalDNAMethylation_20181017.xlsx` available, so that the pilot analysis can be re-run?

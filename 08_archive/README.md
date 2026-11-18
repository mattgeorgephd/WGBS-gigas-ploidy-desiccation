# 08_archive

Everything that was in this repository before the reorganization of 2026-11-18, kept for provenance. The pipeline in `01_pipeline/` does not read from this folder, with one exception: the checks in `06_diagnostics/` compare regenerated outputs against the legacy outputs stored here.

Do not edit files in this folder. If a legacy result needs correcting, the correction belongs in the pipeline.

## Path rule

A file that was at `<path>` at commit `523d734` is now at `08_archive/<path>`. Paths cited in [`issues.md`](../issues.md) and in `07_documentation/audit/` resolve by that rule. Git history follows each file through the move (`git log --follow`).

Exceptions, moved out of the legacy tree because they are raw inputs or reference documents:

| Was | Now |
|---|---|
| `ELISA/docs/20190206_FirstReadingSTRIGG.xls`, `...2ndReading...`, `...3rdReading...` | `03_input_files/elisa/raw/` |
| `ELISA/docs/20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx` | `03_input_files/elisa/raw/20190131_RonitsSamplesforGlobalDNAMeth.xlsx` (apostrophe dropped from the name) |
| `ELISA/docs/MethylFlashGlobalDNAELISA_P-1030.pdf` | `07_documentation/protocols/` |
| `gene_expression/data/qpcr_ct_values/*.csv` | `03_input_files/qpcr/raw/` |
| `gene_expression/data/zr3534_wgbs_info.csv` | `03_input_files/metadata/raw/` |
| `metadata.xlsx` | `03_input_files/metadata/raw/legacy_metadata.xlsx` |
| `gene_expression/data/*.pdf` (two Zymo reports) | `07_documentation/sequencing_reports/` |
| `old_WGBS/WGBS/GFF/*` (11 genome feature tracks) | `03_input_files/genome_features/*.gz` (gzipped) |

Removed, not archived: five Office lock files (`~$*.xlsx`, `~$*.pptx`), one `.DS_Store`, one Jupyter checkpoint.

## Legacy dataset labels

File and folder names in this archive use three labels:

- `ronit` or `RONIT`: the ten heat-treated libraries (Zymo project zr3534), analysed alone.
- `yaamini` or `YAAMINI`: the control libraries (Zymo project zr3644), analysed alone.
- `both`, `BOTH`, or `all-gigas`: the 22-library comparison reported in the manuscript.

`COMBO` and `INTERACTION` both name the set of loci present in the ploidy and heat DML lists. `SHOCK` names the 45 C acute heat shock contrast within the heat-treated libraries. "Desiccation" is the legacy word for the heat treatment.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | The legacy analysis scripts and notebooks. See the table below. |
| `analyses/` | Hand-edited workbooks that sat between analysis steps (CpG counts by coverage, DML counts by threshold, feature overlap counts, the chi-squared sheet). |
| `data/` | The two hand-built workbooks behind Figure 1 and the per-library Bismark alignment reports that were in the repository (a mix of two alignment runs; see issue I-13). |
| `ELISA/` | Global 5-mC ELISA analysis written in 2019 (`GlobalDNAMeth_Polyploids.*`), the 2018 pilot, exploratory workbooks, and a qPCR notebook that was filed here. |
| `gene_expression/` | qPCR scripts and plots. |
| `GOterms/` | topGO result tables and annotation joins, one folder per dataset label and contrast. `BOTH-*` are the manuscript results. |
| `tables/` | DML lists for the manuscript contrasts (`DML-getMethylDiff-*-Cov10-20.csv`), GO tables, and table images. |
| `stats/` | Model tables for Figure 1 as written by the legacy notebook. |
| `figures/`, `plots/` | Figure panels, assembled figures (`figures/final_figures/`, PNG and PowerPoint), exploratory plots. |
| `old_WGBS/` | An earlier copy of the WGBS analysis (2021), duplicated when the repository was restructured the first time. |
| `presentations/` | Research poster. |

## Legacy code, by manuscript item

| Script (`code/`) | What it did | Manuscript item |
|---|---|---|
| `5-mC & 5mCpG plots.Rmd` | Box plots and two-way ANOVA for ELISA %5-mC and Bismark %mCpG, from two hand-built workbooks | Figure 1; Results 3.1 statistics |
| `1_coverage_analysis.R` | CpG counts at 1x to 20x for five sample sets; per-library histograms; whole-data PCA | Figure S1; Results 3.1 |
| `2_WGBS_Methylkit-all-gigas.R` | methylKit differential methylation for ploidy and heat on 22 libraries; overlap set; threshold sweep; PCA; heat maps; BED export | Results 3.2; Figures 2, 3B, 4 |
| `3_WGBS_DML_feature_location-new.ipynb` | bedtools overlaps of DML and CpG sets with 13 genome feature tracks (Python and bash) | Results 3.3; Figure 3A |
| `3_Genome_Feature_Location_plots.R` | Stacked bars drawn from `analyses/DML_locations-all-datasets.xlsx` | Figure 3A |
| `4_WGBS_GOterm-Annotation-BOTH.ipynb` | blastx to UniProt to GO annotation; DML to gene to transcript joins; background transcript list (Python and bash) | Input to enrichment |
| `5_WGBS_Functional-enrichement-analysis-BOTH.R` | topGO Fisher tests; "top ten" term tables; DML within those genes; simplifyEnrichment plots | Results 3.4, 3.5; Table 1; inputs to Figures 5 to 9 |
| `gene-plots.Rmd` | Bar charts drawn from numbers typed into the script | Figures 5 and 9 |
| `Revigo.R`, `../figures/revigo_figure/Revigo-PLOIDY-BP.R` | Plot scripts exported by the REVIGO web server | Figure 6 (Figures 7 and 8 have no script) |
| `2_WGBS_Methylkit.R`, `2_WGBS_Methylkit-heat-only.R`, `*-RONIT.*`, `*-YAAMINI.*`, `1_WGBS_PE_report_analysis.R`, `5_WGBS_Functional-enrichement-analysis.R`, `OLD/` | Earlier versions and single-dataset variants | Not in the manuscript |

None of these scripts runs as committed (absolute paths, saved-workspace dependencies, manual steps); see issues I-31 and I-32. The issues each one carries are listed in `issues.md`.

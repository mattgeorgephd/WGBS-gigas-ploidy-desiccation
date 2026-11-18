# elisa

Raw data of the global DNA methylation assay: MethylFlash Global DNA Methylation (5-mC) ELISA Easy Kit (Epigentek P-1030; manual in `07_documentation/protocols/`), one 96-well plate read on 2019-02-06, 32 ctenidia samples (8 per ploidy by treatment group) from the heat experiment, in duplicate wells.

## Files

| File | Contents |
|---|---|
| `raw/20190206_FirstReadingSTRIGG.xls` | Plate reader export, first reading. Sheet `plate list`: well and absorbance at 450 nm for all 96 wells. This is the reading used for the calculation. |
| `raw/20190206_2ndReadingSTRIGG.xls`, `raw/20190206_3rdReadingSTRIGG.xls` | Second and third readings taken a few minutes later; used only to check that readings were stable. |
| `raw/20190131_RonitsSamplesforGlobalDNAMeth.xlsx` | Sheet `Concentrations`: stock DNA concentration of each sample. Sheet `PlateMap`: planned and actual plate layouts. Sheet `PlateMap_list`: well, type, and sample for all 96 wells, including the flags `no_sample`, `contam`, and `no_ab`. The file was named `20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx` before the reorganization. |

## Plate layout

Standards occupy strips 1 and 2 (negative control, then 0.1, 0.2, 0.5, 1, and 5% 5-mC, and the undiluted positive control). Samples occupy strips 3 to 10 in group order: diploid control, diploid heat, triploid control, triploid heat. Three wells in strip 3 are flagged `contam` and one well in strip 8 `no_ab`; D01 was run in two wells of strip 1 instead, and D02 and T08 have one well each.

## Open points

The mass of DNA loaded per well, the meaning of an alternative input table in a legacy workbook, and the reason for the `contam` flags are not documented; see issue I-12 in [`issues.md`](../../issues.md) and `07_documentation/audit/elisa_qpcr_audit.md`.

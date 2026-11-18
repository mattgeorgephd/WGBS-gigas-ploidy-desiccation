# 05_manuscript

The manuscript's figures and tables under the numbers the manuscript uses, and the table that ties every quoted number to the script that produces it. Written by `01_pipeline/12_manuscript_values.R` from the outputs in `04_output/`; nothing here is drawn or typed by hand.

## Planned contents

| Path | Contents |
|---|---|
| `figures/Figure_<n>.png`, `.pdf` | Complete figures, panels assembled in R (`patchwork`), including significance letters and sample annotations |
| `tables/Table_1.csv` | DML within genes under enriched GO terms |
| `supplement/` | Lists of DML, GO terms, and genes cited in Section 7 of the manuscript |
| `manuscript_values.csv` | One row per number in the manuscript text: section, the sentence it comes from, the manuscript value, the regenerated value, the output file and script that hold it, and a status (`matches`, `differs`, `not reproducible`, `corrected`) |

The manuscript text itself is a Google Doc (linked from the root README). The legacy figures, assembled in PowerPoint, are in `08_archive/figures/final_figures/`.

No files exist here yet; this is the last phase of the plan.

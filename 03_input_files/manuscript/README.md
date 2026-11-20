# manuscript

The two tables that tie the manuscript draft to the pipeline. Both are read by `01_pipeline/13_manuscript_values.R`, which writes [`05_manuscript/`](../../05_manuscript/).

## Files

| File | Contents |
|---|---|
| `manuscript_values_map.csv` | One row per quantitative statement of the draft, with an R expression that recomputes the quantity from the pipeline outputs. 545 rows: Methods, Results 3.1 to 3.5, Table 1, the figure legends, the Discussion and Conclusions, and the section headed 3.5 that follows the Conclusions. Every row was written by one agent and checked by a second one against the draft and against the tables. |
| `manuscript_items_map.csv` | One row per file that step 13 copies into `05_manuscript/` under the number the manuscript uses (figures, Table 1, supplementary lists). |

### Columns of `manuscript_values_map.csv`

| Column | Meaning |
|---|---|
| `id` | Section prefix and running number: `M2.5-003` (Methods 2.5), `R3.2-010` (Results 3.2), `T1-004` (Table 1), `F3-001` (legend of Figure 3), `FS1-001` (Figure S1), `D4.3-002` (Discussion 4.3; `D4.9` is the Conclusions), `D3.5x-001` (the section headed 3.5 that follows the Conclusions) |
| `section` | Section of the draft |
| `manuscript_line` | Line of the text copy of the draft (see below) |
| `quote` | The clause that holds the value, copied exactly from the draft, at most 250 characters. This is what ties a row to the manuscript: search for it in the Google Doc |
| `item` | What the value is, in plain words |
| `manuscript_value`, `manuscript_numeric` | The value as written, and as a number where it is one |
| `kind` | count, percent, statistic, p_value, parameter, design, coordinate, identifier, inequality, ranking, software_version, other |
| `tolerance` | Largest difference from the manuscript value that still counts as a match: 0 for counts, identifiers and statements answered by TRUE or FALSE; half a unit of the last printed digit for rounded values (0.05 for a percentage printed with one decimal, 0.5 for one printed as a whole number) |
| `expression_legacy` | R expression that recomputes the quantity as the draft's analysis did. Empty when the pipeline cannot regenerate the statement (a laboratory detail, for example); the note says why |
| `expression_optional` | Template of an R expression for the same quantity in an optional analysis that the paper did not adopt, where that differs; `A` is the analysis id (for example `corrected`: the upper coverage filter of the Methods, I-05; `mnchisq`: the overdispersion correction; for the GO rows the gene-level analysis of step 10). Step 13 evaluates it in an optional run for every analysis folder and for the four single options, and reports "not run" where the analysis or a table it reads is absent. Until 2026-11-20 the column was `expression_corrected`, with `corrected` where `A` stands now; the rewrite (`07_documentation/toggles_migration/rewrite_maps.py`) changed no other column, and with `A = "corrected"` every expression gives its former value. The template is therefore the corrected computation (as the former column was), not the paper's computation applied to the analysis: applied to the paper's outputs the templates differ from `expression_legacy` in 56 of 263 values (2026-11-20). Changed since: R3.1-001 and R3.1-003 take the column of the analysis's coverage filter (`cpg_in_any_library` in the legacy reading, `cpg_in_any_library_unfiltered` with option `corrected`, as the former corrected expression did); D4.9-014 gives NA ("not comparable") when the two AK9 loci are not both in the analysis's loci in both lists; and the seven templates that read the united table of step 05 (R3.2-029, R3.2-030, R3.2-031, D4.2-004, D4.3-005, D4.9-014, D3.5x-006) read `united_at_dml_<A>.csv`, its rows at the DML, which is tracked where the per-CpG tables are not. For a statement about what a contrast means ("baseline" differences, loci that respond to heat differently by ploidy) it is the test of that meaning with the diagnostic analysis that fits it (`06_diagnostics/dml_sensitivity/`), and the note says which |
| `issue` | The entries of [`issues.md`](../../issues.md) that concern the statement |
| `note` | Why a value differs, why it cannot be regenerated, or what must be known to read a match |

An expression is one line of R that returns one value, with these helpers in scope: `tab("<step folder>", "<file>")` reads a table of `04_output/<step folder>/tables/` (for `expression_optional`, of `04_output/not_adopted/<A>/<step folder>/tables/`, or of the paper folder for a table that does not depend on the analysis), `paper_tab()` a table of the paper's outputs, `diagnostic("<file>")` a table of `06_diagnostics/`, `config` the settings list of `run_config.R` (with the coverage standard of the analysis), and `A` the analysis id (optional expressions only). A number is compared with the manuscript value within the tolerance; a text value must be identical; `TRUE` or `FALSE` is the answer for a statement that is an inequality or a ranking. The header of step 13 has the details.

The file is written with a CSV library (quoted fields, doubled quotes inside fields, CRLF line ends). Edit it with a CSV library or a spreadsheet, not by hand, and read it in R with `read.csv()`: `data.table::fread()` does not undo the doubled quotes.

### Columns of `manuscript_items_map.csv`

| Column | Meaning |
|---|---|
| `item` | The manuscript item (Figure 3, Table 1, a supplementary list) |
| `analysis` | `legacy` (the paper; copied by the paper run), `optional` (a template: `{A}` in `source` and `destination` stands for the analysis id; copied by an optional run for every analysis folder whose source exists, into `05_manuscript/not_adopted/<id>/`), `comparisons` (a table or figure of `04_output/not_adopted/comparisons/`; copied by an optional run), or empty when the item does not depend on the analysis (paper run) |
| `source` | Pipeline output to copy, relative to the repository root |
| `destination` | Path under `05_manuscript/` |
| `rows` | Optional R expression; when given, only the rows of the source table that satisfy it are written |

## The manuscript text is not in this repository

The draft is an unpublished Google Doc ([link in the root README](../../README.md)); the repository is public. The first author decided on 2026-11-19 that the text is not to be committed.

What the tables refer to:

- **The copy.** The rows were built from a plain-text copy of the Google Doc as it stood on **2026-11-18**: 205 lines, 80,244 bytes, md5 `f852956ab244b18fbf05f17c22f83630`. It was read through the Google Drive connector of the audit session, which returns the document as Markdown-like text with one paragraph, heading, table row, or figure legend per line. An export made another way (File, Download, Plain text) breaks lines differently, so its line numbers and checksum will not match.
- **Line ranges of that copy.**

  | Lines | Part of the draft |
  |---|---|
  | 1 to 16 | Title page, abstract placeholder, keywords, Introduction |
  | 17 to 32 | Materials and Methods, sections 2 to 2.6 |
  | 33 to 37 | Results 3.1 |
  | 38 to 46 | Results 3.2 and 3.3 |
  | 47 to 58 | Results 3.4 and 3.5 |
  | 59 to 113 | Discussion (sections 4.1 to 4.8) and Conclusions (4.9) |
  | 114 to 115 | The section headed 3.5 that follows the Conclusions |
  | 116 to 158 | References |
  | 159 to 183 | Table 1 |
  | 185 to 199 | Legends of Figures 1 to 9 |
  | 200 to 204 | Supplemental material: legend of Figure S1 and the three lists it names |

- **What to rely on.** `manuscript_line` is a convenience for whoever holds that copy. The durable reference is the `quote` column. If the Google Doc has been edited since 2026-11-18, a quote that can no longer be found marks a statement that changed, and its row has to be revised.

To extend or rebuild the map with the agent workflow in `07_documentation/agent_workflows/` a text copy is needed; export the Google Doc again and check the line ranges in `scripts/manuscript-values.js` against the new copy.

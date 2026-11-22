# Build specification: settings and toggles of the pipeline

Judge D3 of the scoping pass. Nothing in the repository was changed or run. Code read on branch `reproducibility-audit` between 02:09 and 02:20 UTC on 2026-11-20, at HEAD `671feca` (02:00 UTC), with a clean working tree (`git status --short` empty at 02:09 and 02:20). Line numbers below are those of that commit. The hardening pass may still move them; the builder re-reads every file named here before editing it.

Result in one sentence: build design D2 (one place for the paper, one labelled folder per optional analysis), with eleven repairs and five ideas taken from D1, as specified in sections 3 to 9.

Contents: 1 review of the two proposals; 2 scores; 3 principles and run kinds; 4 settings and switches; 5 names and folders; 6 changes by file and step; 7 order of work; 8 tests; 9 run time and memory; 10 decisions for the author; 11 what was verified and what is assumed.

---

## 1. Review of the two proposals

### 1.1 What D1 (smallest change) gets wrong or leaves out

| # | Defect | Where (code) | Effect |
|---|---|---|---|
| D1-1 | 51 corrected expressions of the values map read spanning tables with `mode == "corrected"`: `cpg_counts_by_coverage.csv` 7, `cpg_counts_by_feature_track.csv` 10, `library_correlation_summary.csv` 2, `dml_counts_by_mode.csv` 29, `annotation_summary.csv` 3. D1 moves the corrected rows into `_corrected` twins (its sections 3 and 7) but its change list (section 8) does not rewrite the map. | `03_input_files/manuscript/manuscript_values_map.csv` (count by `code/map_counts.py`); `13_manuscript_values.R` `evaluate()` and the final `stop()` | each of the 51 returns 0 values, is recorded as an error, and step 13 stops in every run with `--corrected` |
| D1-2 | Step 13 evaluates optional expressions for `corrected` only. Analyses with `_5x`, `_mn_chisq`, `_ctexcl` get no values and no "not run" entry. | D1 section 4, row 13 | the brief's "step 13 reports optional values as not run when absent" holds only for (a) |
| D1-3 | Every optional run re-runs the paper's analysis in the same processes and rewrites the paper's files (`legacy` is always first in the derived `config$modes`, D1 section 1). Paper bytes then rest on the restructured shared loops being free of bugs (step 04 becomes loops over readings and (reading, coverage) pairs). `finish_step()` rewrites `session_info.txt` of every paper folder. | `02_R_functions/logging.R` line 36; `01_pipeline/04_cpg_coverage.R` lines 149, 240, 309, 349 | "impossible to overwrite" is not met; it is only tested afterwards. An optional run also pays for steps 01 to 03, 11 and 12 |
| D1-4 | The variant step is `04b`, after step 04, so the exclusion cannot reach the CpG sets of step 04. The 2022 study and the author's instruction remove the variant CpGs "from the characterised CpGs and from the candidate DML". | D1 open point 2; `07_documentation/session_handoff.md` item 2b (line 95); variant spec `V2-method/spec.md` 3.5 | option (d) would not do what was asked |
| D1-5 | Optional files are not labelled "not adopted": `_corrected`, `_5x`, `_legacy_5x`. Handoff item 3 asks for "names or folders that say 'not adopted'". `dml_ploidy_legacy_5x.csv` starts with the paper's name and matches a glob `dml_ploidy_legacy*`. The 44 committed corrected files stay beside the paper's files. | `session_handoff.md` line 97; D1 section 3 | a reader can mistake an optional file for a paper file |
| D1-6 | Two positions of the tag in step 04: `library_correlation_5x_legacy.csv` but `library_correlation_summary_legacy_5x.csv`. | D1 section 3, examples | no single naming rule |
| D1-7 | Step 10's comparison tables (`go_terms_by_analysis.csv`, `legacy_top10_at_gene_level.csv`) are untagged and hold only the analyses of the run that wrote them, so a later run with fewer analyses removes rows. | `10_go_enrichment_gene_level.R` lines 170-179, 209-212; D1 open point 7 | violates "tables that span analyses do not lose rows" |
| D1-8 | `--diagnostics` turns on (a) and (e), so a diagnostics run refits corrected mode and reruns step 10. | D1 section 2 | an hour or more of fits that the diagnostics do not need (they can read the committed corrected outputs) |
| D1-9 | Count error: `config$modes` is used 25 times in 5 diagnostics scripts, not 31 in 6. | `grep -c 'config\$modes' 06_diagnostics/*.R` at 02:1x UTC | minor |

What D1 gets right and this specification keeps: a separate manifest for the alignments (step 01 downloads every row of `gannet_manifest.csv`, `01_fetch_inputs.R` lines 46-47); the paper's two `calculateDiffMeth()` calls stay literally unchanged and only the MN calls add `overdispersion = "MN", test = "Chisq"`; `config$min_coverage` keeps meaning the paper's 10; the step 11 gene-level outputs get their own `_gene_level` names; a missing table in a paper expression stays an error; a unit test of the settings over every combination of toggles; the toggles printed in every banner.

### 1.2 What D2 (clearest separation) gets wrong or leaves out

| # | Defect | Where (code) | Effect |
|---|---|---|---|
| D2-1 | `finish_step()` writes `session_info.txt` to `output_dir(step)`, the paper folder. D2's change list for `logging.R` names only `start_step()` and `require_mode()`. | `02_R_functions/logging.R` lines 35-41 | every optional run writes into six paper folders; D2's own test 8 would fail |
| D2-2 | Mode-independent outputs of steps 04 and 08 are written to the paper folder whatever the analysis: `library_coverage.csv` (line 165; its column `n_cpg_10x` uses `config$min_coverage`, line 145), `genome_cpg_total.csv` (line 193); step 08 `go_map_id_repair.csv` (85), `gene_go_map.csv.gz` (93), `transcript_table.csv.gz` (99). D2 says only "writes the same file set into `not_adopted/<id>/`". | `04_cpg_coverage.R`, `08_gene_annotation.R` | either a write into a paper folder or duplicated paper tables with a wrong label |
| D2-3 | Step 11 contradicts the run kinds. Section 1 says a paper run never writes under `not_adopted/`; section 5 has the paper's step 11 write its gene-level parts to `not_adopted/legacy/11_semantic_similarity/` "when step 10 output exists". Step 10's legacy output will be committed there, so every default run would write optional files. Step 11 is "not started" in optional runs, so the gene-level part would never run on request. | D2 sections 1, 3, 5 (row 11), 10 | default run not paper-only; (e) incomplete |
| D2-4 | The run kind of `--gene-level-go` alone is ambiguous: "any switch makes an optional run" (section 1) against "steps of the paper run: 01 to 09, 11, 12, 13; 10 only with (e)" (section 3). | D2 sections 1 and 3 | two readings of one command |
| D2-5 | In the paper part of step 13 a missing diagnostic table gives "not run" instead of an error. | D2 section 8; `13_manuscript_values.R` header lines 67-68 ("a broken row cannot pass unnoticed") | 15 paper values could silently stop being checked |
| D2-6 | No place for the alignment list. Adding 22 BAM and 22 index rows to `gannet_manifest.csv` makes step 01 of the default run download 109.1 GB. Free disk is 28 GB (`df -h`, 02:1x UTC); the largest BAM is 6.76 GB. | `01_fetch_inputs.R` lines 46-47; `V1-data/bam_files.csv` | default run breaks or fills the disk |
| D2-7 | The six `WGBS_*` variables of D2 can be left exported in a shell and switch options on silently (the hazard D1 names as R10), and D2 does not address it. | D2 section 2.2 | a "default" run that is not the paper |
| D2-8 | After `config$modes` is removed, a loop `for (mode in config$modes)` that the builder misses runs zero times without an error (`config$modes` is `NULL`). D2 says a missed site "fails loudly". | 25 such loops in 01_pipeline (count below) and 25 in 06_diagnostics | silent loss of outputs |
| D2-9 | The post-test exclusion breaks two lines of step 05 that D2 does not mention: `stopifnot(nrow(d) == nrow(dml_table(fit$diffmeth[[contrast]], ...)))` (line 238) compares with the unexcluded fit, and the threshold sweep (lines 263-271) reads the fit directly, so it would still count excluded loci. | `05_differential_methylation.R` | the step stops, or the sweep table counts loci that the DML lists exclude |
| D2-10 | No table for the number of DML removed per contrast. Adding columns to `dml_counts_by_mode.csv` would change the header of the paper's copy. | D2 section 5 row 05 | lost information, or a paper byte change |
| D2-11 | Step 13 lists only the analysis folders present; the four standard options that were never run do not appear as "not run". | D2 section 8 | weaker than the brief |

What D2 gets right and this specification keeps: the paper and the optional analyses in different folders, with `not_adopted` in every optional path; one analysis per run, named by tokens in a fixed order; the paper's id `legacy`, so every paper file keeps its name; `not_adopted/comparisons/`, rebuilt from every analysis folder on disk; a stamp file per analysis folder; the step 05 cache per fit id, with the exclusion after the test so it never refits; the scripted rewrite of the 263 corrected expressions into templates, with the check that `A = "corrected"` gives today's values; the README table; `--plan`; the frozen paper settings (as a recommendation, section 10).

### 1.3 Points of the brief, proposal by proposal

| Point of the brief | D1 | D2 | This spec |
|---|---|---|---|
| Default run gives the paper only, names and bytes kept | yes, except the one-time drop of optional rows (both) | yes, same | yes, same; one-time changes listed in 7.3 |
| Each option off by default, setting plus switch | yes | yes | yes |
| Optional outputs can never overwrite the paper's | tested only (D1-3) | by folders, but D2-1, D2-2, D2-3 | by folders and a write guard in `paths.R` (4.6) |
| Cannot be mistaken for paper outputs | weak (D1-5) | yes | yes |
| Combinations | suffix per run, paper always rerun | one analysis per run | one analysis per run |
| Steps 06 to 13 for whichever analyses ran | 06-10 yes; 13 corrected only (D1-2) | yes, except D2-3 | yes |
| Step 13 reports absent values as not run | corrected only | folders present only | every folder present plus the four standard options |
| Spanning tables do not lose or mix rows | step 04, 05, 08 yes; step 10 no (D1-7) | yes | yes, rebuild order specified (6.11) |
| Cache valid per analysis | yes | yes, plus stamp | yes |
| One README table | partial | yes | yes |
| Run time and memory of 5x | estimated | estimated | both estimates, pilot required (section 9) |

---

## 2. Scores

Scale 1 (poor) to 5 (good). For "amount and risk of code change", 5 means least.

| Criterion | D1 | D2 | Reason |
|---|---|---|---|
| Safety of the paper's outputs | 2 | 4 | D1 rewrites paper files in every optional run (D1-3). D2 separates by folder but has three leaks (D2-1 to D2-3) |
| Clarity for a user | 2 | 5 | D2: "a path with `not_adopted` is not the paper", one command per choice, `--plan` |
| Amount and risk of code change | 4 | 2 | D1 about 33 code places but undercounted (D1-1, D1-2). D2 about 300 sites in about 30 files, about 90 file moves, 263 expressions rewritten by script |
| Combinations | 3 | 4 | D1: (b) to (d) apply to every reading and the paper always runs. D2: any combination is one analysis; several analyses by several runs |
| Step 13 and cross-analysis tables | 2 | 4 | D1-1, D1-2, D1-7; D2-11 |
| Cache validity | 4 | 5 | both cache per fit id with a model part; D2 adds the stamp |

The extra code of D2 is mechanical (paths, loops) and every missed site either trips the write guard of 4.6 or the static test T0, so its risk can be bought down; the defects of D1 sit in the structure.

---

## 3. Principles and run kinds

- **P1 Paper bytes.** Every paper output depends only on the paper's settings and inputs. No optional run opens a file in a paper folder for writing. Enforced in code (4.6) and tested (T2).
- **P2 Labels.** Every optional output lies under a folder whose path contains `not_adopted`, or under `06_diagnostics/` (option (f), whose folder already says what it is; renumbered (g) on 2026-11-21 and (h) on 2026-11-22, see the note in section 3). A path without either is a paper output.
- **P3 Spanning tables.** The paper's copy holds the paper's rows only. Each analysis folder holds its own copy. `04_output/not_adopted/comparisons/` holds the rebuilt side-by-side table of every analysis folder on disk, paper first.
- **P4 One analysis per run.** Switches (a) to (d) compose one analysis id. (e) and (f) add steps.

> Note of 2026-11-21 (after this specification was implemented): the option `d50` (`--difference 50`, author's decision of 2026-11-21) was added as option (e), so the letters of this specification are now: (a) to (e) compose one analysis id (tokens `corrected`, `cov5`, `mnchisq`, `ctexcl`, `d50`), (f) is the gene-level GO enrichment and (g) the diagnostics. On 2026-11-22 the option `q05` (`--qvalue 0.05`, author's request) was added as option (f): (a) to (f) compose one analysis id (tokens `corrected`, `cov5`, `mnchisq`, `ctexcl`, `d50`, `q05`), (g) is the gene-level GO enrichment and (h) the diagnostics. The text below keeps its original letters; `run_config.R` and the READMEs use the new ones.

Run kinds, decided in `run_config.R` from the settings and the switches:

| Run kind | When | Steps started | Writes |
|---|---|---|---|
| `paper` | id `legacy`, (e) off, (f) empty | 01, 02, 03, 04, 05, 06, 07, 08, 09, 11 (paper part), 12, 13 (paper part) | `04_output/<step>/`, `05_manuscript/` outside `not_adopted/` |
| `optional` with id other than `legacy` | any of (a) to (d) on | variant step if (d) and its list is absent; 04, 05, 06, 07, 08, 09 for the analysis; 10 and 11 (gene-level part) if (e); diagnostics if (f); 13 (optional part) | `04_output/not_adopted/<id>/`, `04_output/not_adopted/comparisons/`, `04_output/not_adopted/ct_variants/`, `05_manuscript/not_adopted/`, `06_diagnostics/` (only with (f)) |
| `optional` with id `legacy` | only (e) and/or (f) | 10 and 11 (gene-level part) if (e), on the paper's DML; diagnostics if (f); 13 (optional part) | `04_output/not_adopted/legacy/`, `comparisons/`, `05_manuscript/not_adopted/`, `06_diagnostics/` |

Step 10 is never part of a paper run. Steps 01 to 03 and 12 are never part of an optional run; started alone in an optional configuration they say so and end with status 0. An optional run reads the paper's committed outputs (step 09 Fisher tables, step 08 gene-GO map) through `paper_path()` and never regenerates them.

---

## 4. Settings and switches

### 4.1 `run_config.R`

Remove the block "Analysis modes" (lines 23-40) and its check (lines 209-213). Add after the inputs block:

```r
  # ---- Optional analyses: NOT ADOPTED, all off by default ----------------------
  # With every setting below at its default the run produces the results of the
  # paper and nothing else. Settings (a) to (d) together define ONE optional
  # analysis per run, named by its id (04_output/not_adopted/<id>/); (e) and (f)
  # add steps. The switches of run_pipeline.R turn the same settings on for one
  # run. A step started alone reads WGBS_ANALYSIS, WGBS_GENE_LEVEL_GO and
  # WGBS_DIAGNOSTICS, which run_pipeline.R sets for its steps.
  optional = list(
    coverage_filter_of_methods = FALSE,         # (a) token "corrected": CpGs at or above the 99.9th percentile of
                                                #     coverage in a library removed; every line of each file read (I-05, I-38)
    coverage_standard          = 10L,           # (b) reads per CpG in each of the 22 libraries: 10 (the paper) or
                                                #     5 (token "cov5"; Venkataraman et al. 2022)
    overdispersion_chisq       = FALSE,         # (c) token "mnchisq": methylKit overdispersion = "MN", test = "Chisq" (I-06);
                                                #     no other test is offered
    exclude_ct_variants        = FALSE,         # (d) token "ctexcl": CpGs at C-to-T variants removed from the CpG sets of
                                                #     step 04 and, after the test, from the tested CpGs and DML (I-08)
    gene_level_go              = FALSE,         # (e) step 10 and the gene-level part of step 11 (I-41)
    diagnostics                = character(0)   # (f) "all", or names of scripts of 06_diagnostics/ without .R
  ),
```

No other setting is added. `min_coverage = 10`, `import_min_coverage = 2`, `upper_coverage_percentile = 99.9`, `dml_qvalue = 0.01`, `dml_difference = 20` keep their names and values and keep meaning the paper's.

### 4.2 Checks block (`run_config.R`, inside `local()`)

1. Types: `coverage_filter_of_methods`, `overdispersion_chisq`, `exclude_ct_variants`, `gene_level_go` are `TRUE` or `FALSE`; `coverage_standard` is 5 or 10 (whole number); `diagnostics` is `character(0)`, `"all"`, or names of existing files `06_diagnostics/<name>.R` listed in `06_diagnostics/diagnostics.csv`.
2. Environment: `WGBS_ANALYSIS` (an id, 5.1), `WGBS_GENE_LEVEL_GO` (`TRUE`/`FALSE`), `WGBS_DIAGNOSTICS` (`none`, `all`, or comma-separated names). An empty variable is an error, as for `WGBS_CORES` (lines 218-223). A set `WGBS_ANALYSIS` replaces (a) to (d) as a whole; it is parsed by `analysis_spec()` and must equal the id rebuilt from its tokens (so `cov5-corrected`, `legacy-cov5`, `Cov5`, `5x` stop with a message listing the valid tokens in order).
3. Frozen paper settings (decision 10.6, recommended): compare `min_coverage`, `import_min_coverage`, `upper_coverage_percentile`, `dml_qvalue`, `dml_difference` with `c(10, 2, 99.9, 0.01, 20)`; a difference stops with "config$<name> is a setting of the paper. An optional analysis is chosen with config$optional or a switch of run_pipeline.R". Reason: 22 paper expressions read `config$min_coverage` and 21 read `config$dml_difference` (counts of D1 and D2, consistent with `code/map_counts.py`).
4. `coverage_levels` must contain 1, 5 and 10 (step 04 counts at the standard of either analysis).
5. Derived values:
   - `cfg$analysis <- analysis_spec(id)`: `list(id, paper = id == "legacy", tokens, filter = "legacy" | "corrected", min_coverage = 10L | 5L, overdispersion = "none" | "MN", ct_excluded, fit_id)`. `fit_id` is the id without the token `ctexcl`, and `"legacy"` when nothing is left (`ctexcl` -> `legacy`; `cov5-mnchisq-ctexcl` -> `cov5-mnchisq`).
   - `cfg$addons <- list(gene_level_go, diagnostics)`.
   - `cfg$run_kind <- if (cfg$analysis$paper && !cfg$addons$gene_level_go && length(cfg$addons$diagnostics) == 0) "paper" else "optional"`.
   - `cfg$modes` is not defined.

`analysis_spec()` must be available in `run_config.R` before `02_R_functions/` is sourced; define it in `02_R_functions/analyses.R` and source that one file at the top of the checks block (`source(here::here("02_R_functions", "analyses.R"), local = TRUE)`), or define it inline in `run_config.R`. The builder picks one; the function exists once.

### 4.3 Switches of `run_pipeline.R`

```
Rscript run_pipeline.R [--steps <selection>] [--coverage-filter] [--coverage 5|10] [--overdispersion]
                       [--exclude-ct-variants] [--gene-level-go] [--diagnostics [all|name,...]]
                       [--plan] [--list] [--help]
```

| Switch | Setting | Token |
|---|---|---|
| `--coverage-filter` | (a) `coverage_filter_of_methods = TRUE` | `corrected` |
| `--coverage 5` | (b) `coverage_standard = 5L` (`--coverage 10` is allowed and means the paper's 10) | `cov5` |
| `--overdispersion` | (c) `overdispersion_chisq = TRUE` | `mnchisq` |
| `--exclude-ct-variants` | (d) `exclude_ct_variants = TRUE` | `ctexcl` |
| `--gene-level-go` | (e) `gene_level_go = TRUE` | none |
| `--diagnostics` or `--diagnostics name1,name2` | (f) `"all"` or the names | none |
| `--plan` | prints the plan and stops | |

Rules:
- Arguments in any order; each switch at most once; anything else stops before a step starts (replaces lines 44-69). A switch can only turn a setting on; the result is the union of `run_config.R` and the switches.
- `run_pipeline.R` refuses to start when `WGBS_ANALYSIS`, `WGBS_GENE_LEVEL_GO` or `WGBS_DIAGNOSTICS` is set in its own environment ("unset it; choose the analysis with run_config.R or a switch"). It then sets the three variables with `Sys.setenv()` before the first `system2()` and unsets them on exit. This closes the hazard of a variable left in a shell (D2-7). `Sys.setenv()` is used, not the `env` argument of `system2()`, so that it works on every platform.
- The plan is printed at the start of every run and alone with `--plan`: run kind, analysis id, its folder, steps to run, steps skipped and why, whether the variant list exists, and the expected time and memory (section 9, with a warning for any `cov5` id). The id and run kind are printed again at the end with the timing table (replaces lines 75-76, 111-116).
- `--steps` selects among the steps of the run kind. A selected step that is not part of it appears in the closing table as "not run: <reason>" (replaces `needs_legacy`, lines 78-97).
- `--list` prints every step with its run kinds, plus the unnumbered variant step and the diagnostics.
- Variant step: with (d) on and `04_output/not_adopted/ct_variants/tables/ct_variant_cpgs.csv.gz` absent, `01_pipeline/optional/ct_variant_loci.R` runs first, after the plan has said that it downloads 22 alignments one at a time (109.1 GB in all). The script pattern `^[0-9]{2}_.*\\.R$` (line 35) stays; the variant script is not numbered and lies in a subfolder, so `--steps` never selects it.
- Diagnostics: run after step 11 and before step 13, in the order of `06_diagnostics/diagnostics.csv`. A script whose needs are not met (for example network) is reported as "not run: needs <x>" and the run goes on.

### 4.4 Examples

```sh
Rscript run_pipeline.R                                                     # the paper
Rscript run_pipeline.R --coverage-filter                                   # -> not_adopted/corrected/
Rscript run_pipeline.R --coverage 5                                        # -> not_adopted/cov5/
Rscript run_pipeline.R --overdispersion                                    # -> not_adopted/mnchisq/
Rscript run_pipeline.R --exclude-ct-variants                               # -> not_adopted/ctexcl/
Rscript run_pipeline.R --gene-level-go                                     # step 10 on the paper's DML -> not_adopted/legacy/
Rscript run_pipeline.R --diagnostics                                       # 06_diagnostics/, every script
Rscript run_pipeline.R --coverage-filter --coverage 5 --exclude-ct-variants   # -> not_adopted/corrected-cov5-ctexcl/
Rscript run_pipeline.R --coverage 5 --overdispersion --exclude-ct-variants    # -> not_adopted/cov5-mnchisq-ctexcl/
Rscript run_pipeline.R --plan --coverage 5 --overdispersion                # print only
WGBS_ANALYSIS=cov5 Rscript 01_pipeline/05_differential_methylation.R       # one step alone
```

### 4.5 Step 05 model and cache

- `fit_mode(a)`: `cfg_a <- modifyList(config, list(min_coverage = a$min_coverage))`; `read_filtered_counts(path, a$filter, cfg_a, chr_levels)` (the helper itself is unchanged; it accepts only `"legacy"` and `"corrected"`, `methylkit_helpers.R` line 30). Log lines use `a$fit_id` and `cfg_a$min_coverage` (line 153 reads `config$min_coverage` today).
- The two tests are written as literal calls in an `if`:
  ```r
  dm_ploidy <- if (a$overdispersion == "none") {
    calculateDiffMeth(meth_ploidy, covariates = data.frame(heat = lib$treatment), mc.cores = config$cores)
  } else {
    calculateDiffMeth(meth_ploidy, covariates = data.frame(heat = lib$treatment),
                      overdispersion = "MN", test = "Chisq", mc.cores = config$cores)
  }
  ```
  and the same for heat. The paper's call keeps exactly its present arguments (lines 156-160); `test` is not passed to it. `do.call()` is not used (it would pass the methylBase object by value into the call).
- Memory (recommended in the same edit, so the paper refits only once): drop `kept[[i]]` as soon as object `i` is built, and build `meth_heat` after the ploidy test. Neither changes a number; T1 proves it.
- Cache file `04_output/05_differential_methylation/cache/methylkit_<a$fit_id>.rds` (paper: `methylkit_legacy.rds`, corrected: `methylkit_corrected.rds`, unchanged names). `fit_fingerprint(a)` keeps its six parts with `mode = a$fit_id` and `coverage_settings` from `cfg_a`, and gains `model = c(overdispersion = a$overdispersion, test = if (a$overdispersion == "MN") "Chisq" else "not passed")`. The cached object keeps `mode = a$fit_id`; `write_mode(fit, a)` takes names from `a$id`.
- An optional run whose fit id is `legacy` (id `ctexcl`) reads the paper's cache; if that cache is stale it refits and rewrites `methylkit_legacy.rds`. That file is a git-ignored intermediate, not an output; the write guard allows the `cache` subfolder.
- One-time cost: `fit_mode()` changes, its parsed code is part of the fingerprint (lines 178-187), so the first run refits the paper (about 18 minutes at `WGBS_CORES=2`) and, when asked for, `corrected`.

### 4.6 Write guards (`02_R_functions/paths.R`)

- `output_dir(step, sub)`: when `config$run_kind == "optional"` and `sub != "cache"`, stop with "an optional run may not write into the paper folder 04_output/<step>/<sub>; use analysis_path()". Every step builds its paths through `output_dir()` or `output_path()` today, so every site the builder misses fails at once instead of writing.
- `paper_path(step, sub, file)`: `here::here("04_output", step, sub, file)`, never creates a folder; for reading paper outputs in optional runs (step 10 reads step 09 Fisher tables and step 08 `gene_go_map.csv.gz`; step 11 gene-level part reads `named_terms_results_3_4.csv`; step 13 `paper_tab()`).
- `manuscript_path(...)`: in an optional run only `manuscript_path("not_adopted", ...)` is allowed.
- `analysis_dir(a, step, sub)` (in `analyses.R`): the paper analysis in a paper run returns `output_dir(step, sub)`; the paper analysis in an optional run (id `legacy`, only (e)) returns `04_output/not_adopted/legacy/<step>/<sub>`; any other analysis `04_output/not_adopted/<id>/<step>/<sub>`. It creates the folder. In a paper run it stops if asked for a non-paper analysis.
- `rebuild_comparison()` stops in a paper run.

---

## 5. Names and folders

### 5.1 Analysis ids

- Tokens, in this fixed order: `corrected`, `cov5`, `mnchisq`, `ctexcl`. The id is the tokens that are on, joined by `-`. No token: `legacy` (the paper).
- 16 ids: `legacy`, `ctexcl`, `mnchisq`, `mnchisq-ctexcl`, `cov5`, `cov5-ctexcl`, `cov5-mnchisq`, `cov5-mnchisq-ctexcl`, `corrected`, `corrected-ctexcl`, `corrected-mnchisq`, `corrected-mnchisq-ctexcl`, `corrected-cov5`, `corrected-cov5-ctexcl`, `corrected-cov5-mnchisq`, `corrected-cov5-mnchisq-ctexcl`.
- No optional id contains `legacy`. Token `cov5`, not `5x`, because step 04 names already carry the level (`library_correlation_5x_cov5.csv` is readable; `..._5x_5x` would not be).
- Order of analyses wherever several are listed (comparison tables, step 13): `legacy` first, then the other ids in byte order (`sort(method = "radix")`). With `legacy` and `corrected` this is today's order.
- The column `mode` of every table keeps its name and holds the id (renaming it would change every paper table that has it). READMEs say "mode: the analysis id".

### 5.2 Folders

```text
04_output/
  01_fetch_inputs/ ... 13_manuscript_values/    THE PAPER (step 10 folder no longer exists)
  not_adopted/
    README.md                                   the table of 5.4
    <id>/                                       one optional analysis
      analysis.csv                              stamp (5.3)
      04_cpg_coverage/ ... 09_go_enrichment_legacy/
      10_go_enrichment_gene_level/              only with (e)
      11_semantic_similarity/                   gene-level part only, with (e)
    legacy/                                     (e) on the paper's DML: 10_... and 11_... only
    ct_variants/tables/                         ct_variant_cpgs.csv.gz, ct_variant_run.csv, ct_variant_log.txt
    comparisons/                                rebuilt side-by-side tables; analyses_included.csv
    13_manuscript_values/tables/                expression_errors_by_analysis.csv
05_manuscript/
  manuscript_values.csv, manuscript_values_summary.csv, figures/legacy/, supplement/legacy/, tables/legacy/   THE PAPER
  not_adopted/
    manuscript_values_by_analysis.csv, manuscript_values_summary_by_analysis.csv
    <id>/figures/, <id>/supplement/             items of an analysis
    comparisons/figures/, comparisons/supplement/
01_pipeline/optional/ct_variant_loci.R          the variant step (unnumbered)
```

### 5.3 Stamp `analysis.csv`

One row: `id, filter, min_coverage, overdispersion, test, ct_excluded, fit_id, fit_fingerprint_md5` (md5 of `serialize()` of the fingerprint, written by step 05), `ct_variant_list_md5` (with `ctexcl`), `git_commit`, `written_utc`, `steps_done`. Step 04 creates it, each later step adds itself to `steps_done`. Steps 05 to 13 stop if the stamp's settings differ from `config$analysis`, or if a step they read from is not in `steps_done`. `comparisons/analyses_included.csv` copies every stamp it used.

### 5.4 File names

Rule: inside `not_adopted/<id>/`, a paper file name that ends in the mode (`_legacy` before the extension) gets `_<id>` instead; a paper name that holds no mode keeps its name (the folder says the analysis). Names that begin with `legacy_` (step 09; `legacy` names the procedure there) keep their names.

| Step | Paper (unchanged) | Optional example (`<id>` = `cov5-mnchisq-ctexcl` unless stated) |
|---|---|---|
| 04 | `tables/cpg_counts_by_coverage.csv` (legacy rows), `library_correlation_10x_legacy.csv`, `pca_scores_1x_legacy.csv`, `figures/figure_S1_legacy.png` | `not_adopted/cov5/04_cpg_coverage/tables/cpg_counts_by_coverage.csv` (mode `cov5`), `library_correlation_5x_cov5.csv`, `pca_scores_1x_cov5.csv`, `figures/figure_S1_cov5.png`. Not written in optional runs: `library_coverage.csv`, `genome_cpg_total.csv` |
| 05 | `dml_ploidy_legacy.csv`, `united_cpg_legacy.csv.gz`, `diffmeth_heat_legacy.csv.gz`, `dml_counts_by_mode.csv` (legacy row), `bed/dml_shared_legacy.bed` | `dml_ploidy_<id>.csv`, `united_cpg_<id>.csv.gz`, `diffmeth_heat_<id>.csv.gz`, `dml_counts_by_threshold_<id>.csv`, `dml_counts_by_mode.csv` (one row), `ct_exclusion_summary_<id>.csv` and `ct_excluded_cpgs_<id>.csv.gz` (only with `ctexcl`), `bed/dml_shared_<id>.bed`; cache `methylkit_cov5-mnchisq.rds` |
| 06 | `figures/figure_2_legacy.png`, `tables/pca_dml_legacy.csv` | `figures/figure_2_<id>.png`, `tables/pca_dml_<id>.csv` |
| 07 | `feature_overlap_by_track_legacy.csv`, `figures/figure_3_legacy.png` | `feature_overlap_by_track_<id>.csv`, `figure_3_<id>.png` |
| 08 | `dml_gene_overlap_ploidy_legacy.csv`, `gene_table_legacy.csv.gz`, `dml_go_annotation_heat_legacy.tsv.gz`, `annotation_summary.csv` | `dml_gene_overlap_ploidy_<id>.csv`, `gene_table_<id>.csv.gz`, `dml_go_annotation_heat_<id>.tsv.gz` (written for every analysis, step 09 needs it), `annotation_summary.csv` (its rows). Not written in optional runs: `go_map_id_repair.csv`, `gene_go_map.csv.gz`, `transcript_table.csv.gz` |
| 09 | `legacy_fisher_ploidy_BP.csv`, `figure_5_legacy.png`, `figure_5_values_legacy.csv`, `table_1_legacy.csv` | `legacy_fisher_ploidy_BP.csv` (in the analysis folder), `figure_5_<id>.png`, `figure_5_values_<id>.csv`, `table_1_<id>.csv` |
| 10 | none (not a paper step) | `not_adopted/legacy/10_go_enrichment_gene_level/tables/go_gene_level_ploidy_legacy.csv` (moved file); `not_adopted/corrected/10_.../figures/figure_go_bias_corrected.png`; per analysis also `go_gene_level_summary.csv`, `go_enriched_terms.csv`, `go_enriched_term_genes.csv`, `dml_probability_by_cpg_count.csv`, `legacy_top10_at_gene_level.csv`, `go_terms_by_analysis.csv` (its rows) |
| 11 | `named_terms_results_3_4.csv` (without the five `gene_level_*` columns), `semantic_summary.csv` and `semantic_terms_dropped.csv` (without gene-level rows), Figures 6 to 8 | `not_adopted/<id>/11_semantic_similarity/tables/named_terms_results_3_4_gene_level.csv` (set, go_id, five `gene_level_*` columns), `semantic_clusters_gene_level.csv`, `semantic_summary_gene_level.csv`, `semantic_terms_dropped_gene_level.csv`, `figures/figure_semantic_gene_level_<set>_<ont>_<id>` |
| 13 | `05_manuscript/manuscript_values.csv` (without `regenerated_corrected`, `corrected_differs`, `expression_corrected`), `manuscript_values_summary.csv` (without `change_in_corrected_mode`), `13_.../tables/expression_errors.csv` (without `expression_corrected`), `manuscript_items_copied.csv` (paper items) | `05_manuscript/not_adopted/manuscript_values_by_analysis.csv`, `..._summary_by_analysis.csv`, `05_manuscript/not_adopted/corrected/figures/Figure_2.png` |
| comparisons | | `not_adopted/comparisons/04_cpg_coverage/tables/cpg_counts_by_coverage.csv` (every analysis), `.../05_differential_methylation/tables/dml_counts_by_mode.csv`, `.../08_gene_annotation/tables/annotation_summary.csv`, `.../10_go_enrichment_gene_level/tables/go_terms_by_analysis.csv` and the other five step 10 cross tables, `figures/figure_go_terms_by_analysis.png` |

### 5.5 README table (root README shortened; `01_pipeline/README.md` replaces its section "Modes"; `04_output/not_adopted/README.md`)

| Analysis | Adopted in the paper? | Command | Outputs | What it changes | Issue |
|---|---|---|---|---|---|
| The paper | yes | `Rscript run_pipeline.R` | `04_output/<step>/`, `05_manuscript/` | at least 10 reads per CpG in each of the 22 libraries, no upper coverage filter, first line of each coverage file dropped, logistic regression without overdispersion correction, DML at q < 0.01 and a difference of more than 20 points | |
| Coverage filter of the Methods | no | `--coverage-filter` | `not_adopted/corrected/` | CpGs at or above the 99.9th percentile of coverage in a library removed; every line read | I-05, I-38 |
| 5x standard | no | `--coverage 5` | `not_adopted/cov5/` | 5 reads per CpG in each library (2,297,236 CpGs tested instead of 480,579) | Venkataraman et al. 2022 |
| Overdispersion correction | no | `--overdispersion` | `not_adopted/mnchisq/` | `overdispersion = "MN"`, `test = "Chisq"`; no other test | I-06 |
| C-to-T variant exclusion | no | `--exclude-ct-variants` | `not_adopted/ctexcl/` | CpGs at C-to-T variants removed from the CpG sets and, after the test, from the tested CpGs and the DML | I-08 |
| Any combination | no | the switches together | `not_adopted/<tokens joined by ->/` | | |
| Gene-level GO enrichment | no | `--gene-level-go` (with any of the above) | `not_adopted/<id>/10_go_enrichment_gene_level/` | step 10 and the gene-level part of step 11 | I-41 |
| Diagnostics | no | `--diagnostics` | `06_diagnostics/` | the checks and controls of the audit | |
| Side by side | no | made by every optional run | `not_adopted/comparisons/` | every analysis folder present, paper first | |

---

## 6. Changes by file and step

Counts are from the code at `671feca` (grep counts named; line counts are estimates).

### 6.1 `02_R_functions/analyses.R` (new, about 200 lines)

`analysis_spec(id)`, `analysis_id(tokens)`, `analysis_order(ids)`, `analysis_dir(a, step, sub)`, `analysis_path(a, step, sub, file)`, `paper_path(step, sub, file)`, `write_stamp(a, step, ...)`, `check_stamp(a, needs)`, `analysis_folders()` (ids of folders under `not_adopted/` with a stamp), `rebuild_comparison(step, file, order_keys = NULL)`, `ct_variant_table()` (reads the list, stops with the command to make it), `require_run_kind(kind, step_does)` (replaces `require_mode()`: message, exit 0 under Rscript, error when interactive), `diagnostic_modes()` (fixed `c("legacy", "corrected")` for the diagnostics, each checked to exist).

### 6.2 `02_R_functions/paths.R`, `logging.R`

- `paths.R`: guards and `paper_path()` of 4.6 (about 20 lines).
- `logging.R`: `start_step()` prints run kind, analysis id, its folder and the add-ons (line 14); `finish_step()` writes `session_info.txt` to `analysis_dir(config$analysis, step)` in optional runs (line 36; D2-1); `require_mode()` removed (lines 18-31).
- `methylkit_helpers.R`, `coverage_io.R`: no code change.

### 6.3 `run_config.R` (about 70 lines) and `run_pipeline.R` (about 100 lines)

As in 4.1 to 4.3. Remove lines 23-40 and 209-213 of `run_config.R`; lines 44-97 and 111-116 of `run_pipeline.R` are rewritten (parser, plan, run kinds, variant step, diagnostics runner, closing table).

### 6.4 Step 04 (about 25 sites)

- `a <- config$analysis; mode <- a$id`. The four loops over `config$modes` (lines 149, 240, 309, 349) run once. Lines 150 and 207 test `a$filter == "legacy"` instead of `mode == "legacy"`. Levels of distributions, correlation and PCA: `c(1, a$min_coverage)` (lines 149, 309). Counts by coverage keep every level of `config$coverage_levels`.
- With `a$ct_excluded`: after the first line is handled, drop the rows of each library whose `(chr, start)` is in `ct_variant_table()`; the unfiltered count (`cpg_in_any_library_unfiltered`) drops them too.
- Paths through `analysis_path()`; `library_coverage.csv` and `genome_cpg_total.csv` written only in a paper run (D2-2). Line 270 prints the rows of `a`.
- Optional run: `rebuild_comparison()` for the six spanning tables.

### 6.5 Step 05 (about 20 sites)

As in 4.5, plus `write_mode(fit, a)`:
- With `a$ct_excluded`: `keep <- !paste(united$chr, united$start) %in% paste(v$chr, v$start)`; drop the same rows from the united table and both test tables; renumber `united_row` after the drop; then apply the thresholds. q-values and p-values stay those of the fit on all CpGs (2022 study; variant spec 1.4 and 3.5).
- Line 238: compare with `dml_table()` before the drop (or drop the check for excluded analyses and replace it with `stopifnot(all(d$united_row <= nrow(united)))`); D2-9.
- Lines 263-271 (threshold sweep): computed from the excluded tables.
- With `ctexcl`, write `ct_exclusion_summary_<id>.csv` (contrast, tested CpGs before and after, DML before and after, removed) and `ct_excluded_cpgs_<id>.csv.gz` (the removed tested CpGs with p, q, difference per contrast and whether each was a DML). Not columns of `dml_counts_by_mode.csv` (D2-10).
- `dml_counts_by_mode.csv`: one row (the analysis) in its folder; `rebuild_comparison()` in optional runs. Stamp with the fingerprint md5.

### 6.6 Steps 06, 07 (about 7 and 11 sites)

Loops at line 79 (06) and 163 (07) run once for `a`; paths through `analysis_path()`. No change to the drawing code. Step 06 clusters every DML of a set with `pheatmap` (line 169): memory grows with the square of the set size; measured in run T3-b.

### 6.7 Step 08 (about 14 sites)

Lines 116-120 removed; loop at 122 runs once; line 177 `if (mode == "legacy")` becomes unconditional, with the file name `dml_go_annotation_<set>_<id>.tsv.gz` (`sprintf("dml_go_annotation_%s_%s.tsv.gz", set, a$id)`, which for the paper is today's name). The three mode-independent tables (lines 85, 93, 99) written only in a paper run. `annotation_summary.csv` per folder plus comparison.

### 6.8 Step 09 (about 17 sites)

Guard lines 98-105 replaced by `check_stamp(a, needs = "08_gene_annotation")`; line 131 `[mode == "legacy"]` becomes `[mode == a$id]`; 14 paths; figure and value-table names `figure_5_<id>`, `figure_9_<id>`, `figure_5_values_<id>.csv`, `figure_9_values_<id>.csv`, `table_1_<id>.csv`. Runs for optional analyses (handoff item 3: "steps 06 to 09 make their tables and figures for it and say so where a list is too short"); its behaviour on 14-locus lists is checked in T3-c.

### 6.9 Step 10 (about 20 sites)

`require_run_kind("optional", ...)` plus `config$addons$gene_level_go`. Loops at lines 80, 134, 145, 188, 223 run once. `gene_go_map.csv.gz` and the step 09 Fisher tables come from `paper_path()` (lines 76, 156-158). Lines 170 and 209 (`"corrected" %in% config$modes`) are replaced: each folder writes its own column or rows, and `rebuild_comparison()` builds `legacy_top10_at_gene_level.csv` (columns merged on `set, ontology, go_id` in analysis order), `go_terms_by_analysis.csv` and its figure, `go_gene_level_summary.csv`, `go_enriched_terms.csv`, `go_enriched_term_genes.csv`, `dml_probability_by_cpg_count.csv` from every folder that has step 10. With `legacy` and `corrected` present these equal today's committed tables byte for byte (T3-a2).

### 6.10 Step 11 (about 10 sites; edited by the hardening pass at 01:38)

- `require_mode()` at line 173 removed. The input check (lines 196-210) splits: paper inputs in a paper run; gene-level inputs (`go_gene_level_<set>_<id>.csv`, `go_enriched_terms.csv` of the analysis folder, the paper's `named_terms_results_3_4.csv`) in an optional run with (e).
- Paper run: sections up to 5 as now, without the merge with step 10 (lines 868-878) and without the five `gene_level_*` columns (line 884); no section 6.
- Optional run with (e): only the gene-level work, written under `analysis_dir(a, "11_semantic_similarity", ...)` with the names of 5.4 (D2-3). Implementation: guard the paper sections with `if (config$run_kind == "paper")` and the gene-level ones with `if (config$run_kind == "optional" && config$addons$gene_level_go)`; shared function definitions stay outside both.
- No paper expression reads a `gene_level` column or `semantic_summary.csv` (checked by `code/map_counts.py`: 0 legacy expressions name step 10 or `gene_level`), so the paper values do not move.

### 6.11 Comparison rebuild order

`rebuild_comparison()` binds the copies of every folder in analysis order (5.1). Two tables are interleaved today and are ordered after binding with a stable sort:
- `methylation_distribution_by_library.csv`: library (order of first appearance in the paper's table, which is the sample-sheet order), `min_coverage` ascending, analysis order, bin (order of first appearance).
- `methylation_distribution_by_group.csv`: grouping (`ploidy`, `treatment`, `all`), group (first appearance), `min_coverage` ascending, analysis order, bin.
Verified on the committed files (02:1x UTC): blocks of 10 bins alternate legacy and corrected within (library, level) and within (group, level). Every other spanning table is mode-major today (2 runs of the `mode` column each). A `cov5` analysis brings level 5, which sorts between 1 and 10.

### 6.12 Step 13 (about 80 lines)

- Paper part (paper run): as today for `expression_legacy`, with `tab()` reading paper folders, `config` the paper's, `diagnostic()` reading `06_diagnostics/`. A missing table stays an error (D2-5). Writes the paper files of 5.4. Items: rows with `analysis == "legacy"`.
- Optional part (optional run): for every id in `union(analysis_folders(), c("corrected", "cov5", "mnchisq", "ctexcl"), config$analysis$id)`:
  - absent folder: every row `not run: analysis not run`;
  - stamp not matching its id: every row `not run: folder stale`;
  - else evaluate `expression_optional` with `A = id`, `tab(step, file)` reading `not_adopted/A/<step>/tables/` (for `A == "legacy"`, falling back to the paper folder), `paper_tab()` reading paper folders, `diagnostic()`, and `config` with `min_coverage = a$min_coverage`. A missing input table gives `not run: <step>/<file> absent`; any other error is recorded in `not_adopted/13_manuscript_values/tables/expression_errors_by_analysis.csv` and stops the step after writing.
  - Output columns: `id, section, item, manuscript_value, analysis, regenerated, differs_from_manuscript, status`.
  - Items: rows with `analysis == "optional"` are templates with `{A}` in source and destination, copied for every folder whose source exists; rows with `analysis == "comparisons"` copied from `not_adopted/comparisons/`.
- Map reading (lines 87-91): `expression_corrected` becomes `expression_optional` in `colClasses`.

### 6.13 Maps (one-time script, `07_documentation/toggles_migration/rewrite_maps.R`)

- `manuscript_values_map.csv`: rename the column; in its 263 expressions replace each file name ending in `_corrected.` (201 sites in 190 rows, D2's count) with `sprintf("..._%s.<ext>", A)` and each literal `"corrected"` (56) with `A`; keep the 2 rows that read paper step 09 tables (D4.3-008, D4.7-005; through `paper_tab()`) and the 7 that read `dml_sensitivity/dml_counts.csv` with `mode == "legacy"`. Every other column byte-identical (T0).
- `manuscript_items_map.csv` (CRLF line endings kept): 16 `corrected` rows and the 4 step-10 bias-figure rows become `optional` templates; the 3 `both` rows and `legacy_top10_at_gene_level.csv` (row 59) become `comparisons`; allowed values `legacy`, `optional`, `comparisons` (line 164).

### 6.14 Variant step `01_pipeline/optional/ct_variant_loci.R` (new; method from `V2-method/spec.md`)

- R only (Rsamtools 2.18.0, GenomicAlignments 1.38.2, Biostrings 2.70.2 are installed; to be added to `renv.lock`). BS-Snper is not run; the rule of the 2022 study (variant spec 1.4): a CpG is excluded if a record at the C of the CpG on the plus strand has REF C and ALT T or "T,G", any FILTER, in the pooled call of all 22 libraries or in any one library.
- Inputs: `03_input_files/gannet_manifest_alignments.csv` (new, read only by this script; 22 BAM and 22 index rows with size and, where gannet offers one, md5; D2-6), the genome of the cache, every CpG of the reference.
- One alignment at a time: download into `<gannet cache>/alignments/`, check size, process, delete before the next. Before each download check free space of at least twice the file (largest 6.76 GB; 28 GB free at 02:1x UTC).
- Rule parameters in `run_config.R` under `ct_variants = list(...)`, set to the 2022 values, not offered as toggles (decision 10.8).
- Outputs (committed once): `04_output/not_adopted/ct_variants/tables/ct_variant_cpgs.csv.gz` (`chr`, `start` = 0-based position of the C as in the coverage files, `pos1` = start + 1, evidence columns of the variant spec 3.4), `ct_variant_run.csv` (settings, size of each alignment read, counts), `ct_variant_log.txt`. Steps 04 and 05 join on `(chr, start)`.

### 6.15 `06_diagnostics/`

- 5 scripts use `config$modes` 25 times (`check_methylkit_equivalence.R` 14, `dml_sensitivity.R` 6, `cohort_separation.R` 2, `library_methylation_level.R` 2, `feature_odds_matched.R` 1): replace with `diagnostic_modes()`, reading corrected outputs through `analysis_path(analysis_spec("corrected"), ...)`. 4 scripts name corrected files or the literal (`coverage_filter_variants.R` 13, `check_methylkit_equivalence.R` 12, `locus_nc047563.R` 6, `dml_sensitivity.R` 1 lines): same. `dml_sensitivity.R` line 430 then no longer meets a mode outside the list.
- New `06_diagnostics/diagnostics.csv`: `script, order, reads_analyses, needs_network, minutes, memory_gb, feeds_paper_values` (`TRUE` for `fastq_read_names` and `table1_group_means`).
- `check_legacy_outputs.R`: the one-line guard the hardening pass requested (handoff line 70).

### 6.16 Other files

- `.gitignore`: `04_output/not_adopted/*cov5*/05_differential_methylation/tables/*.csv.gz` (5x united and test tables, about 107 MB and 71 MB, above the 50 MB rule of `04_output/README.md` line 17).
- READMEs (root, `01_pipeline/`, `02_R_functions/`, `04_output/`, new `04_output/not_adopted/`, `05_manuscript/`, `06_diagnostics/`), `07_documentation/reproducibility_plan.md` section 4, `issues.md` (I-05, I-06, I-08; 2 path references): the table of 5.5 and run times. Docs reference the moved paths in 9 places outside the maps and generated tables (`git grep` at 02:1x UTC). Records under `07_documentation/agent_workflows/` are not edited (author, 00:45 UTC).
- `tests/` (new): `check_settings.R` (T0), `make_expected_paper.R` and `check_paper_outputs.R` (T1), `snapshot_paper.R` (T2), `check_step13_not_run.R` (T4). R only.

---

## 7. Order of work

### 7.1 Preconditions

1. The hardening pass has finished and is integrated, and the F-test removal and the repeat-track rename are committed: all three change bytes of committed files that the tests use as the baseline. Tag that commit `toggles-base` (BASE). `git status` clean.
2. Nothing else running on the machine during the runs of 7.2 steps 6 to 8 (15 GB, 4 cores; the other workflow needs up to 10 GB).

### 7.2 Steps

1. **Baseline.** `tests/baseline/paper_manifest_BASE.csv`: path, bytes, md5, md5 of the decompressed content for `.gz`, for every tracked file under `04_output/`, `05_manuscript/`, `06_diagnostics/` at BASE.
2. **Commit "move" (renames only, `git mv`):** 44 `*_corrected.*` files of 04 to 10 into `04_output/not_adopted/corrected/<step>/<sub>/`; the 19 files of step 10 into `not_adopted/legacy/10_...` (5 legacy files plus `session_info.txt`), `not_adopted/corrected/10_...` (5) and `not_adopted/comparisons/10_...` (8); `semantic_clusters_gene_level.csv` into `not_adopted/legacy/11_semantic_similarity/tables/`; the 22 non-paper files of `05_manuscript/` into `05_manuscript/not_adopted/{corrected,legacy,comparisons}/`. No content change, so git shows renames.
3. **Commit "split" (one-time R script `07_documentation/toggles_migration/split_tables.R`):** for each of the 8 spanning tables, write the paper copy (header plus legacy rows, `write_table()`), the corrected copy into `not_adopted/corrected/<step>/tables/` (header plus corrected rows) and the full BASE table into `not_adopted/comparisons/<step>/tables/`. Project step 11's `named_terms_results_3_4.csv` (drop five columns) and `semantic_summary.csv`, `semantic_terms_dropped.csv` (drop gene-level rows, writing them to the `_gene_level` files of `not_adopted/legacy/11_...`). Write the stamps of `not_adopted/legacy/` and `not_adopted/corrected/`. Run `rewrite_maps.R` (6.13). These files are the expected paper state EXPECTED; step 13's and step 01's outputs are not projected here, they come from the run of step 6.
4. **Commit(s) "code":** `analyses.R`, `paths.R`, `logging.R`; `run_config.R`; `run_pipeline.R`; steps 04 to 13; the variant step and its manifest; diagnostics; `.gitignore`; tests. Docs last.
5. **T0** (seconds), then step 13 alone in both configurations on the committed files (minutes): T-map and T4.
6. **Default run** `Rscript run_pipeline.R` with `WGBS_CORES=2`: refits the paper once (about 18 minutes in step 05). T1. Commit the paper outputs that changed by design (step 13 tables, `input_inventory.csv`) after checking them against 7.3.
7. **Optional runs** in this order, each followed by T2 and committed (5x large tables are git-ignored): T3-c `--overdispersion`; T3-a `--coverage-filter`; T3-a2 `--coverage-filter --gene-level-go`; T3-e `--gene-level-go`; C3 `--coverage-filter --overdispersion --gene-level-go`; variant list (`--exclude-ct-variants` builds it, hours); T3-d `--exclude-ct-variants`; C2 `--overdispersion --exclude-ct-variants --gene-level-go`.
8. **5x:** first the memory pilot (section 9) with `WGBS_CORES=1`, alone on the machine; then T3-b `--coverage 5` and C1 `--coverage-filter --coverage 5 --exclude-ct-variants`. Optionally the author's example `--coverage 5 --overdispersion --exclude-ct-variants`.
9. **T3-f** `--diagnostics`. Then a final default run and T1 again (proves that the optional runs left nothing the paper reads). Independent review.

### 7.3 Paper files whose bytes change once (all expected; each checked in T1)

| File | Change | Check |
|---|---|---|
| 6 spanning tables of step 04, `dml_counts_by_mode.csv`, `annotation_summary.csv` | corrected rows removed | header and every data line equal, in order, to the BASE lines with `mode == "legacy"` (`compare_paper_rows`, ported to R from D2's script) |
| `11_semantic_similarity/tables/named_terms_results_3_4.csv` | five `gene_level_*` columns removed | every remaining field and the header equal BASE as text |
| `semantic_summary.csv`, `semantic_terms_dropped.csv` | gene-level rows removed | remaining lines equal BASE lines in order |
| `05_manuscript/manuscript_values.csv` | 3 columns removed | remaining columns equal BASE as text |
| `05_manuscript/manuscript_values_summary.csv` | `change_in_corrected_mode` removed | same |
| `13_manuscript_values/tables/expression_errors.csv` | `expression_corrected` removed (header) | same |
| `13_manuscript_values/tables/manuscript_items_copied.csv` | only paper items | rows equal BASE rows with `analysis == "legacy"` minus the 3 step-10 items |
| `01_fetch_inputs/tables/input_inventory.csv` | rows of the two maps and the new alignment manifest | every other row equal |

`semantic_clusters_gene_level.csv` leaves the paper folder (moved). Nothing else changes; any other difference fails T1.

---

## 8. Tests

Byte comparisons use md5 of the file and, for `.gz`, of the decompressed content (a gzip header names the operating system, `logging.R` lines 49-50). PDFs and `session_info.txt` are left out of T1, as `04_output/README.md` allows; T2 includes them.

### T0 Static and unit tests (`tests/check_settings.R`, seconds, no data)

- `config$modes` appears nowhere in `01_pipeline/`, `02_R_functions/`, `run_config.R`, `run_pipeline.R`, `06_diagnostics/` (grep; D2-8). `require_mode` appears nowhere.
- `run_config.R` sourced for each of the 16 ids times `WGBS_GENE_LEVEL_GO` TRUE/FALSE times `WGBS_DIAGNOSTICS` none/all (64 combinations): `config$min_coverage == 10`; `analysis$min_coverage` is 5 exactly when `cov5` is a token; `fit_id` as in 4.2; `run_kind == "paper"` only for (`legacy`, FALSE, none); `analysis_path()` of an optional analysis always contains `not_adopted`, of the paper in a paper run never.
- Invalid values stop: ids `cov5-corrected`, `legacy-cov5`, `Cov5`, `5x`, empty; `coverage_standard = 7`; `min_coverage = 5` (frozen); `WGBS_ANALYSIS` set when `run_pipeline.R` starts.
- Guards: in an optional configuration `output_dir("05_differential_methylation", "tables")` and `manuscript_path("figures")` stop; `output_dir(..., "cache")` and `manuscript_path("not_adopted")` do not.
- `run_pipeline.R --plan` for the 11 commands of 4.4 prints the expected step lists (compared with a table in the test).
- Maps: in `manuscript_values_map.csv` every column except the renamed one is byte-identical to BASE (md5 per column); 263 non-empty `expression_optional`; no `_corrected.` and no literal `"corrected"` left; items map keeps CRLF.

### T-map Templates reproduce today's corrected values (minutes, no fit)

After 7.2 step 3, step 13's optional part with only the moved `corrected` and `legacy` folders present: for `A = "corrected"` the 263 values equal, as text, the column `regenerated_corrected` of BASE `05_manuscript/manuscript_values.csv`.

### T1 (i): a default run reproduces every committed legacy output

`tests/check_paper_outputs.R` after `Rscript run_pipeline.R`:
- every file of EXPECTED (BASE manifest minus moved files, with the projected files of 7.2 step 3) has the same bytes after the run; the step 13 and step 01 files of 7.3 satisfy their projection checks against BASE;
- `git status --porcelain 04_output 05_manuscript` shows no added or deleted path and no path under `not_adopted/`;
- the step 05 log says "fitting again" once (the code change) and, on a second run of step 05, "using the cached fit";
- `06_diagnostics/check_legacy_outputs.R` still gives 329 of 331.
Catches: a changed paper path, a changed model call (R3 of D1: the paper call must have its present arguments), a wrong reading of `config$min_coverage`, the dropping of rows or columns beyond 7.3.

### T2 (ii): no optional run changes a legacy output

`tests/snapshot_paper.R before|after`, around every run of 7.2 steps 7 to 9: path, size, md5 and modification time of every file under `04_output/` and `05_manuscript/` outside `not_adopted/` (PDFs and `session_info.txt` included), of `06_diagnostics/fastq_read_names.csv` and `table1_group_means.csv` (read by 15 paper rows), and the list of folders there. Before and after must be identical: no file added, removed, rewritten (mtime) or changed. The cache folder is left out. Catches D2-1, D2-2, D2-3 and any path the guard missed.

### T3 (iii): each toggle alone and three combinations

Every run: all steps of its plan finish with status 0; T2 holds; the stamp of the folder matches; `comparisons/analyses_included.csv` lists every folder; `comparisons/05_.../dml_counts_by_mode.csv` has one row per folder, paper first.

| Run | Command | Expected (checked) |
|---|---|---|
| T3-a | `--coverage-filter` | every moved file under `not_adopted/corrected/` (steps 04 to 09) byte-identical to its committed copy; the 8 tables of `comparisons/` byte-identical to the BASE spanning tables (with `legacy` and `corrected` the only folders with those steps); row 454,555 CpGs, 952 ploidy, 2,273 heat, 159 shared |
| T3-a2 | `--coverage-filter --gene-level-go` | the 5 moved step 10 corrected files and the 8 step 10 comparison files byte-identical to BASE |
| T3-b | `--coverage 5` | `cpg_analysed` 2,297,236 = `cpg_in_all_libraries` of `not_adopted/cov5/04_.../cpg_counts_by_coverage.csv` at level 5; fit `methylkit_cov5.rds`; peak memory and time of steps 05 and 06 recorded |
| T3-c | `--overdispersion` | 480,579 tested; 14 ploidy and 185 heat DML at the positions of `06_diagnostics/dml_sensitivity/dml_A1_ploidy_mn_chisq_legacy.csv` and `dml_A1_heat_mn_chisq_legacy.csv`, every one in the paper's lists; steps 08 and 09 finish on these small lists |
| T3-d | `--exclude-ct-variants` | no refit (log: cached fit `methylkit_legacy.rds`); `cpg_analysed` = 480,579 minus the tested CpGs in the list; DML lists = paper lists minus the listed loci, exactly; p and q of every kept row equal the paper's; `ct_exclusion_summary_ctexcl.csv` adds up; step 04 counts at every level reduced by the listed CpGs present |
| T3-d0 | unit run of step 05 with an empty variant list (scratch copy) | every table equals that of the base analysis except the `mode` value and the names |
| T3-e | `--gene-level-go` | the 5 moved step 10 legacy files and `semantic_clusters_gene_level.csv` byte-identical; no paper step started |
| T3-f | `--diagnostics` | every script finished or "not run: needs <x>"; the two tables read by paper rows unchanged (T2) |
| C1 | `--coverage-filter --coverage 5 --exclude-ct-variants` | id `corrected-cov5-ctexcl`; fit `methylkit_corrected-cov5.rds`, its `cpg_analysed` before exclusion 2,265,759 |
| C2 | `--overdispersion --exclude-ct-variants --gene-level-go` | id `mnchisq-ctexcl`; no refit (shares `methylkit_mnchisq.rds` of T3-c); DML = T3-c lists minus listed loci; step 10 and 11 gene-level files in its folder |
| C3 | `--coverage-filter --overdispersion --gene-level-go` | id `corrected-mnchisq`; new fit `methylkit_corrected-mnchisq.rds`; step 10 and 11 gene-level files in its folder |
| T3-cache | default run after T3-b | `methylkit_legacy.rds` reused (no refit); `methylkit_cov5.rds` untouched |

### T4 (iv): step 13 reports optional values as not run when absent

- In the state after 7.2 step 3 (only `corrected` and `legacy` folders): `manuscript_values_by_analysis.csv` has 263 rows `not run: analysis not run` for each of `cov5`, `mnchisq`, `ctexcl`.
- In a scratch copy of the repository (`git worktree add` outside the working tree, never in place): remove `not_adopted/corrected/10_go_enrichment_gene_level/` and run step 13 with `WGBS_ANALYSIS=corrected`: exactly the 55 rows whose expression reads step 10 become `not run: 10_go_enrichment_gene_level/<file> absent`, the other 208 are unchanged, exit status 0.
- In the same scratch copy, a stamp edited to another `min_coverage`: every row of that id `not run: folder stale`.
- In the same scratch copy, a paper table removed: the paper part of step 13 stops with an error (the paper rule is kept).

---

## 9. Run time and memory

Measured (repository): step 05 at 10x fits both modes in about 37 minutes at `WGBS_CORES=2` (about 18.5 per fit), 18 minutes at 4 cores, about an hour on one core; import about 4 GB; about 6.5 GB at 2 cores and 10 GB at 4 (`05_differential_methylation.R` lines 89-93). Step 04: 12 minutes, 9.5 GB (`01_pipeline/README.md`, quoted by both designs). MN with the chi-squared test was not slower than the paper's model at 10x (6.5 against 9.96 minutes for the heat contrast; D1, from the sensitivity cache).

Scaling for 5x: fit time proportional to tested CpGs, 2,297,236 / 480,579 = 4.78 (legacy reading) and 2,265,759 / 454,555 = 4.98 (corrected); rows held during the import 1.49 times (zr3534_1: 7,520,025 at 5x against 5,054,063 at 10x; D1).

| Run | Time (estimate unless measured) | Memory |
|---|---|---|
| default run, first after the change | step 05 refit about 18 min at 2 cores; whole pipeline as today plus that | as today (step 04 9.5 GB) |
| `--overdispersion` | step 04 12 min; step 05 about 15 to 20 min at 2 cores; steps 06 to 09 under 5 min | as the paper |
| `--coverage-filter` | step 04 12 min; step 05 about 18 min; 06 to 09 under 5 min | as the paper |
| `--exclude-ct-variants`, list present | step 04 12 min; step 05 about 1 min (cached fit) | as the paper |
| variant list, once | 3.6 to 10.9 hours, 109.1 GB downloaded one file at a time (variant scoping figure, quoted by D2; not verified here) | to be measured; disk 2 x 6.76 GB at most |
| `--coverage 5` | step 05 about 85 to 92 min at 2 cores, about 2.4 to 2.5 h on one core; steps 06 to 08 a few minutes each | not measured: estimates 6 to 8 GB import (D1), 8 to 12 GB on one core and 15 GB or more on two (D2). Must be measured |

Pilot before any full 5x run: `WGBS_CORES=1 /usr/bin/time -v Rscript 01_pipeline/05_differential_methylation.R` with `WGBS_ANALYSIS=cov5`, alone on the machine, after a run of step 04 for `cov5`. Stop if the maximum resident set size passes 13 GB. `--plan` prints the warning for every id with `cov5`. Disk per 5x fit: cache about 245 MB, united about 107 MB, each test table about 71 MB (git-ignored).

---

## 10. Decisions that belong to the author

Already decided and not reopened: the C-to-T exclusion is off by default and computed in R only; the correction uses the chi-squared test only; 10x stays the default coverage; the records under `07_documentation/agent_workflows/` are left as they are.

| # | Decision | Recommendation | Reason |
|---|---|---|---|
| 10.1 | Names of settings and switches (`coverage_filter_of_methods`/`--coverage-filter`, `coverage_standard`/`--coverage 5`, `overdispersion_chisq`/`--overdispersion`, `exclude_ct_variants`/`--exclude-ct-variants`, `gene_level_go`/`--gene-level-go`, `diagnostics`/`--diagnostics`) | as listed | each says what it changes; `coverage_standard` uses the author's word "standard" |
| 10.2 | Token of (a): keep `corrected` or rename (for example `methodsfilter`) | keep `corrected` | 44 committed files, 263 expressions and issues.md use it; renaming changes every name |
| 10.3 | Move the 44 corrected files, the 19 step 10 files and 22 files of `05_manuscript/` into `not_adopted/` | move | handoff item 3 asks for names or folders that say "not adopted"; after the move a path alone says whether a file is the paper's |
| 10.4 | Accept the one-time byte changes of 7.3 (14 paper files lose optional rows or columns; the paper rows and values stay byte-identical) | accept | the alternative keeps optional rows inside paper tables, against the standing rule of 2026-11-19 |
| 10.5 | Step 11's gene-level parts leave the paper (five columns of `named_terms_results_3_4.csv`, one row of `semantic_summary.csv`, `semantic_clusters_gene_level.csv`) | yes | they read step 10, which is not adopted; no paper value reads them |
| 10.6 | Freeze the five paper settings in `run_config.R` | yes | 22 paper values read `config$min_coverage` and 21 `config$dml_difference`; an edit there would change the paper silently. Cost: thresholds can no longer be explored by editing; that would be a new toggle |
| 10.7 | Apply the exclusion to the CpG sets of step 04 and to the background of step 07 as well as to the DML after the test | yes | as the 2022 study and the author's instruction of 00:15 UTC |
| 10.8 | Sub-settings of the variant rule (C only or C and G; any FILTER or PASS only; pooled over 22 libraries or per experiment) | the 2022 rule only, fixed in `config$ct_variants`; no toggle | "as in Venkataraman"; each further sub-setting would need its own token so that two lists never share a folder |
| 10.9 | One analysis per optional run, several by several runs | yes | each command has one meaning; `comparisons/` shows them side by side |
| 10.10 | 15 paper values read two diagnostic tables (`fastq_read_names.csv`, `table1_group_means.csv`) that the default run does not regenerate | keep them as committed inputs now; T2 proves `--diagnostics` does not change them; consider promoting `table1_group_means.R` into a paper step later | the default run then depends only on committed files; promoting is a separate change of the paper run |
| 10.11 | 5x per-CpG tables git-ignored (a clone must rerun step 05 to get them) | yes | 107 MB and 71 MB exceed the 50 MB rule and GitHub's 100 MB limit |
| 10.12 | Folder of the tests: `tests/` at the root | yes | they are not diagnostics and must not run with `--diagnostics` |

---

## 11. What was verified, and how; what is assumed

Verified by me (read only, 02:09 to 02:20 UTC, HEAD `671feca`):
- `run_config.R` lines 23-40 (modes) and 209-213 (check); `run_pipeline.R` lines 35, 44-69, 78-97, 111-116; `logging.R` lines 12-41 (`finish_step()` writes to `output_dir(step)`); `paths.R` lines 13-41.
- Step 05 lines 116-163 (`fit_mode`), 175-188 (fingerprint), 195-219 (cache), 222-282 (`write_mode`, check at 238, sweep at 263-271), 285-291; header lines 89-93 (times and memory).
- Step 04 lines 145, 149-150, 165, 193, 207, 240, 264-288, 309, 338-339, 349, 377; step 08 lines 85, 93, 99, 116-122, 177, 204; step 09 lines 98-105, 131; step 10 lines 66-80, 145-158, 170-220; step 11 lines 173, 196-216, 868-947; step 13 lines 78-170; step 01 lines 46-47; `methylkit_helpers.R` lines 29-50 and 101-104.
- `config$modes` uses: 01_pipeline 37 lines (04: 4, 05: 2, 06: 1, 07: 1, 08: 6, 09: 5, 10: 7, 11: 3, 12: 1, 13: 7), `run_pipeline.R` 4, `logging.R` 4, diagnostics 25 in 5 scripts (`grep -c`).
- Values map (`code/map_counts.py`, run with `python3 -I`): 545 rows; 498 legacy and 263 corrected expressions; 51 corrected expressions read spanning tables (counts in D1-1); 55 corrected expressions read step 10; no legacy expression reads step 10 or names `gene_level`; 15 legacy and 14 corrected expressions call `diagnostic()`.
- Committed tables: CpGs in all 22 libraries 2,297,236 and 480,579 (legacy, 5x and 10x), 2,265,759 and 454,555 (corrected); DML 1,038 / 2,351 / 200 (legacy) and 952 / 2,273 / 159 (corrected); MN with chi-squared at 10x 14 ploidy and 185 heat DML (`06_diagnostics/dml_sensitivity/dml_counts.csv` rows 6-7). Row order of the 8 spanning tables (mode-major except the two methylation distribution tables, 10-bin blocks alternating).
- 44 tracked files with `corrected` under `04_output/`; 19 tracked files in step 10; 70 tracked files under `05_manuscript/`; items map 68 lines with CRLF rows for corrected and step 10 items (lines 6-59 listed).
- Diagnostics write only into `06_diagnostics/` (grep of `fwrite`/`write_table`/`saveRDS`/`save_figure`; no write to `04_output`).
- Machine: 15 GB memory, 10 GB available, 4 cores, 28 GB free disk. Alignments: 22 BAM files, 109.1 GB, largest 6.76 GB (`V1-data/bam_files.csv`, a file of the variant scoping agent). Rsamtools 2.18.0, GenomicAlignments 1.38.2, Biostrings 2.70.2, methylKit 1.28.0 installed.
- Handoff item 3 (line 97) and item 2b (line 95) of `07_documentation/session_handoff.md`.

Taken from the proposals and not re-verified: D2's counts of 201 file names in 190 rows and 56 literals in the corrected expressions; D1's row counts of zr3534_1 and the MN timings; the 3.6 to 10.9 hour estimate of the variant step; every 5x time and memory figure (estimates by scaling, no measurement exists). Assumed: gzip output of the same R build on the same OS is byte-stable (the repository already relies on it); the variant rule of `V2-method/spec.md` is the one to build.

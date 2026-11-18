# WGBS-gigas-ploidy-desiccation

DNA methylation in diploid and triploid Pacific oysters (*Crassostrea gigas*) after an aerial heat exposure: global 5-mC by ELISA and locus-level methylation by whole-genome bisulfite sequencing (WGBS).

- **Manuscript (draft):** [*Epigenetic divergence in the thermal stress response of diploid and triploid Pacific oysters*](https://docs.google.com/document/d/17mcGDI-TWmU4vgBXmiXmeofe4qEuFH5inBKBHhG9tzg/edit), M. N. George, R. Jain, S. B. Roberts
- **Lead:** Matt George
- **Status:** The repository is being rebuilt as a reproducible, R-only pipeline on the branch `reproducibility-audit`. The audit and the folder reorganization are done; the pipeline scripts are being added step by step. Until they are complete, the analysis that produced the current manuscript draft is the legacy code in `08_archive/`.

> **Where the work stands.** Read [`issues.md`](issues.md) first: it lists 37 problems found in the audit, three of which bear on the main conclusions (the heat treatment is confounded with sequencing batch in the WGBS data; the "interaction" loci are an overlap of two lists rather than a tested interaction; the GO enrichment counts transcript isoforms as genes). The plan and the decisions needed from the authors are in [`07_documentation/reproducibility_plan.md`](07_documentation/reproducibility_plan.md).

---

## The study in brief

Adult diploid and triploid oysters were either held as controls or exposed to air at 27 C for 24 h ("heat"; older files call this treatment "desiccation"). Ctenidia were sampled.

| Assay | Samples | Design |
|---|---|---|
| Global 5-mC ELISA | 32 animals, 8 per ploidy by treatment group | Controls and heat-treated animals from the same experiment, one plate |
| WGBS | 22 libraries: 5 diploid and 5 triploid heat-treated, 6 diploid and 6 triploid controls | Heat-treated libraries are Zymo project zr3534 (the heat experiment); control libraries are Zymo project zr3644 (animals from a separate experiment). Four of the ten heat-treated animals also received 1 h at 45 C. |
| qPCR | 40 animals, 9 to 12 genes | Same experiment as the ELISA |

Sample names in the heat experiment (D = diploid, T = triploid; suffix `-C` = ctenidia):

| Animals | Treatment |
|---|---|
| D01 to D08, T01 to T08 | Control |
| D09, D10, T09, T10 | Control, then 1 h acute heat shock at 45 C |
| D11 to D18, T11 to T18 | Heat (24 h aerial emersion at 27 C) |
| D19, D20, T19, T20 | Heat, then 1 h acute heat shock at 45 C |

## Repository layout

Numbered folders follow the order of the work: scripts, helpers, inputs, outputs, manuscript products, checks, documentation, archive. Each folder has a `README.md` with a file inventory.

```text
WGBS-gigas-ploidy-desiccation/
├── 01_pipeline/          Numbered R scripts, one per analysis step
├── 02_R_functions/       Helper functions sourced by every step
├── 03_input_files/       Inputs kept in GitHub, and the manifest of inputs fetched from gannet
├── 04_output/            Results, one subfolder per pipeline step
├── 05_manuscript/        Figures and tables as numbered in the manuscript; value-by-value traceability table
├── 06_diagnostics/       Checks of regenerated results against legacy results and manuscript values
├── 07_documentation/     Plan, audit notes, upstream processing records, protocols
├── 08_archive/           Everything from before the reorganization, original paths preserved
├── issues.md             Issue log: what is wrong, the evidence, what needs deciding
├── README.md             This file
├── .gitignore
└── WGBS-gigas-ploidy-desiccation.Rproj
```

| Folder | Contents | README |
|---|---|---|
| `01_pipeline/` | The twelve analysis steps and how a step is written | [01_pipeline/README.md](01_pipeline/README.md) |
| `02_R_functions/` | Shared helpers | [02_R_functions/README.md](02_R_functions/README.md) |
| `03_input_files/` | Sample sheets, ELISA plate reads, qPCR Cq values, Bismark reports, genome feature tracks, annotation tables, and the gannet manifest | [03_input_files/README.md](03_input_files/README.md) |
| `04_output/` | Tables and figures written by the pipeline | [04_output/README.md](04_output/README.md) |
| `05_manuscript/` | Final figures, tables, supplementary lists, and the table tying each manuscript number to its source | [05_manuscript/README.md](05_manuscript/README.md) |
| `06_diagnostics/` | Reproduction checks | [06_diagnostics/README.md](06_diagnostics/README.md) |
| `07_documentation/` | Plan, audit, provenance | [07_documentation/README.md](07_documentation/README.md) |
| `08_archive/` | Legacy code, outputs, figures, slides, with a guide to which legacy script produced which manuscript item | [08_archive/README.md](08_archive/README.md) |

At the root: `run_config.R` is the single control surface (thresholds, analysis modes, cache location). Still to be added as the pipeline arrives: `run_pipeline.R` (orchestrator) and the `renv` lockfile.

### How paths work

Scripts build every path with `here::here()`, which anchors on the repository root (found through the `.Rproj` file), not on the location of the script or the working directory. A script therefore runs unchanged from RStudio or from a terminal, on any operating system, and no script calls `setwd()`. If a numbered folder is renamed, the matching string in `02_R_functions/paths.R` is the only place to update.

## Where the data live

**GitHub first.** Any input that fits in GitHub is in `03_input_files/` (ceiling: 50 MB per file, after gzip where possible). That covers everything except the per-library coverage files.

**gannet for the large files.** The 22 Bismark coverage files (about 440 MB each, 9.1 GB in total) stay at <https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/cov_files/>. They are listed with md5 checksums in [`03_input_files/gannet_manifest.csv`](03_input_files/gannet_manifest.csv), downloaded once by the first pipeline step into a git-ignored cache, and verified.

Other locations, for reference:

| What | Where |
|---|---|
| Alignment outputs (BAM, reports) for all libraries | gannet `/panopea/061021-big/` (June 2021 run, used for the analysis) and `/panopea/030521-ronrosM/` (March 2021 run, superseded) |
| Genome feature tracks, original copies | OSF, doi:[10.17605/OSF.IO/YGCTB](https://doi.org/10.17605/OSF.IO/YGCTB) |
| Reference genome | Roslin v1, RefSeq GCF_902806645.1, with mitochondrial genome NC_001276.1; gannet `/panopea/Cg-roslin/` |
| Raw reads | NCBI SRA BioProject PRJNA678408 for the heat-treated libraries; the accession for the control libraries is not recorded here (issue I-36) |
| Legacy intermediates and saved workspaces | gannet project directory; inventory in `07_documentation/audit/gannet_manifest_2026-11-18.tsv` |

## Running the pipeline

Not yet runnable end to end; the scripts are being added (see the status column in [`01_pipeline/README.md`](01_pipeline/README.md)). When complete:

1. Install R (4.3 or later) and open `WGBS-gigas-ploidy-desiccation.Rproj`, or start R in the repository root.
2. Restore the package library with `renv::restore()`.
3. Run `Rscript run_pipeline.R`, or source the numbered scripts in `01_pipeline/` in order.

The first run downloads 9.1 GB from gannet. Steps 04 and 05 need about 12 GB of memory; all other steps run from files in the repository. A Linux setup script for a clean machine is in `07_documentation/cloud_environment_setup.sh`.

## Conventions

- The treatment is "heat" or "heat-treated" in all new code, outputs, and text. "Desiccation" survives only in `08_archive/` and in this repository's name.
- Legacy dataset labels: `ronit` = the heat-treated animals alone; `yaamini` = the control animals alone; `both` or `all-gigas` = the 22-library comparison reported in the manuscript.
- Diploid is drawn in blue and triploid in orange-red throughout.

## Related material

- Research poster by Ronit Jain: `08_archive/presentations/Triploid-Diploid-Oyster-Poster_RonitJain.pptx`
- Earlier manuscript draft by Ronit Jain: [Google Doc](https://docs.google.com/document/d/1nwY9I3pVzF5Xlfdzb7SdO1QKaNp5qyj0DuPCCa_YXi8/edit?usp=sharing)
- Global DNA methylation analysis by Shelly Trigg (2019), the origin of `08_archive/ELISA/GlobalDNAMeth_Polyploids.*`: [script](https://github.com/shellytrigg/C_gigas/blob/master/Polyploids/GlobalDNAMeth_Polyploids.R), [rendered analysis](https://github.com/shellytrigg/C_gigas/blob/master/Polyploids/GlobalDNAMeth_Polyploids.md)

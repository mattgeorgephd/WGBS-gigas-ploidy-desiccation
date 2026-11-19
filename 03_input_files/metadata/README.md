# metadata

Sample sheets. `wgbs_libraries.csv` is the single description of the 22 WGBS libraries that every pipeline step reads; no script carries its own list of samples or its own group labels.

## Files

| File | Contents |
|---|---|
| `wgbs_libraries.csv` | One row per WGBS library. Built by `build_wgbs_libraries.R`; do not edit by hand. |
| `build_wgbs_libraries.R` | Builds `wgbs_libraries.csv` from the two raw files and checks that they agree. Run with `Rscript 03_input_files/metadata/build_wgbs_libraries.R`. |
| `raw/legacy_metadata.xlsx` | The legacy `metadata.xlsx`. Sheet `list` gives the 22 libraries in the order the legacy methylKit analysis used; sheets `list2`, `list3`, and `all` describe subsets and other *C. gigas* datasets that are not part of this study. Its `path` column holds paths on the original Windows machine and is not used. |
| `raw/zr3534_wgbs_info.csv` | Sequencing sample sheet for the ten heat-treated libraries: animal ID, ploidy, treatment flags, library kit, SRA accessions. |

## Columns of `wgbs_libraries.csv`

| Column | Meaning |
|---|---|
| `library_id` | Sequencing library ID (Zymo project and number), for example `zr3534_1` |
| `legacy_sample_index` | 1 to 22: the sample ID and column order used in the legacy methylKit objects and outputs |
| `zymo_project` | `zr3534` or `zr3644` |
| `cohort` | `heat_experiment` or `control_experiment`: which experiment the animal came from |
| `legacy_label` | `ronit` or `yaamini`: the label legacy files use for the cohort |
| `tissue` | ctenidia for all 22 (flash frozen and stored at -80 C before DNA extraction, not fixed or embedded; stated by the author on 2026-11-19, without distinction between the two experiments) |
| `ploidy` | `diploid` or `triploid` |
| `treatment` | `heat` or `control`. The draft's Methods describe the heat treatment as 24 h of aerial emersion at 27 C; the author states on 2026-11-19 that the 24 h emersion is not confirmed, and that all heat-treated animals count as one treatment group |
| `animal_id` | Animal ID in the heat experiment (D = diploid, T = triploid); empty for the control cohort, whose animal IDs are not recorded in this repository |
| `acute_heat_shock` | `yes` for the four heat-treated animals that also received 1 h at 45 C (D19, D20, T19, T20) |
| `library_kit`, `sra_bioproject`, `sra_experiment` | As recorded for zr3534; empty for zr3644 |
| `cov_file` | Name of the library's coverage file on gannet (see `../gannet_manifest.csv`) |

## Things to know

- Treatment is aligned with cohort and with sequencing project: every heat library is zr3534 and every control library is zr3644. See issue I-01 in [`issues.md`](../../issues.md).
- Animal numbering in the heat experiment: 01 to 08 control, 09 and 10 control plus acute heat shock, 11 to 18 heat, 19 and 20 heat plus acute heat shock. Only heat-treated animals from that experiment were sequenced.
- Missing for the control cohort: animal IDs, source and holding details, accessions (issue I-36).
- Not recorded for any of the 22 animals: sex (author, 2026-11-19). Not known at present: how the triploids were produced and from how many parents (issue I-49).

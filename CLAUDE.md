# CLAUDE.md

## Overview

This repo maps **Northern Leaf Blight (NLB) resistance QTL in maize**. Phenotypes come from
multi-year, multi-location field disease scoring of the **IBM (B73 × Mo17)** population — the
recombinant inbred lines (RILs) and their reciprocal backcrosses to B73 and Mo17 — plus a set
of **NAM** lines. The analysis fits mixed models to get line BLUPs, quantifies heterosis
(mid- and best-parent), and runs linkage QTL mapping (single-QTL MQM scans and two-QTL/epistasis
`scantwo`) per population. The work is organized as a numbered R pipeline in `scripts/`.

## Pipeline — run in order (`scripts/01`–`07`)

Each stage writes files consumed by the next.

1. **`01_make_supp_data.R`** — Assembles raw field-scoring CSVs (NLB 2022 Clayton, 2023 Clayton,
   2024 Clayton, 2024 Illinois) into the multi-sheet supplemental workbook
   `data/nlb_mdh_file_s1.xlsx`. Computes weighted mean disease (`wmd`) and plot x/y coordinates
   (`x_y_coords()` helper). Also reads height/DTA tables.
2. **`02_get_blups.R`** — Reads the workbook and fits mixed models (`lme4` / `blme`, with an
   optional `brms` heterogeneous-variance model) to get line BLUPs for the IBM populations (RIL,
   B73 BC, Mo17 BC) and NAM. Also computes heritability, within-year rep correlations, cross-env
   BLUP correlations, and heterosis (MPH/BPH via `add_heterosis()`). Writes
   `analyses/IBM_NLB_BLUPs.csv` (+ `_RILs`, `_B73BC`, `_Mo17BC` splits), `analyses/NAM_NLB_BLUPs.csv`,
   and `analyses/line_blups_envs.csv`.
3. **`03_plot_blup_distributions.R`** — `ggplot2` / `cowplot` frequency-polygon figures of BLUP
   distributions and MPH/BPH heterosis → `figures/*.pdf`.
4. **`04_make_cross_file.sh`** — Shell driver that calls `scripts/convert_rqtl2_to_rqtl.R` to
   merge genotypes + BLUP phenotypes + the genetic map (`data/ibm302map.csv`) into R/qtl cross
   CSVs in `analyses/` (`RIL_cross.csv`, `B73_cross.csv`, `Mo17_cross.csv`, plus old-data and
   per-environment BLUP crosses).
5. **`05_qtl_analysis.R`** — R/qtl **MQM** scan (`mqmaugment` → `mqmautocofactors` → `mqmscan`)
   with batched permutations for significance thresholds. Writes `analyses/qtl_analyses/<name>.RDS`
   (a list of `scan` + `permutations`); permutation batches go to `analyses/qtl_analyses/tmp/`.
6. **`06_identify_qtl.R`** — Loads the MQM `.RDS`, converts to R/qtl2, and calls main-effect QTL
   peaks with `qtl2::find_peaks` against per-trait permutation thresholds, per population.
   Writes `analyses/main_effect_peaks.csv`.
7. **`07_epistatic_qtl.R`** — R/qtl `scantwo` two-QTL / epistasis scan with batched permutations.
   Unlike 05–06, this is a **CLI**: `Rscript scripts/07_epistatic_qtl.R <population> [permutations]
   [cores]` where `<population>` is `RIL`, `B73_BC`, or `Mo17_BC` (permutations default 100, cores
   default 4). A `presets` list holds the per-population `input_file`/`genotype`/`na.strings`/
   `crosstype`/`phenos`. Both the main scan and the permutations are **resumable**: the main scan
   is cached at `analyses/qtl_analyses/tmp/<name>.scantwo.tmp.RDS`, and permutations run in batches
   of 10 to `tmp/<name>.trait.<pheno>.batch.<j>.scantwo.perm.tmp` (existing files are skipped on
   re-run). Writes `analyses/qtl_analyses/<name>_scantwo.RDS` (a list of `scan` + `permutations`).

## How these scripts are run

- **Run from the project root** (`~/projects/nlb_mdh`). Scripts mix root-relative paths
  (`analyses/...`, `data/...`) with hardcoded absolute paths (`~/projects/nlb_mdh/...`, and a
  few legacy `~/projects/mdh_qtl/...` inputs), so the working directory matters.
- **Interactive, per-population execution (05–06).** Scripts 05 and 06 have their `commandArgs()`
  CLI parsing commented out. The active pattern is to re-assign the config variables — `input_file`,
  `genotype`, `alleles`, `na.strings`, `crosstype`, `phenos`, `output_dir`, `output_name` — in a
  block per population and run that block interactively. These scripts are **not** meant to be
  `source()`-d top to bottom: later blocks overwrite the same variables, so the last block wins.
- **CLI execution (07).** Script 07 has been converted to a proper CLI driven by `commandArgs()`
  and a `presets` list — pass the population (and optionally permutation count and core count) as
  arguments instead of editing the script (see step 7 above).
- A `snakefile.smk` exists but the numbered scripts are the current workflow.

## Conventions & glossary

- **`wmd`** — weighted mean disease (a disease-progress score from repeated field ratings).
  Some models flip it to a resistance scale as `100 - wmd`.
- **Line naming** — `M0####` = an IBM RIL; an `x` denotes a hybrid/cross (e.g. `M0230xMo17`);
  founders are `B73`, `Mo17`, and the F1 `B73xMo17`. `02_get_blups.R` normalizes a couple of
  reciprocal/alternate-parent name variants.
- **Populations & R/qtl cross codes**:
  - RIL — `crosstype = "ril"`, genotypes `A`/`B`, `na.strings = "-"`
  - B73 BC — `crosstype = "bc"`, genotypes `AA`/`AB`, `na.strings = "A-"`
  - Mo17 BC — `crosstype = "bc"`, genotypes `AB`/`BB`, `na.strings = "-B"`
  - NAM — 26 NAM founders crossed to B73/Mo17 (BLUPs only; not QTL-mapped here)
- **Directory layout**:
  - `data/` — raw field CSVs, the supplemental workbook, genetic map (`ibm302map.csv`), genotype files.
  - `analyses/` — derived tables and QTL outputs; `analyses/qtl_analyses/` holds the `.RDS`
    scan objects and `analyses/qtl_analyses/tmp/` the permutation batches.
  - `figures/` — output PDFs.
  - `scripts/` — the code.
  - `nlb_mdh.Rproj` — RStudio project.

## Tooling / dependencies

R, with:
- QTL mapping: `qtl`, `qtl2`, `qtl2convert`; parallel permutations via `parallel` / `snow`.
- Mixed models: `lme4`, `blme`, optional `brms`.
- Data & I/O: `dplyr`, `tidyr`, `stringr`, `lubridate`, `openxlsx`, `readxl`.
- Plotting: `ggplot2`, `cowplot`, `gridGraphics`.

## Scope note

The numbered `01`–`07` pipeline above is the current workflow. The other, unnumbered scripts in
`scripts/` (e.g. `run_rqtl*.R`, `get_*.R`, `convert_rqtl2_to_rqtl.R`, `bayesian_blups.R`,
`genomic_prediction.R`, `gp_helper_functions.R`) are older/superseded or supporting helpers and
are not documented here yet — except `convert_rqtl2_to_rqtl.R`, which the numbered pipeline still
calls from step 04.

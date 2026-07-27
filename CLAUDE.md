# CLAUDE.md

## Overview

This repo maps **Northern Leaf Blight (NLB) resistance QTL in maize**. Phenotypes come from
multi-year, multi-location field disease scoring of the **IBM (B73 × Mo17)** population — the
recombinant inbred lines (RILs) and their reciprocal backcrosses to B73 and Mo17 — plus a set
of **NAM** lines. The analysis fits mixed models to get line BLUPs, quantifies heterosis
(mid- and best-parent), and runs linkage QTL mapping (single-QTL MQM scans and two-QTL/epistasis
`scantwo`) per population. The work is organized as a numbered R pipeline in `scripts/`.

## Pipeline — run in order (`scripts/01`–`11`)

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
8. **`08_identify_epistatic_qtl.R`** — Calls significant epistatic pairs from a `07` `scantwo`
   RDS. CLI: `Rscript scripts/08_identify_epistatic_qtl.R <population> [alpha]` (alpha default
   0.05). Per phenotype, gets every local-maximum chr-pair (`summary(scan, what="int",
   thresholds=0)`) and a permutation p-value for each (`mean(perms$int >= lod.int)`); flags
   `same_chr_close` (same chromosome, peaks <20 cM apart — a linked-marker artifact, not real
   epistasis) and buckets each pair into `sig_level` ∈ `{alpha, 0.20, "ns"}`. Writes/merges
   `analyses/epistatic_peaks.csv` (replacing that population's rows on re-run).
9. **`09_estimate_qtl_effects.R`** — Fits additive (`a`) and dominance (`d`) effects at the
   `06`/`08` peaks with `fitqtl(get.ests=TRUE)`, adding `Qi:Qj` terms for any genuine
   (`sig_level=="0.05"`, not `same_chr_close`) epistatic partner from `08`. CLI:
   `Rscript scripts/09_estimate_qtl_effects.R <population>`. Reports everything on one
   convention — the B73-allele effect on the disease-scale BLUP, positive = B73 raises disease —
   which requires per-population sign handling of R/qtl's internal genotype-code trap (see
   `plans/identify_qtl_and_effects.md`): RIL peaks give `a` directly; a BC's MPH peaks give `d`;
   a BC's BLUP peaks give a confounded `d-a` (B73 BC) or `-a-d` (Mo17 BC) contrast, reported with
   `estimate_type="confounded"` (never as `a`/`d`). Each `Qi:Qj` term also gets its own output row
   (`effect_class="epistatic"`, `chr`/`pos` = locus 1, `chr2`/`pos2` = locus 2, `lod` = the pair's
   `lod.int` from `08`, no sign flip — both loci share the same population/trait so the flip
   squares to +1). Every population additionally borrows a dummy-QTL fit of its own gene-action
   parameter (RIL → BLUP `a`; either BC → MPH `d`) at every QTL position from `06` where it lacks
   an independent peak on that parameter — so every position gets an `a` and both BCs' `d`s, real
   or borrowed. A boolean `borrowed` column marks these (not `lod = NA`, which is now reserved for
   genuinely absent info: borrowed rows still get a real `fitqtl_lod`, just no `lod` from `06`/`08`
   since they weren't an independently significant peak). Writes/merges `analyses/qtl_effects.csv`
   (columns: `cross, trait, chr, pos, chr2, pos2, ci_lo, ci_hi, lod, estimate, estimate_type,
   effect_class, borrowed, fitqtl_lod`).
10. **`10_qtl_gene_action.R`** — No population arg; combines all of `09`'s output. Filters to
    `effect_class=="main"` first (epistatic rows have NA CIs and describe a locus pair, not a
    single QTL, so they're excluded from colocalization). Colocalizes the remaining peaks into a
    `qtl_id` per chromosome via connected components of overlapping confidence intervals (a base-R
    `find_overlaps()` on `[ci_lo, ci_hi]`, doing what `IRanges::findOverlaps` would — written this
    way to avoid a Bioconductor dependency), joins each cluster's `a` with its B73/Mo17 BC `d`
    (real value preferred over borrowed when both exist), and classifies gene action from `d/a`
    (cutoffs 0.2/0.8/1.2 → additive/pd/dominant/od/ud, direction-aware). Every colocalized QTL now
    gets a classification, real or borrowed; `a_sig`/`B73_d_sig`/`Mo17_d_sig` (`= !borrowed`) mark
    which. Writes `analyses/qtl_gene_action.csv`.
11. **`11_genome_wide_effect_scan.R`** — Sliding whole-genome effect profile: at every imputed
    map position, fits a dummy QTL alongside that population's real peaks (dropping any real
    peak within 20 cM of the test position) plus its genuine `08` epistatic partners, and
    records the dummy QTL's effect. CLI: `Rscript scripts/11_genome_wide_effect_scan.R
    <population>`. Writes/merges `analyses/qtl_effects_whole_genome.csv`, one effect column per
    population/phenotype (e.g. `NLB_ril_effect`, `NLB_b73_effect`, `NLB_b73_mph_effect`,
    `NLB_mo17_effect`, `NLB_mo17_mph_effect`) — plot effect vs. LOD along the genome. Positions
    are rounded to 4 decimals so repeated per-population runs join cleanly on `(chr, pos)`
    despite `write.csv`/`read.csv` precision loss on re-read.
12. **`12_epistatic_effect_table.R`** — No population arg; builds a genotype-class effect table
    for the significant epistatic pairs (`sig_level=="0.05" & !same_chr_close` in
    `epistatic_peaks.csv`, currently 6 pairs: 3 RIL + 3 B73_BC). For each pair, calls raw R/qtl
    `effectplot()` genotype-class means/SEs (no fitqtl / main-effect adjustment — faithful to the
    original published table this distills), centers on the population mean, and reports
    B73-allele-count genotype cells (`RIL` 0/2, `B73_BC` 1/2, `Mo17_BC` 0/1 — per-population
    label→count map handles R/qtl's internal-code trap, see `plans/epistatic_effect_table.md`).
    Writes `analyses/epistatic_effects_long.csv` (tidy, one row per genotype cell) and
    `analyses/epistatic_effects_wide.csv` (publication layout, one row per pair, nine genotype
    columns `2_2...0_0` formatted `"mean_dev +/- se"`, ASCII `+/-` rather than `±` since this
    environment has no UTF-8 locale). Also saves one interaction plot per pair to
    `figures/epistasis_<cross>_chr<chr1>-<pos1>_chr<chr2>-<pos2>.pdf`, styled to match the published
    `mdh_qtl/figures/epi_qtl_1.pdf`: a single figure faceted into RIL / B73 BC / Mo17 BC panels
    (each pair's genotype-class mean ± SE is recomputed in all three populations via `pair_effects()`,
    not just the one where it was significant), default grey theme, size-20 text, x-labels rotated
    45°, second locus genotype on x, first locus genotype as color, y = mean deviation from the
    population mean.
13. **`13_circos_plots.R`** — No population argument; builds three whole-genome `circlize` circos
    PDFs from `06`'s peaks, `09`'s thresholds/scans, and `11`'s effect grid, for the five analysis
    slots in published order (RIL, B73 BC, B73 MPH, Mo17 BC, Mo17 MPH). Because the MQM LOD scans
    (`analyses/qtl_analyses/{RIL,B73_BC,Mo17_BC}.RDS`, on their own pseudomarker grids) and the
    effect table (`analyses/qtl_effects_whole_genome.csv`, on a step-2.5 cM grid) don't share rows,
    `build_lod_qtl()` first `approx()`-interpolates each slot's LOD onto the effect grid's `(chr,
    pos)` so LOD and effect line up 1:1; significance uses `get_thr()` (copied from `06`) against
    each RDS's `$permutations`. Writes `figures/nlb_circos_lod.pdf` (five LOD tracks, red where
    significant, dashed peak markers from `analyses/main_effect_peaks.csv`),
    `figures/nlb_circos_effect.pdf` (same five tracks as signed effect size, masked to black/red by
    the same significance grid, dashed zero baseline), and `figures/nlb_circos_paper.pdf` (both
    stacked in one page with A/B labels and a legend).

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

## Effect sign convention (scripts 09, 11)

All additive/dominance estimates are reported as **the substitution effect of the B73 allele on
`NLB_WMD_BLUP`, disease scale as-is** (positive = B73 allele raises disease). Getting this right
requires per-population sign handling because R/qtl's raw `fitqtl` coefficient is on an internal
genotype code that isn't the same for every cross:
- **RIL**: internal code runs B73=-1, Mo17=+1, so `a = -raw_estimate`.
- **B73 BC**: R/qtl's internal "AB" genotype is the true B73/Mo17 heterozygote, so
  `raw_estimate = het - hom` already — no flip. BLUP peaks give the confounded contrast `d - a`;
  MPH peaks give `d` directly.
- **Mo17 BC**: R/qtl's internal "AA" is the true B73/Mo17 heterozygote (the reverse of B73 BC —
  an internal-code trap), so `raw_estimate = hom - het`. BLUP peaks are still reported as-is
  (already `-a - d`, no flip); MPH peaks must be negated to report `d = het - hom`.

See `plans/identify_qtl_and_effects.md` for the full derivation and worked reference values
(e.g. RIL `9@100` → `a ≈ +1.9`, matching effectplot genotype-class means).

## Tooling / dependencies

R, with:
- QTL mapping: `qtl`, `qtl2`, `qtl2convert`; parallel permutations via `parallel` / `snow`.
- Mixed models: `lme4`, `blme`, optional `brms`.
- Data & I/O: `dplyr`, `tidyr`, `stringr`, `lubridate`, `openxlsx`, `readxl`.
- Plotting: `ggplot2`, `cowplot`, `gridGraphics`.

## Scope note

The numbered `01`–`11` pipeline above is the current workflow. The other, unnumbered scripts in
`scripts/` (e.g. `run_rqtl*.R`, `get_*.R`, `convert_rqtl2_to_rqtl.R`, `bayesian_blups.R`,
`genomic_prediction.R`, `gp_helper_functions.R`) are older/superseded or supporting helpers and
are not documented here yet — except `convert_rqtl2_to_rqtl.R`, which the numbered pipeline still
calls from step 04. `scripts/get_effects/` holds the legacy `mdh_qtl` notebooks that `08`–`11`
were distilled from; they reference files/crosses from that other project and won't run here.

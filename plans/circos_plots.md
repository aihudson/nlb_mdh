# Plan: Circos plots for NLB QTL (LOD + effect, plus a combined publication figure)

## Context

The `mdh_qtl` project (the SLB/GLS predecessor of this NLB analysis) produced circular
genome-wide QTL summaries with `circlize`: a **LOD circos** (`mdh_qtl/figures/slb_circos_lod.pdf`)
and an **effect circos** (`mdh_qtl/figures/slb_effect_circos.pdf`), each drawing five stacked tracks
— one per mapping analysis (RIL, B73 BC, B73 MPH, Mo17 BC, Mo17 MPH) — around a cM-based ideogram of
the 10 maize chromosomes. Those scripts live (byte-identical to the `mdh_qtl` originals) at
`scripts/get_effects/circos_plot.R` and `scripts/get_effects/effect_circos_plot.R`. They were never
finished for `nlb_mdh`: they read `mdh_qtl`-era inputs (`data/RIL_rqtl.RDS`, SLB/GLS trait names,
`analyses/qtl_effects_whole_genome.csv` with `pos..cM.` columns,
`analyses/rqtl_combined_qtl_fitqtl.csv`), they contain a lot of interactive scratch code, and they
save to the graphics device manually (no `pdf()` call).

We want to reuse that plotting logic to produce, for the single NLB trait, three PDFs:

1. `figures/nlb_circos_lod.pdf` — five-track LOD circos.
2. `figures/nlb_circos_effect.pdf` — five-track effect circos.
3. `figures/nlb_circos_paper.pdf` — **one publication-quality figure with the LOD circos on top and
   the effect circos on the bottom** (matching the composed `mdh_qtl/figures/*_circos_paper.pdf`
   layout, but assembled in R rather than by hand in Inkscape).

The work is: adapt the two `make_*_circos()` functions to the current numbered-pipeline data
products, wire them to real `pdf()` output, and add the stacked combined figure. All inputs already
exist — no re-running of the QTL pipeline is required.

## Deliverable

One new pipeline script: `scripts/13_circos_plots.R`. No population argument (it reads all three
populations and both traits at once, like scripts 10 and 12). Run from the project root.

## Inputs (current pipeline, replacing the `mdh_qtl` inputs)

The five analysis "slots" keep the published order **RIL, B73 BC, B73 MPH, Mo17 BC, Mo17 MPH**:

| Slot | LOD source | LOD column | Effect column (`qtl_effects_whole_genome.csv`) |
|---|---|---|---|
| RIL | `analyses/qtl_analyses/RIL.RDS` `$scan` (single scanone) | `LOD NLB_WMD_BLUP` | `NLB_ril_effect` |
| B73 BC | `analyses/qtl_analyses/B73_BC.RDS` `$scan[[1]]` | `LOD NLB_WMD_BLUP` | `NLB_b73_effect` |
| B73 MPH | `analyses/qtl_analyses/B73_BC.RDS` `$scan[[2]]` | `LOD NLB_WMD_BLUP_MPH` | `NLB_b73_mph_effect` |
| Mo17 BC | `analyses/qtl_analyses/Mo17_BC.RDS` `$scan[[1]]` | `LOD NLB_WMD_BLUP` | `NLB_mo17_effect` |
| Mo17 MPH | `analyses/qtl_analyses/Mo17_BC.RDS` `$scan[[2]]` | `LOD NLB_WMD_BLUP_MPH` | `NLB_mo17_mph_effect` |

- `$scan` for RIL is a plain `scanone`; for the BCs it is an `mqmmulti` (list of two `scanone`s,
  element 1 = BLUP, element 2 = MPH), exactly as script 06 normalizes it. **Verify these structures
  before coding** (the interrupted check: `readRDS(...)$scan` class/columns and
  `$permutations` null pattern per file).
- **Significance thresholds:** reuse script 06's `get_thr()` (`scripts/06_identify_qtl.R:44-53`) on
  each RDS's `$permutations`, keyed by the `"LOD <trait>"` column name. Do **not** use the old
  circos's `mqmprocesspermutation()` path — the current `$permutations` objects are consumed
  directly by `summary(pp)[[1]]` in `get_thr()`.
- **Effect values:** `analyses/qtl_effects_whole_genome.csv` (columns `chr`, `pos`, and the five
  `NLB_*_effect` columns above). Note the key column is `pos` here, not the `pos..cM.` the old
  script expects.
- **QTL peak markers:** `analyses/main_effect_peaks.csv` (columns `trait`, `chr`, `pos`, `cross`, …)
  replaces `rqtl_combined_qtl_fitqtl.csv`. Take the distinct `(chr, pos)` peaks across all three
  populations as the combined peak set (there is only one underlying trait, NLB). Strip the
  `"LOD "` prefix from `trait` as script 11 does (`scripts/11_genome_wide_effect_scan.R:101`).
- **Ideogram:** build `IcM_bed` (chr, start=0, end=max cM per chr) from the common grid below — same
  shape as the old `IcM_bed`, so the `circos.par`/`circos.initializeWithIdeogram` block is reused
  verbatim (start.degree 70, gap.degree `c(rep(1,9),40)` for the top gap after chr 10).

## Key adaptation: put LOD and effect on one common grid

The old `make_effect_circos()` masks effect values by LOD significance using **positional row
alignment** — `effect_df2[which(is.na(lod_qtl$ril_lod)), ...] <- NA` — which only works if the LOD
table and the effect table share identical rows. In `nlb_mdh` they do not: the LOD scans are on the
MQM step-1 pseudomarker grid (and the three crosses may not even share an identical grid after
`jittermap`), while `qtl_effects_whole_genome.csv` is on a step-2.5 cM grid.

**Resolve both mismatches in one step:** interpolate every population's LOD curve onto the effect
CSV's `(chr, pos)` grid with `approx()` per chromosome, producing a `lod_qtl` table whose rows line
up 1:1 with the effect table. Then:

- the LOD circos plots this interpolated `lod_qtl` (2.5 cM resolution — fine for a genome overview;
  exact peaks are drawn separately from `main_effect_peaks.csv`);
- the effect circos reuses the old significance-masking logic **verbatim**, because rows now align.

This single change is what makes the two reused `make_*_circos()` functions work unmodified in
their core drawing loops.

## Script structure (`scripts/13_circos_plots.R`)

1. `library(circlize); library(qtl); library(dplyr)` (drop `gap`; the old scripts loaded it but
   never used it).
2. **Loaders** (small helpers): read the three RDS files; extract the five `(scan, threshold)` pairs
   in table order using `get_thr()` copied from script 06; read the effect CSV and
   `main_effect_peaks.csv`.
3. **`build_lod_qtl(effect_grid, scans, thresholds)`** — the new grid-alignment helper: for each of
   the five slots, `approx()` the slot's `LOD` onto `effect_grid`'s `(chr, pos)`, then split into
   the `*_lod` (below-threshold) / `*_lod_sig` (above-threshold) column pair exactly as the old
   `mutate(..., ifelse(...))` blocks do. Output columns match what `make_lod_circos()` /
   `make_effect_circos()` index by position (`ril_lod`, `b73_lod`, …, plus the `_sig` variants).
4. **`make_lod_circos(lod_qtl, combined_qtl)`** — adapted from
   `scripts/get_effects/circos_plot.R:213-310`. Keep the five `circos.genomicTrackPlotRegion` blocks
   (black below-threshold line; red line+points for significant; dashed red `type="h"` peak marker),
   the `circos.text` track labels, `track.height=0.1`, and the `ylim=c(0,ymax)` from the max
   significant LOD. Only the data-assembly at the top changes (now takes the pre-built `lod_qtl`).
5. **`make_effect_circos(lod_qtl, effect_df, combined_qtl)`** — adapted from
   `scripts/get_effects/effect_circos_plot.R:187-321`. Keep the `no2` background/overlay column
   duplication, the five symmetric `ylim=c(-ymax,ymax)` tracks, the dashed zero baseline
   (`circos.segments(x0=0,y0=0,x1=CELL_META$xlim,y1=0,lty=2)`), and the peak marker. Update the
   column selection to `pos` (not `pos..cM.`) and the `NLB_*_effect` names; the significance-masking
   block is reused as-is now that rows align.
6. **Render the three PDFs** (the piece the old scripts lacked):
   - `pdf("figures/nlb_circos_lod.pdf", ...)` → `make_lod_circos(...)` → `circos.clear()` → `dev.off()`.
   - `pdf("figures/nlb_circos_effect.pdf", ...)` → `make_effect_circos(...)` → `circos.clear()` → `dev.off()`.
   - **Combined figure** `figures/nlb_circos_paper.pdf`: open a tall portrait `pdf()` (e.g.
     `width=7, height=13`), set `layout(matrix(c(1,2), nrow=2))` (or `par(mfrow=c(2,1))`), call
     `make_lod_circos()` then `circos.clear()`, then `make_effect_circos()` then `circos.clear()`,
     with `par(mar=...)` room for a panel letter (`mtext("A"...)`, `mtext("B"...)`) and a small
     shared legend (red = significant LOD / effect, black = below threshold, dashed red = QTL peak,
     dashed grey = zero effect). circlize draws each circle into the current `mfrow`/`layout` cell,
     so LOD lands on top and effect on the bottom.
7. Add a scope note to `CLAUDE.md`'s pipeline section documenting `13_circos_plots.R` (mirroring the
   style of the 10/12 entries).

## Notes / decisions already made

- **Five tracks, one trait.** NLB has exactly the five analyses the old figure used, so the layout
  transfers directly; there is no GLS/DTA/PH/EH looping to port.
- **Sign convention** is already baked into `qtl_effects_whole_genome.csv` by script 11
  (`sign_mult`), so the effect circos plots those values as-is — positive = B73 allele raises
  disease, consistent with the rest of the pipeline.
- **No new data computation.** This script only reads existing `analyses/` outputs and draws.

## Verification

1. `Rscript scripts/13_circos_plots.R` from `~/projects/nlb_mdh` completes without error and writes
   the three PDFs into `figures/`.
2. Open each PDF: confirm 10 chromosome sectors with the top gap, five labelled tracks (RIL, B73 BC,
   B73 MPH, Mo17 BC, Mo17 MPH), red highlighting where LOD exceeds threshold, and dashed peak
   markers at the `main_effect_peaks.csv` positions (spot-check against known peaks, e.g. RIL chr9
   ~100 cM and chr6 ~245 cM).
3. On the effect circos, confirm the dashed zero baseline and that significant (red) segments
   coincide with the LOD circos's significant regions for the same track — a direct check that the
   grid interpolation aligned LOD and effect correctly.
4. `nlb_circos_paper.pdf` shows the LOD circle on top and the effect circle on the bottom in a
   single page, with A/B panel labels and a legend.

# Plan: Per-chromosome LOD + effect figures (chr3_slb_frreal.pdf style)

## Context

The published `mdh_qtl` project has a per-chromosome QTL figure,
`~/projects/mdh_qtl/figures/chr3_slb_frreal.pdf`: a 3-row × 2-column ggplot facet
grid for a single chromosome, rows = populations (RIL, B73 BC, Mo17 BC), columns =
**Effect** (left) and **LOD** (right). Each panel is a black line vs genetic
position (`Pos (IcM)`), with red dotted vertical lines at QTL peak positions, a
dotted horizontal zero/baseline, and — on the LOD panels — a solid horizontal
significance threshold.

We want the same figure for the **NLB** analysis in this repo. Every input the
original notebook built in-memory is already persisted by the numbered pipeline, and
`scripts/13_circos_plots.R` already loads all of it and aligns LOD to the effect
grid. This is a new, file-driven plotting script that reuses script 13's data
assembly and swaps the circos tracks for Cartesian ggplot facets.

Decisions (confirmed with user):
- **3 rows** only: RIL, B73 BC, Mo17 BC (drop the B73 MPH / Mo17 MPH slots).
- **Only chromosomes that carry a QTL** (one figure per chromosome appearing in
  `main_effect_peaks.csv`).
- **One PDF per chromosome**, named `figures/chr<N>_nlb.pdf` (not "frreal").

## Reference: how the original was drawn

`mdh_qtl/scripts/mqm_single_qtl_plot_frreal.Rmd` (chr3 block ~lines 381–402):
`ggplot(aes(x=Pos, y=value)) + geom_line()`, `facet_wrap(vars(Cross, measurement),
scales="free_y", nrow=3, strip.position="left")`, red dotted `geom_vline`s at peak
cM, dotted `geom_hline` at the baseline (0 for Effect, off-screen for LOD), solid
`geom_hline` at the per-cross LOD threshold, `theme(panel.grid=element_blank(),
axis.line=element_line("black"), text=element_text(size=20), axis.title.y=element_blank())`,
and `ggh4x::facetted_pos_scales()` to fix Effect y ∈ [-0.33,0.33] and LOD y ∈ [0,6].
It depended on an in-memory workspace and hard-coded thresholds/positions — we
replace all of that with the pipeline's saved files.

## Data sources (all already on disk)

| Track | Source | Notes |
|---|---|---|
| Effect line | `analyses/qtl_effects_whole_genome.csv` | cols `chr, pos, NLB_ril_effect, NLB_b73_effect, NLB_mo17_effect` (+ the two `_mph` we ignore); step-2.5 cM grid |
| LOD line | `analyses/qtl_analyses/{RIL,B73_BC,Mo17_BC}.RDS` | `$scan` LOD, interpolated onto the effect grid |
| LOD threshold | same RDS `$permutations` via `get_thr()` | one solid hline per cross/row |
| Peak verticals | `analyses/main_effect_peaks.csv` | red dotted `geom_vline`s; also selects which chromosomes to render |

## Implementation — new file `scripts/14_per_chromosome_qtl_plots.R`

Run from project root: `Rscript scripts/14_per_chromosome_qtl_plots.R` (no
population arg, like scripts 10/12/13).

### 1. Reuse script 13's loader verbatim, restricted to 3 slots
Copy from `scripts/13_circos_plots.R:13-91`:
- `get_thr()` (lines 14–19)
- RDS loads + `scans` list + `thresholds` list (lines 21–43), **keeping only**
  `ril`, `b73`, `mo17` entries. Set `slot_names <- c("ril","b73","mo17")`,
  `slot_labels <- c("RIL","B73 BC","Mo17 BC")`,
  `effect_cols <- c("NLB_ril_effect","NLB_b73_effect","NLB_mo17_effect")`.
- `effect_grid` + `main_effect_peaks` reads (lines 45–52).
- `build_lod_qtl()` (lines 64–91) — it loops over `slot_names`, so with 3 slots it
  returns `chr, start(=pos), end`, plus `<slot>_lod` / `<slot>_lod_sig` columns.
  Recover the **raw** interpolated LOD per slot with
  `dplyr::coalesce(<slot>_lod, <slot>_lod_sig)` (below-threshold values live in
  `_lod`, at/above-threshold in `_lod_sig`; exactly one is non-NA per row).

### 2. Build one tidy long data frame for plotting
For each of the 3 slots produce two measurement layers keyed by `(chr, pos, Cross,
measurement, value)`:
- **Effect**: `value = effect_grid[[effect_cols[slot]]]`, `measurement="Effect"`.
- **LOD**: `value = raw interpolated LOD`, `measurement="LOD"`.

`Cross` = factor with levels `c("RIL","B73 BC","Mo17 BC")` (row order).
Strip the `chr` prefix so `chr` is numeric (or keep the effect_grid's numeric `chr`).
Attach the x label column as `pos` in cM.

Also build small helper frames:
- `thr_df`: one row per `Cross` with the LOD `yintercept` from `thresholds`
  (used only in LOD panels — join on `measurement=="LOD"`).
- `peaks`: `main_effect_peaks %>% distinct(chr, pos)` (strip `"LOD "` from `trait`
  as script 13 line 49 does) for the red verticals.

### 3. Render per chromosome (loop over chromosomes with a peak)
`for (ch in sort(unique(peaks$chr)))`:
Filter the long frame, `thr_df`, and peak positions to `ch`.

Match the reference layout. Two viable assemblies — **use the two-plot + cowplot
route** to avoid adding the `ggh4x` dependency (cowplot is already a project dep,
see `scripts/03_plot_blup_distributions.R`):

- `p_eff`: `ggplot(effect rows) + geom_line() +
  facet_grid(Cross ~ ., switch="y") + geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(data=peaks_ch, aes(xintercept=pos), color="red", linetype="dotted",
  linewidth=1) + ylim(symmetric effect limit)`.
- `p_lod`: same but LOD rows, `geom_hline(data=thr_df_ch, aes(yintercept=yintercept),
  linewidth=1)` (solid per-row threshold) and the same red verticals; `ylim(0, lodmax)`.
- Shared theme: `theme_bw()`-free — `theme(panel.grid=element_blank(),
  axis.line=element_line("black"), text=element_text(size=20))`, strips on the left
  (`strip.placement="outside"`, blank `axis.title.y`) so the row label (population)
  and column meaning read like the reference. Give each plot an x title `"Pos (IcM)"`.
- `cowplot::plot_grid(p_eff, p_lod, ncol=2, labels=NULL)` → matches the Effect-left /
  LOD-right, 3-row grid.

`facet_grid(Cross ~ .)` keeps a **single shared y-scale down each column**
(all 3 Effect panels share limits, all 3 LOD panels share limits) — this reproduces
the `facetted_pos_scales` behavior of the original without `ggh4x`.

*(Alternative if an exact single-`facet_wrap` match is preferred: one ggplot with
`facet_wrap(vars(Cross, measurement), nrow=3, scales="free_y", strip.position="left")`
plus `ggh4x::facetted_pos_scales(y=list(measurement=="Effect" ~ scale_y_continuous(...),
measurement=="LOD" ~ scale_y_continuous(...)))`. Only take this route if adding the
`ggh4x` dep is acceptable.)*

### 4. Auto y-limits (per figure, over the 3 crosses on that chromosome)
- Effect: `lim <- max(abs(effect values on ch)) * 1.1`; use `c(-lim, lim)`.
- LOD: `c(0, ceiling(max(max LOD on ch, max threshold) * 1.05))`.
This generalizes the reference's hand-tuned per-figure limits.

### 5. Save
`ggsave(sprintf("figures/chr%s_nlb.pdf", ch), combined, width=8, height=10)`
(or `pdf()`/`print()`/`dev.off()`), matching the reference's 8×10 in.

### 6. (Optional) Document as pipeline step 14
Add a `14_per_chromosome_qtl_plots.R` bullet to the pipeline list in `CLAUDE.md` so
the numbered workflow stays complete (no population arg; consumes 06/11 outputs + the
RDS scans; writes `figures/chr<N>_nlb.pdf` per QTL chromosome).

## Verification

1. From `~/projects/nlb_mdh`, run `Rscript scripts/14_per_chromosome_qtl_plots.R`;
   confirm it exits 0 with no missing-file/threshold errors.
2. Confirm one `figures/chr<N>_nlb.pdf` exists for each distinct chromosome in
   `analyses/main_effect_peaks.csv` (and none for QTL-free chromosomes).
3. Open one output and check against `~/projects/mdh_qtl/figures/chr3_slb_frreal.pdf`:
   3 rows (RIL/B73 BC/Mo17 BC) × 2 columns (Effect left, LOD right), black lines, red
   dotted verticals at the peak cM positions, dotted zero baseline on Effect panels,
   solid threshold line on LOD panels, size-20 text, no grid lines.
4. Sanity-check numbers: the LOD peaks and threshold crossings should match the
   red-highlighted significant regions in `figures/nlb_circos_lod.pdf`, and peak
   x-positions should match `main_effect_peaks.csv` rows for that chromosome.

## Files

- **New**: `scripts/14_per_chromosome_qtl_plots.R`
- **Read only**: `analyses/qtl_effects_whole_genome.csv`,
  `analyses/main_effect_peaks.csv`,
  `analyses/qtl_analyses/{RIL,B73_BC,Mo17_BC}.RDS`
- **Reuse as template**: `scripts/13_circos_plots.R` (data loader, `get_thr`,
  `build_lod_qtl`), `scripts/03_plot_blup_distributions.R` (ggplot2/cowplot style)
- **Output**: `figures/chr<N>_nlb.pdf` (one per QTL chromosome)
- **Optional edit**: `CLAUDE.md` (add step-14 entry)

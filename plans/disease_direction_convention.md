# Disease-direction convention (resistance-forward) — audit & remediation

> **Permanent record.** This document is the single source of truth for which direction ("higher =
> more resistant" vs. "higher = more disease") every output of the numbered pipeline (`scripts/01`–`16`)
> is on, why, and what still needs to be brought into line. It is self-contained: a fresh agent can
> execute the remediation from this file alone. Last audited **2026-08-12** against the
> `analyses/` CSVs dated Aug 4–11 2026.

---

## 1. The convention (decided)

- **Raw field scores** (`wmd` = weighted mean disease, and the per-date `NLB_*` columns) are on the
  **disease / susceptibility scale**: 0 = healthy, 100 = dead, **higher = worse**. This is how NLB
  has traditionally been reported, and it is the scale of the supplemental workbook File S1.
- **Everything downstream of `02_get_blups.R`** is on the **resistance scale**, because `02` flips
  the phenotype with `wmd = 100 - wmd` **before** fitting any BLUP:
  - `scripts/02_get_blups.R:46` — `mutate(wmd = 100-wmd)` (IBM lines)
  - `scripts/02_get_blups.R:152` — `mutate(wmd = 100 - wmd)` (NAM lines)
  So all BLUPs, heterosis, cross-file phenotypes, QTL effects, gene-action, and their figures are
  **higher = more resistant**.
- **Why the flip exists:** heterosis *percentages* `(hybrid − parent)/parent` divide by a
  near-zero disease value on the disease scale, producing wild ("wonky") percentages. On the
  resistance scale the denominator is ~70–100, so the percentages behave. The flip was originally
  motivated by heterosis but it propagates through the *entire* pipeline.
- **Manuscript-wide rule:** report **resistance-forward** (higher = more resistant) everywhere
  except the raw File S1 scores, which stay disease-severity and are labeled as such.
- **Direction-agnostic quantities** (no disease/resistance meaning): all LOD scores, percent
  variance explained, correlations, and heritabilities. Positions/CIs are cM, not effects.

### Empirical proof of direction (so no one has to re-derive it)

From the current `analyses/IBM_NLB_BLUPs.csv` (`NLB_WMD_BLUP` column):

| Line / population | `NLB_WMD_BLUP` | Meaning |
|---|---|---|
| **B73** (founder) | **80.96** | susceptible parent → **low** |
| **Mo17** (founder) | **98.12** | resistant parent → **high** |
| B73×Mo17 (F1) | 92.98 | between parents, above mid-parent (89.54) |
| RIL population | mean 88.26 (range 70.6–101.3) | — |
| B73 BC population | mean 88.98 | — |
| **Mo17 BC population** | **mean 93.60** (highest) | backcross to resistant parent → most resistant |

The susceptible parent (B73) sits at the **low** end and the resistant parent (Mo17) at the
**high** end. Values live near 70–100 because they are `100 − wmd` (raw disease ~0–30). This is the
anchor: **higher `NLB_WMD_BLUP` = more resistant**, on every figure/table derived from this column.

QTL-effect sign anchors (from `plans/identify_qtl_and_effects.md`, unchanged numbers):
`RIL 9@100 a ≈ +1.9` because B73/B73 = 89.79 > Mo17/Mo17 = 85.94 (B73 allele **raises resistance**
here); `RIL 6@245 a ≈ −1.7` because B73/B73 = 86.27 < Mo17/Mo17 = 89.71 (B73 allele **lowers
resistance = raises disease** here — B73 is *not* protective at chr6).

---

## 2. Per-output direction table (audit result)

| Script → output | Direction | Note |
|---|---|---|
| 01 `data/nlb_mdh_file_s1.xlsx` (raw `wmd`, per-date scores) | **Disease** (higher = worse) | Raw data; no flip. Label as "disease severity". |
| 02 `IBM_NLB_BLUPs*.csv`, `NAM_NLB_BLUPs.csv`, `line_blups_envs.csv` | **Resistance** | `100 - wmd` at 02:46 (IBM), 02:152 (NAM). |
| 02 heterosis `*_MPH/BPH(_PCT)` | **Resistance** | `add_heterosis()` 02:70–105 on the flipped BLUP; positive = hybrid more resistant; BPH `max` = most-resistant parent. |
| 02 `data/old_nlb_data.csv` | **Mixed** | `NLB_WMD_BLUP_new` = resistance; old AUDPC cols = disease. |
| 03 `blup_distributions.pdf`, `mph_bph_distributions.pdf` | **Resistance** | Plots the CSV as-is (no re-flip). Confirmed by founder BLUPs (§1). Axis currently labeled only "NLB". |
| 04 cross files (`RIL/B73/Mo17/env_blups_cross.csv`) | **Resistance** | `convert_rqtl2_to_rqtl.R` passes phenotypes through untouched; `old_nlb_data_cross.csv` = mixed. |
| 05/06/07/08 scans & peaks (`main_effect_peaks.csv`, `epistatic_peaks.csv`) | **Direction-agnostic** | LOD / positions / p-values only; no signed effects written. |
| 09 `qtl_effects.csv` `estimate` | **Resistance** | Fit on the resistance BLUP; positive = B73 allele raises resistance. **Code comment says "disease" — WRONG (fix in §3).** |
| 10 `qtl_gene_action.csv` (`a`, `*_d`, `d/a`, `*_action`) | **Resistance** | Inherits 09 verbatim. The `ud`/`od` split is sign-dependent, so scale matters; additive/pd/dominant use `abs()` only. |
| 11 `qtl_effects_whole_genome.csv` `NLB_*_effect` | **Resistance** | Same `sign_mult` convention as 09. **Comment says "disease" — WRONG (fix in §3).** |
| 12 `epistatic_effects_long/wide.csv`, `epistasis_*.pdf` | **Resistance** | Genotype-class deviation from the resistance-scale trait mean (12:61,75); y-axis titled "Mean". |
| 13 `nlb_circos_effect.pdf`, panel B of `nlb_circos_paper.pdf` | **Resistance** | Effect track from `qtl_effects_whole_genome.csv`, no re-flip; neutral legend. LOD track direction-agnostic. |
| 14 `chr*_nlb.pdf` Effect panel | **Resistance** | Same effect columns; y-axis blanked. LOD panel direction-agnostic. |
| 15 `qtl_overlap_old_study.csv` | **Direction-agnostic** | Overlap uses positions only. The transcribed old-study `a` column is old **disease**-scale (AUDPC), carried for reference, unused in logic. |
| 16 S1/S2/S3 (correlations), S5 (H²) | **Direction-agnostic** | Correlation sign / heritability are flip-invariant. |
| 16 S6/S7/S8 (heterosis), S9 (QTL effects), S10 (gene action) | **Resistance** | Captions **already correct** ("on the resistance scale" / "positive = B73 allele increases resistance"). |
| 16 S11 (old-vs-new correlation) | **Direction-agnostic**; negative **by design** | Old = disease (AUDPC), new = resistance (100−WMD); caption already explains the negative sign. |

---

## 3. Defects to fix (documentation only — numbers are already correct)

The `sign_mult()` **mechanics are correct** (RIL flip; B73-BC/Mo17-BC internal-genotype-code trap).
Only the *scale wording* is wrong: it calls the resistance-scale BLUP "disease scale" and says
"positive = raises disease." Swap "disease" → "resistance" and "raises disease" → "raises
resistance" (noting: lower resistance = more disease). Keep all per-population mechanics verbatim.

1. **`CLAUDE.md`**
   - Pipeline step 9 description (~line 63): "the B73-allele effect on the disease-scale BLUP,
     positive = B73 raises disease" → resistance-scale, positive = raises resistance.
   - Section **"## Effect sign convention (scripts 09, 11)"** (starts ~line 221): retitle the scale
     as resistance and note `NLB_WMD_BLUP = 100 − wmd`. Keep the RIL / B73-BC / Mo17-BC bullets.
   - Add one line near the Overview/Conventions pointing to this file as the project-wide direction
     convention.
2. **`scripts/09_estimate_qtl_effects.R`** — header comment (~lines 32–34): "disease-BLUP
   convention" → "resistance-BLUP convention"; per-population bullets unchanged.
3. **`scripts/11_genome_wide_effect_scan.R`** — header comment (~lines 26–27): same swap.
4. **`plans/identify_qtl_and_effects.md`** — the substantive fix:
   - ~line 151 "disease scale as-is (positive = B73 allele raises the BLUP)" → resistance scale.
   - **Correct the inverted chr6 interpretation (~lines 44–48):** B73/B73 = 86.27 < Mo17/Mo17 =
     89.71 are *resistance* means, so at chr6 the B73 allele **lowers resistance = raises disease**;
     B73 is **not** protective there. Rewrite the "B73 lowers disease / protective at chr6" prose.
     Cross-check the chr9 example (B73/B73 89.79 > Mo17/Mo17 85.94 → B73 raises resistance).
   - Add a pointer to this file.

**Not done (deliberate, per chosen scope):** the misleading column name `NLB_WMD_BLUP` (reads
"disease" but holds resistance) is **not** renamed; instead this doc + CLAUDE.md state it is a
`100 − wmd` resistance quantity.

---

## 4. Figure labels + re-render (figures are stale — this is not cosmetic)

**Staleness:** current `analyses/` CSVs are Aug 4–11 2026, but saved figures predate them and still
show the *pre-flip* direction:
- `figures/blup_distributions.pdf` — **Apr 5 2026** (stale; shows Mo17 BC at the *low* end, the old
  direction). Current data has Mo17 BC **highest** (93.60).
- `figures/mph_bph_distributions.pdf` — **Apr 28 2026** (stale).
- `figures/chr*_nlb.pdf`, `figures/epistasis_*.pdf`, `figures/nlb_circos_*.pdf` — **Jul 29 2026**
  (trail the Aug 4/11 CSVs slightly).
- `figures/blup_distributions_paper.pdf`, `figures/mph_bph_distributions_paper.pdf` — **Feb 2025**,
  produced by **broken/dead template code** in `03` (references an undefined `blues` object and
  non-NLB traits). Leave or delete; do **not** try to fix.

So applying labels requires re-rendering from the current CSVs, which also corrects the stale
direction. This reads current `analyses/` CSVs — it does **not** re-run the modeling in
02/05/07/09/etc.

Label edits (then re-run the script):
- **`scripts/03_plot_blup_distributions.R`** — `labs(x=...)` at lines 65/112/152. BLUP panel:
  `"NLB resistance BLUP (100 - WMD)"`; MPH/BPH panels: `"NLB resistance MPH (%)"` / `"... BPH (%)"`.
  Optional "→ more resistant" cue.
- **`scripts/12_epistatic_effect_table.R`** — y-axis (~line 138) `"Mean"` →
  `"Deviation from mean resistance (100 - WMD)"`.
- **`scripts/13_circos_plots.R`** — annotate the effect track/legend (~lines 189–191):
  "B73-allele effect on resistance (+ = more resistant)". LOD panel stays magnitude.
- **`scripts/14_per_chromosome_qtl_plots.R`** — Effect-panel y-axis is blanked
  (`axis.title.y = element_blank()`, ~line 111) / title "Effect" (~line 130). Add a shared y-title:
  "B73-allele effect on resistance (+ = more resistant)".
- **`scripts/16_supplement_tables.R`** — captions **already correct**; verify only. Optionally add
  a File S1 disease-scale note to the front-matter. The docx is current (Aug 11). ⚠ A Word lock
  file `manuscript/~$nlb_mdh_supplement.docx` was present — close Word before re-running 16.

Re-render order: `Rscript scripts/03_...`, then `12`, `13`, `14`. Overwrites
`blup_distributions.pdf`, `mph_bph_distributions.pdf`, `epistasis_*.pdf`, `nlb_circos_*.pdf`,
`chr*_nlb.pdf`.

---

## 5. Known pre-existing issues (out of scope here, recorded so they aren't lost)

- `03` `*_paper.pdf` blocks: dead/broken template code (undefined `blues`, non-NLB traits) — never
  produces valid NLB output.
- `04_make_cross_file.sh` references `analyses/IBM_NLB_BLUPS_*.csv` (uppercase `BLUPS`) while `02`
  writes `IBM_NLB_BLUPs_*.csv` (lowercase `BLUPs`) — a case mismatch that would break on a
  case-sensitive filesystem.
- `01` height/DTA reads at the tail have an unterminated string; the file won't source cleanly past
  that point (unrelated to direction).

---

## 6. Verification checklist

1. **Grep** — `grep -rn "raises disease\|disease-scale\|disease scale" CLAUDE.md scripts/ plans/`
   returns only legitimate *raw*-disease references (File S1 / old-study AUDPC), not the BLUP/effect
   convention.
2. **Cross-doc agreement** — CLAUDE.md's effect-sign section, `09`/`11` comments, `plans/
   identify_qtl_and_effects.md`, and `16` captions all say resistance-forward (positive = raises
   resistance); chr6 now reads "B73 raises disease / lowers resistance."
3. **Figure re-render** — after running 03/12/13/14, `figures/blup_distributions.pdf` shows **Mo17
   BC at the higher (more-resistant) end** (mean 93.60), and every axis/legend states the direction.
4. **Reference values unchanged** — `qtl_effects.csv` still has `RIL 9@100 a ≈ +1.9`, `RIL 6@245
   a ≈ −1.7`; only their interpretation ("+ = more resistant") is corrected.
5. **Founder anchor** — `IBM_NLB_BLUPs.csv`: B73 = 80.96 < Mo17 = 98.12 (susceptible below
   resistant), the standing proof that higher = more resistant.

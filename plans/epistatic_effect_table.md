# Plan: Epistatic-effect genotype-class table for significant QTL pairs

> Self-contained implementation spec. Repo root: `~/projects/nlb_mdh` (run everything from there).
> A neighboring project `~/projects/mdh_qtl` holds the *original* (older) analysis this is distilled
> from; its files are referenced for provenance but are **not** run here.

## Context / goal

Produce a table of **epistatic QTL effects** for the significant interacting QTL pairs in this
NLB (Northern Leaf Blight, IBM B73×Mo17) mapping project. For each significant pair, report the
**phenotypic mean of every genotype-class combination** at the two loci — expressed relative to the
population mean, ± its standard error — plus the pair's full-model LOD, interaction LOD, and which
analysis it was significant in. This mirrors a previously published table whose layout is:

```
Trait  Chr1 Chr2  Pos1  Pos2   2_2  2_1  1_2  1_1  1_0  0_1  0_0   LOD_full  LOD_int  Sig
```
where each cell `x_y` = (# B73 alleles at QTL1)_(# B73 alleles at QTL2), holding `mean ± SE`
centered on the population mean; `Sig` names the population/analysis where the pair was significant.

### How the old table was actually produced (verified)
It was **hand-assembled** in `~/projects/mdh_qtl/scripts/rqtl_effectplot.Rmd` (see lines 104–117,
252–261, 339–348) from **raw** R/qtl `effectplot()` genotype-class means: `$Means` centered by
subtracting the population mean, `$SEs` used verbatim, with **no** adjustment for other main-effect
QTL. The `scripts/get_effects/` folder in *this* repo is just a copy of those mdh_qtl scripts. There
is **no** single reusable table-generator — the reusable core is the `effectplot` means/SE pattern,
which the new script automates and generalizes. `effectplot` genotype-class means are exactly the
"phenotypic mean of each genotype" the caption describes.

### Decisions locked with the user
- **Pairs = significant only:** `sig_level == "0.05" & !same_chr_close` from
  `analyses/epistatic_peaks.csv`. This yields **exactly 6 pairs** (see list below). (Marginal 0.20
  pairs, all MPH-trait pairs, and all Mo17 BC significant pairs are excluded — the latter two are
  same-chromosome-close linked-marker artifacts, not real epistasis.)
- **Cell means = raw `effectplot` means**, centered on the population mean, ± class SE. No fitqtl /
  main-effect-QTL adjustment (faithful to the old published table).
- **Two outputs:** a tidy long CSV (source of truth) and a wide publication-style CSV.

The pair list plus `lod.int` and `lod.full` are already columns in `epistatic_peaks.csv`, so **no
scantwo recomputation is needed** — carry those through.

### The exact 6 pairs (from `analyses/epistatic_peaks.csv`, sig_level=="0.05", same_chr_close==FALSE)
| cross | trait | chr1 | pos1 | chr2 | pos2 | lod.int | lod.full |
|-------|-------|------|------|------|------|---------|----------|
| RIL | NLB_WMD_BLUP | 1 | 891.600213 | 3 | 412.600112 | 6.693 | 6.871 |
| RIL | NLB_WMD_BLUP | 3 | 389.600107 | 5 | 284.500071 | 5.715 | 5.779 |
| RIL | NLB_WMD_BLUP | 7 | 253.500051 | 8 | 153.400022 | 4.430 | 5.000 |
| B73_BC | NLB_WMD_BLUP | 3 | 37.900013 | 3 | 249.700070 | 5.100 | 5.503 |
| B73_BC | NLB_WMD_BLUP | 3 | 138.100029 | 8 | 567.500106 | 4.703 | 5.132 |
| B73_BC | NLB_WMD_BLUP | 4 | 447.200086 | 9 | 219.800047 | 4.516 | 5.646 |

(The script must not hardcode these — it derives them by filtering the CSV — but they are the
expected result and the verification target.)

## Deliverable: `scripts/12_epistatic_effect_table.R`

New numbered-pipeline script. **No population CLI argument** — it reads all significant pairs from
`epistatic_peaks.csv` and loops over the populations present, loading each cross once. Follow the
existing numbered-pipeline conventions: run from repo root, root-relative paths, and **duplicate**
the `presets`/`read_cross` config inline (the numbered scripts do not `source()` each other).

### Config to copy verbatim from `scripts/09_estimate_qtl_effects.R`

`presets` (09 lines 11–18 — the `ga_trait` field is unused here, keep or drop):
```r
presets <- list(
  RIL = list(input_file = "analyses/RIL_cross.csv", genotype = c("A", "B"),
             na.strings = "-", crosstype = "ril"),
  B73_BC = list(input_file = "analyses/B73_cross.csv", genotype = c("AA", "AB"),
                na.strings = "A-", crosstype = "bc"),
  Mo17_BC = list(input_file = "analyses/Mo17_cross.csv", genotype = c("AB", "BB"),
                 na.strings = "-B", crosstype = "bc")
)
```
Note the `epistatic_peaks.csv` `cross` values are exactly `RIL`, `B73_BC`, `Mo17_BC` — same keys.

`read_cross()` (09 lines 53–64, copy verbatim):
```r
read_cross <- function(input_file, genotype, na.strings, crosstype) {
  ril <- crosstype == "ril"
  if (ril) crosstype <- "bc"
  cross <- read.cross(format = "csv", file = input_file, genotype = genotype,
                       alleles = c("A", "B"), na.strings = na.strings, crosstype = crosstype)
  cross <- jittermap(cross)
  if (ril) cross <- convert2riself(cross)
  cross$pheno <- cross$pheno %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), as.numeric))
  cross
}
```
Load each cross with imputation so `effectplot` can resolve arbitrary `chr@pos` (09 line 187):
```r
cross <- sim.geno(read_cross(preset$input_file, preset$genotype, preset$na.strings,
                             preset$crosstype), step = 2.5)
```
Requires `library(qtl)` and `library(dplyr)`.

### Genotype-class → B73-allele-count mapping (the R/qtl internal-code trap)
`effectplot` orders the two genotype classes by R/qtl internal code (1 then 2), and
`read.cross(genotype = c(x, y))` maps CSV allele-string `x`→code 1, `y`→code 2. Per the sign/label
logic documented in `09` lines 26–40 and `CLAUDE.md` "Effect sign convention", the `geno1`/`geno2`
label vectors passed to `effectplot` (in code order) and their B73-allele counts are:

| population | genotype=c(...) | geno labels (code order) | B73 counts (code order) |
|-----------|------------------|--------------------------|--------------------------|
| RIL | c("A","B") | c("B73/B73","Mo17/Mo17") | c(2, 0) |
| B73_BC | c("AA","AB") | c("B73/B73","B73/Mo17") | c(2, 1) |
| Mo17_BC | c("AB","BB") | c("B73/Mo17","Mo17/Mo17") | c(1, 0) |

The Mo17_BC row is the "trap": internal code 1 is the **het** (count 1), not B73/B73. (No Mo17_BC
pairs survive the significance filter, but keep the mapping for correctness/generality.) These match
the `geno1`/`geno2` args used in `~/projects/mdh_qtl/scripts/rqtl_effectplot.Rmd` lines 80–100.

### Core primitive (the reused `effectplot` pattern)
For a pair on a loaded `cross`, trait `trait`, with label vectors `g1lab`/`g2lab` and count vectors
`c1`/`c2` for that population:
```r
ep <- effectplot(cross, pheno.col = match(trait, colnames(cross$pheno)),
                 mname1 = sprintf("%s@%s", chr1, pos1),
                 mname2 = sprintf("%s@%s", chr2, pos2),
                 geno1 = g1lab, geno2 = g2lab, draw = FALSE)
# ep$Means and ep$SEs are 2x2: rows = locus1 classes (code order), cols = locus2 classes.
gm <- mean(cross$pheno[[trait]], na.rm = TRUE)   # population-mean centering
# cell (i,j): geno1 = c1[i], geno2 = c2[j], mean_dev = ep$Means[i,j] - gm, se = ep$SEs[i,j]
```
`effectplot` averages over the `sim.geno` imputations, so `draw = FALSE` just suppresses plotting.

### Algorithm
1. `library(qtl); library(dplyr)`. Read `analyses/epistatic_peaks.csv`
   (`read.csv(..., stringsAsFactors = FALSE)`). Filter `sig_level == "0.05" & !same_chr_close`.
2. For each population in `unique(filtered$cross)`:
   - Load `cross` via `sim.geno(read_cross(...), step = 2.5)` using that population's preset.
   - For each `trait` in that population's rows, `gm <- mean(cross$pheno[[trait]], na.rm = TRUE)`.
   - For each pair row: call `effectplot` as above; emit one tidy row per genotype cell (`i,j`).
3. Assemble the long table; derive the wide table by pivoting on `paste(geno1, geno2, sep = "_")`.

### Outputs
- **`analyses/epistatic_effects_long.csv`** (source of truth), columns:
  `cross, trait, chr1, pos1, chr2, pos2, geno1, geno2, mean_dev, se, lod_full, lod_int, sig_level`
  (`geno1`/`geno2` = integer B73-allele counts at locus1/locus2; `lod_full`/`lod_int` copied from
  the pair's `lod.full`/`lod.int`).
- **`analyses/epistatic_effects_wide.csv`** (publication layout), one row per pair:
  `cross, trait, chr1, chr2, pos1, pos2`, then the **nine** genotype cells in this fixed order
  `2_2, 2_1, 2_0, 1_2, 1_1, 1_0, 0_2, 0_1, 0_0` (full union; only the 4 cells relevant to each
  population are filled, the rest `NA`/`"N/A"`), each cell formatted as `"mean_dev ± se"` rounded to
  2 dp, then `lod_full, lod_int, sig_level`. Use `enc2utf8`/`"±"` for the ± glyph, or `"+/-"`
  if ASCII is preferred.

For the significant set, expect these filled cells: **RIL** pairs fill `2_2, 2_0, 0_2, 0_0`;
**B73_BC** pairs fill `2_2, 2_1, 1_2, 1_1`.

## Files
- **Create:** `scripts/12_epistatic_effect_table.R` (the only new code).
- **Update (optional but recommended):** `CLAUDE.md` — add a step-12 entry to the pipeline list and
  note the two new `analyses/` outputs, mirroring the style of the existing steps 8–11.
- **Also copy this plan** into the repo at `plans/epistatic_effect_table.md` (the repo already keeps
  design docs there, e.g. `plans/identify_qtl_and_effects.md`) once implementation begins.

## Verification (end-to-end)
1. `Rscript scripts/12_epistatic_effect_table.R` — runs clean, writes both CSVs.
2. Row counts: exactly **6 pairs** (RIL 3, B73_BC 3, no Mo17_BC, no MPH trait). Long =
   6 × 4 cells = **24 rows**; wide = **6 rows**. Confirm no `same_chr_close` pair and no
   `NLB_WMD_BLUP_MPH` row leaked in.
3. Spot-check RIL `1@891.6 × 3@412.6` interactively: the four cell values (`2_2/2_0/0_2/0_0`) must
   equal `effectplot(ril_cross, pheno.col, "1@891.600213", "3@412.600112", geno1=c("B73/B73",
   "Mo17/Mo17"), geno2=c("B73/B73","Mo17/Mo17"), draw=FALSE)$Means - mean(RIL NLB_WMD_BLUP)`, and the
   SEs equal `$SEs`. Its `lod_full`/`lod_int` must equal `6.871`/`6.693` from `epistatic_peaks.csv`.
4. Confirm B73_BC pairs populate only `2_*`/`1_*` cells (never `0_*`), verifying the code→B73-count
   map and the internal-code handling.

## Key file references
- `scripts/09_estimate_qtl_effects.R` — `presets` (11–18), `read_cross` (53–64), `sim.geno` load
  (187), sign/internal-code notes (26–40).
- `scripts/08_identify_epistatic_qtl.R` — produced `epistatic_peaks.csv`; defines `sig_level`
  buckets and `same_chr_close`.
- `analyses/epistatic_peaks.csv` — columns `cross, trait, chr1, pos1, chr2, pos2, lod.int, lod.full,
  int_p, sig_level, same_chr_close`; the pair source.
- `~/projects/mdh_qtl/scripts/rqtl_effectplot.Rmd` — original `effectplot` means/SE + centering
  pattern (lines 80–117, 339–348). `scripts/get_effects/rqtl_effectplot.Rmd` is a copy.
- `CLAUDE.md` — pipeline overview, effect-sign convention, per-population cross codes.

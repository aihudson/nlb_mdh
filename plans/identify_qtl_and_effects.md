# Plan: Identify significant QTL and estimate their effects (scripts 08–11)

> **Self-contained plan.** Everything a fresh session needs is in this file: the design (scripts 08–11), the decided effect convention, the concrete reference values discovered from the data, and copy-paste R verification snippets. Run all commands from the project root `~/projects/nlb_mdh`. A copy of this plan also lives at `~/.claude/plans/i-want-to-identify-spicy-crescent.md`; consider copying it into the repo (e.g. `plans/`) if you want it version-controlled.

## Context

The NLB QTL pipeline currently stops after peak/scan generation:
- `06_identify_qtl.R` already calls significant **main-effect** QTL from the MQM scans → `analyses/main_effect_peaks.csv`.
- `07_epistatic_qtl.R` produces `scantwo` scans + permutations (`analyses/qtl_analyses/<pop>_scantwo.RDS`) but **nothing consumes them** — significant epistatic pairs are never called.
- No script estimates **effect sizes** (additive `a`, dominance `d`) at the identified QTL.

The goal is to finish the pipeline: (1) identify significant epistatic QTL, and (2) estimate QTL effects — additive from the BLUP phenotype, dominance from the MPH phenotype — then classify gene action (`d/a`), and produce a genome-wide effect profile. This reconstructs (and cleans up) the analysis done in the older `mdh_qtl` project, whose code survives as the undocumented, semi-manual notebooks in `scripts/get_effects/`.

### What the data shows (verified during planning)
- Phenotype columns are identical across `RIL_cross.csv` / `B73_cross.csv` / `Mo17_cross.csv`: col **1 = `NLB_WMD_BLUP`** (additive), col **5 = `NLB_WMD_BLUP_MPH`** (dominance). Only one disease trait — far simpler than the legacy 5-trait analysis.
- Main-effect peaks already called (6 total): RIL 6@245 & 9@100; B73_BC 6@245 (BLUP) & 1@480 (MPH); Mo17_BC 5@535 (BLUP) & 1@760 (MPH).
- **Epistasis is real but needs artifact-filtering.** Against the interaction-LOD permutation threshold: RIL BLUP **c1×c3** (int LOD 6.69, p≈0.000) and **c3×c5** (p≈0.008); B73_BC BLUP **c3×c8** (p≈0.042) and **c4×c9** (p≈0.05) are genuine inter-chromosomal interactions. The MPH hits and all same-chromosome hits are pairs 1–2 cM apart (linked-marker artifacts). So `08` must flag/exclude same-chromosome, closely-spaced pairs and also surface near-significant (α=0.20) pairs so you can decide whether to run more permutations.

## Approach

Add four CLI-style numbered scripts that follow the **`07_epistatic_qtl.R` skeleton** (usage comments → `commandArgs()` with `ifelse` defaults → `presets` list keyed by `RIL`/`B73_BC`/`Mo17_BC` → shared cross-reading/RIL-conversion/numeric-pheno block → `paste0` paths → outputs under `analyses/`). Reuse the effect-estimation mechanics distilled from the legacy code rather than reinventing them.

**Effect-estimation mechanics to reuse** (all from `qtl`): `sim.geno(cross, step=2.5)` → `makeqtl(cross, chr, pos)` → `fitqtl(cross, pheno.col, qtl, get.ests=TRUE, formula=...)` → effect estimate = `summary(out.qtl)$ests[-1,1]`; per-QTL drop-one LOD = `summary(out.qtl)$result.drop[,"LOD"]`. The cleanest reference is `run_fitqtl()` in [rqtl_mqm_effect.Rmd](scripts/get_effects/rqtl_mqm_effect.Rmd#L47); the epistasis-aware, formula-building version is `run_fitqtl()`/`mqm_effect_scan()` in [qtl_effect_all.Rmd](scripts/get_effects/qtl_effect_all.Rmd#L60) and [get_qtl_effects_frreal.R](scripts/get_effects/get_qtl_effects_frreal.R#L1); the partial NLB port to build on is [get_effects_mqm.R](scripts/get_effects_mqm.R). Note the legacy `negative=TRUE` sign flip is **not** used here (see convention below).

### Effect sign & allele convention (decided)
All effects are the **substitution effect of the B73 allele on `NLB_WMD_BLUP`, disease scale as-is** (positive = B73 allele raises the BLUP; no sign flip). Parameterization: B73/B73 = m+a, B73/Mo17 = m+d, Mo17/Mo17 = m−a.
- **Additive `a` — RIL only.** RIL gives both homozygous classes, so `a = (μ[B73/B73] − μ[Mo17/Mo17]) / 2` (verified at 9@100: (89.79−85.94)/2 = +1.93). For a QTL significant **only in a BC**, obtain its `a` by fitting at that map position **in the RIL cross** (the legacy `qtl_remain_effect`/`get_qtl_effects` approach).
- **Dominance `d` — BC MPH only.** The MPH phenotype directly estimates `d` (heterozygote deviation from mid-parent).
- **BC BLUP is confounded, not additive.** The two-class BC BLUP contrast equals `d − a` (B73 BC) or `−a − d` (Mo17 BC) (verified: B73 BC 6@245 = +2.07; Mo17 BC 5@535 = −2.11). Report it only if useful, explicitly labeled as the confounded contrast — never as `a`.
- **Internal-code trap (must handle).** R/qtl relabels both BCs internally as `AA`/`AB`; for **Mo17 BC** internal `AA` = the real **B73/Mo17 het** and internal `AB` = **Mo17/Mo17**. So `09`/`11` must map internal genotype codes → real B73/Mo17 genotypes **per population** before assigning sign, and **verify every reported effect's sign against `effectplot` genotype-class means** (as done during planning). `fitqtl`'s `get.ests` sign is only trusted after this per-population mapping/verification; effects are negated as needed so all are on the B73-allele/disease-BLUP convention.

### `08_identify_epistatic_qtl.R` — call significant epistatic QTL (CLI: `<population> [alpha]`)
- Read `analyses/qtl_analyses/<pop>_scantwo.RDS` (`list(scan, permutations)`; `permutations` keyed by pheno name — see [07_epistatic_qtl.R](scripts/07_epistatic_qtl.R#L98)).
- Per phenotype: interaction-LOD thresholds via `summary(perms, alpha=c(0.05,0.20))[,"int"]`; top interacting pairs via `summary(scan, what="int", thresholds=0)`; per-pair permutation p = `mean(perms$int >= lod.int)`.
- Add a **`same_chr_close`** flag: `chr1==chr2 & abs(pos1-pos2) < min_sep_cM` (default ~20 cM) → the linked-marker artifact; keep them in the output but mark them non-genuine.
- Add **`sig_level`** ∈ {`0.05`,`0.20`,`ns`} so near-significant pairs are visible.
- Follows the legacy `summary(what="int", ...)` idea from [get_peaks_scantwo.R](scripts/get_peaks_scantwo.R) but adapted to 07's RDS shape.
- Output: `analyses/epistatic_peaks.csv` — cols `cross, trait, chr1, pos1, chr2, pos2, lod.int, lod.full, int_p, sig_level, same_chr_close`.

### `09_estimate_qtl_effects.R` — additive & dominance effects at peaks (CLI: `<population>`)
- Read `analyses/main_effect_peaks.csv` (strip the `"LOD "` prefix from `trait` to get the pheno column name) and, for genuine significant interactions, `analyses/epistatic_peaks.csv`.
- Reuse `read_cross()`/RIL-conversion from [06_identify_qtl.R](scripts/06_identify_qtl.R#L17); `sim.geno(cross, step=2.5)`.
- For each phenotype present for that population: `makeqtl` from that pheno's peaks, build the `fitqtl` formula adding `Qi:Qj` terms for any significant epistatic pairs (formula-builder from [qtl_effect_all.Rmd](scripts/get_effects/qtl_effect_all.Rmd#L69)), `fitqtl(get.ests=TRUE)`, extract estimate + drop-one LOD.
- Apply the **B73-allele / disease-BLUP convention** above: map internal codes → real genotypes per population, negate the estimate where needed so `+` = B73 allele raises BLUP, and cross-check the sign against `effectplot` marginal means. Label each estimate by what it is: **RIL BLUP → `a`**; **BC MPH → `d`**; **BC BLUP → confounded contrast (`d−a` / `−a−d`)**, not additive.
- Also fit the **RIL** cross at any BC-only peak positions to supply their `a` (needed by `10`).
- Output: `analyses/qtl_effects.csv` — peaks with `cross, trait, chr, pos, ci_lo, ci_hi, lod, estimate, estimate_type` (`a`/`d`/`confounded`), `fitqtl_lod`.

### `10_qtl_gene_action.R` — colocalize across populations + classify gene action (no population arg)
- **Colocalization (what you did before, now automated):** the legacy [overlapping_qtl.Rmd](scripts/get_effects/overlapping_qtl.Rmd) computed CI overlaps per chromosome with `IRanges::findOverlaps` on `[ci_lo, ci_hi]`, **but then assigned the shared QTL id integers by hand** (hardcoded per-chromosome vectors) and hand-annotated prior-literature matches. This script replaces the manual id assignment with automatic connected-component grouping of overlapping CIs (still `IRanges::findOverlaps`), yielding a `qtl_id` per colocalized cluster. With only ~6 peaks and 1 trait this is small and fully automatable.
- Join additive `a` (RIL, `estimate_type=="a"`) with dominance `d` (each BC's MPH, `estimate_type=="d"`) by `qtl_id` — all already on the B73-allele/disease-BLUP convention from `09`. Compute `B73 d/a`, `Mo17 d/a`; classify gene action with the legacy `case_when` cutoffs 0.2 / 0.8 / 1.2 (additive / partial-dominant / dominant / over- / under-dominant), direction-aware, from [qtl_effect_all.Rmd](scripts/get_effects/qtl_effect_all.Rmd#L319).
- Output: `analyses/qtl_gene_action.csv` — combined table `qtl_id, chr, pos, a, B73_d, B73_d/a, B73_action, Mo17_d, Mo17_d/a, Mo17_action, sig_in`.

### `11_genome_wide_effect_scan.R` — sliding effect profile (CLI: `<population>`)
- Reuse the sliding "dummy-QTL" `fitqtl` scan: `sim.geno(cross, step=2.5)`, iterate the imputed map (`attr(cross$geno[[k]]$draws, "map")`), at each position add a dummy QTL to a `makeqtl` built from the population's real peaks (dropping any real peak within a window of the test position), `fitqtl(get.ests=TRUE)`, record the dummy QTL's effect. Adapt `mqm_effect_scan()` from [get_effects_mqm.R](scripts/get_effects_mqm.R#L1) to read `analyses/main_effect_peaks.csv` + `analyses/epistatic_peaks.csv` instead of the stale `peaks_mqm.csv`/`peaks_scantwo.csv`.
- Output: `analyses/qtl_effects_whole_genome.csv` — one effect column per population/phenotype (e.g. `NLB_ril_effect`, `NLB_b73_effect`, `NLB_b73_mph_effect`, …), for plotting effect vs LOD along the genome.

## Files
- **New:** `scripts/08_identify_epistatic_qtl.R`, `scripts/09_estimate_qtl_effects.R`, `scripts/10_qtl_gene_action.R`, `scripts/11_genome_wide_effect_scan.R`.
- **Reference/reuse (not modified):** `scripts/07_epistatic_qtl.R` (skeleton), `scripts/06_identify_qtl.R` (`read_cross`), `scripts/get_effects/rqtl_mqm_effect.Rmd`/`qtl_effect_all.Rmd`/`get_qtl_effects_frreal.R`, `scripts/get_effects_mqm.R`, `scripts/get_peaks_scantwo.R`.
- **Update after implementation:** add scripts 08–11 to the pipeline list in [CLAUDE.md](CLAUDE.md).
- New packages beyond current deps: `IRanges` (or `GenomicRanges`) for `10`; `tidyr` for reshaping. `qtl` already provides `fitqtl`/`makeqtl`/`sim.geno`/`scantwo`.

## Verification
Run from the project root (all inputs already exist on disk):
1. `Rscript scripts/08_identify_epistatic_qtl.R RIL` (then `B73_BC`, `Mo17_BC`) → inspect `analyses/epistatic_peaks.csv`: RIL should flag c1×c3 (`int_p`≈0) and c3×c5 as genuine `sig_level=0.05`; same-chr close pairs flagged `same_chr_close=TRUE`.
2. `Rscript scripts/09_estimate_qtl_effects.R RIL` (+ BC pops) → `analyses/qtl_effects.csv`: peak rows with finite `estimate`, correct `estimate_type` (`a`/`d`/`confounded`), and `fitqtl_lod`. Confirm the B73-allele/disease-BLUP convention by spot-checking against `effectplot` marginal means — e.g. RIL 9@100 should give `a ≈ +1.93`; Mo17 BC mapping must treat internal `AA` as the B73/Mo17 het.
3. `Rscript scripts/10_qtl_gene_action.R` → `analyses/qtl_gene_action.csv`: colocalized `qtl_id`s with `a`, `d`, `d/a`, and an `action` label in {additive, pd, dominant, od, ud}.
4. `Rscript scripts/11_genome_wide_effect_scan.R RIL` (+ BC pops) → `analyses/qtl_effects_whole_genome.csv`: one effect value per map position; `plot(pos, effect, type="l")` should peak near the significant QTL and cross ~0 elsewhere.

## Open follow-ups
- Mo17 MPH epistasis is only near-significant (p≈0.075); `08`'s α=0.20 output will surface it so you can decide whether to run more `07` permutations.

---

# Appendix — facts, reference values, and test snippets (verified during planning)

Everything below was confirmed by running R against the on-disk data on 2026-07-26. A new session can trust these without re-deriving, and use the snippets as regression checks after implementing each script.

## A. Data facts a cold session needs
- **Run from** `~/projects/nlb_mdh`. `Rscript` at `/usr/local/bin/Rscript`.
- **Phenotype columns (identical order in `RIL_cross.csv`, `B73_cross.csv`, `Mo17_cross.csv`):** `1=NLB_WMD_BLUP`, `2=Line`, `3=NLB_WMD_BLUP_BPH`, `4=NLB_WMD_BLUP_BPH_PCT`, `5=NLB_WMD_BLUP_MPH`, `6=NLB_WMD_BLUP_MPH_PCT`. Only one disease trait. Presets use `phenos=1` (RIL) and `phenos=c(1,5)` (BC = BLUP + MPH).
- **Per-population cross config** (also in `07_epistatic_qtl.R` `presets`):
  - RIL: `input_file="analyses/RIL_cross.csv"`, `genotype=c("A","B")`, `na.strings="-"`, `crosstype="ril"` (read as `bc` then `convert2riself`). A=B73, B=Mo17.
  - B73_BC: `analyses/B73_cross.csv`, `genotype=c("AA","AB")`, `na.strings="A-"`, `crosstype="bc"`. AA=B73/B73, AB=B73/Mo17.
  - Mo17_BC: `analyses/Mo17_cross.csv`, `genotype=c("AB","BB")`, `na.strings="-B"`, `crosstype="bc"`. AB=B73/Mo17, BB=Mo17/Mo17. **Trap:** R/qtl internally relabels these as `AA`/`AB`, so internal `AA`=real B73/Mo17 het, internal `AB`=real Mo17/Mo17.
  - `alleles=c("A","B")` always; `output_dir="analyses/qtl_analyses/"`.
- **RDS structures:**
  - Main-effect (`05` → `analyses/qtl_analyses/<pop>.RDS`): `list(scan, permutations)`; `permutations` indexed by **numeric pheno column** (NULL gaps); each perm's `colnames()[1]` is `"LOD <trait>"`.
  - Epistasis (`07` → `analyses/qtl_analyses/<pop>_scantwo.RDS`): `list(scan, permutations)`; `scan` is a `scantwo` object (multi-pheno for BC: `dim(scan$lod)=1339x1339x2`); `permutations` keyed by **pheno name** (`"NLB_WMD_BLUP"`, `"NLB_WMD_BLUP_MPH"`), each a `scantwoperm` with `$full/$fv1/$int/$add/$av1`.
- **`analyses/main_effect_peaks.csv` (current, 6 peaks)** — columns `"",lodindex,trait,chr,pos,lod,ci_lo,ci_hi,cross`; `trait` carries a `"LOD "` prefix (strip it):
  | cross | trait | chr | pos | lod | ci_lo | ci_hi |
  |---|---|---|---|---|---|---|
  | RIL | LOD NLB_WMD_BLUP | 6 | 245 | 4.30 | 235.2 | 255 |
  | RIL | LOD NLB_WMD_BLUP | 9 | 100 | 6.32 | 90 | 112.6 |
  | B73_BC | LOD NLB_WMD_BLUP | 6 | 245 | 3.65 | 235 | 284.0 |
  | B73_BC | LOD NLB_WMD_BLUP_MPH | 1 | 480 | 4.17 | 466.9 | 495 |
  | Mo17_BC | LOD NLB_WMD_BLUP | 5 | 535 | 4.28 | 515 | 549.3 |
  | Mo17_BC | LOD NLB_WMD_BLUP_MPH | 1 | 760 | 3.33 | 735.2 | 917.9 |

## B. Epistasis reference values (expected output of `08`)
Interaction-LOD permutation thresholds and observed genome-wide max interaction LOD (positions maximizing interaction), per population/phenotype:

| pop | pheno | int_thr 5% | int_thr 20% | max obs int LOD | perm p | genuine? |
|---|---|---|---|---|---|---|
| RIL | NLB_WMD_BLUP | 4.43 | 3.70 | 6.69 | 0.000 | **yes** — c1×c3 (892×413) p≈0.000; c3×c5 (390×285) p≈0.008; c7×c8 (254×153) p≈0.05 |
| B73_BC | NLB_WMD_BLUP | 4.49 | 3.66 | 5.10 | 0.017 | **yes** — c3×c8 (138×568) p≈0.042; c4×c9 (447×220) p≈0.05 (c3×c3 is same-chr, 38 vs 250 cM) |
| B73_BC | NLB_WMD_BLUP_MPH | 6.71 | 6.35 | 7.67 | 0.000 | **artifact** — top pairs c1×c1/c9×c9/c5×c5 with pos1≈pos2 (1–2 cM apart) |
| Mo17_BC | NLB_WMD_BLUP | 6.30 | 5.24 | 6.88 | 0.033 | **artifact** — c1×c1/c9×c9/c6×c6, all same-chr ~1–2 cM apart |
| Mo17_BC | NLB_WMD_BLUP_MPH | 9.56 | 8.87 | 9.39 | 0.075 | near-sig (α=0.20 only); c9×c9 same-chr artifact |

Rule for `08`: flag `same_chr_close = chr1==chr2 & abs(pos1-pos2) < ~20 cM`; only inter-chromosomal (or well-separated) pairs below the p threshold are "genuine". Genuine significant epistasis = **RIL c1×c3, c3×c5** and **B73_BC BLUP c3×c8, c4×c9**.

## C. Effect-convention reference values (expected sanity checks for `09`)
Convention: **effect of the B73 allele on `NLB_WMD_BLUP`, disease scale, no flip** (positive = B73 raises BLUP). Genotype-class means from `effectplot` (after `sim.geno(step=2)`):
- RIL 9@100: B73/B73=89.79, Mo17/Mo17=85.94 → **a = +1.93** (clean additive).
- RIL 6@245: B73/B73 > Mo17/Mo17 (small positive a).
- B73_BC 6@245 (BLUP): B73/B73=87.77, B73/Mo17=89.84 → contrast +2.07 = **d − a** (confounded).
- B73_BC 1@480 (MPH): "hom"=3.94, het=5.27 → positive **d**.
- Mo17_BC 5@535 (BLUP): B73/Mo17=94.46, Mo17/Mo17=92.35 → contrast −2.11 = **−a − d** (confounded; remember internal `AA`=het).

## D. Copy-paste verification snippets

**D1. Re-check epistasis (regression test for `08`):**
```r
suppressMessages(library(qtl))
for (pop in c("RIL","B73_BC","Mo17_BC")) {
  o <- readRDS(paste0("analyses/qtl_analyses/", pop, "_scantwo.RDS"))
  for (nm in names(o$permutations)) {
    idx   <- match(nm, names(o$permutations))
    perms <- o$permutations[[nm]]
    s     <- tryCatch(subset(o$scan, lodcolumn = idx), error=function(e) o$scan)
    thr   <- as.data.frame(summary(perms, alpha=c(0.05,0.20)))
    top   <- summary(s, thresholds=c(0,0,0,0,0), what="int")
    top   <- top[order(-top$lod.int),]; t3 <- head(top,3)
    t3$int_p <- sapply(t3$lod.int, function(x) mean(as.numeric(perms$int) >= x))
    cat("\n==", pop, nm, " int_thr5%=", round(thr[1,3],2), "\n", sep="")
    print(t3[,c("chr1","chr2","pos1","pos2","lod.int","int_p","lod.full")], digits=3, row.names=FALSE)
  }
}
```
Expect the Section B pairs/p-values.

**D2. Re-check effect direction (regression test for `09`):**
```r
suppressMessages(library(qtl))
rc <- function(f,g,na,ct){ ril<-ct=="ril"; if(ril) ct<-"bc"
  x<-read.cross("csv","",f,genotype=g,alleles=c("A","B"),na.strings=na,crosstype=ct)
  x<-jittermap(x); if(ril) x<-convert2riself(x)
  x$pheno<-as.data.frame(lapply(x$pheno,function(v) as.numeric(as.character(v)))); x }
ril <-sim.geno(rc("analyses/RIL_cross.csv", c("A","B"),  "-", "ril"),step=2)
b73 <-sim.geno(rc("analyses/B73_cross.csv", c("AA","AB"),"A-","bc"), step=2)
mo17<-sim.geno(rc("analyses/Mo17_cross.csv",c("AB","BB"),"-B","bc"), step=2)
effectplot(ril, "NLB_WMD_BLUP",     mname1="9@100", draw=FALSE)$Means   # B73/B73 89.79 > Mo17/Mo17 85.94 -> a=+1.93
effectplot(b73, "NLB_WMD_BLUP",     mname1="6@245", draw=FALSE)$Means   # AA 87.77, AB 89.84 -> d-a
effectplot(b73, "NLB_WMD_BLUP_MPH", mname1="1@480", draw=FALSE)$Means   # d
effectplot(mo17,"NLB_WMD_BLUP",     mname1="5@535", draw=FALSE)$Means   # internal AA=het 94.46, AB=Mo17hom 92.35
```

**D3. End-to-end (after implementing 08–11):**
```sh
Rscript scripts/08_identify_epistatic_qtl.R RIL
Rscript scripts/08_identify_epistatic_qtl.R B73_BC
Rscript scripts/08_identify_epistatic_qtl.R Mo17_BC   # -> analyses/epistatic_peaks.csv
Rscript scripts/09_estimate_qtl_effects.R RIL
Rscript scripts/09_estimate_qtl_effects.R B73_BC
Rscript scripts/09_estimate_qtl_effects.R Mo17_BC     # -> analyses/qtl_effects.csv
Rscript scripts/10_qtl_gene_action.R                  # -> analyses/qtl_gene_action.csv
Rscript scripts/11_genome_wide_effect_scan.R RIL      # (+ BC pops) -> analyses/qtl_effects_whole_genome.csv
```
Pass criteria: `epistatic_peaks.csv` marks RIL c1×c3 & c3×c5 genuine and flags same-chr pairs; `qtl_effects.csv` gives RIL 9@100 `a≈+1.93` with `estimate_type=="a"`, BC MPH rows `estimate_type=="d"`, BC BLUP rows `estimate_type=="confounded"`; `qtl_gene_action.csv` has an `action` label per colocalized `qtl_id`; `qtl_effects_whole_genome.csv` effect profile peaks near significant QTL.

## E. Reference-code index (legacy `scripts/get_effects/`, from `~/projects/mdh_qtl`)
- `rqtl_mqm_effect.Rmd` `run_fitqtl()` — simplest `makeqtl`→`fitqtl(get.ests=TRUE)`→`summary$ests[-1,1]` (main effects only).
- `qtl_effect_all.Rmd` / `get_qtl_effects_frreal.R` — epistasis-aware `fitqtl` with programmatic `y ~ Q1+Q2+...+Qi:Qj` formula; drop-one LOD from `summary$result.drop[,"LOD"]`; d/a gene-action `case_when` at 0.2/0.8/1.2.
- `get_effects_mqm.R` (repo root, partial NLB port) — sliding "dummy-QTL" whole-genome effect scan; **stale below the RIL/B73 blocks** (references undefined `lod`, `combined_qtl_ril`, `epi_qtl`) — fix when adapting for `11`.
- `overlapping_qtl.Rmd` — CI-overlap grouping via `IRanges::findOverlaps`, but QTL ids were assigned **by hand**; `10` automates this.
- `get_peaks_scantwo.R` — legacy epistasis peak caller (`summary(what="int", perms=, alphas=c(1,0,0.05,0,0))`); `08` adapts it to `07`'s RDS shape.
- **Missing helpers:** `run_rqtl2_mqm()` / `get_effect_mqm()` (sourced as `scripts/get_effect_mqm.R`) live only in `~/projects/mdh_qtl` — do not depend on them; the self-contained logic above is sufficient.
</content>
</invoke>

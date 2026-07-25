# Plan: Make `07_epistatic_qtl.R` CLI-runnable and run B73/Mo17 backcrosses

## Context

`scripts/07_epistatic_qtl.R` is the `scantwo` two-QTL / epistasis stage of the NLB QTL pipeline.
Today it has two problems:

1. **It only runs interactively for RILs.** The config (input file, genotypes, phenotypes, etc.)
   is a single hardcoded block, so running the B73 and Mo17 backcrosses means hand-editing the
   script each time.
2. **Its permutation-combining section is broken.** The bottom of the script is copy-pasted from
   the MQM script (05) and references variables that don't exist in a `scantwo` context
   (`mqm_temp_file`, `augmentedcross`, `mqmprocesspermutation`). The final output line also reuses
   `{output_name}.RDS`, which would clobber the **MQM** output `analyses/qtl_analyses/RIL.RDS`
   written by 05.

Because of (2), the completed RIL run left its temp files —
`analyses/qtl_analyses/tmp/RIL.scantwo.tmp.RDS` (main scan) plus 10 batches × 10 =
**100 permutations already computed** (`RIL.trait.NLB_WMD_BLUP.batch.{1..10}.scantwo.perm.tmp`) —
but **never produced a final combined output**. RIL therefore only needs *finalizing*, not
recomputing.

Goal: rewrite 07 to take a population name on the command line, fix the `scantwo` permutation
combine/save, finalize the existing RIL result, and enable running B73 BC and Mo17 BC (100
permutations each, matching RIL).

## Approach

Rewrite `scripts/07_epistatic_qtl.R` with a **population-preset CLI**.

### 1. CLI + presets
Parse `commandArgs(trailingOnly = TRUE)`:
- `args[1]` = population name (`RIL`, `B73_BC`, or `Mo17_BC`) → selects a preset.
- `args[2]` = permutations, optional, **default 100**.
- `args[3]` = cores, optional, **default 4** (see §6).

Hold the three known configs (from 05/06) as a named list, keyed by population name:

| pop | input_file | genotype | na.strings | crosstype | phenos |
|-----|-----------|----------|-----------|-----------|--------|
| RIL | analyses/RIL_cross.csv | A,B | `-` | ril | 1 |
| B73_BC | analyses/B73_cross.csv | AA,AB | `A-` | bc | 1,5 |
| Mo17_BC | analyses/Mo17_cross.csv | AB,BB | `-B` | bc | 1,5 |

(`alleles = c("A","B")` for all.) `phenos` are **numeric column indices** — col 1 =
`NLB_WMD_BLUP`, col 5 = `NLB_WMD_BLUP_MPH`. Fixed: `output_dir = "analyses/qtl_analyses/"`,
`output_name = <pop>`. Keep example `Rscript` invocations as comments at the top (mirrors
`scripts/run_rqtl_scantwo.R`).

### 2. Cross loading (unchanged logic)
Keep the existing `ril`→`convert2riself` handling and `read.cross` call, and the
`mutate(across(everything(), as.numeric))` phenotype coercion (current lines 18–39).

### 3. Main scan — make resumable
Guard the expensive main scan so re-runs don't recompute it:
```r
scantwo_temp_file <- paste0(output_dir, "tmp/", output_name, ".scantwo.tmp.RDS")
if (file.exists(scantwo_temp_file)) {
  out2 <- readRDS(scantwo_temp_file)
} else {
  out2 <- scantwo(cross, n.cluster = cores, pheno.col = phenos, verbose = TRUE)
  saveRDS(out2, scantwo_temp_file)
}
```

### 4. Permutation batches — match existing filenames
Keep the `file.exists()`-guarded batch loop, but derive the filename from the **trait name** via
`colnames(cross$pheno)[i]` (as 05 does), and pass the **numeric** `pheno.col = i`:
```r
for (i in phenos) {
  pheno_name <- colnames(cross$pheno)[i]
  for (j in 1:(permutations/10)) {
    batch_filename <- paste0(output_dir,"tmp/",output_name,".trait.",pheno_name,".batch.",j,".scantwo.perm.tmp")
    if (!file.exists(batch_filename)) {
      results <- scantwo(cross, n.cluster=cores, n.perm=10, pheno.col=i, verbose=TRUE)
      saveRDS(results, batch_filename)
    }
  }
}
```
For RIL (`phenos = 1`) this reproduces exactly the existing
`RIL.trait.NLB_WMD_BLUP.batch.*.scantwo.perm.tmp` names, so the 100 done perms are **skipped**.

### 5. Fix the combine + save (the broken section)
Replace the MQM copy-paste with `scantwo`-correct combining (`c()` concatenates `scantwoperm`
objects — see `scripts/combine_perms.R`):
```r
perm_results <- list()
for (i in phenos) {
  pheno_name <- colnames(cross$pheno)[i]
  perm_list <- lapply(1:(permutations/10), function(j)
    readRDS(paste0(output_dir,"tmp/",output_name,".trait.",pheno_name,".batch.",j,".scantwo.perm.tmp")))
  perm_results[[pheno_name]] <- Reduce(c, perm_list)
}
output <- list(scan = out2, permutations = perm_results)
saveRDS(output, paste0(output_dir, output_name, "_scantwo.RDS"))
```
Note the **`_scantwo.RDS`** suffix — distinct from the MQM `{output_name}.RDS`, so 05's output is
not overwritten. The `list(scan=, permutations=)` shape mirrors 05's output convention.

### 6. Cores
Make cores an **optional `args[3]` with a modest default** so a background run leaves the machine
usable for other work — do **not** use `detectCores() - 1` (it hogs the box). Default
`cores <- 4`; override per-invocation, e.g. `Rscript scripts/07_epistatic_qtl.R B73_BC 100 6`.

## Files
- **Rewrite:** `scripts/07_epistatic_qtl.R` (only file changed).
- Reference patterns reused: `scripts/05_qtl_analysis.R` (batch/temp naming, phenotype coercion),
  `scripts/combine_perms.R` (`Reduce(c,...)`), `scripts/run_rqtl_scantwo.R` (`commandArgs` style).

## Verification
1. **Finalize RIL (fast, no recompute):**
   `Rscript scripts/07_epistatic_qtl.R RIL`
   Expect: main scan loaded from the existing `.scantwo.tmp.RDS`, all 10 perm batches skipped
   (`file.exists`), and a new `analyses/qtl_analyses/RIL_scantwo.RDS` written. This exercises the
   whole rewritten path in minutes and confirms MQM `RIL.RDS` is untouched.
2. **Check output structure:** in R, `x <- readRDS("analyses/qtl_analyses/RIL_scantwo.RDS")`;
   confirm `class(x$scan)` is `scantwo` and `summary(x$permutations[["NLB_WMD_BLUP"]])` returns
   per-model thresholds (full/fv1/int/add/av1). This is the shape the legacy
   `scripts/get_peaks_scantwo.R` consumes.
3. **BC quick smoke test (before the full run):** run each BC population with a tiny perm count
   first — `Rscript scripts/07_epistatic_qtl.R B73_BC 10` and `... Mo17_BC 10`. This exercises the
   full path for the BC configs end-to-end: cross loads with the BC genotype/na.strings, the main
   `scantwo` runs and saves its `.scantwo.tmp.RDS`, one perm batch per trait (cols 1 and 5)
   is written, and `{pop}_scantwo.RDS` is produced with two entries in `$permutations`. Confirm
   `summary()` on each works. Because runs are resumable, this work is **not wasted** — the full
   run reuses the saved main scan and this first batch.
4. **BC full runs (deferred — long-running, user will launch):**
   `Rscript scripts/07_epistatic_qtl.R B73_BC` and `... Mo17_BC` (100 perms).
   Each is a fresh main scan + 100 perms over two traits; expect ~a day+ each (RIL's 100 perms
   took ~20 h). Both are resumable — safe to stop/restart; existing batches are skipped.

## Notes / decisions
- Population-preset CLI (not full positional args).
- 100 permutations for BC (matches RIL).
- Script rewrite only — no runs launched as part of implementing this plan. Commands above are
  ready to launch when you choose.

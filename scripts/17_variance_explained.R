require(qtl)
require(dplyr)

# sim.geno imputes missing/pseudomarker genotypes by Monte Carlo, so the %var
# wiggles run-to-run; pin a seed and use plenty of draws for a stable,
# reproducible reported number (the numbered pipeline elsewhere uses the
# n.draws default, which is fine for scans but too noisy for a headline stat).
set.seed(1)
n_draws <- 256

# 17_variance_explained.R -- proportion of genetic variance explained by the
# identified NLB QTL, per population, mirroring the sister SLB/GLS paper's
# paragraph ("... explained 42.71% of the genetic variance for SLB ... in the
# RILs [additive] ... 24.49% of the MPH ... in the B73 backcross [dominance]").
#
# Method (recovered from mdh_qtl/scripts/qtl_effect_all.Rmd, helper run_fitqtl):
# the percentage is the full-model %var reported by R/qtl fitqtl,
# summary(out.qtl)$result.full["Model", "%var"]. The model fits the *union* of
# identified (colocalized) QTL for the trait in EVERY population -- RIL on the
# line BLUP (its variance ~ additive genetic variance) and each backcross on
# the mid-parent-heterosis (MPH) BLUP (its variance ~ dominance genetic
# variance) -- with genuine epistatic pairs added as Qi:Qj interaction terms.
# The additive-vs-dominance reading comes entirely from which trait is fit; no
# division by heritability is involved (the fitted trait is already a line-mean
# genetic value).
#
# We report BOTH the with-epistasis %var (faithful to mdh_qtl's model) and the
# main-effects-only %var, since the published sentence's exact %var call was
# interactive and not committed, leaving which variant it used ambiguous.
#
# No arguments; loops all three populations. Run from the project root:
#   Rscript scripts/17_variance_explained.R

# per-population config: cross file + genotype encoding (see 09_estimate_qtl_effects.R
# presets and CLAUDE.md "Populations & R/qtl cross codes") and the gene-action
# trait whose variance the fit's %var is a proportion of.
presets <- list(
  RIL = list(input_file = "analyses/RIL_cross.csv", genotype = c("A", "B"),
             na.strings = "-", crosstype = "ril",
             ga_trait = "NLB_WMD_BLUP", variance_type = "additive"),
  B73_BC = list(input_file = "analyses/B73_cross.csv", genotype = c("AA", "AB"),
                na.strings = "A-", crosstype = "bc",
                ga_trait = "NLB_WMD_BLUP_MPH", variance_type = "dominance"),
  Mo17_BC = list(input_file = "analyses/Mo17_cross.csv", genotype = c("AB", "BB"),
                 na.strings = "-B", crosstype = "bc",
                 ga_trait = "NLB_WMD_BLUP_MPH", variance_type = "dominance")
)

gene_action_file <- "analyses/qtl_gene_action.csv"
epistatic_peaks_file <- "analyses/epistatic_peaks.csv"
genuine_alpha <- "0.05"
output_file <- "analyses/qtl_variance_explained.csv"

# copied from 09_estimate_qtl_effects.R: read an R/qtl cross, coercing the
# ril crosstype to bc for read.cross then convert2riself, and numeric-ify pheno.
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

# fit the union QTL set (main effects) plus, when include_epi, the genuine
# epistatic pairs as interaction-only Qi:Qj terms (the epi loci enter makeqtl
# but not the main-effect part of the formula -- same construction as
# mdh_qtl run_fitqtl and 09_estimate_qtl_effects.R::fit_main_effects). Returns
# the full-model %var, or NA if the fit fails (e.g. a rank-deficient model).
pve_full_model <- function(cross, trait, qtl_df, epi_df, include_epi) {
  main_n <- nrow(qtl_df)
  chrs <- qtl_df$chr
  poss <- qtl_df$pos
  qnam <- paste0("Q", seq_len(main_n))

  epinam <- ""
  if (include_epi && nrow(epi_df) > 0) {
    for (j in seq_len(nrow(epi_df))) {
      idx1 <- main_n + (2 * j - 1)
      idx2 <- main_n + (2 * j)
      chrs <- c(chrs, epi_df$chr1[j], epi_df$chr2[j])
      poss <- c(poss, epi_df$pos1[j], epi_df$pos2[j])
      epinam <- paste0(epinam, "+Q", idx1, ":Q", idx2)
    }
  }

  pheno.col <- match(trait, colnames(cross$pheno))
  out <- tryCatch({
    qtlobj <- makeqtl(cross, chr = chrs, pos = poss)
    formula <- as.formula(paste0("y ~ ", paste(qnam, collapse = "+"), epinam))
    fitqtl(cross, pheno.col = pheno.col, qtl = qtlobj, formula = formula, get.ests = FALSE)
  }, error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  summary(out)$result.full["Model", "%var"]
}

# union QTL set: one representative (chr, pos) per colocalization cluster from
# script 10, the direct NLB analog of mdh_qtl's per-trait deduplicated QTL list.
gene_action <- read.csv(gene_action_file, stringsAsFactors = FALSE)
union_qtl <- gene_action %>%
  distinct(qtl_id, .keep_all = TRUE) %>%
  transmute(chr = as.numeric(chr), pos = as.numeric(pos)) %>%
  arrange(chr, pos)

# genuine epistatic pairs (union across populations, deduped on locus pair) --
# mdh_qtl filtered epi_qtl by trait only, carrying the same pairs into every
# population's fit.
genuine_epi <- data.frame()
if (file.exists(epistatic_peaks_file)) {
  genuine_epi <- read.csv(epistatic_peaks_file, stringsAsFactors = FALSE) %>%
    filter(sig_level == genuine_alpha, !same_chr_close) %>%
    distinct(chr1, pos1, chr2, pos2) %>%
    mutate(across(everything(), as.numeric))
}

results <- data.frame()
for (population in names(presets)) {
  preset <- presets[[population]]
  cross <- sim.geno(read_cross(preset$input_file, preset$genotype,
                               preset$na.strings, preset$crosstype),
                    step = 2.5, n.draws = n_draws)

  pct_epi <- pve_full_model(cross, preset$ga_trait, union_qtl, genuine_epi, include_epi = TRUE)
  pct_main <- pve_full_model(cross, preset$ga_trait, union_qtl, genuine_epi, include_epi = FALSE)

  results <- bind_rows(results, data.frame(
    population = population,
    trait = preset$ga_trait,
    variance_type = preset$variance_type,
    n_qtl = nrow(union_qtl),
    n_epi_pairs = nrow(genuine_epi),
    pct_var_with_epistasis = round(pct_epi, 2),
    pct_var_main_only = round(pct_main, 2)
  ))
}

write.csv(results, output_file, row.names = FALSE)

# paragraph-style console summary mirroring the SLB/GLS sentence
cat("\nProportion of genetic variance explained by the identified NLB QTL\n")
cat("(full-model %var from fitqtl; RIL = additive, backcrosses = dominance on MPH)\n\n")
print(results, row.names = FALSE)

ril <- results[results$population == "RIL", ]
b73 <- results[results$population == "B73_BC", ]
mo17 <- results[results$population == "Mo17_BC", ]
cat(sprintf(
  paste0("\nMain-effects-only model: the identified QTL explained %.2f%% of the ",
         "additive\ngenetic variance in the RILs, %.2f%% of the MPH (dominance ",
         "variance) in the\nB73 backcross, and %.2f%% in the Mo17 backcross.\n"),
  ril$pct_var_main_only, b73$pct_var_main_only, mo17$pct_var_main_only))
cat(sprintf(
  paste0("With epistatic interaction terms included, the comparable values are ",
         "%.2f%%,\n%.2f%%, and %.2f%%.\n"),
  ril$pct_var_with_epistasis, b73$pct_var_with_epistasis, mo17$pct_var_with_epistasis))

require(qtl)
require(dplyr)

# Rscript scripts/09_estimate_qtl_effects.R RIL
# Rscript scripts/09_estimate_qtl_effects.R B73_BC
# Rscript scripts/09_estimate_qtl_effects.R Mo17_BC

args <- commandArgs(trailingOnly = TRUE)
population <- args[1]

presets <- list(
  RIL = list(input_file = "analyses/RIL_cross.csv", genotype = c("A", "B"),
             na.strings = "-", crosstype = "ril"),
  B73_BC = list(input_file = "analyses/B73_cross.csv", genotype = c("AA", "AB"),
                na.strings = "A-", crosstype = "bc"),
  Mo17_BC = list(input_file = "analyses/Mo17_cross.csv", genotype = c("AB", "BB"),
                 na.strings = "-B", crosstype = "bc")
)

genuine_alpha <- "0.05"
main_effect_peaks_file <- "analyses/main_effect_peaks.csv"
epistatic_peaks_file <- "analyses/epistatic_peaks.csv"
output_file <- "analyses/qtl_effects.csv"

# Sign of the raw fitqtl Q-term estimate, mapped onto the B73-allele /
# disease-BLUP convention (see plans/identify_qtl_and_effects.md, "Effect
# sign & allele convention"):
#  - RIL: the additive-only genotype code runs B73=-1, Mo17=+1, so
#    a = -raw_est = (B73/B73 - Mo17/Mo17) / 2.
#  - B73_BC: R/qtl's internal "AB" is the true B73/Mo17 het, so
#    raw_est = het - hom already: BLUP -> d-a, MPH -> d. No flip.
#  - Mo17_BC: internal "AA" is the true B73/Mo17 het (the internal-code
#    trap), so raw_est = hom - het. For BLUP that's already -a-d as-is
#    (matches the confounded-contrast convention, no flip). For MPH it must
#    be flipped to report d = het - hom.
sign_mult <- function(population, trait) {
  if (population == "RIL") return(-1)
  if (population == "Mo17_BC" && trait == "NLB_WMD_BLUP_MPH") return(-1)
  1
}

estimate_type_for <- function(population, trait) {
  if (population == "RIL") return("a")
  if (trait == "NLB_WMD_BLUP_MPH") return("d")
  "confounded"
}

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

# fit one phenotype's main peaks together, adding Qi:Qj terms for any
# genuine epistatic partners so the main-peak estimates/LODs are adjusted
# for background epistasis; returns one row per main peak, same order.
fit_main_effects <- function(cross, trait, main_peaks, epi_peaks) {
  main_n <- nrow(main_peaks)
  chrs <- main_peaks$chr
  poss <- main_peaks$pos
  qnam <- paste0("Q", seq_len(main_n))

  epinam <- ""
  epi_n <- nrow(epi_peaks)
  if (epi_n > 0) {
    for (j in seq_len(epi_n)) {
      idx1 <- main_n + (2 * j - 1)
      idx2 <- main_n + (2 * j)
      chrs <- c(chrs, epi_peaks$chr1[j], epi_peaks$chr2[j])
      poss <- c(poss, epi_peaks$pos1[j], epi_peaks$pos2[j])
      epinam <- paste0(epinam, "+Q", idx1, ":Q", idx2)
    }
  }

  qtlobj <- makeqtl(cross, chr = chrs, pos = poss)
  formula <- as.formula(paste0("y ~ ", paste(qnam, collapse = "+"), epinam))
  pheno.col <- match(trait, colnames(cross$pheno))
  out.qtl <- fitqtl(cross, pheno.col = pheno.col, qtl = qtlobj, get.ests = TRUE, formula = formula)
  s <- summary(out.qtl)

  # result.drop is NULL when there's exactly one QTL term (main_n==1, no
  # epistasis partners) -- dropping the only term is the same comparison as
  # the full-vs-null model, so fall back to the full-model LOD in that case.
  if (is.null(s$result.drop)) {
    fitqtl_lod <- rep(s$result.full["Model", "LOD"], main_n)
  } else {
    fitqtl_lod <- as.numeric(s$result.drop[, "LOD"])[seq_len(main_n)]
  }

  data.frame(
    raw_estimate = as.numeric(s$ests[-1, 1])[seq_len(main_n)],
    fitqtl_lod = fitqtl_lod
  )
}

preset <- presets[[population]]
cross <- sim.geno(read_cross(preset$input_file, preset$genotype, preset$na.strings, preset$crosstype), step = 2.5)

main_effect_peaks <- read.csv(main_effect_peaks_file, stringsAsFactors = FALSE)
main_effect_peaks$trait <- sub("^LOD ", "", main_effect_peaks$trait)

genuine_epi <- data.frame()
if (file.exists(epistatic_peaks_file)) {
  epistatic_peaks <- read.csv(epistatic_peaks_file, stringsAsFactors = FALSE)
  genuine_epi <- epistatic_peaks %>% filter(sig_level == genuine_alpha, !same_chr_close)
}

peaks <- main_effect_peaks %>% filter(cross == population)

results <- data.frame()
for (trait in unique(peaks$trait)) {
  trait_peaks <- peaks %>% filter(trait == !!trait)
  trait_epi <- genuine_epi %>% filter(cross == population, trait == !!trait)
  fit <- fit_main_effects(cross, trait, trait_peaks, trait_epi)
  trait_peaks$raw_estimate <- fit$raw_estimate
  trait_peaks$fitqtl_lod <- fit$fitqtl_lod
  trait_peaks$estimate_type <- estimate_type_for(population, trait)
  trait_peaks$estimate <- trait_peaks$raw_estimate * sign_mult(population, trait)
  results <- rbind(results, trait_peaks)
}

# supplemental RIL "a" fits at this population's BC-only peak positions (no
# RIL peak already covers that location), needed by 10_qtl_gene_action.R
supplemental <- data.frame()
if (population != "RIL") {
  ril_peaks <- main_effect_peaks %>% filter(cross == "RIL")
  ril_preset <- presets[["RIL"]]
  ril_cross <- sim.geno(read_cross(ril_preset$input_file, ril_preset$genotype, ril_preset$na.strings, ril_preset$crosstype), step = 2.5)
  ril_genuine_epi <- genuine_epi %>% filter(cross == "RIL", trait == "NLB_WMD_BLUP")

  already_covered <- mapply(function(c, p) any(ril_peaks$chr == c & abs(ril_peaks$pos - p) <= 5),
                             peaks$chr, peaks$pos)
  bc_only <- peaks[!already_covered, ]
  bc_only <- bc_only[!duplicated(bc_only[, c("chr", "pos")]), ]

  if (nrow(bc_only) > 0) {
    fit <- fit_main_effects(ril_cross, "NLB_WMD_BLUP", bc_only, ril_genuine_epi)
    supplemental <- data.frame(
      cross = "RIL",
      trait = "NLB_WMD_BLUP",
      chr = bc_only$chr,
      pos = bc_only$pos,
      # NA marks this as a borrowed fit at a BC peak's position, not an
      # independently significant RIL peak (kept for 10's colocalization
      # join, but excluded from its "sig_in" significance tally)
      lod = NA,
      ci_lo = bc_only$ci_lo,
      ci_hi = bc_only$ci_hi,
      raw_estimate = fit$raw_estimate,
      fitqtl_lod = fit$fitqtl_lod,
      estimate_type = "a"
    )
    supplemental$estimate <- supplemental$raw_estimate * sign_mult("RIL", "NLB_WMD_BLUP")
  }
}

results <- bind_rows(results, supplemental) %>%
  select(cross, trait, chr, pos, ci_lo, ci_hi, lod, estimate, estimate_type, fitqtl_lod)

# merge into the combined output, replacing any existing rows with the same
# (cross, trait, chr, pos) identity so re-running a population is idempotent
if (file.exists(output_file)) {
  existing <- read.csv(output_file, stringsAsFactors = FALSE)
  key <- function(df) paste(df$cross, df$trait, df$chr, df$pos)
  existing <- existing[!key(existing) %in% key(results), ]
  results <- rbind(existing, results)
}

write.csv(results, output_file, row.names = FALSE)

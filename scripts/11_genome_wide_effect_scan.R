require(qtl)
require(dplyr)

# Rscript scripts/11_genome_wide_effect_scan.R RIL
# Rscript scripts/11_genome_wide_effect_scan.R B73_BC
# Rscript scripts/11_genome_wide_effect_scan.R Mo17_BC

args <- commandArgs(trailingOnly = TRUE)
population <- args[1]

min_sep_cM <- 20
genuine_alpha <- "0.05"
main_effect_peaks_file <- "analyses/main_effect_peaks.csv"
epistatic_peaks_file <- "analyses/epistatic_peaks.csv"
output_file <- "analyses/qtl_effects_whole_genome.csv"

presets <- list(
  RIL = list(input_file = "analyses/RIL_cross.csv", genotype = c("A", "B"),
             na.strings = "-", crosstype = "ril", phenos = 1),
  B73_BC = list(input_file = "analyses/B73_cross.csv", genotype = c("AA", "AB"),
                na.strings = "A-", crosstype = "bc", phenos = c(1, 5)),
  Mo17_BC = list(input_file = "analyses/Mo17_cross.csv", genotype = c("AB", "BB"),
                 na.strings = "-B", crosstype = "bc", phenos = c(1, 5))
)

# same raw-fitqtl-estimate -> reported-estimate sign convention as 09: the
# B73-allele effect on the resistance-scale BLUP (NLB_WMD_BLUP = 100 - wmd,
# positive = B73 raises resistance)
# (see plans/identify_qtl_and_effects.md, "Effect sign & allele convention")
sign_mult <- function(population, trait) {
  if (population == "RIL") return(-1)
  if (population == "Mo17_BC" && trait == "NLB_WMD_BLUP_MPH") return(-1)
  1
}

column_suffix <- function(population, trait) {
  base <- switch(population, RIL = "ril", B73_BC = "b73", Mo17_BC = "mo17")
  if (grepl("_MPH$", trait)) base <- paste0(base, "_mph")
  base
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

# fit main peaks (minus any within min_sep_cM of the dummy position on the
# same chromosome) + genuine epistatic partners (same exclusion) + one dummy
# QTL at (chr, pos); returns the dummy QTL's raw estimate
dummy_qtl_effect <- function(cross, trait, chr, pos, main_peaks, epi_peaks) {
  near_dummy <- function(c, p) c == chr & abs(p - pos) < min_sep_cM

  keep_main <- main_peaks[!near_dummy(main_peaks$chr, main_peaks$pos), , drop = FALSE]
  keep_epi <- epi_peaks
  if (nrow(keep_epi) > 0) {
    drop_epi <- near_dummy(keep_epi$chr1, keep_epi$pos1) | near_dummy(keep_epi$chr2, keep_epi$pos2)
    keep_epi <- keep_epi[!drop_epi, , drop = FALSE]
  }

  main_n <- nrow(keep_main)
  dummy_idx <- main_n + 1
  chrs <- c(keep_main$chr, chr)
  poss <- c(keep_main$pos, pos)
  qnam <- paste0("Q", seq_len(main_n + 1))

  epinam <- ""
  epi_n <- nrow(keep_epi)
  if (epi_n > 0) {
    for (j in seq_len(epi_n)) {
      idx1 <- main_n + 1 + (2 * j - 1)
      idx2 <- main_n + 1 + (2 * j)
      chrs <- c(chrs, keep_epi$chr1[j], keep_epi$chr2[j])
      poss <- c(poss, keep_epi$pos1[j], keep_epi$pos2[j])
      epinam <- paste0(epinam, "+Q", idx1, ":Q", idx2)
    }
  }

  qtlobj <- makeqtl(cross, chr = chrs, pos = poss)
  formula <- as.formula(paste0("y ~ ", paste(qnam, collapse = "+"), epinam))
  pheno.col <- match(trait, colnames(cross$pheno))
  out.qtl <- tryCatch(
    fitqtl(cross, pheno.col = pheno.col, qtl = qtlobj, get.ests = TRUE, formula = formula),
    error = function(e) NULL
  )
  if (is.null(out.qtl)) return(NA_real_)

  ests <- summary(out.qtl)$ests[-1, 1]
  as.numeric(ests[dummy_idx])
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

# genomic-position grid, from this cross's own imputed map
grid <- data.frame()
for (k in seq_along(cross$geno)) {
  chr_name <- names(cross$geno)[k]
  map <- attr(cross$geno[[k]]$draws, "map")
  grid <- rbind(grid, data.frame(chr = as.numeric(chr_name), pos = as.numeric(map)))
}
# round so the join key survives a write.csv/read.csv roundtrip (which
# truncates to ~7 significant digits) on later re-runs for other populations
grid$pos <- round(grid$pos, 4)

for (trait_col in preset$phenos) {
  trait <- colnames(cross$pheno)[trait_col]
  trait_peaks <- main_effect_peaks %>% filter(cross == population, trait == !!trait)
  trait_epi <- genuine_epi %>% filter(cross == population, trait == !!trait)

  effect <- mapply(function(c, p) dummy_qtl_effect(cross, trait, c, p, trait_peaks, trait_epi),
                    grid$chr, grid$pos)
  grid[[paste0("NLB_", column_suffix(population, trait), "_effect")]] <-
    as.numeric(effect) * sign_mult(population, trait)
}

# merge into the combined output (one row per chr/pos across all crosses),
# replacing this population's effect columns if present from a prior run
if (file.exists(output_file)) {
  existing <- read.csv(output_file, stringsAsFactors = FALSE)
  new_cols <- setdiff(colnames(grid), c("chr", "pos"))
  existing <- existing[, !colnames(existing) %in% new_cols, drop = FALSE]
  grid <- full_join(existing, grid, by = c("chr", "pos"))
}

grid <- grid %>% arrange(chr, pos)
write.csv(grid, output_file, row.names = FALSE)

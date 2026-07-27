require(qtl)
require(dplyr)

# Rscript scripts/12_epistatic_effect_table.R
# No population arg -- reads all significant epistatic pairs from
# epistatic_peaks.csv and loops over the populations present, loading each
# cross once.

presets <- list(
  RIL = list(input_file = "analyses/RIL_cross.csv", genotype = c("A", "B"),
             na.strings = "-", crosstype = "ril"),
  B73_BC = list(input_file = "analyses/B73_cross.csv", genotype = c("AA", "AB"),
                na.strings = "A-", crosstype = "bc"),
  Mo17_BC = list(input_file = "analyses/Mo17_cross.csv", genotype = c("AB", "BB"),
                 na.strings = "-B", crosstype = "bc")
)

epistatic_peaks_file <- "analyses/epistatic_peaks.csv"
long_output_file <- "analyses/epistatic_effects_long.csv"
wide_output_file <- "analyses/epistatic_effects_wide.csv"

# Genotype-class label vectors (in R/qtl internal code order 1, 2) and their
# B73-allele counts, per population -- see CLAUDE.md "Effect sign convention"
# and plans/epistatic_effect_table.md for the internal-code trap in Mo17_BC.
geno_map <- list(
  RIL = list(labels = c("B73/B73", "Mo17/Mo17"), counts = c(2, 0)),
  B73_BC = list(labels = c("B73/B73", "B73/Mo17"), counts = c(2, 1)),
  Mo17_BC = list(labels = c("B73/Mo17", "Mo17/Mo17"), counts = c(1, 0))
)

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

epistatic_peaks <- read.csv(epistatic_peaks_file, stringsAsFactors = FALSE)
sig_pairs <- epistatic_peaks %>% filter(sig_level == "0.05", !same_chr_close)

long_rows <- data.frame()
for (population in unique(sig_pairs$cross)) {
  preset <- presets[[population]]
  gm_map <- geno_map[[population]]
  cross <- sim.geno(read_cross(preset$input_file, preset$genotype, preset$na.strings,
                                preset$crosstype), step = 2.5)

  pop_pairs <- sig_pairs %>% filter(cross == population)
  for (trait in unique(pop_pairs$trait)) {
    gm <- mean(cross$pheno[[trait]], na.rm = TRUE)
    trait_pairs <- pop_pairs %>% filter(trait == !!trait)

    for (i in seq_len(nrow(trait_pairs))) {
      pair <- trait_pairs[i, ]
      ep <- effectplot(cross, pheno.col = match(trait, colnames(cross$pheno)),
                        mname1 = sprintf("%s@%s", pair$chr1, pair$pos1),
                        mname2 = sprintf("%s@%s", pair$chr2, pair$pos2),
                        geno1 = gm_map$labels, geno2 = gm_map$labels, draw = FALSE)

      for (r in seq_along(gm_map$counts)) {
        for (c in seq_along(gm_map$counts)) {
          long_rows <- bind_rows(long_rows, data.frame(
            cross = population, trait = trait,
            chr1 = pair$chr1, pos1 = pair$pos1,
            chr2 = pair$chr2, pos2 = pair$pos2,
            geno1 = gm_map$counts[r], geno2 = gm_map$counts[c],
            mean_dev = ep$Means[r, c] - gm,
            se = ep$SEs[r, c],
            lod_full = pair$lod.full, lod_int = pair$lod.int,
            sig_level = pair$sig_level
          ))
        }
      }
    }
  }
}

write.csv(long_rows, long_output_file, row.names = FALSE)

# wide publication-style table: nine genotype cells in fixed order, formatted
# "mean_dev +/- se" rounded to 2dp; only cells populated for a given population
# are filled, the rest are "N/A"
cell_order <- c("2_2", "2_1", "2_0", "1_2", "1_1", "1_0", "0_2", "0_1", "0_0")

wide_rows <- long_rows %>%
  distinct(cross, trait, chr1, pos1, chr2, pos2, lod_full, lod_int, sig_level)

wide_cells <- lapply(seq_len(nrow(wide_rows)), function(i) {
  pair <- wide_rows[i, ]
  cells <- long_rows %>%
    filter(cross == pair$cross, trait == pair$trait,
           chr1 == pair$chr1, pos1 == pair$pos1,
           chr2 == pair$chr2, pos2 == pair$pos2)
  vals <- setNames(rep("N/A", length(cell_order)), cell_order)
  for (j in seq_len(nrow(cells))) {
    key <- paste(cells$geno1[j], cells$geno2[j], sep = "_")
    vals[key] <- sprintf("%.2f +/- %.2f", round(cells$mean_dev[j], 2), round(cells$se[j], 2))
  }
  as.data.frame(as.list(vals), stringsAsFactors = FALSE, check.names = FALSE)
})
wide_cells <- bind_rows(wide_cells)

wide_out <- bind_cols(
  wide_rows %>% select(cross, trait, chr1, chr2, pos1, pos2),
  wide_cells,
  wide_rows %>% select(lod_full, lod_int, sig_level)
)

write.csv(wide_out, wide_output_file, row.names = FALSE)

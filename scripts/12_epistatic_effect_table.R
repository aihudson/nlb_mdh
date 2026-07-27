require(qtl)
require(dplyr)
require(ggplot2)

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

# Load every population's cross once (all three are needed for the faceted
# figures, which show each significant pair's effect across RIL / B73 BC /
# Mo17 BC even where the pair is only significant in one of them).
crosses <- lapply(names(presets), function(pop) {
  p <- presets[[pop]]
  sim.geno(read_cross(p$input_file, p$genotype, p$na.strings, p$crosstype), step = 2.5)
})
names(crosses) <- names(presets)

# Genotype-class means/SEs for one pair in one population, as B73-allele-count
# cells centered on that population's trait mean. Returns NULL if effectplot
# can't be computed (e.g. a trait that isn't meaningful in this population).
pair_effects <- function(cross, gm_map, trait, chr1, pos1, chr2, pos2) {
  gm <- mean(cross$pheno[[trait]], na.rm = TRUE)
  ep <- tryCatch(
    effectplot(cross, pheno.col = match(trait, colnames(cross$pheno)),
               mname1 = sprintf("%s@%s", chr1, pos1),
               mname2 = sprintf("%s@%s", chr2, pos2),
               geno1 = gm_map$labels, geno2 = gm_map$labels, draw = FALSE),
    error = function(e) NULL)
  if (is.null(ep)) return(NULL)
  cells <- data.frame()
  for (r in seq_along(gm_map$counts)) {
    for (c in seq_along(gm_map$counts)) {
      cells <- bind_rows(cells, data.frame(
        geno1 = gm_map$counts[r], geno2 = gm_map$counts[c],
        geno1_label = gm_map$labels[r], geno2_label = gm_map$labels[c],
        mean_dev = ep$Means[r, c] - gm,
        se = ep$SEs[r, c]
      ))
    }
  }
  cells
}

# --- data table: genotype-class effects for each pair in its significant
# population only (drives the long/wide CSV outputs) ---
long_rows <- data.frame()
for (population in unique(sig_pairs$cross)) {
  cross <- crosses[[population]]
  gm_map <- geno_map[[population]]
  pop_pairs <- sig_pairs %>% filter(cross == population)
  for (i in seq_len(nrow(pop_pairs))) {
    pair <- pop_pairs[i, ]
    cells <- pair_effects(cross, gm_map, pair$trait,
                          pair$chr1, pair$pos1, pair$chr2, pair$pos2)
    long_rows <- bind_rows(long_rows, data.frame(
      cross = population, trait = pair$trait,
      chr1 = pair$chr1, pos1 = pair$pos1,
      chr2 = pair$chr2, pos2 = pair$pos2,
      geno1 = cells$geno1, geno2 = cells$geno2,
      mean_dev = cells$mean_dev, se = cells$se,
      lod_full = pair$lod.full, lod_int = pair$lod.int,
      sig_level = pair$sig_level
    ))
  }
}

write.csv(long_rows, long_output_file, row.names = FALSE)

# --- figures: one faceted plot per pair, matching mdh_qtl/figures/epi_qtl_1.pdf
# (panels RIL / B73 BC / Mo17 BC, default grey theme, size-20 text, x-labels
# rotated 45 deg; x = second locus genotype, colour = first locus genotype,
# y = mean deviation from the population mean). Each pair's effect is computed
# in all three populations regardless of where it was significant. ---
geno_levels <- c("B73/B73", "B73/Mo17", "Mo17/Mo17")
pop_labels <- c(RIL = "RIL", B73_BC = "B73 BC", Mo17_BC = "Mo17 BC")

unique_pairs <- sig_pairs %>% distinct(cross, trait, chr1, pos1, chr2, pos2)
for (i in seq_len(nrow(unique_pairs))) {
  pair <- unique_pairs[i, ]
  fig_df <- data.frame()
  for (population in names(crosses)) {
    cells <- pair_effects(crosses[[population]], geno_map[[population]], pair$trait,
                          pair$chr1, pair$pos1, pair$chr2, pair$pos2)
    if (is.null(cells)) next
    cells$Population <- pop_labels[[population]]
    fig_df <- bind_rows(fig_df, cells)
  }
  fig_df$Population <- factor(fig_df$Population, levels = pop_labels)
  fig_df$geno1_label <- factor(fig_df$geno1_label, levels = geno_levels)
  fig_df$geno2_label <- factor(fig_df$geno2_label, levels = geno_levels)

  locus1_lab <- sprintf("Chr %s %.1f IcM", pair$chr1, pair$pos1)
  locus2_lab <- sprintf("Chr %s %.1f IcM", pair$chr2, pair$pos2)
  plot <- ggplot(fig_df, aes(x = geno2_label, y = mean_dev,
                             color = geno1_label, group = geno1_label)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = mean_dev - se, ymax = mean_dev + se), width = 0.2) +
    labs(x = locus2_lab, y = "Mean", color = locus1_lab) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_wrap(vars(Population), ncol = 3)

  fig_file <- sprintf("figures/epistasis_%s_chr%s-%.1f_chr%s-%.1f.pdf",
                       pair$cross, pair$chr1, pair$pos1, pair$chr2, pair$pos2)
  ggsave(fig_file, plot, width = 10, height = 6)
}

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

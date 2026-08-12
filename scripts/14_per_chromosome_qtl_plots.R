library(qtl)
library(dplyr)
library(ggplot2)
library(cowplot)

# Rscript scripts/14_per_chromosome_qtl_plots.R
# Per-chromosome LOD + effect figures, one PDF per chromosome carrying a QTL,
# matching the style of mdh_qtl/figures/chr3_slb_frreal.pdf but for the 3
# NLB populations (RIL, B73 BC, Mo17 BC).

slot_names <- c("ril", "b73", "mo17")
slot_labels <- c("RIL", "B73 BC", "Mo17 BC")
effect_cols <- c("NLB_ril_effect", "NLB_b73_effect", "NLB_mo17_effect")

# copied from scripts/06_identify_qtl.R:44-53 (also scripts/13_circos_plots.R:14-19)
get_thr <- function(perms){
  perms <- perms[!vapply(perms, is.null, logical(1))]
  thr <- vapply(perms, function(pp) summary(pp)[[1]], numeric(1))
  names(thr) <- vapply(perms, function(pp) colnames(pp)[1], character(1))
  thr
}

ril_rds <- readRDS("analyses/qtl_analyses/RIL.RDS")
b73_rds <- readRDS("analyses/qtl_analyses/B73_BC.RDS")
mo17_rds <- readRDS("analyses/qtl_analyses/Mo17_BC.RDS")

ril_thr <- get_thr(ril_rds$permutations)
b73_thr <- get_thr(b73_rds$permutations)
mo17_thr <- get_thr(mo17_rds$permutations)

scans <- list(
  ril  = list(scan = as.data.frame(ril_rds$scan),       lodcol = "LOD NLB_WMD_BLUP"),
  b73  = list(scan = as.data.frame(b73_rds$scan[[1]]),  lodcol = "LOD NLB_WMD_BLUP"),
  mo17 = list(scan = as.data.frame(mo17_rds$scan[[1]]), lodcol = "LOD NLB_WMD_BLUP")
)

thresholds <- list(
  ril  = ril_thr[["LOD NLB_WMD_BLUP"]],
  b73  = b73_thr[["LOD NLB_WMD_BLUP"]],
  mo17 = mo17_thr[["LOD NLB_WMD_BLUP"]]
)

effect_grid <- read.csv("analyses/qtl_effects_whole_genome.csv", stringsAsFactors = FALSE) %>%
  arrange(chr, pos)

main_effect_peaks <- read.csv("analyses/main_effect_peaks.csv", stringsAsFactors = FALSE)
main_effect_peaks$trait <- sub("^LOD ", "", main_effect_peaks$trait)
peaks <- main_effect_peaks %>%
  transmute(chr = as.numeric(chr), pos = pos) %>%
  distinct(chr, pos)

# copied from scripts/13_circos_plots.R:64-91, restricted to the 3 NLB slots;
# interpolates every population's LOD curve onto the effect CSV's (chr, pos)
# grid so LOD and effect rows line up 1:1
build_lod_qtl <- function(effect_grid, scans, thresholds){
  raw <- list()
  for (nm in slot_names){
    slot <- scans[[nm]]
    scan_chr <- as.numeric(as.character(slot$scan$chr))
    scan_pos <- slot$scan[["pos (cM)"]]
    scan_lod <- slot$scan[[slot$lodcol]]
    vals <- rep(NA_real_, nrow(effect_grid))
    for (ch in unique(effect_grid$chr)){
      grid_idx <- which(effect_grid$chr == ch)
      scan_idx <- which(scan_chr == ch)
      vals[grid_idx] <- approx(scan_pos[scan_idx], scan_lod[scan_idx],
                                xout = effect_grid$pos[grid_idx], rule = 2)$y
    }
    raw[[nm]] <- vals
  }

  lod_qtl <- data.frame(chr = effect_grid$chr, pos = effect_grid$pos)
  for (nm in slot_names){
    lod_qtl[[paste0(nm, "_lod")]] <- ifelse(raw[[nm]] < thresholds[[nm]], raw[[nm]], NA)
  }
  for (nm in slot_names){
    lod_qtl[[paste0(nm, "_lod_sig")]] <- ifelse(raw[[nm]] >= thresholds[[nm]], raw[[nm]], NA)
  }
  lod_qtl
}

lod_qtl <- build_lod_qtl(effect_grid, scans, thresholds)

# tidy long frame: two measurement layers (Effect, LOD) per slot, aligned 1:1
# with effect_grid/lod_qtl row order
long_list <- list()
for (i in seq_along(slot_names)){
  nm <- slot_names[i]
  lod_raw <- coalesce(lod_qtl[[paste0(nm, "_lod")]], lod_qtl[[paste0(nm, "_lod_sig")]])
  long_list[[paste0(nm, "_eff")]] <- data.frame(
    chr = effect_grid$chr, pos = effect_grid$pos,
    Cross = slot_labels[i], measurement = "Effect", value = effect_grid[[effect_cols[i]]]
  )
  long_list[[paste0(nm, "_lod")]] <- data.frame(
    chr = effect_grid$chr, pos = effect_grid$pos,
    Cross = slot_labels[i], measurement = "LOD", value = lod_raw
  )
}
long_df <- bind_rows(long_list)
long_df$Cross <- factor(long_df$Cross, levels = slot_labels)

thr_df <- data.frame(
  Cross = factor(slot_labels, levels = slot_labels),
  yintercept = c(thresholds$ril, thresholds$b73, thresholds$mo17)
)

shared_theme <- theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = "black"),
  text = element_text(size = 20),
  axis.title.y = element_blank(),
  strip.placement = "outside"
)

for (ch in sort(unique(peaks$chr))){
  eff_df <- long_df %>% filter(chr == ch, measurement == "Effect")
  lod_df <- long_df %>% filter(chr == ch, measurement == "LOD")
  peaks_ch <- peaks %>% filter(chr == ch)

  eff_lim <- max(abs(eff_df$value), na.rm = TRUE) * 1.1
  lod_lim <- ceiling(max(c(lod_df$value, thr_df$yintercept), na.rm = TRUE) * 1.05)

  p_eff <- ggplot(eff_df, aes(x = pos, y = value)) +
    geom_line() +
    facet_grid(Cross ~ ., switch = "y") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(data = peaks_ch, aes(xintercept = pos),
               color = "red", linetype = "dotted", linewidth = 1) +
    ylim(-eff_lim, eff_lim) +
    labs(x = "Pos (cM)", y = NULL, title = "Effect") +
    shared_theme

  p_lod <- ggplot(lod_df, aes(x = pos, y = value)) +
    geom_line() +
    facet_grid(Cross ~ ., switch = "y") +
    geom_hline(data = thr_df, aes(yintercept = yintercept), linewidth = 1) +
    geom_vline(data = peaks_ch, aes(xintercept = pos),
               color = "red", linetype = "dotted", linewidth = 1) +
    ylim(0, lod_lim) +
    labs(x = "Pos (cM)", y = NULL, title = "LOD") +
    shared_theme

  combined <- plot_grid(p_eff, p_lod, ncol = 2)
  # the Effect (left) panel blanks its y-title in shared_theme; add it here on
  # the assembled grid so the effect direction is stated
  combined <- ggdraw(combined) +
    draw_label("B73-allele effect on resistance (+ = more resistant)",
               x = 0, y = 0.5, vjust = 1.5, angle = 90, size = 12)
  ggsave(sprintf("figures/chr%s_nlb.pdf", ch), combined, width = 8, height = 10)
}

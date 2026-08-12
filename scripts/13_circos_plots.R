library(circlize)
library(qtl)
library(dplyr)

# Rscript scripts/13_circos_plots.R

# five analysis slots, published order: RIL, B73 BC, B73 MPH, Mo17 BC, Mo17 MPH
slot_names <- c("ril", "b73", "b73_mph", "mo17", "mo17_mph")
slot_labels <- c("RIL", "B73 BC", "B73 MPH", "Mo17 BC", "Mo17 MPH")
effect_cols <- c("NLB_ril_effect", "NLB_b73_effect", "NLB_b73_mph_effect",
                  "NLB_mo17_effect", "NLB_mo17_mph_effect")

# copied from scripts/06_identify_qtl.R:44-53
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
  ril      = list(scan = as.data.frame(ril_rds$scan),       lodcol = "LOD NLB_WMD_BLUP"),
  b73      = list(scan = as.data.frame(b73_rds$scan[[1]]),  lodcol = "LOD NLB_WMD_BLUP"),
  b73_mph  = list(scan = as.data.frame(b73_rds$scan[[2]]),  lodcol = "LOD NLB_WMD_BLUP_MPH"),
  mo17     = list(scan = as.data.frame(mo17_rds$scan[[1]]), lodcol = "LOD NLB_WMD_BLUP"),
  mo17_mph = list(scan = as.data.frame(mo17_rds$scan[[2]]), lodcol = "LOD NLB_WMD_BLUP_MPH")
)

thresholds <- list(
  ril      = ril_thr[["LOD NLB_WMD_BLUP"]],
  b73      = b73_thr[["LOD NLB_WMD_BLUP"]],
  b73_mph  = b73_thr[["LOD NLB_WMD_BLUP_MPH"]],
  mo17     = mo17_thr[["LOD NLB_WMD_BLUP"]],
  mo17_mph = mo17_thr[["LOD NLB_WMD_BLUP_MPH"]]
)

effect_grid <- read.csv("analyses/qtl_effects_whole_genome.csv", stringsAsFactors = FALSE) %>%
  arrange(chr, pos)

main_effect_peaks <- read.csv("analyses/main_effect_peaks.csv", stringsAsFactors = FALSE)
main_effect_peaks$trait <- sub("^LOD ", "", main_effect_peaks$trait)
combined_qtl <- main_effect_peaks %>%
  distinct(chr, pos) %>%
  transmute(chr = paste0("chr", chr), start = pos, peak = TRUE)

IcM_bed <- effect_grid %>%
  group_by(chr) %>%
  summarize(end = max(pos), .groups = "drop") %>%
  mutate(chr = paste0("chr", chr), start = 0) %>%
  select(chr, start, end) %>%
  as.data.frame()

# interpolate every population's LOD curve onto the effect CSV's (chr, pos)
# grid so LOD and effect rows line up 1:1 (the LOD scans are on the MQM
# pseudomarker grid; qtl_effects_whole_genome.csv is on a step-2.5 cM grid)
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

  lod_qtl <- data.frame(chr = paste0("chr", effect_grid$chr),
                         start = effect_grid$pos,
                         end = effect_grid$pos + 1)
  for (nm in slot_names){
    lod_qtl[[paste0(nm, "_lod")]] <- ifelse(raw[[nm]] < thresholds[[nm]], raw[[nm]], NA)
  }
  for (nm in slot_names){
    lod_qtl[[paste0(nm, "_lod_sig")]] <- ifelse(raw[[nm]] >= thresholds[[nm]], raw[[nm]], NA)
  }
  lod_qtl
}

# adapted from scripts/get_effects/circos_plot.R:213-310
make_lod_circos <- function(lod_qtl, combined_qtl){
  lod_qtl <- merge(lod_qtl, combined_qtl, by = c("chr", "start"), all = TRUE)
  lod_qtl$end <- lod_qtl$start + 1
  ymax <- lod_qtl %>% select(ends_with("_lod_sig")) %>% unlist() %>% max(na.rm = TRUE) %>% ceiling()
  lod_qtl$peak[which(lod_qtl$peak == TRUE)] <- ymax

  circos.par("start.degree" = 70)
  circos.par("gap.degree" = c(rep(1, 9), 40))
  par(cex = 1)
  circos.initializeWithIdeogram(IcM_bed, plotType = c("axis", "labels"))

  invisible(lapply(1:5, function(i){
    circos.genomicTrackPlotRegion(lod_qtl, ylim = c(0, ymax), track.height = 0.1,
      panel.fun = function(region, value, ...){
        circos.genomicLines(region, value, numeric.column = i, col = "black")
        circos.genomicLines(region, value, numeric.column = 5 + i, col = "red",
                             type = "o", cex = 0.20, pch = 16, pt.col = "red")
        circos.genomicLines(region, value, numeric.column = 11, col = "red",
                             type = "h", lty = 3, ylim = c(0, ymax))
      })
    circos.text(x = 1, y = -1, adj = c(degree(317.5), degree(0)),
                labels = slot_labels[i], font = 2, cex = 0.5)
  }))
}

# adapted from scripts/get_effects/effect_circos_plot.R:187-321
make_effect_circos <- function(lod_qtl, effect_df, combined_qtl){
  lod_cols <- paste0(slot_names, "_lod")
  sig_cols <- paste0(slot_names, "_lod_sig")

  effect_df2 <- data.frame(chr = lod_qtl$chr, start = lod_qtl$start, end = lod_qtl$end,
                            effect_df[, effect_cols])
  effect_df2 <- cbind(effect_df2,
                       setNames(effect_df2[, effect_cols], paste0(effect_cols, "no2")))

  for (i in 1:5){
    effect_df2[is.na(lod_qtl[[lod_cols[i]]]), effect_cols[i]] <- NA
    effect_df2[is.na(lod_qtl[[sig_cols[i]]]), paste0(effect_cols[i], "no2")] <- NA
  }

  ymax <- effect_df2 %>% select(all_of(c(effect_cols, paste0(effect_cols, "no2")))) %>%
    unlist() %>% abs() %>% max(na.rm = TRUE) * 1.1

  effect_df2 <- merge(effect_df2, combined_qtl, by = c("chr", "start"), all = TRUE)
  effect_df2$end <- effect_df2$start + 1
  effect_df2$peak[which(effect_df2$peak == TRUE)] <- ymax

  circos.par("start.degree" = 70)
  circos.par("gap.degree" = c(rep(1, 9), 40))
  par(cex = 1)
  circos.initializeWithIdeogram(IcM_bed, plotType = c("axis", "labels"))

  invisible(lapply(1:5, function(i){
    circos.genomicTrackPlotRegion(effect_df2, ylim = c(-ymax, ymax), track.height = 0.1,
      panel.fun = function(region, value, ...){
        circos.genomicLines(region, value, numeric.column = i, col = "black")
        circos.genomicLines(region, value, numeric.column = 5 + i, col = "red",
                             type = "o", cex = 0.20, pch = 16, pt.col = "red")
        circos.genomicLines(region, value, numeric.column = 11, col = "red",
                             type = "h", lty = 3, ylim = c(-ymax, ymax))
        circos.segments(x0 = 0, y0 = 0, x1 = CELL_META$xlim, y1 = 0, lty = 2)
      })
    circos.text(x = 1, y = 0, adj = c(degree(317.5), degree(0)),
                labels = slot_labels[i], font = 2, cex = 0.5)
  }))
}

lod_qtl <- build_lod_qtl(effect_grid, scans, thresholds)

pdf("figures/nlb_circos_lod.pdf", width = 7, height = 7)
make_lod_circos(lod_qtl, combined_qtl)
circos.clear()
dev.off()

pdf("figures/nlb_circos_effect.pdf", width = 7, height = 7)
make_effect_circos(lod_qtl, effect_grid, combined_qtl)
mtext("B73-allele effect on resistance (+ = more resistant)",
      side = 1, line = 0, cex = 0.9)
circos.clear()
dev.off()

pdf("figures/nlb_circos_paper.pdf", width = 7, height = 13)
layout(matrix(c(1, 2, 3), ncol = 1), heights = c(1, 1, 0.12))

par(mar = c(1, 1, 2, 1))
make_lod_circos(lod_qtl, combined_qtl)
mtext("A", side = 3, line = 0, adj = 0, cex = 1.3, font = 2)
circos.clear()

par(mar = c(1, 1, 2, 1))
make_effect_circos(lod_qtl, effect_grid, combined_qtl)
mtext("B", side = 3, line = 0, adj = 0, cex = 1.3, font = 2)
mtext("B73-allele effect on resistance (+ = more resistant)",
      side = 1, line = 0, cex = 0.8)
circos.clear()

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", ncol = 2, bty = "n", cex = 0.9,
       legend = c("Significant LOD / effect", "Below threshold",
                  "QTL peak", "Zero effect"),
       col = c("red", "black", "red", "grey40"),
       lty = c(1, 1, 3, 2), pch = c(16, NA, NA, NA))
dev.off()

#!/usr/bin/env Rscript
# 15_compare_old_study_qtl.R
#
# Compare the NLB QTL identified in this project (analyses/main_effect_peaks.csv)
# against the QTL from the earlier published IBM B73xMo17 study (its Table 3),
# and report which old QTL overlap this project's peaks.
#
# No population argument; run from the project root:
#   Rscript scripts/15_compare_old_study_qtl.R
#
# --- Why anchor via flanking markers, not raw positions --------------------
# Both studies sit on the IBM2 genetic map and share named markers, so they ARE
# comparable. BUT the absolute Imu position numbers are offset between the two
# map builds: e.g. the old AU06WMD chr1 QTL is reported at Imu interval 534-546,
# yet its flanking markers bnlg1598 / umc1396 sit at 484.6 / 500.5 in this
# project's data/ibm302map.csv (an interval cannot lie outside its own flanking
# markers on the same map). So we place each old QTL onto THIS project's map by
# looking up its two flanking markers in ibm302map.csv, and compare those anchored
# intervals against this project's peak support intervals.
#
# Trait scope: WMD is the primary (same trait mapped in this project). Old
# incubation-period (IP) QTL are carried through and flagged as a secondary
# cross-trait observation. DTA is encoded for reference but excluded from the
# overlap comparison.

NEAR_CM <- 20  # a non-overlapping old QTL within this many cM of a peak = "near_miss"
              # (matches the min_sep_cM = 20 linked-locus threshold used in script 09)

# --- 1. Old published Table 3, hand-transcribed from the PDF ----------------
# imu_lo/imu_hi are the OLD study's raw 2-LOD interval (reference only; NOT used
# for overlap because of the map-build offset described above). flank_lo/flank_hi
# are the reported flanking markers, used to anchor onto this project's map.
# (hda.03 in the PDF is marker hda103 in ibm302map.csv.)
old <- read.csv(text = "
env,measure,chr,bin,imu_lo,imu_hi,lod,r2,a,flank_lo,flank_hi
BLUP,WMD,2,2.00-2.01,0,21.8,4.1,6.7,2.28,isu053a,isu144a
BLUP,WMD,4,4.08,449.2,456.9,3.12,4.3,-1.23,umc1775,mmp3
AU06,WMD,1,1.06,534,546,7.26,10.2,4.2,bnlg1598,umc1396
AU06,WMD,2,2.00-2.01,0,23.8,3.45,4.6,2.79,isu053a,isu144a
AU06,WMD,3,3.05,312.5,318.2,3.78,4.9,2.92,rz296b,mmc0022
AU06,WMD,3,3.09,727.2,732.5,3.98,5.3,-2.94,csu845,sho89
AU06,WMD,4,4.07-4.08,437.5,449.4,5.56,7.6,-3.51,asg33,umc1667
AU07,WMD,2,2.00-2.01,0,27.4,4.55,9.3,3.99,isu053a,isu144a
AU07,WMD,4,4.07-4.08,438.5,498.7,2.8,4.5,-2.6,asg33,ufg23
CL07,WMD,4,4.08,449.2,455.9,3.03,4.9,-1.06,umc1775,mmp3
BLUP,IP,2,2.02,148.1,153.1,3.81,4.8,-0.42,bnlg2277,umc1262
BLUP,IP,4,4.05,277.3,283.3,3,3.8,0.38,umc1303,mmp125
BLUP,IP,6,6.05,299.2,325.9,3.29,3.9,-0.38,npi560,umc1020
BLUP,IP,6,6.07,480.5,502.9,3.06,3.9,0.39,umc1490,umc2165
BLUP,IP,8,8.05,374.5,379.2,4.32,5.3,-0.47,bnlg1651,hda103
BLUP,IP,8,8.07,464,500.7,3.23,3.9,-0.37,rz538a,lim301
AU06,IP,2,2.02,148.1,153.1,3.14,4.9,-1.29,bnlg2277,umc1262
AU06,IP,4,4.06,353.8,367,4.55,6.9,1.55,umc2027,zm1
AU06,IP,9,9.02,105.1,128.6,3.71,5.8,1.5,umc1170,umc1636
AU06,IP,9,9.04,290,300,3.74,5.7,-1.4,bnlg1159b,ufg70
AU07,IP,2,2.01,33.4,47.4,4.16,6.6,-1.15,isu144a,umc1165
AU07,IP,2,2.02-2.03,164.6,179.4,4.41,6.7,-1.17,umc1422,mmc0231
AU07,IP,8,8.08,544.3,562.5,3.42,5.1,-1,umc1933,npi107
CL07,IP,2,2.02,65.7,116.3,3.22,4.1,-0.37,bnlg1017,eks1
CL07,IP,6,6.05,314,320.7,4.32,5.8,-0.44,npi252,bnlg1702
BLUP,DTA,2,2.02,149,154,3.81,4.5,0.41,bnlg1327,umc1261
BLUP,DTA,4,4.09,565,579.8,3.89,4.7,-0.45,umc52,php10025
BLUP,DTA,8,8.05,358,367,8.43,10.4,0.7,psb107a,bnlg666
BLUP,DTA,9,9.02,138.1,147.5,4.46,5.3,-0.46,umc1636,bnlg1401
", stringsAsFactors = FALSE, strip.white = TRUE)

# --- 2. Anchor each old QTL onto this project's map via flanking markers -----
map <- read.csv("data/ibm302map.csv", header = FALSE,
                stringsAsFactors = FALSE)   # V1 marker, V2 chr, V3 position
names(map) <- c("marker", "chr", "pos")
marker_pos <- setNames(map$pos, map$marker)

lo_pos <- marker_pos[old$flank_lo]
hi_pos <- marker_pos[old$flank_hi]

missing <- unique(c(old$flank_lo[is.na(lo_pos)], old$flank_hi[is.na(hi_pos)]))
if (length(missing)) {
  warning("Flanking markers not found in ibm302map.csv (anchor set NA): ",
          paste(missing, collapse = ", "))
}

old$anchor_lo <- pmin(lo_pos, hi_pos)
old$anchor_hi <- pmax(lo_pos, hi_pos)

# --- 3. This project's peaks -------------------------------------------------
peaks <- read.csv("analyses/main_effect_peaks.csv", stringsAsFactors = FALSE)
peaks <- peaks[, c("cross", "trait", "chr", "pos", "lod", "ci_lo", "ci_hi")]
names(peaks) <- c("new_cross", "new_trait", "chr", "new_pos", "new_lod",
                  "new_ci_lo", "new_ci_hi")
peaks$chr <- as.integer(peaks$chr)

# --- 4. Overlap test: WMD (primary) + IP (secondary), on shared chromosomes --
cmp <- old[old$measure %in% c("WMD", "IP"), ]

# All same-chromosome (new peak, old QTL) pairs.
pairs <- merge(peaks, cmp, by = "chr", suffixes = c("", "_old"))

pairs$interval_overlap <- with(pairs, anchor_lo <= new_ci_hi & anchor_hi >= new_ci_lo)
pairs$gap_cM <- with(pairs, ifelse(interval_overlap, 0,
                                   pmax(anchor_lo - new_ci_hi, new_ci_lo - anchor_hi)))
pairs$same_measure <- pairs$measure == "WMD"   # this project maps WMD

classify <- function(overlap, gap, measure) {
  if (isTRUE(overlap) && measure == "WMD") return("recovered_WMD")
  if (isTRUE(overlap) && measure == "IP")  return("cross_trait_IP")
  if (!overlap && gap <= NEAR_CM && measure == "WMD") return("near_miss_WMD")
  if (!overlap && gap <= NEAR_CM && measure == "IP")  return("near_miss_IP")
  "same_chr_no_overlap"
}
pairs$overlap_call <- mapply(classify, pairs$interval_overlap, pairs$gap_cM, pairs$measure)

# New peaks with no old WMD/IP QTL anywhere on their chromosome.
matched_chr <- unique(pairs[, c("new_cross", "new_pos", "chr")])
peaks_nomatch <- peaks[!interaction(peaks$new_cross, peaks$new_pos, peaks$chr) %in%
                         interaction(matched_chr$new_cross, matched_chr$new_pos, matched_chr$chr), ]

# Old WMD/IP QTL whose chromosome carries no new peak (i.e. not recovered here).
old_unmatched <- cmp[!cmp$chr %in% peaks$chr, ]

# --- 5. Assemble tidy output -------------------------------------------------
out_cols <- c("new_cross", "new_trait", "chr", "new_pos", "new_lod",
              "new_ci_lo", "new_ci_hi",
              "env", "measure", "bin", "imu_lo", "imu_hi", "lod", "r2", "a",
              "flank_lo", "flank_hi", "anchor_lo", "anchor_hi",
              "interval_overlap", "gap_cM", "same_measure", "overlap_call")

pair_rows <- pairs[, out_cols]

blank_old <- function(df, call) {
  data.frame(new_cross = df$new_cross, new_trait = df$new_trait, chr = df$chr,
             new_pos = df$new_pos, new_lod = df$new_lod,
             new_ci_lo = df$new_ci_lo, new_ci_hi = df$new_ci_hi,
             env = NA, measure = NA, bin = NA, imu_lo = NA, imu_hi = NA,
             lod = NA, r2 = NA, a = NA, flank_lo = NA, flank_hi = NA,
             anchor_lo = NA, anchor_hi = NA, interval_overlap = NA, gap_cM = NA,
             same_measure = NA, overlap_call = call, stringsAsFactors = FALSE)
}
blank_new <- function(df, call) {
  data.frame(new_cross = NA, new_trait = NA, chr = df$chr, new_pos = NA,
             new_lod = NA, new_ci_lo = NA, new_ci_hi = NA,
             env = df$env, measure = df$measure, bin = df$bin,
             imu_lo = df$imu_lo, imu_hi = df$imu_hi, lod = df$lod, r2 = df$r2,
             a = df$a, flank_lo = df$flank_lo, flank_hi = df$flank_hi,
             anchor_lo = df$anchor_lo, anchor_hi = df$anchor_hi,
             interval_overlap = NA, gap_cM = NA, same_measure = df$measure == "WMD",
             overlap_call = call, stringsAsFactors = FALSE)
}

result <- rbind(
  pair_rows,
  if (nrow(peaks_nomatch)) blank_old(peaks_nomatch, "novel_no_old_on_chr"),
  if (nrow(old_unmatched))  blank_new(old_unmatched, "not_recovered")
)
result <- result[order(result$chr, result$new_pos, result$measure), ]

dir.create("analyses", showWarnings = FALSE)
write.csv(result, "analyses/qtl_overlap_old_study.csv", row.names = FALSE)

# --- 6. Concise stdout summary ----------------------------------------------
cat("\n=== Old-study QTL overlap with this project's peaks ===\n")
cat("(old QTL anchored onto this project's IBM map via shared flanking markers)\n\n")

new_peaks <- unique(peaks[, c("new_cross", "chr", "new_pos", "new_trait")])
new_peaks <- new_peaks[order(new_peaks$chr, new_peaks$new_pos), ]
cat("-- Per project peak --\n")
for (i in seq_len(nrow(new_peaks))) {
  p <- new_peaks[i, ]
  sub <- pairs[pairs$new_cross == p$new_cross & pairs$chr == p$chr &
                 pairs$new_pos == p$new_pos, ]
  ip_calls <- c("cross_trait_IP", "near_miss_IP")
  best <- if (any(sub$overlap_call == "recovered_WMD")) "RECOVERED old WMD QTL" else
          if (any(sub$overlap_call == "near_miss_WMD"))  "near an old WMD QTL" else
          if (any(sub$overlap_call %in% ip_calls)) "near an old IP QTL (cross-trait)" else
          "NOVEL (no old WMD/IP QTL here)"
  hit_calls <- c("recovered_WMD", "near_miss_WMD", ip_calls)
  hits <- sub[sub$overlap_call %in% hit_calls, ]
  detail <- if (nrow(hits)) paste0(" [", paste(unique(paste0(hits$env, hits$measure,
              " ", hits$flank_lo, "-", hits$flank_hi)), collapse = "; "), "]") else ""
  cat(sprintf("  chr%s @ %g cM (%s, %s): %s%s\n", p$chr, p$new_pos, p$new_cross,
              p$new_trait, best, detail))
}

cat("\n-- Old WMD QTL not recovered (no project peak on that chromosome / interval) --\n")
old_wmd <- old[old$measure == "WMD", ]
for (i in seq_len(nrow(old_wmd))) {
  q <- old_wmd[i, ]
  recovered <- any(pairs$measure == "WMD" & pairs$chr == q$chr &
                     pairs$interval_overlap &
                     pairs$flank_lo == q$flank_lo & pairs$flank_hi == q$flank_hi)
  if (!recovered) {
    cat(sprintf("  chr%s bin %s (%sWMD, %s-%s, anchored %g-%g cM)\n", q$chr, q$bin,
                q$env, q$flank_lo, q$flank_hi, q$anchor_lo, q$anchor_hi))
  }
}

cat("\nWrote analyses/qtl_overlap_old_study.csv\n")

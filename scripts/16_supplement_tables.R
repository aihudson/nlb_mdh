# 16_supplement_tables.R
# ---------------------------------------------------------------------------
# Regenerates the NLB supplemental-tables Word document
# (manuscript/nlb_mdh_supplement.docx), renumbered to mirror the sister
# SLB/GLS supplement (slb_gls_ms/slb_gls_heterosis_supp_figures_tables.docx).
#
# Table map (NLB label -> content):
#   S1  within-year between-rep correlations of WMD           (from 02_get_blups.R:256-299)
#   S2  between-scoring-date correlations of raw disease      (NEW; date columns from data/nlb_mdh_file_s1.xlsx)
#   S3  cross-environment correlations of WMD BLUPs           (from 02_get_blups.R:301-359)
#   ( S4 inter-trait correlations - N/A, NLB is single-trait; number intentionally skipped )
#   S5  broad-sense heritability per population               (from 02_get_blups.R:242-254)
#   S6  MPH/BPH of the B73 x Mo17 reference hybrid            (IBM_NLB_BLUPs.csv + emmeans sig)
#   S7  MPH/BPH of NAM diverse-line testcrosses              (NAM_NLB_BLUPs.csv)
#   S8  mean MPH/BPH of the backcross populations            (IBM_NLB_BLUPs_{B73BC,Mo17BC}.csv)
#   S9  main-effect + epistatic-pair-member QTL              (main_effect_peaks.csv + qtl_effects.csv)
#   S10 gene-action (d/a) classification of QTL              (qtl_gene_action.csv)
#   S11 correlation with the earlier published IBM NLB study (data/old_nlb_data.csv)
#
# Run from the project root:  Rscript scripts/16_supplement_tables.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(lme4)
  library(emmeans)
  library(officer)
  library(flextable)
})

out_docx <- "manuscript/nlb_mdh_supplement.docx"
xlsx     <- "~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx"

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
# Pearson correlation that degrades to NA when a population lacks enough data
# (e.g. the 2024_IL trial carried only backcross hybrids, no pure RILs)
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok]))
  c(r = unname(ct$estimate), p = ct$p.value)
}

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, "**", ifelse(p < 0.05, "*", "")))
fmt_r    <- function(r) sprintf("%.2f", r)
# correlation cell: plain r when the whole table is significant at P<1e-4
# (matching the SLB blanket caption), otherwise r with significance stars
cor_cell <- function(r, p, all_sig)
  ifelse(is.na(r), "NA", paste0(sprintf("%.2f", r), if (!all_sig) stars(p) else ""))
fmt_pct  <- function(x, p = NA) ifelse(is.na(x), "NA", paste0(sprintf("%.2f", x * 100), stars(p)))
fmt_num  <- function(x, d = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmt_p    <- function(p) ifelse(p < 2.2e-16, "< 2.2e-16", formatC(p, format = "e", digits = 2))

# flextable styled to match the SLB/GLS supplement (Calibri, bold header, header rule)
style_ft <- function(df, widths = NULL) {
  ft <- flextable(df)
  ft <- font(ft, fontname = "Calibri", part = "all")
  ft <- fontsize(ft, size = 11, part = "all")
  ft <- bold(ft, part = "header")
  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header", border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "header", border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "body", border = fp_border(color = "black", width = 1))
  ft <- align(ft, part = "all", align = "left")
  ft <- padding(ft, padding = 2, part = "all")
  ft <- set_table_properties(ft, layout = "autofit", width = 1)
  ft
}

doc <- read_docx()
add_caption <- function(doc, txt)
  body_add_fpar(doc, fpar(ftext(txt, prop = fp_text(bold = TRUE, font.size = 11, font.family = "Calibri"))))
add_note <- function(doc, txt)
  body_add_fpar(doc, fpar(ftext(txt, prop = fp_text(font.size = 9, font.family = "Calibri"))))
add_table <- function(doc, cap, df, notes = character(0)) {
  doc <- add_caption(doc, cap)
  doc <- body_add_flextable(doc, style_ft(df))
  for (n in notes) doc <- add_note(doc, n)
  doc <- body_add_par(doc, "", style = "Normal")
  doc
}

envs <- c("2022_CL", "2023_CL", "2024_CL", "2024_IL")

# ---------------------------------------------------------------------------
# load / rebuild phenotype data (mirrors 02_get_blups.R:8-46)
# ---------------------------------------------------------------------------
# normalize Line names the way 02_get_blups.R:31-36 does, so population filters work
# on every sheet (notably 2024_IL, where RILs are recorded as selfs "M0###xM0###")
normalize_line <- function(l) {
  l[l == "Mo17xM0230"] <- "M0230xMo17"
  l[l == "Mo17xM0318"] <- "M0318xMo17"
  gsub("M0(\\d+x)M0(\\d+)", "M0\\2", l)
}
sheets <- lapply(1:4, function(i) { s <- read.xlsx(xlsx, sheet = i); s$Line <- normalize_line(s$Line); s })
d1 <- sheets[[1]] %>% mutate(wmd = NLB_0728)
data <- do.call(rbind, lapply(list(d1, sheets[[2]], sheets[[3]], sheets[[4]]),
                              function(x) dplyr::select(x, Row, Line, Loc, Year, Rep, wmd)))
data$Year <- as.factor(data$Year); data$Rep <- as.factor(data$Rep); data$Loc <- as.factor(data$Loc)
data$YearLoc <- paste(data$Year, data$Loc, sep = "_")
data$Line[data$Line == "Mo17xM0230"] <- "M0230xMo17"
data$Line[data$Line == "Mo17xM0318"] <- "M0318xMo17"
data <- data %>% mutate(Line = gsub("M0(\\d+x)M0(\\d+)", "M0\\2", Line))

data_ibm <- data %>%
  filter(Line %in% c("B73", "Mo17", "B73xMo17") | grepl("M0", Line)) %>%
  mutate(wmd = 100 - wmd)

# population filters (as in 02_get_blups.R)
is_ril  <- function(l) grepl("M0", l) & !grepl("x", l)
is_b73  <- function(l) grepl("B73", l) & grepl("M0", l)
is_mo17 <- function(l) grepl("Mo17", l) & grepl("M0", l)

# ---------------------------------------------------------------------------
# S1  within-year between-rep correlations of WMD
# ---------------------------------------------------------------------------
rep_cor <- function(df) {
  w <- df %>% dplyr::select(YearLoc, Line, Rep, wmd) %>%
    pivot_wider(names_from = Rep, values_from = wmd)
  sapply(envs, function(e) {
    s <- filter(w, YearLoc == e)
    safe_cor(s$`1`, s$`2`)
  })
}
c_ril  <- rep_cor(filter(data, is_ril(Line)))
c_b73  <- rep_cor(filter(data, is_b73(Line)))
c_mo17 <- rep_cor(filter(data, is_mo17(Line)))

s1_p <- c(c_ril["p", ], c_b73["p", ], c_mo17["p", ])
s1_all_sig <- all(s1_p < 1e-4, na.rm = TRUE)
s1 <- data.frame(Environment = envs,
                 RIL = cor_cell(c_ril["r", ], c_ril["p", ], s1_all_sig),
                 `B73 BC` = cor_cell(c_b73["r", ], c_b73["p", ], s1_all_sig),
                 `Mo17 BC` = cor_cell(c_mo17["r", ], c_mo17["p", ], s1_all_sig),
                 check.names = FALSE)

# ---------------------------------------------------------------------------
# S2  between-scoring-date correlations (NEW)
# ---------------------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))
# per-environment date matrices (columns = timepoints in order)
dates_2023 <- transform(sheets[[2]], t1 = num(NLB_0717), t2 = num(NLB_0725), t3 = num(NLB_0801))
dates_2024cl <- transform(sheets[[3]], t1 = num(NLB_0717), t2 = num(NLB_0728))
il <- sheets[[4]]
dates_2024il <- transform(il,
                          t1 = num(NLB_0724),
                          t2 = ifelse(!is.na(num(NLB_0801)), num(NLB_0801), num(NLB_0802)),
                          t3 = ifelse(!is.na(num(NLB_0808)), num(NLB_0808), num(NLB_0809)))

date_pairs <- list(
  list(env = "2023_CL", df = dates_2023, pairs = list(c("t1", "t2"), c("t1", "t3"), c("t2", "t3"))),
  list(env = "2024_CL", df = dates_2024cl, pairs = list(c("t1", "t2"))),
  list(env = "2024_IL", df = dates_2024il, pairs = list(c("t1", "t2"), c("t1", "t3"), c("t2", "t3"))))
pair_lab <- c("t1t2" = "Date 1 & 2", "t1t3" = "Date 1 & 3", "t2t3" = "Date 2 & 3")

s2_rows <- list(); s2_p <- c()
for (dp in date_pairs) {
  first <- TRUE
  for (pr in dp$pairs) {
    rp <- sapply(list(RIL = is_ril, `B73 BC` = is_b73, `Mo17 BC` = is_mo17), function(f) {
      s <- dp$df[f(dp$df$Line), ]
      safe_cor(s[[pr[1]]], s[[pr[2]]])
    })
    s2_p <- c(s2_p, rp["p", ])
    s2_rows[[length(s2_rows) + 1]] <- list(
      env = ifelse(first, dp$env, ""), dates = unname(pair_lab[paste0(pr[1], pr[2])]),
      r = rp["r", ], p = rp["p", ])
    first <- FALSE
  }
}
s2_all_sig <- all(s2_p < 1e-4, na.rm = TRUE)
s2 <- do.call(rbind, lapply(s2_rows, function(x) data.frame(
  Environment = x$env, Dates = x$dates,
  RIL = cor_cell(x$r["RIL"], x$p["RIL"], s2_all_sig),
  `B73 BC` = cor_cell(x$r["B73 BC"], x$p["B73 BC"], s2_all_sig),
  `Mo17 BC` = cor_cell(x$r["Mo17 BC"], x$p["Mo17 BC"], s2_all_sig),
  check.names = FALSE, row.names = NULL)))

# ---------------------------------------------------------------------------
# S3  cross-environment correlations of WMD BLUPs
# ---------------------------------------------------------------------------
env_blups <- lapply(envs, function(e) {
  s <- filter(data_ibm, YearLoc == e)
  m <- lmer(wmd ~ (1 | Line) + Rep, data = s)
  b <- ranef(m)$Line; b$Line <- rownames(b)
  b$blup <- b[["(Intercept)"]] + summary(m)$coefficients[1]
  b$YearLoc <- e
  b[, c("Line", "blup", "YearLoc")]
})
wide <- do.call(rbind, env_blups) %>% pivot_wider(names_from = YearLoc, values_from = blup)
cor_pop <- function(sub) {
  m <- as.matrix(sub[, envs]); r <- cor(m, use = "pairwise.complete.obs")
  p <- outer(envs, envs, Vectorize(function(a, b)
    if (a == b) 0 else safe_cor(m[, a], m[, b])["p"]))
  dimnames(p) <- list(envs, envs); list(r = r, p = p)
}
r_ril  <- cor_pop(filter(wide, is_ril(Line)))
r_b73  <- cor_pop(filter(wide, is_b73(Line)))
r_mo17 <- cor_pop(filter(wide, is_mo17(Line)))
pairs6 <- combn(envs, 2, simplify = FALSE)
s3_p <- unlist(lapply(pairs6, function(pr)
  c(r_ril$p[pr[1], pr[2]], r_b73$p[pr[1], pr[2]], r_mo17$p[pr[1], pr[2]])))
s3_all_sig <- all(s3_p < 1e-4, na.rm = TRUE)
s3 <- do.call(rbind, lapply(pairs6, function(pr) data.frame(
  Environments = paste(pr[1], "&", pr[2]),
  RIL = cor_cell(r_ril$r[pr[1], pr[2]], r_ril$p[pr[1], pr[2]], s3_all_sig),
  `B73 BC` = cor_cell(r_b73$r[pr[1], pr[2]], r_b73$p[pr[1], pr[2]], s3_all_sig),
  `Mo17 BC` = cor_cell(r_mo17$r[pr[1], pr[2]], r_mo17$p[pr[1], pr[2]], s3_all_sig),
  check.names = FALSE)))

# ---------------------------------------------------------------------------
# S5  broad-sense heritability (mirrors 02_get_blups.R:242-254)
# ---------------------------------------------------------------------------
h2 <- function(sub) {
  m <- lmer(wmd ~ (1 | Line) + YearLoc + (1 | YearLoc:Line) + YearLoc:Rep, REML = TRUE, data = sub)
  v <- as.data.frame(VarCorr(m))$vcov            # order: YearLoc:Line, Line, Residual
  VarGY <- v[1]; VarG <- v[2]; VarEps <- v[3]
  VarG / (VarG + VarGY / 4 + VarEps / (4 * 2))
}
s5 <- data.frame(
  Trait = "NLB WMD",
  `RIL Heritability`     = fmt_num(h2(filter(data_ibm, is_ril(Line))), 2),
  `B73 BC Heritability`  = fmt_num(h2(filter(data_ibm, is_b73(Line))), 2),
  `Mo17 BC Heritability` = fmt_num(h2(filter(data_ibm, is_mo17(Line))), 2),
  check.names = FALSE)

# ---------------------------------------------------------------------------
# S6  MPH/BPH of the B73 x Mo17 reference hybrid
# ---------------------------------------------------------------------------
ibm_blups <- read.csv("analyses/IBM_NLB_BLUPs.csv")
ref <- ibm_blups[ibm_blups$Line == "B73xMo17", ]
sub <- droplevels(subset(data_ibm, Line %in% c("B73", "Mo17", "B73xMo17") & !is.na(wmd)))
sub$Line <- factor(sub$Line)
m6 <- lm(wmd ~ Line + YearLoc + YearLoc:Rep, data = sub)
emm <- emmeans(m6, "Line"); es <- summary(emm); lv <- as.character(es$Line)
mph_c <- as.numeric(lv == "B73xMo17") - 0.5 * as.numeric(lv == "B73") - 0.5 * as.numeric(lv == "Mo17")
best  <- if (es$emmean[lv == "B73"] > es$emmean[lv == "Mo17"]) "B73" else "Mo17"
bph_c <- as.numeric(lv == "B73xMo17") - as.numeric(lv == best)
ct6 <- summary(contrast(emm, method = list(MPH = mph_c, BPH = bph_c)))
s6 <- data.frame(
  Trait = "NLB WMD",
  `MPH (%)` = fmt_pct(ref$NLB_WMD_BLUP_MPH_PCT, ct6$p.value[ct6$contrast == "MPH"]),
  `BPH (%)` = fmt_pct(ref$NLB_WMD_BLUP_BPH_PCT, ct6$p.value[ct6$contrast == "BPH"]),
  check.names = FALSE)

# ---------------------------------------------------------------------------
# S7  MPH/BPH of NAM diverse-line testcrosses
# ---------------------------------------------------------------------------
nam <- read.csv("analyses/NAM_NLB_BLUPs.csv")
get_cell <- function(line, kind) {                     # kind = "MPH" | "BPH"
  r <- nam[nam$Line == line, ]
  if (nrow(r) == 0) return("NA")
  fmt_pct(r[[paste0("NLB_WMD_BLUP_", kind, "_PCT")]], r[[paste0("NLB_WMD_BLUP_", kind, "_p")]])
}
parents <- c("B73", "B97", "CML52", "CML69", "CML103", "CML228", "CML247", "CML277",
             "CML322", "CML333", "HP301", "Il14H", "Ki3", "Ki11", "Ky21", "M37W",
             "M162W", "Mo18W", "Ms71", "NC350", "NC358", "Oh43", "Oh7B", "P39",
             "Tx303", "Tzi8", "Mo17")
parents <- parents[parents %in% sub("^[^x]*x", "", nam$Line[grepl("x", nam$Line)]) |
                     parents %in% c("B73", "Mo17")]
s7 <- do.call(rbind, lapply(parents, function(p) data.frame(
  `Diverse Parent` = p,
  `B73 Cross MPH (%)`  = get_cell(paste0("B73x", p),  "MPH"),
  `Mo17 Cross MPH (%)` = get_cell(paste0("Mo17x", p), "MPH"),
  `B73 Cross BPH (%)`  = get_cell(paste0("B73x", p),  "BPH"),
  `Mo17 Cross BPH (%)` = get_cell(paste0("Mo17x", p), "BPH"),
  check.names = FALSE)))

# ---------------------------------------------------------------------------
# S8  mean MPH/BPH of the backcross populations
# ---------------------------------------------------------------------------
b73bc  <- read.csv("analyses/IBM_NLB_BLUPs_B73BC.csv")
mo17bc <- read.csv("analyses/IBM_NLB_BLUPs_Mo17BC.csv")
mmm <- function(x) { x <- x[!is.na(x)] * 100
  sprintf("%.2f (%.2f, %.2f)", mean(x), min(x), max(x)) }
s8 <- data.frame(
  Trait = "NLB WMD",
  `B73 BC MPH (%)`  = mmm(b73bc$NLB_WMD_BLUP_MPH_PCT),
  `Mo17 BC MPH (%)` = mmm(mo17bc$NLB_WMD_BLUP_MPH_PCT),
  `B73 BC BPH (%)`  = mmm(b73bc$NLB_WMD_BLUP_BPH_PCT),
  `Mo17 BC BPH (%)` = mmm(mo17bc$NLB_WMD_BLUP_BPH_PCT),
  check.names = FALSE)

# ---------------------------------------------------------------------------
# S9  main-effect + epistatic-pair-member QTL
# ---------------------------------------------------------------------------
eff <- read.csv("analyses/qtl_effects.csv")
ga  <- read.csv("analyses/qtl_gene_action.csv")
# chromosome-specific QTL id (i.j) from the colocalized peaks in qtl_gene_action
ga <- ga %>% arrange(chr, pos) %>% group_by(chr) %>%
  mutate(qid = paste0(chr, ".", row_number())) %>% ungroup()
qid_of <- function(chr, pos) {
  hit <- ga$qid[ga$chr == chr & abs(ga$pos - pos) < 1e-6]
  if (length(hit)) hit[1] else paste0(chr, ".?")
}
sig_lab <- c("RIL|NLB_WMD_BLUP" = "RIL",
             "B73_BC|NLB_WMD_BLUP" = "B73 BC", "B73_BC|NLB_WMD_BLUP_MPH" = "B73 MPH",
             "Mo17_BC|NLB_WMD_BLUP" = "Mo17 BC", "Mo17_BC|NLB_WMD_BLUP_MPH" = "Mo17 MPH")

main <- eff %>% filter(effect_class == "main", !borrowed) %>% arrange(chr, pos)
s9_main <- data.frame(
  QTL = mapply(qid_of, main$chr, main$pos),
  Trait = "NLB WMD", Chr = main$chr, `Pos (IcM)` = fmt_num(main$pos, 1),
  `CI Low` = fmt_num(main$ci_lo, 1), `CI High` = fmt_num(main$ci_hi, 1),
  Effect = fmt_num(main$estimate, 2), LOD = fmt_num(main$lod, 2),
  `PVE Joint (%)` = fmt_num(main$fitqtl_pct_var, 2),
  `PVE Single (%)` = fmt_num(main$single_qtl_pct_var, 2),
  Sig = sig_lab[paste(main$cross, main$trait, sep = "|")], check.names = FALSE)

epi <- eff %>% filter(effect_class == "epistatic") %>%
  mutate(pop = ifelse(cross == "RIL", "RIL", "B73 BC"), Ek = row_number())
s9_epi <- do.call(rbind, lapply(seq_len(nrow(epi)), function(i) {
  r <- epi[i, ]; lab <- paste0(r$pop, " Epistasis")
  data.frame(
    QTL = c(paste0(r$chr, ".E", r$Ek, ".1"), paste0(r$chr2, ".E", r$Ek, ".2")),
    Trait = "NLB WMD", Chr = c(r$chr, r$chr2),
    `Pos (IcM)` = fmt_num(c(r$pos, r$pos2), 1),
    `CI Low` = "N/A", `CI High` = "N/A",
    Effect = fmt_num(rep(r$estimate, 2), 2), LOD = fmt_num(rep(r$lod, 2), 2),
    `PVE Joint (%)` = fmt_num(rep(r$fitqtl_pct_var, 2), 2),
    `PVE Single (%)` = fmt_num(rep(r$single_qtl_pct_var, 2), 2),
    Sig = lab, check.names = FALSE)
}))
s9 <- rbind(s9_main, s9_epi)

# ---------------------------------------------------------------------------
# S10  gene-action (d/a) classification of QTL
# ---------------------------------------------------------------------------
classes <- c("OD", "D +", "PD +", "Add", "PD -", "D -", "UD")
dir_class <- function(action, ratio) {
  if (is.na(action)) return(NA_character_)
  if (action == "additive") return("Add")
  if (action == "od") return("OD")
  if (action == "ud") return("UD")
  s <- if (!is.na(ratio) && ratio >= 0) "+" else "-"
  if (action == "pd") return(paste("PD", s))
  if (action == "dominant") return(paste("D", s))
  action
}
b73_cl  <- mapply(dir_class, ga$B73_action,  ga$`B73_d.a`)
mo17_cl <- mapply(dir_class, ga$Mo17_action, ga$`Mo17_d.a`)
count_cl <- function(v) sapply(classes, function(k) sum(v == k, na.rm = TRUE))
s10a <- data.frame(`Mode of action` = classes,
                   `B73 BC` = count_cl(b73_cl), `Mo17 BC` = count_cl(mo17_cl),
                   check.names = FALSE, row.names = NULL)
# B73 BC (rows) x Mo17 BC (cols) comparison of mode of action
tab <- table(factor(b73_cl, classes), factor(mo17_cl, classes))
s10b <- data.frame(`B73 BC \\ Mo17 BC` = classes, check.names = FALSE)
for (k in classes) s10b[[k]] <- as.integer(tab[, k])

# ---------------------------------------------------------------------------
# S11  correlation with the earlier published IBM NLB study
# ---------------------------------------------------------------------------
old <- read.csv("data/old_nlb_data.csv")
old$AUDPCNLBLSMEAN <- as.numeric(old$AUDPCNLBLSMEAN)
old$NLBAUDPCBLUP   <- as.numeric(old$NLBAUDPCBLUP)
cor_row <- function(x, lab) {
  ct <- cor.test(x, old$NLB_WMD_BLUP_new)
  data.frame(Comparison = lab, r = fmt_num(unname(ct$estimate), 3),
             `p-value` = fmt_p(ct$p.value),
             `95% CI` = sprintf("%.3f to %.3f", ct$conf.int[1], ct$conf.int[2]),
             check.names = FALSE)
}
s11 <- rbind(cor_row(old$AUDPCNLBLSMEAN, "AUDPC LS-mean vs new WMD BLUP"),
             cor_row(old$NLBAUDPCBLUP,   "AUDPC BLUP vs new WMD BLUP"))

# ---------------------------------------------------------------------------
# assemble the document
# ---------------------------------------------------------------------------
sig_note <- function(all_sig) if (all_sig)
  "All correlations significant at P < 0.0001." else
  "** significant at P < 0.0001, * significant at P < 0.05."

doc <- add_caption(doc, "Supplementary material")
doc <- add_note(doc, paste("Disease phenotype is weighted mean disease (WMD); BLUPs are on the",
                           "flipped resistance scale (100 - WMD), so a positive heterosis or a",
                           "positive B73-allele effect indicates greater resistance."))
doc <- body_add_par(doc, "", style = "Normal")

doc <- add_table(doc,
  "Table S1: Pearson correlation coefficients between replicates within years for NLB WMD, in the recombinant inbred lines (RIL), B73 backcross (B73 BC), and Mo17 backcross (Mo17 BC) populations.",
  s1, notes = sig_note(s1_all_sig))

doc <- add_table(doc,
  paste("Table S2: Pearson correlation coefficients for NLB disease scores between scoring dates,",
        "in the RIL, B73 BC, and Mo17 BC populations. For 2023_CL, Date 1 was 07/17, Date 2 07/25, and Date 3 08/01.",
        "For 2024_CL, Date 1 was 07/17 and Date 2 07/28. For 2024_IL, Date 1 was 07/24, Date 2 was 08/01 or 08/02,",
        "and Date 3 was 08/08 or 08/09 (depending on plot). 2022_CL was scored on a single date and is omitted."),
  s2, notes = sig_note(s2_all_sig))

doc <- add_table(doc,
  "Table S3: Pearson correlation coefficients of genotypic NLB WMD BLUPs between environments, in the RIL, B73 BC, and Mo17 BC populations.",
  s3, notes = sig_note(s3_all_sig))

doc <- add_table(doc,
  "Table S5: Broad-sense heritability estimates for NLB WMD in each population.",
  s5)

doc <- add_table(doc,
  "Table S6: MPH and BPH values of NLB WMD in the B73 x Mo17 reference hybrid grown in the same trials as the RIL and backcross populations.",
  s6, notes = "Point estimate of the percent deviation from the mid-parent (MPH) or best-parent (BPH) value, on the resistance scale. ** significant at P < 0.0001, * significant at P < 0.05.")

doc <- add_table(doc,
  paste("Table S7: MPH and BPH values of NLB WMD in diverse lines (the NAM diverse lines as well as",
        "B73 and Mo17) crossed to B73 and Mo17. The line under \"Diverse Parent\" is the male parent,",
        "so the results for B73 and Mo17 are reciprocal crosses."),
  s7, notes = "Percent deviation from the mid-parent (MPH) or best-parent (BPH) value, on the resistance scale. NA where the cross or a parent BLUP was unavailable. ** significant at P < 0.0001, * significant at P < 0.05.")

doc <- add_table(doc,
  "Table S8: Mean MPH and BPH values of genotypic NLB WMD BLUPs for the B73 backcross and Mo17 backcross populations.",
  s8, notes = "As percentage of the mid-parent (MPH) or best-parent (BPH) value, on the resistance scale. Mean of all lines with minimum and maximum values in parentheses.")

doc <- add_table(doc,
  paste("Table S9: All main-effect QTL identified as significant based on permutation-derived",
        "thresholds in across-environment BLUPs, plus QTL that are members of pairs showing a",
        "significant epistatic interaction at the 5% level."),
  s9, notes = c(
    "QTL: rows with the same i.j name are co-located based on overlapping confidence intervals (i = chromosome, j = a chromosome-specific identifier). Rows named i.Ek.j are the two members of epistatic pair Ek.",
    "Effect: in the RIL analysis the additive effect a; in the mid-parent-heterosis (MPH) analyses the dominance effect d; in the backcross (BC) BLUP analyses a confounded combination of a and d. Positive values indicate the B73 allele increases resistance (a) or that the heterozygote is more resistant than the mid-parent (d).",
    "PVE Joint: percent variance explained by this QTL alone within the full fitted model for its population/trait (all of that population's other significant QTL and genuine epistatic pairs included), i.e. its unique contribution controlling for the rest of the model.",
    "PVE Single: percent variance explained by this QTL fit on its own against a no-QTL null model, ignoring all other QTL; differs from PVE Joint whenever the population has more than one mapped QTL for that trait.",
    "For epistatic pair members, Effect, LOD, PVE Joint, and PVE Single are the interaction effect/LOD/PVE of the pair (shared by both members); confidence intervals are not given as these QTL were not significant as main effects.",
    "Sig: the analysis in which the QTL was significant at the 5% level."))

doc <- add_table(doc,
  "Table S10: Numbers of NLB WMD QTL estimated to have each ratio of dominance to additive effects in the B73 BC and Mo17 BC populations.",
  s10a, notes = "OD overdominant, D dominant, PD partially dominant, Add additive, UD underdominant. \"+\" indicates the resistance-increasing allele is dominant, \"-\" the resistance-decreasing allele. Counts include borrowed (non-significant) gene-action estimates.")

doc <- add_table(doc,
  "Table S10 (continued): Comparison of the estimated mode of action of each colocalized NLB WMD QTL between the B73 BC (rows) and Mo17 BC (columns) populations.",
  s10b)

doc <- add_table(doc,
  "Table S11: Correlation of the new NLB WMD BLUPs with the earlier published IBM NLB study. Negative because the old scale is disease severity (AUDPC) and the new scale is resistance (100 - WMD).",
  s11)

print(doc, target = out_docx)
cat("Wrote", out_docx, "\n")
cat("S1 all sig:", s1_all_sig, "| S2 all sig:", s2_all_sig, "| S3 all sig:", s3_all_sig, "\n")
cat("Heritability (RIL/B73/Mo17):", unlist(s5[2:4]), "\n")
cat("S9 rows:", nrow(s9), "| S7 parents:", nrow(s7), "\n")

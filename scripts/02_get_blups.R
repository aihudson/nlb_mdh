library(dplyr)
library(readxl)
library(openxlsx)
library(lme4)
library(brms)
library(tidyr)
# nlb
data1 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=1)
data2 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=2)
data3 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=3)
data4 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=4)

data1 <- data1 %>%
  mutate(wmd=NLB_0728)

data_list <- list(data1,data2,data3,data4)
# merge
for(i in 1:length(data_list)){
  data_list[[i]] <- data_list[[i]] %>%
    dplyr::select(Row,Line,Loc,Year,Rep,wmd)
}

data <- do.call(rbind,data_list)
data$Year <- as.factor(data$Year)
data$Rep <- as.factor(data$Rep)
data$Loc <- as.factor(data$Loc)
data$YearLoc <- paste(data$Year, data$Loc, sep = "_")

# there are two lines where a different maternal parent was used due to seed quantity issues
# for the purposes of analyses here we are treating them as equivalent, but in the supplemental data the correct parental genotypes are noted
data$Line[which(data$Line == "Mo17xM0230")] <- "M0230xMo17"
data$Line[which(data$Line == "Mo17xM0318")] <- "M0318xMo17"

# change names for IL lines
data <- data %>%
  mutate(Line = gsub("M0(\\d+x)M0(\\d+)", "M0\\2", Line))

### get BLUPs for IBM lines

data_ibm <- data %>%
  filter(Line %in% c("B73", "Mo17", "B73xMo17") | grepl("M0", Line))

# reverse direction of wmd

data_ibm <- data_ibm %>%
  mutate(wmd = 100-wmd)

model_ibm <- lmer(wmd ~ + (1|Line) + YearLoc + (1|YearLoc:Line) + YearLoc:Rep, data = data_ibm)

model_ranef_ibm <- ranef(model_ibm)
line_blups_ibm <- model_ranef_ibm$Line
colnames(line_blups_ibm) <- c("NLB_WMD_BLUP")
line_blups_ibm <- line_blups_ibm + summary(model_ibm)$coefficients[1]

# brms heterogeneous variance

# brms_model <- brm(
#   formula = bf(wmd ~ YearLoc  + (1 | Line) + (1 | Line:YearLoc), sigma ~ 0 + YearLoc),
#   data = data_ibm,
#   family = gaussian(),
#   chains = 4, 
#   iter = 500
# )
# re <- ranef(brms_model)
# re_line <- re_line[,1,1]
# summary(re_line)
# get heterosis for ibm lines

line_blups_ibm$Line <- rownames(line_blups_ibm)
add_heterosis <- function(df, trait, parent){
  hybrid_index <- grep("*x", df$Line)
  
  hybrids <- df[hybrid_index,]
  
  hybrids$Parent1 <- sub("x.*", "", hybrids$Line)
  # split on the FIRST x only, so a founder whose name contains an x (e.g. NAM's Tx303)
  # is kept whole as Parent2 rather than truncated (greedy ".*x" would leave "303")
  hybrids$Parent2 <- sub("^[^x]*x", "", hybrids$Line)
  hybrids$Parent1_BLUP <- df[[trait]][match(hybrids$Parent1, df$Line)]
  hybrids$Parent2_BLUP <- df[[trait]][match(hybrids$Parent2, df$Line)]
  
  trait_mph <- paste(trait, "_MPH", sep = "") 
  trait_mph_pct <- paste(trait, "_MPH_PCT", sep = "") 
  
  hybrids <- hybrids %>%
    rowwise() %>%
    dplyr::mutate(MP = mean(c(Parent1_BLUP, Parent2_BLUP))) %>%
    dplyr::mutate(!!trait_mph := !!as.name(trait) - MP,
                  !!trait_mph_pct := (!!as.name(trait) - MP)/MP) %>%
    select(Line, !!trait_mph, !!trait_mph_pct, Parent1_BLUP, Parent2_BLUP, NLB_WMD_BLUP)
  
  trait_bph <- paste(trait, "_BPH", sep = "") 
  trait_bph_pct <- paste(trait, "_BPH_PCT", sep = "") 
  
  hybrids <- hybrids %>%
    rowwise() %>%
    dplyr::mutate(BP = max(c(Parent1_BLUP, Parent2_BLUP))) %>%
    dplyr::mutate(!!trait_bph := !!as.name(trait) - BP,
                  !!trait_bph_pct := (!!as.name(trait) - BP)/BP) %>%
    select(Line, !!trait_bph, !!trait_bph_pct, !!trait_mph, !!trait_mph_pct)
  
  df <- full_join(df, hybrids)
  
  return(df)
}

line_blups_ibm2 <- add_heterosis(line_blups_ibm, "NLB_WMD_BLUP")

write.csv(line_blups_ibm2, file = "analyses/IBM_NLB_BLUPs.csv")

# write files for rils, b73 bc, and mo17 bc

line_blups_rils <- line_blups_ibm2 %>%
  filter(grepl("M0", Line)) %>%
  filter(!grepl("x", Line))

line_blups_b73bc <- line_blups_ibm2 %>%
  filter(grepl("B73", Line)) %>%
  filter(grepl("x", Line)) %>%
  filter(!grepl("Mo17", Line))

line_blups_mo17bc <- line_blups_ibm2 %>%
  filter(grepl("Mo17", Line)) %>%
  filter(grepl("x", Line)) %>%
  filter(!grepl("B73", Line))

write.csv(line_blups_rils, file = "analyses/IBM_NLB_BLUPs_RILs.csv")
write.csv(line_blups_b73bc, file = "analyses/IBM_NLB_BLUPs_B73BC.csv")
write.csv(line_blups_mo17bc, file = "analyses/IBM_NLB_BLUPs_Mo17BC.csv")

### get BLUPs for NAM lines
NAM_parents <- c( "B73", "B97", "CML52", "CML69", "CML103", 
                  "CML228", "CML247", "CML277", "CML322", 
                  "CML333", "HP301", "Il14H", "Ki3", "Ki11", "Ky21",
                  "M37W", "M162W", "Mo18W", "Ms71", "NC350", "NC358",
                  "Oh43", "Oh7B", "P39", "Tx303", "Tzi8")

mothers <- c("B73", "Mo17")

line_filter <- c(NAM_parents, mothers, "B73xMo17", as.vector(outer(mothers, NAM_parents, paste, sep="x")))

data_nam <- data %>%
  filter(Line %in% line_filter) 

# filter out B73 and Mo17 rows from 2024, which only had the IBM experiment
data_nam <- data_nam %>%
  filter(Year != "2024")

# reverse direction of wmd to the resistance scale, matching the IBM path so that
# best-parent heterosis (max) means the most resistant parent
data_nam <- data_nam %>%
  mutate(wmd = 100 - wmd)

library(blme)
library(emmeans)

control <- lmerControl(optimizer = "bobyqa", 
                       optCtrl = list(maxfun = 2e5))

(model_nam_blme <- blmer(wmd ~ (1|Line) + YearLoc + (1|YearLoc:Line), data = subset(data_nam, !is.na(wmd)),
              cov.prior = wishart, REML = FALSE, control = control))


(model_nam <- lmer(wmd ~ (1|Line) + YearLoc + (1|YearLoc:Line) + YearLoc:Rep, data = subset(data_nam, !is.na(wmd)), control = control))

model_ranef_nam_blme <- ranef(model_nam_blme)
line_blups_nam_blme <- model_ranef_nam_blme$Line

model_ranef_nam_lmer <- ranef(model_nam)
line_blups_nam_lmer <- model_ranef_nam_lmer$Line

colnames(line_blups_nam_blme) <- c("NLB_WMD_BLUP")
line_blups_nam_blme <- line_blups_nam_blme+ summary(model_nam_blme)$coefficients[1]

# heterosis point estimates (MPH/BPH) from the NAM BLUPs, via the same helper used for the IBM lines
line_blups_nam_blme$Line <- rownames(line_blups_nam_blme)
nam_heterosis <- add_heterosis(line_blups_nam_blme, "NLB_WMD_BLUP")

# significance of MPH/BPH per cross: fixed-effects model + emmeans contrasts
# (Y = mu + YearLoc + Rep(YearLoc) + Line + Line:YearLoc, all fixed; NAM here is single-location
#  2022+2023, so YearLoc stands in for the paper's Year)
nam_lm <- lm(wmd ~ Line + YearLoc + YearLoc:Rep + Line:YearLoc, data = subset(data_nam, !is.na(wmd)))
EMM <- emmeans(nam_lm, "Line")
emm_s <- summary(EMM)

founders <- setdiff(NAM_parents, mothers)   # loop parent1 over founders only (skip B73xB73 etc.)
B73_vec  <- as.numeric(emm_s$Line == "B73")
Mo17_vec <- as.numeric(emm_s$Line == "Mo17")

mph_list <- list(); bph_list <- list()
for (m in mothers) {
  mvec <- if (m == "B73") B73_vec else Mo17_vec
  for (p in founders) {
    hyb <- paste0(m, "x", p)
    # need the hybrid and both parents present as fitted lines; otherwise the contrast would be
    # a degenerate partial one (e.g. missing-parent -> spurious p), so leave it NA like the BLUPs
    if (!any(emm_s$Line == hyb) || !any(emm_s$Line == p) || !any(emm_s$Line == m)) next
    hvec <- as.numeric(emm_s$Line == hyb)
    pvec <- as.numeric(emm_s$Line == p)
    mph_list[[paste(hyb, "MPH")]] <- hvec - 0.5*mvec - 0.5*pvec
    # best = most resistant parent (higher emmean on the resistance scale); na.rm so a single
    # non-estimable line elsewhere doesn't NA out the whole sum
    best <- ifelse(sum(emm_s$emmean*mvec, na.rm=TRUE) > sum(emm_s$emmean*pvec, na.rm=TRUE), m, p)
    bph_list[[paste(hyb, "BPH")]] <- hvec - as.numeric(emm_s$Line == best)
  }
}

# reference B73 x Mo17 F1 (mother x mother, so not generated by the founder loop above)
if (all(c("B73xMo17", "B73", "Mo17") %in% emm_s$Line)) {
  hvec <- as.numeric(emm_s$Line == "B73xMo17")
  mph_list[["B73xMo17 MPH"]] <- hvec - 0.5*B73_vec - 0.5*Mo17_vec
  best <- ifelse(sum(emm_s$emmean*B73_vec, na.rm=TRUE) > sum(emm_s$emmean*Mo17_vec, na.rm=TRUE), "B73", "Mo17")
  bph_list[["B73xMo17 BPH"]] <- hvec - as.numeric(emm_s$Line == best)
}
mph_c <- summary(contrast(EMM, method = mph_list))
bph_c <- summary(contrast(EMM, method = bph_list))
sig_mph <- data.frame(Line = sub(" MPH$", "", mph_c$contrast), NLB_WMD_BLUP_MPH_p = mph_c$p.value)
sig_bph <- data.frame(Line = sub(" BPH$", "", bph_c$contrast), NLB_WMD_BLUP_BPH_p = bph_c$p.value)

nam_out <- nam_heterosis %>%
  left_join(sig_mph, by = "Line") %>%
  left_join(sig_bph, by = "Line")
write.csv(nam_out, file = "analyses/NAM_NLB_BLUPs.csv", row.names = FALSE)

### get heritability for IBM RILs, IBMxB73, and IBMxMo17
data_ibm_ril <- data_ibm %>%
  filter(grepl("M0", Line) & !grepl("x", Line))

model_ibm_ril <- lmer(wmd ~ (1|Line) + YearLoc + (1|YearLoc:Line) + YearLoc:Rep, REML = TRUE, data = data_ibm_ril)

data_ibm_b73 <- data_ibm %>%
  filter(grepl("B73", Line) & grepl("M0", Line))

model_ibm_b73 <- lmer(wmd ~ (1|Line) + YearLoc + (1|YearLoc:Line) + YearLoc:Rep, REML = TRUE, data = data_ibm_b73)

data_ibm_mo17 <- data_ibm %>%
  filter(grepl("Mo17", Line) & grepl("M0", Line))

(model_ibm_mo17 <- lmer(wmd ~ (1|Line) + YearLoc + (1|YearLoc:Line) + YearLoc:Rep, REML = TRUE, data = data_ibm_mo17))


# H^2 = VarG/(VarG + VarG:YL/(n(YL) + VarEps/((n(YL)*n(r))))
models <- list(model_ibm_ril,model_ibm_b73,model_ibm_mo17)
for(i in 1:length(models)){
  model <- models[[i]]
  varcomps <- as.data.frame(print(VarCorr(model),comp="Variance"))$vcov
  VarG <- varcomps[2]
  VarGY <- varcomps[1]
  VarEps <- varcomps[3]
  nY <- 4
  nR <- 2
  H2 <- VarG/(VarG + (VarGY/nY) + (VarEps/(nY*nR)))
  print(H2)
}

# within year rep correlations for IBM RILs, IBMxB73, and IBMxMo17
# for each YearLoc, calculate correlation of wmd for Rep 1 and 2
# RILs 
data_ibm_ril <- data %>%
  filter(grepl("M0", Line) & !grepl("x", Line))
# reshape data to have columns for Rep 1 and Rep 2, while keeping YearLoc and Line as identifiers
data_ibm_ril_wide <- data_ibm_ril %>%
  select(YearLoc, Line, Rep, wmd) %>%
  pivot_wider(names_from = Rep, values_from = wmd)

for(env in unique(data_ibm_ril_wide$YearLoc)){
  data_sub <- data_ibm_ril_wide %>%
    filter(YearLoc == env)
  cor_test <- cor.test(data_sub$`1`, data_sub$`2`)
  print(paste(env, "r =", cor_test$estimate, "p =", cor_test$p.value))
}

# B73 crosses
data_ibm_b73 <- data %>%
  filter(grepl("B73", Line) & grepl("M0", Line))
data_ibm_b73_wide <- data_ibm_b73 %>%
  select(YearLoc, Line, Rep, wmd) %>%
  pivot_wider(names_from = Rep, values_from = wmd)

for(env in unique(data_ibm_b73_wide$YearLoc)){
  data_sub <- data_ibm_b73_wide %>%
    filter(YearLoc == env)
  cor_test <- cor.test(data_sub$`1`, data_sub$`2`)
  print(paste(env, "r =", cor_test$estimate, "p =", cor_test$p.value))
}

# Mo17 crosses
data_ibm_mo17 <- data %>%
  filter(grepl("Mo17", Line) & grepl("M0", Line))
data_ibm_mo17_wide <- data_ibm_mo17 %>%
  select(YearLoc, Line, Rep, wmd) %>% 
  pivot_wider(names_from = Rep, values_from = wmd)

for(env in unique(data_ibm_mo17_wide$YearLoc)){
  data_sub <- data_ibm_mo17_wide %>%
    filter(YearLoc == env)
  cor_test <- cor.test(data_sub$`1`, data_sub$`2`)
  print(paste(env, "r =", cor_test$estimate, "p =", cor_test$p.value))
}

# correlations of raw disease scores between scoring dates for IBM RILs, IBMxB73,
# and IBMxMo17 (2022 CL was scored on a single date and so is excluded). Uses the
# per-date columns on the original sheets (data2/3/4), which the merged `data`
# above collapses into wmd. Normalize line names the same way as above so the
# population filters work on the 2024 IL sheet, where RILs are recorded as selfs
# ("M0###xM0###").
normalize_line <- function(l){
  l[l == "Mo17xM0230"] <- "M0230xMo17"
  l[l == "Mo17xM0318"] <- "M0318xMo17"
  gsub("M0(\\d+x)M0(\\d+)", "M0\\2", l)
}
data2$Line <- normalize_line(data2$Line)
data3$Line <- normalize_line(data3$Line)
data4$Line <- normalize_line(data4$Line)

# 2024 IL date 2 was 08/01 or 08/02 and date 3 was 08/08 or 08/09 depending on the
# plot; collapse each to a single timepoint
data4$NLB_date2 <- ifelse(!is.na(as.numeric(data4$NLB_0801)), as.numeric(data4$NLB_0801), as.numeric(data4$NLB_0802))
data4$NLB_date3 <- ifelse(!is.na(as.numeric(data4$NLB_0808)), as.numeric(data4$NLB_0808), as.numeric(data4$NLB_0809))

date_sets <- list(
  "2023_CL" = list(df = data2, dates = c("NLB_0717", "NLB_0725", "NLB_0801")),
  "2024_CL" = list(df = data3, dates = c("NLB_0717", "NLB_0728")),
  "2024_IL" = list(df = data4, dates = c("NLB_0724", "NLB_date2", "NLB_date3")))
pop_filters <- list(
  RILs   = function(l) grepl("M0", l) & !grepl("x", l),
  B73BC  = function(l) grepl("B73", l) & grepl("M0", l),
  Mo17BC = function(l) grepl("Mo17", l) & grepl("M0", l))

for(env in names(date_sets)){
  dt <- date_sets[[env]]$dates
  for(pop in names(pop_filters)){
    sub <- date_sets[[env]]$df[pop_filters[[pop]](date_sets[[env]]$df$Line), ]
    for(i in 1:(length(dt) - 1)){
      for(j in (i + 1):length(dt)){
        x <- suppressWarnings(as.numeric(sub[[dt[i]]]))
        y <- suppressWarnings(as.numeric(sub[[dt[j]]]))
        ok <- is.finite(x) & is.finite(y)
        if(sum(ok) < 3) next
        cor_test <- cor.test(x[ok], y[ok])
        print(paste(env, pop, "Date", i, "&", j, "r =", cor_test$estimate, "p =", cor_test$p.value))
      }
    }
  }
}

# correlations of BLUPs between LocYears for IBM RILs, IBMxB73, and IBMxMo17
# first get BLUPs for each LocYear
for(env in unique(data_ibm$YearLoc)){
  data_sub <- data_ibm %>%
    filter(YearLoc == env)
  model_sub <- lmer(wmd ~ (1|Line) + Rep, data = data_sub)
  model_ranef_sub <- ranef(model_sub)
  line_blups_sub <- model_ranef_sub$Line
  colnames(line_blups_sub) <- "NLB_WMD_BLUP"
  line_blups_sub$Line <- rownames(line_blups_sub)
  line_blups_sub$YearLoc <- env
  row.names(line_blups_sub) <- NULL
  line_blups_sub$NLB_WMD_BLUP <- line_blups_sub$NLB_WMD_BLUP + summary(model_sub)$coefficients[1]
  if(env == unique(data_ibm$YearLoc)[1]){
    line_blups_all <- line_blups_sub
  } else {
    line_blups_all <- rbind(line_blups_all, line_blups_sub)
  }
}

# pivot wider to have columns for each YearLoc
line_blups_all_wide <- line_blups_all %>%
  pivot_wider(names_from = YearLoc, values_from = NLB_WMD_BLUP)

# calculate correlations between YearLocs for RILs, B73 crosses, and Mo17 crosses
line_blups_all_wide_rils <- line_blups_all_wide %>%
  filter(grepl("M0", Line) & !grepl("x", Line))
line_blups_all_wide_b73 <- line_blups_all_wide %>%
  filter(grepl("B73", Line) & grepl("M0", Line))
line_blups_all_wide_mo17 <- line_blups_all_wide %>%
  filter(grepl("Mo17", Line) & grepl("M0", Line))

# correlation matrices
cor_rils <- cor(line_blups_all_wide_rils[,-1], use = "pairwise.complete.obs")
cor_b73 <- cor(line_blups_all_wide_b73[,-1], use = "pairwise.complete.obs")
cor_mo17 <- cor(line_blups_all_wide_mo17[,-1], use = "pairwise.complete.obs")

# matrices of correlation p-values (pairwise complete obs, same as cor() above)
cor_pmat <- function(df){
  n <- ncol(df)
  pmat <- matrix(NA, n, n, dimnames = list(colnames(df), colnames(df)))
  for(i in 1:n){
    for(j in 1:n){
      if(i == j){
        pmat[i, j] <- 0
      } else {
        pmat[i, j] <- cor.test(df[[i]], df[[j]])$p.value
      }
    }
  }
  return(pmat)
}
cor_rils_p <- cor_pmat(line_blups_all_wide_rils[,-1])
cor_b73_p <- cor_pmat(line_blups_all_wide_b73[,-1])
cor_mo17_p <- cor_pmat(line_blups_all_wide_mo17[,-1])

print("RILs cross-env BLUP correlations"); print(round(cor_rils, 3)); print("p-values"); print(signif(cor_rils_p, 3))
print("B73 BC cross-env BLUP correlations"); print(round(cor_b73, 3)); print("p-values"); print(signif(cor_b73_p, 3))
print("Mo17 BC cross-env BLUP correlations"); print(round(cor_mo17, 3)); print("p-values"); print(signif(cor_mo17_p, 3))

write.csv(line_blups_all_wide_rils, file = "analyses/line_blups_envs.csv")

# correlation with old NLB data
olddata <- read_xls("~/projects/nlb_mdh/data/IBM NLB DATA.xls",sheet=7)
olddata$`Lines QTLCART ORDER` <- gsub("M005", "M0005", olddata$`Lines QTLCART ORDER`)

# add new blups
olddata_index <- match(olddata$`Lines QTLCART ORDER`, line_blups_ibm2$Line)
olddata$NLB_WMD_BLUP_new <- line_blups_ibm2$NLB_WMD_BLUP[olddata_index]
olddata$AUDPCNLBLSMEAN <- as.numeric(olddata$AUDPCNLBLSMEAN)
olddata$NLBAUDPCBLUP <- as.numeric(olddata$NLBAUDPCBLUP)
cor.test(olddata$AUDPCNLBLSMEAN, olddata$NLB_WMD_BLUP_new)
cor.test(olddata$NLBAUDPCBLUP, olddata$NLB_WMD_BLUP_new)
write.csv(olddata, file = "data/old_nlb_data.csv")

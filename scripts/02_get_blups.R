library(dplyr)
library(readxl)
library(openxlsx)
library(lme4)
library(brms)
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
  hybrids$Parent2 <- sub(".*x", "", hybrids$Line)
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
                  "CML333", "Hp301", "Il14H", "Ki3", "Ki11", "Ky21", 
                  "M37W", "M162W", "Mo18W", "MS71", "NC350", "NC358", 
                  "Oh43", "Oh7B", "P39", "Tx303", "Tzi8")

mothers <- c("B73", "Mo17")

line_filter <- c(NAM_parents, mothers, as.vector(outer(mothers, NAM_parents, paste, sep="x")))

data_nam <- data %>%
  filter(Line %in% line_filter) 

# filter out B73 and Mo17 rows from 2024, which only had the IBM experiment
data_nam <- data_nam %>%
  filter(Year != "2024")

library(blme)

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
write.csv(line_blups_nam_blme, file = "analyses/NAM_NLB_BLUPs.csv")

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
  cor_val <- cor(data_sub$`1`, data_sub$`2`, use = "complete.obs")
  print(paste(env, cor_val))
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
  cor_val <- cor(data_sub$`1`, data_sub$`2`, use = "complete.obs")
  print(paste(env, cor_val))
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
  cor_val <- cor(data_sub$`1`, data_sub$`2`, use = "complete.obs")
  print(paste(env, cor_val))
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

write.csv(line_blups_all_wide_rils, file = "analyses/line_blups_envs.csv")

# correlation with old NLB data
olddata <- read_xls("~/projects/nlb_mdh/data/IBM NLB DATA.xls",sheet=7)
olddata$`Lines QTLCART ORDER` <- gsub("M005", "M0005", olddata$`Lines QTLCART ORDER`)

# add new blups
olddata_index <- match(olddata$`Lines QTLCART ORDER`, line_blups_ibm2$Line)
olddata$NLB_WMD_BLUP_new <- line_blups_ibm2$NLB_WMD_BLUP[olddata_index]
olddata$AUDPCNLBLSMEAN <- as.numeric(olddata$AUDPCNLBLSMEAN)
olddata$NLBAUDPCBLUP <- as.numeric(olddata$NLBAUDPCBLUP)
cor(olddata$AUDPCNLBLSMEAN, olddata$NLB_WMD_BLUP_new, use = "complete.obs")
cor(olddata$NLBAUDPCBLUP, olddata$NLB_WMD_BLUP_new, use = "complete.obs")
write.csv(olddata, file = "data/old_nlb_data.csv")

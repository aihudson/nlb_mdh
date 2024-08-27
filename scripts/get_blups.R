library(dplyr)
library(readxl)
library(openxlsx)
library(lme4)

data1 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=1)
data2 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=2)


data1 <- data1 %>%
  mutate(wmd=NLB_0728)

data_list <- list(data1,data2)
# merge
for(i in 1:length(data_list)){
  data_list[[i]] <- data_list[[i]] %>%
    select(Row,Line,Loc,Year,Rep,wmd)
}

data <- do.call(rbind,data_list)
data$Year <- as.factor(data$Year)
data$Rep <- as.factor(data$Rep)

model <- lmer(wmd ~ (1|Line) + Year + (1|Year:Line) + Year:Rep, data = data)

library(dplyr)
library(readxl)
library(openxlsx)
library(brms)

data1 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=1)
data2 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=2)
data3 <- read.xlsx("~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx",sheet=3)


data1 <- data1 %>%
  mutate(wmd=NLB_0728)

data_list <- list(data1,data2,data3)
# merge
for(i in 1:length(data_list)){
  data_list[[i]] <- data_list[[i]] %>%
    select(Row,Line,Loc,Year,Rep,Block,x,y,wmd)
}

data <- do.call(rbind,data_list)
data$Year <- as.factor(data$Year)
data$Rep <- as.factor(data$Rep)
data$Block <- as.factor(data$Block)
data <- data %>%
  mutate(Population = case_when(
    grepl("xB73", Line) ~ "B73_BC",
    grepl("xMo17", Line) ~ "Mo17_BC",
    !grepl("x",Line) ~ "Inbred")) %>%
    filter(!is.na(Population))

model_simple <- brm(
  wmd ~ (1|Line) + Year + Year:Rep + Year:Rep:Block,
  data = data,
  family = gaussian())

model_simple <- brm(
  wmd ~ (1|gr(Line, by = Population)) + Year + Year:Rep,
  data = data,
  family = gaussian())

fit2 <- brm(
  wmd ~ (1|gr(Line, by = Population)) + (1|gr(Line:Year, by = Population)) + Year + Year:Rep,
  data = data,
  family = gaussian())

ranef <- ranef(fit2)
ranef <- ranef$Line
ranef <- as.data.frame(ranef)
ranef$Line <- rownames(ranef)
rownames(ranef) <- NULL

write.csv(ranef,"analyses/nlb_blups_1.csv",row.names = FALSE,quote = FALSE)

east <- north <- 1:10
Grid <- expand.grid(east, north)
K <- nrow(Grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1

# try using for just 2023 block 1

data2 <- data %>%
  filter(Year==2023,
         Block=="Block1") %>%
  filter(!is.na(Population)) %>%
  filter(!is.na(wmd))

K <- nrow(data2)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(data.frame(data2$x,data2$y)))
W <- array(0, c(K, K))
W[distance == 1] <- 1

fit_car_1 <- brm(
  wmd ~ (1|gr(Line, by = Population)),
                 data = data2,
                 family = gaussian())

fit_car_2 <- update(fit_car_1, formula. = ~ . + car(W, type = "icar"), newdata = data2,
                    data2 = list(W = W)) 

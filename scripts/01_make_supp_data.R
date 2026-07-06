library(dplyr)
library(readxl)
library(openxlsx)
library(stringr)
library(lubridate)


wb <- createWorkbook()

### NLB 2022
data <- read.csv("~/projects/mdh_qtl/data/NLB_2022_raw.csv")
colnames(data) <- c("delete", "Row", "Line", "delete", "Rep", "NLB_0728", "NLB_0728_note","Block")
data <- data %>%
  select(-delete)
data <- data %>%
  mutate(Line=gsub(" X ","x",Line))
data <- data %>%
  mutate(Loc="CL") %>%
  mutate(Year = 2022)
data <- data %>%
  relocate(Row,Line,Loc,Year,Rep,Block,NLB_0728,NLB_0728_note)
data <- data %>%
  arrange(Row)

# adding block and x y coordinates

data <- data %>%
  mutate(x_dim = case_when(
    Block == 1 ~ 8,
    Block == 2 ~ 28,
    Block == 3 ~ 8
  ))

data <- x_y_coords(data)

addWorksheet(wb, sheetName = "NLB 2022")
writeDataTable(wb, sheet = 1, x = data, colNames = TRUE, rowNames = FALSE)

### NLB 2023

data <- read.csv("~/projects/mdh_qtl/data/nlb_2023_raw.csv")
colnames(data) <- c("delete", "Row", "Line", "delete","delete","delete","Rep",
                    "NLB_0717", "NLB_0725", "NLB_0801", "NLB_0717_note","NLB_0725_note", 
                    "NLB_0801_note","delete","delete","NLB_0731_pbk_note", "NLB_0731_pbk", 
                    "delete", "delete")
data <- data %>%
  select(-delete)
data <- data %>%
  mutate(Loc="CL") %>%
  mutate(Year = 2023) %>%
  rowwise() %>%
  mutate(wmd = ((mean(c(NLB_0717,NLB_0725)))*8+mean(c(NLB_0725,NLB_0801))*7)/15)

data <- data %>%
  mutate(Line=gsub("[ ]", "", Line)) %>%
  mutate(Line=gsub("X", "x", Line)) %>%
  mutate(Line=gsub("F1", "", Line))

data <- data %>%
  mutate(
    Block = case_when(
      Row > 11906 & Row < 12312 ~ "Block1",
      Row > 12311 & Row < 15101 ~ "Block2",
      Row > 15100 & Row < 15926 ~ "Block3",
      Row > 15925 & Row < 16046 ~ "Block4"
  ))

data <- data %>%
  relocate(Row,Line,Loc,Year,Rep,Block,NLB_0717,NLB_0717_note,NLB_0725,NLB_0725_note,
           NLB_0801,NLB_0801_note,NLB_0731_pbk,NLB_0731_pbk_note,wmd)
data <- data %>%
  arrange(Row)
data <- data %>%
  filter(Row < 15875)

# adding block and x y coordinates

data$Row <- as.numeric(data$Row)

data <- data %>%
  mutate(
  x_dim = case_when(
    Block == "Block1" ~ 11,
    Block == "Block2" ~ 25,
    Block == "Block3" ~ 25
  ))

data <- x_y_coords(data)

data$Line <- gsub("IBM", "M0", data$Line)

addWorksheet(wb, sheetName = "NLB 2023")
writeDataTable(wb, sheet = 2, x = data, colNames = TRUE, rowNames = FALSE)

### NLB 2024

data <- read.csv("~/projects/mdh_qtl/data/2024-07-28-11-29-03_NLB_table.csv")
data <- data %>% 
  select(Row,line,Nlb.7.17.ah,Nlb.7.28.24.ah)
colnames(data) <- c("Row", "Line", "NLB_0717", "NLB_0728")

data <- data %>%
  mutate(Loc="CL") %>%
  mutate(Year = 2024) %>%
  mutate(NLB_0717_note = gsub("[^a-zA-Z]", "", NLB_0717)) %>%
  mutate(NLB_0728_note = gsub("[^a-zA-Z]", "", NLB_0728)) %>%
  mutate(NLB_0717 = (gsub("[a-zA-Z]", "", NLB_0717))) %>%
  mutate(NLB_0717 = as.numeric(gsub(":", "", NLB_0717))) %>%
  mutate(NLB_0728 = (gsub("[a-zA-Z]", "", NLB_0728))) %>%
  mutate(NLB_0728 = as.numeric(gsub(":", "", NLB_0728))) %>%
  rowwise() %>%
  mutate(wmd = ((mean(c(NLB_0717,NLB_0728)))))
  
data <- data %>%
  mutate(Rep = (ifelse(Row<7501,1,2))) %>%
  filter(Row < 6943 | (Row > 7500 & Row < 8443)) %>%
  filter(Line != "BBP")

data <- data %>%
  mutate(Line=gsub("[ ]", "", Line)) %>%
  mutate(Line=gsub("X", "x", Line)) %>%
  mutate(Line=gsub("F1", "", Line)) 

# adding block and x y coordinates

data <- data %>%
  mutate(Block = case_when(
    Row < 6211 ~ "Block1",
    Row > 6210 & Row < 7597 ~ "Block2",
    Row > 7596 & Row < 8444 ~ "Block3"
  ),
  x_dim = case_when(
    Block == "Block1" ~ 9,
    Block == "Block2" ~ 25,
    Block == "Block3" ~ 25
  ))

data <- x_y_coords(data)

data <- data %>%
  relocate(Row,Line,Loc,Year,Rep,Block,NLB_0717,NLB_0717_note,NLB_0728,NLB_0728_note,wmd)
data <- data %>%
  arrange(Row)

data$Line <- gsub("IBM", "M0", data$Line)

addWorksheet(wb, sheetName = "NLB 2024 CL")
writeDataTable(wb, sheet = 3, x = data, colNames = TRUE, rowNames = FALSE)

### NLB 2024 Illinois

data1 <- read.csv("data/2024-07-24-01-22-39_24NP_ForTablet_table.csv")

data1 <- data1 %>% 
  select(Row,Range,Rep,StandCount,Seed.Name, DiseasePCTRating1)
colnames(data1) <- c("Row", "Range", "Rep", "Stand", "Line", "NLB_0724")

data2 <- read.csv("data/2024-08-01-07-11-09_24NP_plots_table.csv")
data2 <- data2 %>% 
  select(Row,Range,Rep,Plot.Name,NLB.080124.ah)
colnames(data2) <- c("Row", "Range", "Rep", "Line", "NLB_0801")

data3 <- read.csv("data/2024-08-02-10-30-19_24NP_plots_table.csv")
data3 <- data3 %>% 
  select(Row,Range,Rep,Plot.Name,NLB.080224.ah)
colnames(data3) <- c("Row", "Range", "Rep", "Line", "NLB_0802")

data4 <- read.csv("data/2024-08-30-11-56-28_24NP_ForTablet_database.csv")
data4 <- data4 %>% 
  select(Row,Range,Rep,Plot.Name,value,timestamp)

data5 <- data4 %>%
  filter(grepl("08-09",timestamp))
data4 <- data4 %>%
  filter(grepl("08-08",timestamp))
colnames(data4) <- c("Row", "Range", "Rep", "Line", "NLB_0808")
colnames(data5) <- c("Row", "Range", "Rep", "Line", "NLB_0809")

data <- merge(data1,data2,by=c("Line","Row","Range","Rep"), all=TRUE)
data <- merge(data,data3,by=c("Line","Row","Range","Rep"), all=TRUE)
data <- merge(data,data4,by=c("Line","Row","Range","Rep"), all=TRUE)
data <- merge(data,data5,by=c("Line","Row","Range","Rep"), all=TRUE)

data <- data %>% 
  select(-c(NA.x,NA.y))

data <- data %>%
  mutate(Loc="IL") %>%
  mutate(Year = 2024) %>%
  mutate(NLB_0801_note = gsub("[^a-zA-Z]", "", NLB_0801)) %>%
  mutate(NLB_0802_note = gsub("[^a-zA-Z]", "", NLB_0802)) %>%
  mutate(NLB_0801 = (gsub("[a-zA-Z]", "", NLB_0801))) %>%
  mutate(NLB_0801 = as.numeric(gsub(":", "", NLB_0801))) %>%
  mutate(NLB_0802 = (gsub("[a-zA-Z]", "", NLB_0802))) %>%
  mutate(NLB_0802 = as.numeric(gsub(":", "", NLB_0802)))

data$NLB_0802[3] <- NA

nlb_cols <- data %>%
  select(contains("NLB")) %>%
  select(!contains("note"))
data$na <- apply(nlb_cols, MARGIN = 1, function(x) sum(is.na(x)))

data <- data %>%
  mutate(NLB_2_switch = case_when(
    !is.na(NLB_0801) ~ 0,
    !is.na(NLB_0802) ~ 1
  ),
  NLB_3_switch = case_when(
    !is.na(NLB_0808) ~ 0,
    !is.na(NLB_0809) ~ 1
  ))

data <- data %>%
  rowwise() %>%
  mutate(NLB_2 = sum(c(NLB_0801,NLB_0802),na.rm=TRUE),
         NLB_3 = sum(c(NLB_0808,NLB_0809),na.rm=TRUE))

data <- data %>%
  rowwise() %>%
  mutate(
    wmd=((mean(c(NLB_0724,NLB_2)))*(8+NLB_2_switch) + (mean(c(NLB_2,NLB_3)))*(7-NLB_2_switch+NLB_3_switch))/(15+NLB_3_switch)
  )

data$wmd[which(data$na>2)] <- NA

data <- data %>%
  select(-c(na,NLB_2_switch,NLB_2,NLB_3,NLB_3_switch))

data <- data %>%
  rename(x = Row,
         y = Range)

data$Row <- c(1:length(data$Line))
data$Block <- 1

data <- data %>%
  relocate(Row,Line,Loc,Year,Rep,NLB_0724,NLB_0801,NLB_0801_note,NLB_0802, NLB_0802_note,
           NLB_0808,NLB_0809,wmd,Stand,x,y)

data$Line <- gsub("_", "", data$Line)
data$Line <- gsub("IBM", "M0", data$Line)

data <- data %>%
  mutate(M0 = sub(".*x","",Line)) %>%
  mutate(parent = sub("x.*", "", Line)) %>%
  mutate(Line2 = paste(M0,"x",parent,sep=""))

index <- grepl("M0", data$M0)
data$Line[index] <- data$Line2[index]

data <- data %>% 
  select(-c(M0,parent,Line2))

addWorksheet(wb, sheetName = "NLB 2024 IL")
writeDataTable(wb, sheet = 4, x = data, colNames = TRUE, rowNames = FALSE)

### saving workbook
saveWorkbook(wb, "~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx", overwrite = TRUE)

# height
# 2022 CL
# 2023 CL
# rep 1
height_2023_cl_r1 <- read.csv("~/projects/nlb_mdh/data/2023-08-11-11-04-12_CL23NLBR1_table.csv")
height_2023_cl_r2 <- read.csv("~/projects/nlb_mdh/data/2023-08-16-10-24-12_CL23NLBR2_table.csv")
# 2024 CL
# 2024 IL

# DTA
# 2022 CL
# 2023 CL
# rep 1
dta_2023_cl_r1 <- read.csv("~/projects/nlb_mdh/data/2023-08-04-08-34-47_CL23NLBR1_table.csv")
# rep 2
dta_2023_cl_r2 <- read.csv("~/projects/nlb_mdh/data/2023-08-04-08-32-44_CL23NLBR2_table.csv)
# 2024 CL
# 2024 IL



### functions
x_y_coords <- function(dataframe){
  dataframe <- dataframe %>%
    group_by(Block) %>%
    mutate(x = Row - (floor((Row-1)/x_dim))*x_dim,
           x = x - min(x) + 1,
           y = ceiling(Row/x_dim),
           y = y - min(y) + 1)
  return(dataframe)
}

dataframe <- data.frame(Block = rep("1",12),
                        Row = 1:12,
                        x_dim = 4)

dataframe <- x_y_coords(dataframe)

dataframe <- data.frame(Block = c(rep("1",12), rep("2",12)),
                        Row = c(1:24),
                        x_dim = 4)

dataframe <- x_y_coords(dataframe)

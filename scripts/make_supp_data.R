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

addWorksheet(wb, sheetName = "NLB 2024 CL")
writeDataTable(wb, sheet = 3, x = data, colNames = TRUE, rowNames = FALSE)

### saving workbook
saveWorkbook(wb, "~/projects/nlb_mdh/data/nlb_mdh_file_s1.xlsx", overwrite = TRUE)


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

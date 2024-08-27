library(BGLR)
library(qtl)
library(qtl2)


### function to run genomic prediction

ril_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_ril_slb_2022.yaml")
bc_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_bc.yaml")

b73_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_b73_hybrids_gp.yaml")
mo17_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_mo17_hybrids_gp.yaml")


Y <- read.csv("analyses/nlb_blups_1.csv")
colnames(Y)[1] <- "wmd"

cross_train <- ril_slb

cross_test <- b73_slb

Y <- Y
reps <- 1
output_dir <- "models/model_ril_"
fold <- 10

metrics_ril <- run_gp(cross_train=ril_slb,Y=Y,reps=10,fold=5,output_dir="models/model_ril_")


metrics_b73 <- run_gp(cross_train=b73_slb,Y=Y,reps=10,fold=5,output_dir="models/model_ril_")

metrics_mo17 <- run_gp(cross_train=mo17_slb,Y=Y,reps=10,fold=5,output_dir="models/model_ril_")

metrics_b73_based_on_ril <- run_gp(cross_train=ril_slb,cross_test=b73_slb,Y=Y,reps=1,fold=1,output_dir="models/model_ril_",
                                   Dom=FALSE,sep_bc=FALSE)

metrics_mo17_based_on_ril <- run_gp(cross_train=ril_slb,cross_test=mo17_slb,Y=Y,reps=1,fold=1,output_dir="models/model_ril_",
                                   Dom=FALSE,sep_bc=FALSE)


metrics_b73_w_ril <- run_gp(cross_train=ril_slb,cross_test=b73_slb,Y=Y,reps=10,output_dir="models/model_ril_", 
                            fold = 5, sep_bc=FALSE)

metrics_b73_w_ril_da <- run_gp(cross_train=ril_slb,cross_test=b73_slb,Y=Y,reps=10,output_dir="models/model_ril_", 
                            fold = 5, Dom=TRUE, sep_bc=FALSE)

metrics_mo17_w_ril <- run_gp(cross_train=ril_slb,cross_test=mo17_slb,Y=Y,reps=10,output_dir="models/model_ril_", 
                            fold = 5,sep_bc=FALSE)

metrics_mo17_w_ril_da <- run_gp(cross_train=ril_slb,cross_test=mo17_slb,Y=Y,reps=10,output_dir="models/model_ril_", 
                             fold = 5, Dom = TRUE)

metrics_all_pops <- run_gp(cross_train=ril_slb,cross_test=bc_slb,Y=Y,reps=10,output_dir="models/model_ril_", 
                             fold = 5, sep_bc=TRUE)

metrics_all_pops_da <- run_gp(cross_train=ril_slb,cross_test=bc_slb,Y=Y,reps=10,output_dir="models/model_ril_", 
                           fold = 5, Dom = TRUE, sep_bc=TRUE)

cross_test=bc_slb


### function to make ETA

get_marker_matrix <- function(cross,matrix_type){
  X <- do.call(cbind,cross$geno)
  if(matrix_type=="additive"){
    X[X==0] <- NA
    X[X==1] <- 0
    X[X==3] <- 1
  } else if(matrix_type=="dominant"){
    X[X==0] <- 0
    X[X==2] <- 0
    X[X==1] <- 1
  }
  X[is.na(X)] <- 0
  return(X)
}

match_phenos_to_matrix <- function(Y,X){
  y <- Y %>%
    filter(Line %in% rownames(X))
  index <- match(y$Line,rownames(X))
  X <- X[index,]
  return(list(y,X))
}

make_folds <- function(y,fold){
  lines <- unique(y$Line)
  tst.pop <- rep(1:fold,length=length(lines))
  # set.seed(1234)
  tst.pop <- sample(tst.pop,length(lines),FALSE)
  tst.pop <- as.data.frame(tst.pop)
  tst.pop$Line <- unique(y$Line)
  index <- match(y$Line,tst.pop$Line)
  tst.pop <- tst.pop$tst.pop[index]
  y.trn <- y
  y.trn$wmd[tst.pop==1] <- NA
  return(y.trn)
}

make_eta <- function(matrices,models){
  ETA <- list()
  ETA_names <- c()
  for(i in 1:length(matrices)){
    ETA[[i]] <- c(list(X=matrices[[i]], model = models[[i]], saveEffects=FALSE))
    ETA_names <- c(ETA_names, paste("PED", i, sep =""))
  }
  names(ETA) <- ETA_names
  return(ETA)
}

ETA_2 <- make_eta(list(x,d),list("BRR","BRR"))

ETA_1 <- list(PED1=list(X=x, model="BRR",saveEffects=FALSE),
              PED2=list(X=d, model="BRR",saveEffects=FALSE))

run_gp <- function(cross_train,cross_test=NULL,Y,output_dir,reps, fold, Dom=FALSE, test_type="bc", sep_bc=FALSE){
  metrics <- data.frame(pa=NA,rmse=NA)
  if(sep_bc==TRUE){
    metrics <- data.frame(pa=NA,rmse=NA,pa_b73=NA,rmse_b73=NA, pa_mo17=NA, rmse_mo17=NA)
  }
  for(i in 1:reps){
    X <- do.call(cbind,cross_train$geno)
    X[X==0] <- NA
    X[X==1] <- 0
    X[X==3] <- 1
    
    rownames(X) <- gsub(" X ", "x", rownames(X))
    
    y <- Y %>%
      filter(Line %in% rownames(X))
    index <- match(y$Line,rownames(X))
    X <- X[index,]
    X[is.na(X)] <- 0
    
    
    
    if(is.null(cross_test)){
      lines <- unique(y$Line)
      x=scale(X)/sqrt(ncol(X))
      tst.pop <- rep(1:fold,length=length(lines))
      # set.seed(1234)
      tst.pop <- sample(tst.pop,length(lines),FALSE)
      tst.pop <- as.data.frame(tst.pop)
      tst.pop$Line <- unique(y$Line)
      index <- match(y$Line,tst.pop$Line)
      tst.pop <- tst.pop$tst.pop[index]
      
      
      
      
      y.trn <- y
      y.trn$wmd[tst.pop==1] <- NA
      
      ## Ridge Regression
      ETA_1 <- list(PED1=list(X=x, model="BRR",saveEffects=FALSE))
      fit <- BGLR(y.trn$wmd, ETA=ETA_1, thin = 5, 
                  verbose = FALSE, saveAt = paste0(output_dir,"_",i))
    
      pa <- cor(fit$yHat[tst.pop==1],y$wmd[tst.pop==1],use = "complete.obs")
      rmse <- mean((fit$yHat[tst.pop==1]-y$wmd[tst.pop==1])^2,na.rm=TRUE)
      df <- data.frame(pa=pa,rmse=rmse)
      metrics <- rbind(metrics,df)
    } else{
      y.trn <- y
      
      X_test <- do.call(cbind,cross_test$geno)
      
      if(test_type=="bc"){
        X_test[X_test==0] <- NA
        X_test[X_test==1] <- 0
        X_test[X_test==3] <- 1
      } else{
        X_test[X_test==0] <- NA
        X_test[X_test==1] <- 0
        X_test[X_test==3] <- 1
      }
      
      rownames(X_test) <- gsub(" X ", "x", rownames(X_test))
      
      if(Dom==TRUE){
        D <- X
        D[D==0] <- 0
        D[D==2] <- 0
        D[D==1] <- 1
        D_test <- X_test
        D_test[D_test==0] <- 0
        D_test[D_test==2] <- 0
        D_test[D_test==1] <- 1
      }
      
      
      y2 <- Y %>%
        filter(Line %in% rownames(X_test))
      index.tst <- match(y2$Line,rownames(X_test))
      X_test <- X_test[index.tst,]
      
      if(Dom==TRUE){
        D_test <- D_test[index.tst,]
      }
      
      lines <- unique(y2$Line)
      
      tst.pop <- rep(1:fold,length=length(lines))
      # set.seed(1234)
      tst.pop <- sample(tst.pop,length(lines),FALSE)
      tst.pop <- as.data.frame(tst.pop)
      tst.pop$Line <- unique(y2$Line)
      index <- match(y2$Line,tst.pop$Line)
      tst.pop <- tst.pop$tst.pop[index]
      
      
      
      
      y2.trn <- y2
      y2.tst <- y2
      
      if(fold>1){
        y2.trn$wmd[tst.pop==1] <- NA
        y.trn.wmd <- c(y.trn$wmd, y2.trn$wmd)
      } else{
        y.trn.wmd <- c(y.trn$wmd, rep(NA,length(y2.tst$wmd)))
      }
      y.tst.wmd <- c(y.trn$wmd,y2.tst$wmd)
      
      x=rbind(X,X_test)
      x=scale(x)/sqrt(ncol(x))
      x[is.na(x)] <- 0
      
      
      if(Dom==TRUE){
        d=rbind(D,D_test)
        d=scale(d)/sqrt(ncol(d))
        d[is.na(d)] <- 0
      }
      ## Ridge Regression
      if(Dom==TRUE){
        ETA_1 <- list(PED1=list(X=x, model="BRR",saveEffects=FALSE),
                      PED2=list(X=d, model="BRR",saveEffects=FALSE))
      } else{
        ETA_1 <- list(PED1=list(X=x, model="BRR",saveEffects=FALSE))
      }
      
      fit <- BGLR(y.trn.wmd, ETA=ETA_1, thin = 5, 
                  verbose = FALSE, saveAt = paste0(output_dir,"_",i))
      if(fold > 1){
        tst.pop <- c(rep(10,length(y$wmd)),tst.pop)
        pa <- cor(fit$yHat[tst.pop==1],y.tst.wmd[tst.pop==1],use = "complete.obs")
        rmse <- mean((fit$yHat[tst.pop==1]-y.tst.wmd[tst.pop==1])^2,na.rm=TRUE)
        if(sep_bc==TRUE){
          b73_index <- grepl("B73", rownames(x))
          mo17_index <- grepl("Mo17", rownames(x))
          pa_b73 <- cor(fit$yHat[tst.pop==1 & b73_index==TRUE],y.tst.wmd[tst.pop==1 & b73_index==TRUE],use = "complete.obs")
          rmse_b73 <- mean((fit$yHat[tst.pop==1 & b73_index==TRUE]-y.tst.wmd[tst.pop==1 & b73_index==TRUE])^2,na.rm=TRUE)
          pa_mo17 <- cor(fit$yHat[tst.pop==1 & mo17_index==TRUE],y.tst.wmd[tst.pop==1 & mo17_index==TRUE],use = "complete.obs")
          rmse_mo17 <- mean((fit$yHat[tst.pop==1 & mo17_index==TRUE]-y.tst.wmd[tst.pop==1 & mo17_index==TRUE])^2,na.rm=TRUE)
        }
      } else{
        pa <- cor(fit$yHat[is.na(y.trn.wmd)],y.tst$wmd,use = "complete.obs")
        rmse <- mean((fit$yHat[is.na(y.trn.wmd)]-y.tst$wmd)^2,na.rm=TRUE)
      }
      if(sep_bc==FALSE){
        df <- data.frame(pa=pa,rmse=rmse)
        metrics <- rbind(metrics,df)
      } else{
        df <- data.frame(pa=pa,rmse=rmse,pa_b73=pa_b73,rmse_b73=rmse_b73,pa_mo17=pa_mo17,rmse_mo17=rmse_mo17)
        metrics <- rbind(metrics,df)
      }
      
    }
  
    
  }
  return(metrics)
}

### RILs

ril_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_ril_slb_2022.yaml")
ril_X <- do.call(cbind,ril_slb$geno)
ril_X[ril_X==0] <- NA
ril_X[ril_X==1] <- 0
ril_x=scale(ril_X)/sqrt(ncol(ril_X))

Y <- read.csv("analyses/nlb_blups_1.csv")
colnames(Y)[1] <- "wmd"

ril_y <- Y %>%
  filter(Line %in% rownames(ril_x))
ril_index <- match(ril_y$Line,rownames(ril_x))
ril_x <- ril_x[ril_index,]
ril_x[is.na(ril_x)] <- 0


ril_lines <- unique(ril_y$Line)
ril_tst.pop <- rep(1:10,length=length(ril_lines))
# set.seed(1234)
ril_tst.pop <- sample(ril_tst.pop,length(ril_lines),FALSE)
ril_tst.pop <- as.data.frame(ril_tst.pop)
ril_tst.pop$Line <- unique(ril_y$Line)
ril_index <- match(ril_y$Line,ril_tst.pop$Line)
ril_tst.pop <- ril_tst.pop$tst.pop[ril_index]


ril_pa <- rep(NA) #vector to collect the predictive ability
ril_rmse <- rep(NA) #vector to collect the RMSE

ril_y.trn <- y
ril_y.trn$wmd[ril_tst.pop==1] <- NA

## Ridge Regression
ETA_1 <- list(PED1=list(X=ril_x, model="BRR",saveEffects=FALSE))
fit <- BGLR(ril_y.trn$wmd, ETA=ETA_1, thin = 5, 
            verbose = FALSE, saveAt = paste0("models/model_ril_"))
ril_pa[1] <- cor(fit$yHat[tst.pop==1],y$wmd[tst.pop==1],use = "complete.obs")
ril_rmse[1] <- mean((fit$yHat[tst.pop==1]-y$wmd[tst.pop==1])^2,na.rm=TRUE)

### B73 BC

b73_rqtl_cross <- readRDS("~/projects/mdh_qtl/data/b73_rqtl_cross.RDS")
b73_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_b73_hybrids_generic.yaml")
X <- do.call(cbind,b73_slb$geno)
rownames(X) <- gsub(" X ", "x", rownames(X))
X[X==0] <- NA
X[X==1] <- 0
x=scale(X)/sqrt(ncol(X))

y <- Y %>%
  filter(Line %in% rownames(x))
index <- match(y$Line,rownames(x))
x <- x[index,]
x[is.na(x)] <- 0


lines <- unique(y$Line)
tst.pop <- rep(1:10,length=length(lines))
# set.seed(1234)
tst.pop <- sample(tst.pop,length(lines),FALSE)
tst.pop <- as.data.frame(tst.pop)
tst.pop$Line <- unique(y$Line)
index <- match(y$Line,tst.pop$Line)
tst.pop <- tst.pop$tst.pop[index]


pa <- rep(NA,2) #vector to collect the predictive ability
rmse <- rep(NA,2) #vector to collect the RMSE

y.trn <- y
y.trn$wmd[tst.pop==1] <- NA

## Ridge Regression
ETA_1 <- list(PED1=list(X=x, model="BRR",saveEffects=FALSE))
fit <- BGLR(y.trn$wmd, ETA=ETA_1, thin = 5, 
            verbose = FALSE, saveAt = paste0("models/model_1"))
pa[1] <- cor(fit$yHat[tst.pop==1],y$wmd[tst.pop==1],use = "complete.obs")
rmse[1] <- mean((fit$yHat[tst.pop==1]-y$wmd[tst.pop==1])^2,na.rm=TRUE)

### Mo17 BC

mo17_rqtl_cross <- readRDS("~/projects/mdh_qtl/data/mo17_rqtl_cross.RDS")
mo17_slb <- read_cross2(file = "~/projects/mdh_qtl/data/control_file_mo17_hybrids_generic.yaml")
X <- do.call(cbind,mo17_slb$geno)
rownames(X) <- gsub(" X ", "x", rownames(X))
X[X==0] <- NA
X[X==1] <- 0
x=scale(X)/sqrt(ncol(X))

y <- Y %>%
  filter(Line %in% rownames(x))
index <- match(y$Line,rownames(x))
x <- x[index,]
x[is.na(x)] <- 0


lines <- unique(y$Line)
tst.pop <- rep(1:10,length=length(lines))

pa <- rep(NA,100) #vector to collect the predictive ability
rmse <- rep(NA,100) #vector to collect the RMSE

for(i in 1:100){
  # set.seed(1234)
  tst.pop <- sample(tst.pop,length(lines),FALSE)
  tst.pop <- as.data.frame(tst.pop)
  tst.pop$Line <- unique(y$Line)
  index <- match(y$Line,tst.pop$Line)
  tst.pop <- tst.pop$tst.pop[index]
  
  
  
  
  y.trn <- y
  y.trn$wmd[tst.pop==1] <- NA
  
  ## Ridge Regression
  ETA_1 <- list(PED1=list(X=x, model="BRR",saveEffects=FALSE))
  fit <- BGLR(y.trn$wmd, ETA=ETA_1, thin = 5, 
              verbose = FALSE, saveAt = paste0("models/model_1"))
  pa[i] <- cor(fit$yHat[tst.pop==1],y$wmd[tst.pop==1],use = "complete.obs")
  rmse[i] <- mean((fit$yHat[tst.pop==1]-y$wmd[tst.pop==1])^2,na.rm=TRUE)
}






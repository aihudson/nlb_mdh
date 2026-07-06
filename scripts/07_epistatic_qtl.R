require(qtl)
require(parallel)
require(snow)
require(dplyr)
# setwd("~/projects/mdh_qtl")

input_file <- "analyses/RIL_cross.csv"
genotype <- c("A", "B")
alleles <- c("A", "B")
na.strings <- "-"
phenos <- 1
output_dir <- "analyses/qtl_analyses/"
output_name <- "RIL"
crosstype <- "ril"
phenos <- c("NLB_WMD_BLUP")
permutations <- 100

ril=FALSE
if(crosstype=="ril"){
  ril=TRUE
  crosstype="bc"
}

cross <- read.cross(format = "csv",
                    file = input_file,
                    genotype = genotype,
                    alleles = alleles,
                    na.strings = na.strings,
                    crosstype = crosstype)

cross <- jittermap(cross)
if(ril==TRUE){
  cross <- convert2riself(cross)
}

# to make sure all phenotypes are being treated as numeric
cross$pheno <- cross$pheno %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), as.numeric))


cores = 4
# permutations

results_list <- list()


# traits <- c("SLB_WMD_BLUE_MPH", "GLS_WMD_BLUE_MPH", "DTA_BLUE_MPH", "PH_BLUE_MPH", "EH_BLUE_MPH")
pheno.cols <- match(phenos, colnames(cross$pheno))
out2 <- scantwo(cross, n.cluster=cores, pheno.col = pheno.cols, verbose = TRUE)

scantwo_temp_file <- paste(output_dir, "tmp/", output_name, ".scantwo.tmp.RDS", sep = "")
saveRDS(out2, file = scantwo_temp_file)

# permutations

batchnumber <- permutations / 10

for(i in phenos){
  pheno.col <- match(i, colnames(cross$pheno))
  for(j in 1:batchnumber){
    batch_filename <- paste(output_dir, "tmp/", output_name,".trait.", i, ".batch.", j, ".scantwo.perm.tmp", sep = "")
    if(!file.exists(batch_filename)){
      results <- scantwo(cross, n.cluster=cores,n.perm=10, pheno.col = i, verbose = TRUE)
      saveRDS(results, file = batch_filename)
      print(paste(i, " batch ", j, " is done", sep=""))
    }
  }
}

perm_list <- list()

# read in temp files
mqm_auto <- readRDS(mqm_temp_file)

perm_results <- list()

for(i in phenos){
  pheno_name <- colnames(augmentedcross$pheno)[i]
  for(j in 1:batchnumber){
    batch_filename <- paste(output_dir, "tmp/", output_name,".trait.", pheno_name, ".batch.", j, ".perm.tmp", sep = "")
    perm_list[[j]] <- readRDS(batch_filename)
  }
  processed_perms <- lapply(perm_list, mqmprocesspermutation)
  combined_perms <- do.call(c, processed_perms)
  perm_results[[i]] <- combined_perms
}


output <- list(scan = mqm_auto, permutations = perm_results)

output_file <- paste(output_dir, output_name, ".RDS", sep="")

saveRDS(output, file = output_file)


require(qtl)
require(parallel)
require(snow)
require(dplyr)


# args <- commandArgs(trailingOnly = TRUE)
# input_file <- args[1]
input_file <- "analyses/RIL_cross.csv"
input_file <- "analyses/old_nlb_data_cross.csv"
input_file <- "analyses/env_blups_cross.csv"
# input_file <- "analyses/B73_cross.csv"
# input_file <- "data/RIL_all_phenos.csv"
# genotype <- args[2]
# genotype <- unlist(strsplit(genotype, ","))
# genotype <- c("AA", "AB")
genotype <- c("A", "B")
# alleles <- args[3]
# alleles <- unlist(strsplit(alleles, ","))
alleles <- c("A", "B")
# na.strings <- args[4]
# na.strings <- unlist(strsplit(na.strings, ","))
# na.strings <- "A-"
na.strings <- "-"
# crosstype <- args[5]
# crosstype <- "bc"
crosstype <- "ril"
phenos <- args[6]
phenos <- "1,5"
phenos <- as.numeric(unlist(strsplit(phenos, ",")))
permutations <- args[7]
# permutations <- 1000
output_dir <- args[8]
# output_dir <- "analyses/qtl_analyses/"
# output_name <- "RIL"
# output_name <- "B73_BC"

# for RILs
input_file <- "analyses/RIL_cross.csv"
genotype <- c("A", "B")
alleles <- c("A", "B")
na.strings <- "-"
phenos <- 1
permutations <- 1000
output_dir <- "analyses/qtl_analyses/"
output_name <- "RIL"

# for B73 BC
input_file <- "analyses/B73_cross.csv"
genotype <- c("AA", "AB")
alleles <- c("A", "B")
na.strings <- "A-"
crosstype <- "bc"
phenos <- "1,5"
phenos <- as.numeric(unlist(strsplit(phenos, ",")))
permutations <- 1000
output_dir <- "analyses/qtl_analyses/"
output_name <- "B73_BC"

# for Mo17 BC
input_file <- "analyses/Mo17_cross.csv"
genotype <- c("AB", "BB")
alleles <- c("A", "B")
na.strings <- "-B"
crosstype <- "bc"
phenos <- "1,5"
phenos <- as.numeric(unlist(strsplit(phenos, ",")))
permutations <- 1000
output_dir <- "analyses/qtl_analyses/"
output_name <- "Mo17_BC"

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


augmentedcross <- mqmaugment(cross, minprob=0.1)

# automatic backwords selection
autocofactors <- mqmautocofactors(augmentedcross, 50)

cores <- detectCores()-3

mqm_auto <- mqmscan(augmentedcross, autocofactors, pheno.col = phenos, multicore = TRUE, n.clusters = cores)

mqm_temp_file <- paste(output_dir, "tmp/", output_name, ".mqm.tmp.RDS", sep = "")
saveRDS(mqm_auto, file = mqm_temp_file)

# permutations


batchnumber <- permutations / 10

for(i in phenos){
  pheno_name <- colnames(augmentedcross$pheno)[i]
  for(j in 1:batchnumber){
    batch_filename <- paste(output_dir, "tmp/", output_name,".trait.", pheno_name, ".batch.", j, ".perm.tmp", sep = "")
    if(!file.exists(batch_filename)){
      results <- mqmpermutation(augmentedcross, scanfunction=mqmscan, cofactors=autocofactors,
                                n.cluster=cores, n.perm=10, pheno.col = i)
      saveRDS(results, file = batch_filename)
      print(paste(pheno_name, " batch ", j, " is done", sep=""))
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

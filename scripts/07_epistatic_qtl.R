require(qtl)
require(parallel)
require(snow)
require(dplyr)

# Rscript scripts/07_epistatic_qtl.R RIL
# Rscript scripts/07_epistatic_qtl.R B73_BC
# Rscript scripts/07_epistatic_qtl.R Mo17_BC
# Rscript scripts/07_epistatic_qtl.R B73_BC 100 6

args <- commandArgs(trailingOnly = TRUE)
population <- args[1]
permutations <- ifelse(length(args) >= 2, as.numeric(args[2]), 100)
cores <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

presets <- list(
  RIL = list(
    input_file = "analyses/RIL_cross.csv",
    genotype = c("A", "B"),
    na.strings = "-",
    crosstype = "ril",
    phenos = 1
  ),
  B73_BC = list(
    input_file = "analyses/B73_cross.csv",
    genotype = c("AA", "AB"),
    na.strings = "A-",
    crosstype = "bc",
    phenos = c(1, 5)
  ),
  Mo17_BC = list(
    input_file = "analyses/Mo17_cross.csv",
    genotype = c("AB", "BB"),
    na.strings = "-B",
    crosstype = "bc",
    phenos = c(1, 5)
  )
)

preset <- presets[[population]]

input_file <- preset$input_file
genotype <- preset$genotype
alleles <- c("A", "B")
na.strings <- preset$na.strings
crosstype <- preset$crosstype
phenos <- preset$phenos
output_dir <- "analyses/qtl_analyses/"
output_name <- population

ril <- FALSE
if (crosstype == "ril") {
  ril <- TRUE
  crosstype <- "bc"
}

cross <- read.cross(format = "csv",
                    file = input_file,
                    genotype = genotype,
                    alleles = alleles,
                    na.strings = na.strings,
                    crosstype = crosstype)

cross <- jittermap(cross)
if (ril == TRUE) {
  cross <- convert2riself(cross)
}

# to make sure all phenotypes are being treated as numeric
cross$pheno <- cross$pheno %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), as.numeric))

# main scan (resumable)
scantwo_temp_file <- paste0(output_dir, "tmp/", output_name, ".scantwo.tmp.RDS")
if (file.exists(scantwo_temp_file)) {
  out2 <- readRDS(scantwo_temp_file)
} else {
  out2 <- scantwo(cross, n.cluster = cores, pheno.col = phenos, verbose = TRUE)
  saveRDS(out2, file = scantwo_temp_file)
}

# permutations (resumable, batched by trait)
batchnumber <- permutations / 10

for (i in phenos) {
  pheno_name <- colnames(cross$pheno)[i]
  for (j in 1:batchnumber) {
    batch_filename <- paste0(output_dir, "tmp/", output_name, ".trait.", pheno_name, ".batch.", j, ".scantwo.perm.tmp")
    if (!file.exists(batch_filename)) {
      results <- scantwo(cross, n.cluster = cores, n.perm = 10, pheno.col = i, verbose = TRUE)
      saveRDS(results, file = batch_filename)
      print(paste(pheno_name, " batch ", j, " is done", sep = ""))
    }
  }
}

# combine permutation batches and save final output
perm_results <- list()
for (i in phenos) {
  pheno_name <- colnames(cross$pheno)[i]
  perm_list <- lapply(1:batchnumber, function(j)
    readRDS(paste0(output_dir, "tmp/", output_name, ".trait.", pheno_name, ".batch.", j, ".scantwo.perm.tmp")))
  perm_results[[pheno_name]] <- Reduce(c, perm_list)
}

output <- list(scan = out2, permutations = perm_results)

output_file <- paste0(output_dir, output_name, "_scantwo.RDS")

saveRDS(output, file = output_file)

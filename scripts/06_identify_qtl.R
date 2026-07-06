library(qtl)
library(qtl2)
library(qtl2convert)
library(dplyr)







#output_file <- paste(output_dir, output_name, ".RDS", sep="")

#output <- readRDS(output_file)


read_cross <- function(input_file, genotype, alleles, na.strings, crosstype){
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
  return(cross)
}

get_map <- function(cross){
  cross2 <- convert2cross2(cross)
  map <- insert_pseudomarkers(cross2$gmap, step=1)
  return(map)
}

get_thr <- function(perms){
  # perms is output$permutations: a list indexed by phenotype column, with a
  # NULL entry for any column that wasn't scanned. Return a named vector of
  # per-trait thresholds keyed by the "LOD <trait>" lod column name so each
  # scanned phenotype gets its own significance threshold.
  perms <- perms[!vapply(perms, is.null, logical(1))]
  thr <- vapply(perms, function(pp) summary(pp)[[1]], numeric(1))
  names(thr) <- vapply(perms, function(pp) colnames(pp)[1], character(1))
  thr
}

# mqm_auto <- output[[1]]

get_main_effect_peaks <- function(mqm_auto, map, thr, peakdrop=1.8, drop=1.5){
  # A single-phenotype mqmscan returns one scanone; scanning multiple
  # phenotypes returns an mqmmulti (a list of scanone objects) that find_peaks
  # can't align. Normalize to a list and run find_peaks on each element.
  scans <- if(inherits(mqm_auto, "mqmmulti")) mqm_auto else list(mqm_auto)
  blups_peaks <- do.call(rbind, lapply(scans, function(s){
    # match this phenotype's scan to its own permutation threshold by the
    # "LOD <trait>" column name
    lodcol <- grep("^LOD ", colnames(s), value = TRUE)[1]
    if(is.na(lodcol)){
      stop("no 'LOD <trait>' column found in scan; cannot match a threshold")
    }
    if(!lodcol %in% names(thr)){
      stop("no permutation threshold for '", lodcol, "'; available: ",
           paste(names(thr), collapse = ", "))
    }
    trait_thr <- thr[[lodcol]]
    find_peaks(scan1_output = s, map = map, threshold = trait_thr, peakdrop = peakdrop, drop = drop, expand2markers = FALSE)
  }))
  blups_peaks <- blups_peaks %>%
    filter(grepl("LOD ", lodcolumn))
  blups_peaks <- blups_peaks %>%
    rename(trait = "lodcolumn")
}

# blups_peaks <- get_main_effect_peaks(mqm_auto, map, thr)

get_main_effect_peaks_wrapper <- function(output_dir=output_dir,output_name=output_name,genotype,alleles,na.strings,crosstype){
  output_file <- paste(output_dir, output_name, ".RDS", sep="")
  output <- readRDS(output_file)
  cross <- read_cross(input_file=input_file, genotype = genotype, alleles = alleles, na.strings = na.strings, crosstype = crosstype)
  map <- get_map(cross)
  thr <- get_thr(output$permutations)
  blups_peaks <- get_main_effect_peaks(output$scan, map, thr)
  blups_peaks$cross <- output_name
  return(blups_peaks)
}

# initialize data frame
peaks_combined <- data.frame()
# for RILs
input_file <- "analyses/RIL_cross.csv"
genotype <- c("A", "B")
alleles <- c("A", "B")
na.strings <- "-"
phenos <- 1
output_dir <- "analyses/qtl_analyses/"
output_name <- "RIL"
crosstype <- "ril"

blups_peaks <- get_main_effect_peaks_wrapper(output_dir=output_dir, output_name=output_name, genotype = genotype, alleles = alleles, na.strings = na.strings, crosstype = crosstype)

peaks_combined <- rbind(peaks_combined, blups_peaks)

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

blups_peaks <- get_main_effect_peaks_wrapper(output_dir, output_name,genotype = genotype,alleles = alleles,na.strings = na.strings,crosstype = crosstype)

peaks_combined <- rbind(peaks_combined, blups_peaks)

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

blups_peaks <- get_main_effect_peaks_wrapper(output_dir, output_name,genotype = genotype,alleles = alleles,na.strings = na.strings,crosstype = crosstype)

peaks_combined <- rbind(peaks_combined, blups_peaks)

peaks_file <- "analyses/main_effect_peaks.csv"
write.csv(peaks_combined, peaks_file)





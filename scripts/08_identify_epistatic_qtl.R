require(qtl)

# Rscript scripts/08_identify_epistatic_qtl.R RIL
# Rscript scripts/08_identify_epistatic_qtl.R B73_BC
# Rscript scripts/08_identify_epistatic_qtl.R Mo17_BC
# Rscript scripts/08_identify_epistatic_qtl.R RIL 0.05

args <- commandArgs(trailingOnly = TRUE)
population <- args[1]
alpha <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.05)
near_alpha <- 0.20
min_sep_cM <- 20

qtl_analyses_dir <- "analyses/qtl_analyses/"
output_file <- "analyses/epistatic_peaks.csv"

input <- readRDS(paste0(qtl_analyses_dir, population, "_scantwo.RDS"))
scan <- input$scan
permutations <- input$permutations

peaks_combined <- data.frame()

for (pheno_name in names(permutations)) {
  perms <- permutations[[pheno_name]]
  idx <- match(pheno_name, names(permutations))
  s <- tryCatch(subset(scan, lodcolumn = idx), error = function(e) scan)

  top <- summary(s, thresholds = c(0, 0, 0, 0, 0), what = "int")
  if (nrow(top) == 0) next

  int_p <- sapply(top$lod.int, function(x) mean(as.numeric(perms$int) >= x))
  sig_level <- ifelse(int_p <= alpha, sprintf("%.2f", alpha),
                ifelse(int_p <= near_alpha, sprintf("%.2f", near_alpha), "ns"))
  same_chr_close <- (as.character(top$chr1) == as.character(top$chr2)) &
    (abs(top$pos1 - top$pos2) < min_sep_cM)

  peaks <- data.frame(
    cross = population,
    trait = pheno_name,
    chr1 = top$chr1,
    pos1 = top$pos1,
    chr2 = top$chr2,
    pos2 = top$pos2,
    lod.int = top$lod.int,
    lod.full = top$lod.full,
    int_p = int_p,
    sig_level = sig_level,
    same_chr_close = same_chr_close
  )
  peaks <- peaks[order(-peaks$lod.int), ]
  peaks_combined <- rbind(peaks_combined, peaks)
}

# this script runs per population (like 07); merge into the combined output,
# replacing any existing rows for this population
if (file.exists(output_file)) {
  existing <- read.csv(output_file, stringsAsFactors = FALSE)
  existing <- existing[existing$cross != population, ]
  peaks_combined <- rbind(existing, peaks_combined)
}

write.csv(peaks_combined, output_file, row.names = FALSE)

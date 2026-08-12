require(qtl)
require(dplyr)

# Rscript scripts/09_estimate_qtl_effects.R RIL
# Rscript scripts/09_estimate_qtl_effects.R B73_BC
# Rscript scripts/09_estimate_qtl_effects.R Mo17_BC

# sim.geno imputes missing/pseudomarker genotypes by Monte Carlo, so every
# estimate/LOD/%var here wiggles run-to-run; pin a seed and use plenty of
# draws for stable, reproducible output (mirrors 17_variance_explained.R).
set.seed(1)
n_draws <- 256

args <- commandArgs(trailingOnly = TRUE)
population <- args[1]

presets <- list(
  RIL = list(input_file = "analyses/RIL_cross.csv", genotype = c("A", "B"),
             na.strings = "-", crosstype = "ril", ga_trait = "NLB_WMD_BLUP"),
  B73_BC = list(input_file = "analyses/B73_cross.csv", genotype = c("AA", "AB"),
                na.strings = "A-", crosstype = "bc", ga_trait = "NLB_WMD_BLUP_MPH"),
  Mo17_BC = list(input_file = "analyses/Mo17_cross.csv", genotype = c("AB", "BB"),
                 na.strings = "-B", crosstype = "bc", ga_trait = "NLB_WMD_BLUP_MPH")
)

min_sep_cM <- 20
genuine_alpha <- "0.05"
main_effect_peaks_file <- "analyses/main_effect_peaks.csv"
epistatic_peaks_file <- "analyses/epistatic_peaks.csv"
output_file <- "analyses/qtl_effects.csv"

# Sign of the raw fitqtl Q-term estimate, mapped onto the B73-allele /
# resistance-BLUP convention (NLB_WMD_BLUP = 100 - wmd, so positive = B73
# allele raises resistance; see plans/identify_qtl_and_effects.md, "Effect
# sign & allele convention"):
#  - RIL: the additive-only genotype code runs B73=-1, Mo17=+1, so
#    a = -raw_est = (B73/B73 - Mo17/Mo17) / 2.
#  - B73_BC: R/qtl's internal "AB" is the true B73/Mo17 het, so
#    raw_est = het - hom already: BLUP -> d-a, MPH -> d. No flip.
#  - Mo17_BC: internal "AA" is the true B73/Mo17 het (the internal-code
#    trap), so raw_est = hom - het. For BLUP that's already -a-d as-is
#    (matches the confounded-contrast convention, no flip). For MPH it must
#    be flipped to report d = het - hom.
#  - Epistatic (aa) terms: both loci in a scantwo pair share the same
#    population/trait, so the flip is sign_mult(pop,trait)^2 = +1 always
#    (see Revision 2026-07-26, problem #1) -- callers of sign_mult() for
#    epistatic rows should not use this function; they apply no flip.
sign_mult <- function(population, trait) {
  if (population == "RIL") return(-1)
  if (population == "Mo17_BC" && trait == "NLB_WMD_BLUP_MPH") return(-1)
  1
}

estimate_type_for <- function(population, trait) {
  if (population == "RIL") return("a")
  if (trait == "NLB_WMD_BLUP_MPH") return("d")
  "confounded"
}

read_cross <- function(input_file, genotype, na.strings, crosstype) {
  ril <- crosstype == "ril"
  if (ril) crosstype <- "bc"
  cross <- read.cross(format = "csv", file = input_file, genotype = genotype,
                       alleles = c("A", "B"), na.strings = na.strings, crosstype = crosstype)
  cross <- jittermap(cross)
  if (ril) cross <- convert2riself(cross)
  cross$pheno <- cross$pheno %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), as.numeric))
  cross
}

# %var from a lone QTL (or, for an epistatic pair, that pair alone) fit
# against the null (no-QTL) model -- the marginal contrast to the joint
# model's drop-one fitqtl_pct_var: "how much this locus explains by itself,"
# not "its unique contribution controlling for the population's other QTL."
single_locus_pct_var <- function(cross, trait, chr, pos) {
  qnam <- paste0("Q", seq_along(chr))
  epinam <- if (length(chr) == 2) paste0("+", qnam[1], ":", qnam[2]) else ""
  qtlobj <- makeqtl(cross, chr = chr, pos = pos)
  formula <- as.formula(paste0("y ~ ", paste(qnam, collapse = "+"), epinam))
  pheno.col <- match(trait, colnames(cross$pheno))
  out <- tryCatch(
    fitqtl(cross, pheno.col = pheno.col, qtl = qtlobj, get.ests = FALSE, formula = formula),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_real_)
  summary(out)$result.full["Model", "%var"]
}

# fit one phenotype's main peaks together, adding Qi:Qj terms for any
# genuine epistatic partners so the main-peak estimates/LODs are adjusted
# for background epistasis; returns list(main=<one row per main peak>,
# epi=<one row per epistatic pair, with its own estimate + drop-one LOD>)
fit_main_effects <- function(cross, trait, main_peaks, epi_peaks) {
  main_n <- nrow(main_peaks)
  chrs <- main_peaks$chr
  poss <- main_peaks$pos
  qnam <- paste0("Q", seq_len(main_n))

  epinam <- ""
  epi_n <- nrow(epi_peaks)
  if (epi_n > 0) {
    for (j in seq_len(epi_n)) {
      idx1 <- main_n + (2 * j - 1)
      idx2 <- main_n + (2 * j)
      chrs <- c(chrs, epi_peaks$chr1[j], epi_peaks$chr2[j])
      poss <- c(poss, epi_peaks$pos1[j], epi_peaks$pos2[j])
      epinam <- paste0(epinam, "+Q", idx1, ":Q", idx2)
    }
  }

  qtlobj <- makeqtl(cross, chr = chrs, pos = poss)
  formula <- as.formula(paste0("y ~ ", paste(qnam, collapse = "+"), epinam))
  pheno.col <- match(trait, colnames(cross$pheno))
  out.qtl <- fitqtl(cross, pheno.col = pheno.col, qtl = qtlobj, get.ests = TRUE, formula = formula)
  s <- summary(out.qtl)

  # only Qi:Qj interaction terms are in the model for epi partners (no main
  # effect terms for those loci), so ests[-1,1]/result.drop have exactly
  # main_n + epi_n rows; identify the trailing epi_n as interaction rows by
  # ":" in the term name -- robust to R/qtl's "Q3:Q4" vs "chr@pos:chr@pos".
  ests <- as.numeric(s$ests[-1, 1])
  is_epi <- grepl(":", rownames(s$ests)[-1])

  # result.drop is NULL when there's exactly one QTL term (main_n==1, no
  # epistasis partners) -- dropping the only term is the same comparison as
  # the full-vs-null model, so fall back to the full-model LOD/%var in that
  # case.
  if (is.null(s$result.drop)) {
    main_lod <- rep(s$result.full["Model", "LOD"], main_n)
    main_pct_var <- rep(s$result.full["Model", "%var"], main_n)
    epi_lod <- numeric(0)
    epi_pct_var <- numeric(0)
  } else {
    drop_is_epi <- grepl(":", rownames(s$result.drop))
    main_lod <- as.numeric(s$result.drop[!drop_is_epi, "LOD"])[seq_len(main_n)]
    main_pct_var <- as.numeric(s$result.drop[!drop_is_epi, "%var"])[seq_len(main_n)]
    epi_lod <- as.numeric(s$result.drop[drop_is_epi, "LOD"])
    epi_pct_var <- as.numeric(s$result.drop[drop_is_epi, "%var"])
  }

  # marginal (single-locus) %var: each QTL/pair fit alone against the null
  # model, ignoring the population's other mapped QTL -- contrast to the
  # joint/adjusted fitqtl_pct_var above.
  main_single_pct_var <- vapply(seq_len(main_n), function(i) {
    single_locus_pct_var(cross, trait, main_peaks$chr[i], main_peaks$pos[i])
  }, numeric(1))
  epi_single_pct_var <- numeric(0)
  if (epi_n > 0) {
    epi_single_pct_var <- vapply(seq_len(epi_n), function(j) {
      single_locus_pct_var(cross, trait, c(epi_peaks$chr1[j], epi_peaks$chr2[j]),
                            c(epi_peaks$pos1[j], epi_peaks$pos2[j]))
    }, numeric(1))
  }

  main <- data.frame(
    raw_estimate = ests[!is_epi][seq_len(main_n)],
    fitqtl_lod = main_lod,
    fitqtl_pct_var = main_pct_var,
    single_qtl_pct_var = main_single_pct_var
  )

  epi <- data.frame()
  if (epi_n > 0) {
    epi <- data.frame(
      chr1 = epi_peaks$chr1, pos1 = epi_peaks$pos1,
      chr2 = epi_peaks$chr2, pos2 = epi_peaks$pos2,
      lod = epi_peaks$lod.int,
      raw_estimate = ests[is_epi],
      fitqtl_lod = epi_lod,
      fitqtl_pct_var = epi_pct_var,
      single_qtl_pct_var = epi_single_pct_var
    )
  }

  list(main = main, epi = epi)
}

# fit the population's real ga_trait peaks + genuine epi partners plus one
# dummy QTL at (chr, pos), dropping any real peak/epi locus within
# min_sep_cM of the dummy; returns the dummy's raw estimate + drop-one LOD.
# Same "dummy-QTL" mechanics as 11_genome_wide_effect_scan.R::dummy_qtl_effect.
dummy_qtl_effect <- function(cross, trait, chr, pos, main_peaks, epi_peaks) {
  near_dummy <- function(c, p) c == chr & abs(p - pos) < min_sep_cM

  keep_main <- main_peaks[!near_dummy(main_peaks$chr, main_peaks$pos), , drop = FALSE]
  keep_epi <- epi_peaks
  if (nrow(keep_epi) > 0) {
    drop_epi <- near_dummy(keep_epi$chr1, keep_epi$pos1) | near_dummy(keep_epi$chr2, keep_epi$pos2)
    keep_epi <- keep_epi[!drop_epi, , drop = FALSE]
  }

  main_n <- nrow(keep_main)
  dummy_idx <- main_n + 1
  chrs <- c(keep_main$chr, chr)
  poss <- c(keep_main$pos, pos)
  qnam <- paste0("Q", seq_len(main_n + 1))

  epinam <- ""
  epi_n <- nrow(keep_epi)
  if (epi_n > 0) {
    for (j in seq_len(epi_n)) {
      idx1 <- main_n + 1 + (2 * j - 1)
      idx2 <- main_n + 1 + (2 * j)
      chrs <- c(chrs, keep_epi$chr1[j], keep_epi$chr2[j])
      poss <- c(poss, keep_epi$pos1[j], keep_epi$pos2[j])
      epinam <- paste0(epinam, "+Q", idx1, ":Q", idx2)
    }
  }

  qtlobj <- makeqtl(cross, chr = chrs, pos = poss)
  formula <- as.formula(paste0("y ~ ", paste(qnam, collapse = "+"), epinam))
  pheno.col <- match(trait, colnames(cross$pheno))
  out.qtl <- tryCatch(
    fitqtl(cross, pheno.col = pheno.col, qtl = qtlobj, get.ests = TRUE, formula = formula),
    error = function(e) NULL
  )
  if (is.null(out.qtl)) {
    return(list(raw_estimate = NA_real_, fitqtl_lod = NA_real_,
                fitqtl_pct_var = NA_real_, single_qtl_pct_var = NA_real_))
  }

  s <- summary(out.qtl)
  raw_estimate <- as.numeric(s$ests[-1, 1])[dummy_idx]

  if (is.null(s$result.drop)) {
    fitqtl_lod <- s$result.full["Model", "LOD"]
    fitqtl_pct_var <- s$result.full["Model", "%var"]
  } else {
    drop_is_epi <- grepl(":", rownames(s$result.drop))
    fitqtl_lod <- as.numeric(s$result.drop[!drop_is_epi, "LOD"])[dummy_idx]
    fitqtl_pct_var <- as.numeric(s$result.drop[!drop_is_epi, "%var"])[dummy_idx]
  }

  # the dummy QTL alone, not the kept real peaks -- marginal contrast to
  # fitqtl_pct_var above.
  single_qtl_pct_var <- single_locus_pct_var(cross, trait, chr, pos)

  list(raw_estimate = raw_estimate, fitqtl_lod = fitqtl_lod,
       fitqtl_pct_var = fitqtl_pct_var, single_qtl_pct_var = single_qtl_pct_var)
}

preset <- presets[[population]]
cross <- sim.geno(read_cross(preset$input_file, preset$genotype, preset$na.strings, preset$crosstype),
                   step = 2.5, n.draws = n_draws)

main_effect_peaks <- read.csv(main_effect_peaks_file, stringsAsFactors = FALSE)
main_effect_peaks$trait <- sub("^LOD ", "", main_effect_peaks$trait)

genuine_epi <- data.frame()
if (file.exists(epistatic_peaks_file)) {
  epistatic_peaks <- read.csv(epistatic_peaks_file, stringsAsFactors = FALSE)
  genuine_epi <- epistatic_peaks %>% filter(sig_level == genuine_alpha, !same_chr_close)
}

peaks <- main_effect_peaks %>% filter(cross == population)

# real peak fits (main effect + epistatic interaction rows), one pass per
# trait this population has independently significant peaks on
main_rows <- data.frame()
epi_rows <- data.frame()
for (trait in unique(peaks$trait)) {
  trait_peaks <- peaks %>% filter(trait == !!trait)
  trait_epi <- genuine_epi %>% filter(cross == population, trait == !!trait)
  fit <- fit_main_effects(cross, trait, trait_peaks, trait_epi)

  trait_peaks$raw_estimate <- fit$main$raw_estimate
  trait_peaks$fitqtl_lod <- fit$main$fitqtl_lod
  trait_peaks$fitqtl_pct_var <- fit$main$fitqtl_pct_var
  trait_peaks$single_qtl_pct_var <- fit$main$single_qtl_pct_var
  trait_peaks$estimate_type <- estimate_type_for(population, trait)
  trait_peaks$estimate <- trait_peaks$raw_estimate * sign_mult(population, trait)
  trait_peaks$effect_class <- "main"
  trait_peaks$borrowed <- FALSE
  trait_peaks$chr2 <- NA_integer_
  trait_peaks$pos2 <- NA_real_
  main_rows <- bind_rows(main_rows, trait_peaks)

  if (nrow(fit$epi) > 0) {
    # both loci share this population/trait, so the aa-term flip is
    # sign_mult(pop,trait)^2 == +1 always -- no sign flip applied here.
    epi <- fit$epi %>%
      transmute(cross = population, trait = trait,
                chr = chr1, pos = pos1, chr2 = chr2, pos2 = pos2,
                ci_lo = NA_real_, ci_hi = NA_real_,
                lod = lod,
                estimate = raw_estimate,
                estimate_type = "epistatic", effect_class = "epistatic", borrowed = FALSE,
                fitqtl_lod = fitqtl_lod, fitqtl_pct_var = fitqtl_pct_var,
                single_qtl_pct_var = single_qtl_pct_var)
    epi_rows <- bind_rows(epi_rows, epi)
  }
}

# symmetric borrowing: at every QTL position (union across all populations'
# main-effect peaks) this population lacks an independent peak for on its
# own gene-action parameter, borrow a dummy-QTL fit of that parameter --
# RIL borrows "a" (from NLB_WMD_BLUP), both BCs borrow "d" (from
# NLB_WMD_BLUP_MPH) -- so every QTL ends up with an a/d set from every
# population, not just RIL.
union_positions <- main_effect_peaks %>%
  distinct(chr, pos, .keep_all = TRUE) %>%
  select(chr, pos, ci_lo, ci_hi)

own_ga_peaks <- peaks %>% filter(trait == preset$ga_trait)
own_ga_epi <- genuine_epi %>% filter(cross == population, trait == preset$ga_trait)

already_covered <- mapply(function(c, p) any(own_ga_peaks$chr == c & abs(own_ga_peaks$pos - p) <= 5),
                           union_positions$chr, union_positions$pos)
borrow_positions <- union_positions[!already_covered, ]

borrowed_rows <- data.frame()
if (nrow(borrow_positions) > 0) {
  fits <- lapply(seq_len(nrow(borrow_positions)), function(i) {
    dummy_qtl_effect(cross, preset$ga_trait, borrow_positions$chr[i], borrow_positions$pos[i],
                      own_ga_peaks, own_ga_epi)
  })
  borrowed_rows <- data.frame(
    cross = population,
    trait = preset$ga_trait,
    chr = borrow_positions$chr,
    pos = borrow_positions$pos,
    chr2 = NA_integer_,
    pos2 = NA_real_,
    ci_lo = borrow_positions$ci_lo,
    ci_hi = borrow_positions$ci_hi,
    lod = NA_real_,
    raw_estimate = sapply(fits, function(x) x$raw_estimate),
    fitqtl_lod = sapply(fits, function(x) x$fitqtl_lod),
    fitqtl_pct_var = sapply(fits, function(x) x$fitqtl_pct_var),
    single_qtl_pct_var = sapply(fits, function(x) x$single_qtl_pct_var),
    estimate_type = ifelse(population == "RIL", "a", "d"),
    effect_class = "main",
    borrowed = TRUE
  )
  borrowed_rows$estimate <- borrowed_rows$raw_estimate * sign_mult(population, preset$ga_trait)
}

results <- bind_rows(main_rows, epi_rows, borrowed_rows) %>%
  select(cross, trait, chr, pos, chr2, pos2, ci_lo, ci_hi, lod, estimate,
         estimate_type, effect_class, borrowed, fitqtl_lod, fitqtl_pct_var,
         single_qtl_pct_var)

# merge into the combined output, replacing any existing rows with the same
# (cross, trait, chr, pos, chr2, pos2) identity so re-running a population
# is idempotent; tolerate an old file written before chr2/pos2 or the
# pct_var columns existed.
if (file.exists(output_file)) {
  existing <- read.csv(output_file, stringsAsFactors = FALSE)
  if (!"chr2" %in% colnames(existing)) existing$chr2 <- NA_integer_
  if (!"pos2" %in% colnames(existing)) existing$pos2 <- NA_real_
  if (!"fitqtl_pct_var" %in% colnames(existing)) existing$fitqtl_pct_var <- NA_real_
  if (!"single_qtl_pct_var" %in% colnames(existing)) existing$single_qtl_pct_var <- NA_real_
  key <- function(df) {
    paste(df$cross, df$trait, df$chr, df$pos,
          ifelse(is.na(df$chr2), "", df$chr2), ifelse(is.na(df$pos2), "", df$pos2))
  }
  existing <- existing[!key(existing) %in% key(results), ]
  results <- bind_rows(existing, results)
}

write.csv(results, output_file, row.names = FALSE)

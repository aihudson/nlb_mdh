require(dplyr)

# Rscript scripts/10_qtl_gene_action.R

qtl_effects_file <- "analyses/qtl_effects.csv"
output_file <- "analyses/qtl_gene_action.csv"

qtl_effects <- read.csv(qtl_effects_file, stringsAsFactors = FALSE)

# epistatic rows carry NA CIs and describe a locus pair, not a single QTL --
# exclude them from colocalization so they don't pollute the connected-
# component grouping; a/d/confounded main-effect rows all keep real CIs.
main_effects <- qtl_effects %>% filter(effect_class == "main")

# find_overlaps: pairwise overlap of [ci_lo, ci_hi] intervals, returned as a
# data.frame(queryHits, subjectHits) -- the same shape as
# IRanges::findOverlaps(IRanges(start, end)), but in base R so this script
# has no Bioconductor dependency
find_overlaps <- function(lo, hi) {
  n <- length(lo)
  hits <- expand.grid(queryHits = seq_len(n), subjectHits = seq_len(n))
  hits <- hits[hits$queryHits <= hits$subjectHits, ]
  overlap <- lo[hits$queryHits] <= hi[hits$subjectHits] & hi[hits$queryHits] >= lo[hits$subjectHits]
  hits[overlap, ]
}

# colocalize: connected components of overlapping [ci_lo, ci_hi] within each
# chromosome, replacing the legacy hand-assigned qtl ids in overlapping_qtl.Rmd
main_effects$qtl_id <- NA_integer_
next_id <- 1L
for (this_chr in sort(unique(main_effects$chr))) {
  idx <- which(main_effects$chr == this_chr)
  ov <- find_overlaps(main_effects$ci_lo[idx], main_effects$ci_hi[idx])

  n <- length(idx)
  parent <- seq_len(n)
  find_root <- function(x) {
    while (parent[x] != x) x <- parent[x]
    x
  }
  for (k in seq_len(nrow(ov))) {
    ra <- find_root(ov$queryHits[k])
    rb <- find_root(ov$subjectHits[k])
    if (ra != rb) parent[ra] <- rb
  }
  comp <- vapply(seq_len(n), find_root, integer(1))
  local_id <- match(comp, unique(comp))

  main_effects$qtl_id[idx] <- next_id - 1L + local_id
  next_id <- next_id + length(unique(comp))
}

# per-cluster chr/pos: mean position of its member peaks
clusters <- main_effects %>%
  group_by(qtl_id) %>%
  summarise(chr = chr[1], pos = mean(pos), .groups = "drop")

# a/d values now come from real peaks *or* borrowed dummy-QTL fits (09's
# symmetric borrowing), so every qtl_id gets a value; when both a real and a
# borrowed fit exist for the same qtl_id, prefer the real one.
a_vals <- main_effects %>%
  filter(estimate_type == "a") %>%
  group_by(qtl_id) %>%
  arrange(borrowed, .by_group = TRUE) %>%
  summarise(a = estimate[1], a_sig = !borrowed[1], .groups = "drop")

b73_d <- main_effects %>%
  filter(estimate_type == "d", cross == "B73_BC") %>%
  group_by(qtl_id) %>%
  arrange(borrowed, .by_group = TRUE) %>%
  summarise(B73_d = estimate[1], B73_d_sig = !borrowed[1], .groups = "drop")

mo17_d <- main_effects %>%
  filter(estimate_type == "d", cross == "Mo17_BC") %>%
  group_by(qtl_id) %>%
  arrange(borrowed, .by_group = TRUE) %>%
  summarise(Mo17_d = estimate[1], Mo17_d_sig = !borrowed[1], .groups = "drop")

# which cross/trait combos were independently significant (excludes
# borrowed/supplemental fits) at each qtl_id
sig_in <- main_effects %>%
  filter(!borrowed) %>%
  mutate(tag = ifelse(grepl("_MPH$", trait), paste0(cross, "_MPH"), cross)) %>%
  group_by(qtl_id) %>%
  summarise(sig_in = paste(sort(unique(tag)), collapse = ","), .groups = "drop")

classify_action <- function(d_a) {
  case_when(
    is.na(d_a) ~ NA_character_,
    abs(d_a) > 1.2 & d_a < 0 ~ "ud",
    abs(d_a) < 0.2 ~ "additive",
    abs(d_a) < 0.8 ~ "pd",
    abs(d_a) <= 1.2 ~ "dominant",
    TRUE ~ "od"
  )
}

combined <- clusters %>%
  left_join(a_vals, by = "qtl_id") %>%
  left_join(b73_d, by = "qtl_id") %>%
  left_join(mo17_d, by = "qtl_id") %>%
  left_join(sig_in, by = "qtl_id") %>%
  mutate(
    `B73_d/a` = B73_d / abs(a),
    `Mo17_d/a` = Mo17_d / abs(a),
    B73_action = classify_action(`B73_d/a`),
    Mo17_action = classify_action(`Mo17_d/a`)
  ) %>%
  select(qtl_id, chr, pos, a, a_sig, B73_d, `B73_d/a`, B73_action, B73_d_sig,
         Mo17_d, `Mo17_d/a`, Mo17_action, Mo17_d_sig, sig_in) %>%
  arrange(qtl_id)

write.csv(combined, output_file, row.names = FALSE)

# ============================================================================
# 04b_key_cytokines.R
# [REVISED SPEC] "Identify top weighted cytokines based on LV1 (top 3) and on
# LV2 (top 3)."
#
# Selection is by the magnitude of the LOOCV mean loading (stage 04), per
# genotype, from that genotype's own model — top 3 on LV1 and top 3 on LV2, then
# unioned (so a genotype yields at most 6, fewer when an analyte tops both LVs).
# Magnitude, not signed value, because an LV separates groups in both directions:
# a strongly negative loading is as informative as a strongly positive one.
#
# IMPORTANT CAVEAT, carried into every figure built on this table: a PLS-DA
# loading is a weighting, not a test. Where the model itself fails its
# permutation test (stage 04, p_permute >= 0.05), the ranking describes which
# analytes the model leaned on, and nothing more — the group differences come
# from the per-analyte tests in stage 05. This is the Garza et al. structure:
# weight the profile with PLS-DA, then test individual analytes separately.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
if (!exists("combined_xlsx")) source("config.R")

processed_dir <- file.path(results_dir, "processed")
load_tbl <- read.csv(file.path(stats_dir, "plsda_loadings.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
perm_tbl <- read.csv(file.path(stats_dir, "plsda_permutation.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)

rows <- list()
for (g in genotypes) {
  d <- load_tbl[load_tbl$genotype == g, , drop = FALSE]
  if (!nrow(d)) next
  for (lv in c("LV1", "LV2")) {
    m  <- d[[paste0(lv, "_loocv_mean")]]
    s  <- d[[paste0(lv, "_loocv_sd")]]
    o  <- order(abs(m), decreasing = TRUE)[seq_len(min(key_n_per_lv, nrow(d)))]
    rows[[paste(g, lv)]] <- data.frame(
      genotype = g, LV = lv, rank = seq_along(o),
      analyte = d$analyte[o],
      loading = m[o], loading_sd = s[o],
      # a loading whose LOOCV error bar crosses zero is not stable across
      # cross-validation folds -- reported so an unstable pick is visible
      stable  = abs(m[o]) > s[o],
      is_myeloid = d$analyte[o] %in% cytokines_myeloid,
      stringsAsFactors = FALSE)
  }
}
key_tbl <- do.call(rbind, rows)
rownames(key_tbl) <- NULL

# model validity attached to every row, so no one reads the table without it
key_tbl <- merge(key_tbl,
                 perm_tbl[, c("subgroup", "loocv_accuracy", "p_permute")],
                 by.x = "genotype", by.y = "subgroup", all.x = TRUE)
key_tbl$model_valid <- key_tbl$p_permute < alpha
key_tbl <- key_tbl[order(match(key_tbl$genotype, genotypes), key_tbl$LV, key_tbl$rank), ]
write.csv(key_tbl, file.path(stats_dir, "key_cytokines.csv"), row.names = FALSE)

key_cytokines <- lapply(genotypes, function(g) {
  a <- key_tbl$analyte[key_tbl$genotype == g]
  cytokines_all[cytokines_all %in% unique(a)]   # canonical order, deduped
})
names(key_cytokines) <- genotypes
saveRDS(key_cytokines, file.path(processed_dir, "key_cytokines.rds"))

for (g in genotypes) {
  pv <- perm_tbl$p_permute[perm_tbl$subgroup == g]
  cat(sprintf("[04b] %-6s key cytokines (%d): %s\n", g, length(key_cytokines[[g]]),
              paste(key_cytokines[[g]], collapse = ", ")))
  for (lv in c("LV1", "LV2")) {
    d <- key_tbl[key_tbl$genotype == g & key_tbl$LV == lv, ]
    cat(sprintf("[04b]   %s top %d: %s\n", lv, nrow(d),
                paste(sprintf("%s (%s%s)", d$analyte, fmt(d$loading),
                              ifelse(d$stable, "", ", unstable")), collapse = "; ")))
  }
  if (length(pv) && pv >= alpha)
    cat(sprintf("[04b]   !! %s PLS-DA fails permutation (p_permute = %s): this ranking is\n",
                g, fmt(pv)),
        "         a descriptive weighting, NOT evidence of a group difference.\n")
}
cat("[04b] key cytokines ->", file.path(stats_dir, "key_cytokines.csv"), "\n")

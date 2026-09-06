# ============================================================================
# 05b_cytokine_mixed_effects.R
# [CORRECTION, initialComments.md line 14] The Wood-lab alternative: a
# multiple/mixed-effects model that tests the flicker effect while CONTROLLING
# for sex (rather than only splitting by sex). Per analyte, within genotype:
#     y ~ Flicker + Sex + (1 | litter)   [lmer; falls back to y ~ Flicker + Sex
#     when the litter variance is singular -- see lib/stats_helpers.R]
# We report the sex-controlled flicker effect, the sex main effect, and the
# Flicker x Sex interaction (i.e. whether the flicker effect itself differs by
# sex), plus sex-averaged Tukey pairwise contrasts. This complements stage 05
# (which splits by sex); together they answer "is the effect real controlling
# for sex, and does it differ between sexes?".
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")

processed_dir     <- file.path(results_dir, "processed")
animals           <- readRDS(file.path(processed_dir, "animal_table.rds"))
analysis_analytes <- readRDS(file.path(processed_dir, "analysis_analytes.rds"))
key_cytokines     <- readRDS(file.path(processed_dir, "key_cytokines.rds"))

om <- list(); ph <- list()
for (g in genotypes) {
  d <- animals[animals$Genotype == g, , drop = FALSE]
  r <- compare_flicker_mixed_many(d, analysis_analytes, subgroup = g)
  om[[g]] <- r$omnibus; ph[[g]] <- r$posthoc
}
omnibus <- do.call(rbind, om)
posthoc <- do.call(rbind, ph[!vapply(ph, is.null, logical(1))])

# DIAGNOSTIC ONLY (see 05): FDR across analytes per genotype on the sex-
# controlled flicker effect. Inference uses the emmeans-Tukey pairwise p.
omnibus$p_flicker_fdr <- NA_real_
for (g in unique(omnibus$subgroup)) {
  i <- omnibus$subgroup == g & !is.na(omnibus$p_flicker)
  omnibus$p_flicker_fdr[i] <- p.adjust(omnibus$p_flicker[i], method = "BH")
}
write_stats(omnibus, "cytokine_mixed_omnibus")
write_stats(posthoc, "cytokine_mixed_posthoc")

# how many analytes used the mixed model vs the fixed-effects fallback
cat("[05b] model used per analyte:\n")
print(table(omnibus$subgroup, sub(":.*", "", omnibus$model)))

# clear-signal summary: sex-controlled flicker effect, sex effect, interaction
cat("\n[05b] sex-CONTROLLED flicker effect (p<0.05); * = a key cytokine:\n")
sig <- omnibus[!is.na(omnibus$p_flicker) & omnibus$p_flicker < alpha,
               c("subgroup","variable","model","F_flicker","p_flicker","p_flicker_fdr")]
if (nrow(sig)) {
  sig$key <- ifelse(mapply(function(sg, v) v %in% key_cytokines[[sg]],
                           sig$subgroup, sig$variable), "*", "")
  print(sig, row.names = FALSE)
} else cat("   (none)\n")

cat("\n[05b] analytes with a sex effect or Flicker x Sex interaction (p<0.05):\n")
sx <- omnibus[(!is.na(omnibus$p_sex) & omnibus$p_sex < alpha) |
              (!is.na(omnibus$p_flickerXsex) & omnibus$p_flickerXsex < alpha),
              c("subgroup","variable","p_sex","p_flickerXsex")]
if (nrow(sx)) print(sx, row.names = FALSE) else cat("   (none)\n")
cat("\n[05b] mixed-effects (sex-controlled) cytokine stats done ->", stats_dir, "\n")

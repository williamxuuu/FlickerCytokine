# ============================================================================
# 07_correlations.R
# cytokine<->amyloid (Fig 4) and cytokine<->morphology (Fig 5) correlations.
# [CORRECTION] Per initialComments.md and the revised spec:
#   - statistic is SPEARMAN (rho + p); the trend line is Pearson/linear for
#     display only (in the scatter helper);
#   - the p-value drives which pairs are shown, NOT a rho cutoff (in biology a
#     rho of ~.46 is a perfectly reportable relationship, so no rho gate);
#   - [REVISED SPEC] FDR is applied over the PRE-SELECTED HEATMAP PANEL, by row
#     (per target across the retained analytes) -- not across all 32 analytes and
#     not blanket across the whole matrix. A key-cytokine-restricted variant is
#     written alongside as a sensitivity table;
#   - scatterplots label males vs females (different symbols) and now carry the
#     FDR-adjusted q as well as the raw p.
# The original showed RAW Pearson/Spearman p on the plots and buried FDR in CSVs.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")

processed_dir     <- file.path(results_dir, "processed")
animals           <- readRDS(file.path(processed_dir, "animal_table.rds"))
analysis_analytes <- readRDS(file.path(processed_dir, "analysis_analytes.rds"))
key_cytokines     <- readRDS(file.path(processed_dir, "key_cytokines.rds"))

# Scatter dirs are cleared before writing: a previous run's panels for analytes
# that have since been pruned would otherwise sit alongside the current ones and
# get pasted into a figure by mistake.
for (tag in c("cytokine_amyloid", "cytokine_morphology")) {
  d <- file.path(corr_dir, tag)
  if (dir.exists(d)) unlink(list.files(d, pattern = "\\.pdf$", full.names = TRUE))
}

save_scatters <- function(df, pairs, colors, sub, tag) {
  if (is.null(pairs) || !nrow(pairs)) return(invisible())
  d <- file.path(corr_dir, tag); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(nrow(pairs))) {
    p <- scatter_sex(df, xcol = pairs$target[i], ycol = pairs$analyte[i],
                     subgroup = sub, colors = colors, q = pairs$p_adj[i])
    ggplot2::ggsave(file.path(d, sprintf("%s_%s_vs_%s.pdf", sub,
                    make.names(pairs$analyte[i]), make.names(pairs$target[i]))),
                    p, width = 5, height = 4.2)
  }
  cat(sprintf("[07]   %d scatterplots -> %s/\n", nrow(pairs), tag))
}

# ----------------------------------------------------------------------------
# 1. cytokine <-> amyloid, within 5xFAD (Fig 4).
# ----------------------------------------------------------------------------
fad <- animals[animals$Genotype == "5XFAD", , drop = FALSE]
amy_targets <- intersect(c("Abeta 40", "Abeta 42"), colnames(fad))

cyt_amy <- correlate(fad, analysis_analytes, amy_targets, subgroup = "5XFAD",
                     fdr_scope = "by_row")                       # PRIMARY
cyt_amy_key <- correlate(fad, analysis_analytes, amy_targets, subgroup = "5XFAD",
                         fdr_scope = "key",
                         restrict = key_cytokines[["5XFAD"]])    # sensitivity
write_stats(cyt_amy,     "corr_cytokine_amyloid_byrow")
write_stats(cyt_amy_key, "corr_cytokine_amyloid_keyset")

# scatter every pair passing raw p (p drives inclusion); the corrected q is
# annotated on each panel so the figure can never outrun the correction.
sig_amy <- cyt_amy[!is.na(cyt_amy$p_raw) & cyt_amy$p_raw < corr_p_threshold, ]
save_scatters(fad, sig_amy, flicker_colors_fad, "5XFAD", "cytokine_amyloid")

# ----------------------------------------------------------------------------
# 2. cytokine <-> morphology, within genotype (Fig 5).
# ----------------------------------------------------------------------------
morph_here <- intersect(morph_vars, colnames(animals))
cm_all <- list(); cm_key <- list()
for (g in genotypes) {
  sub <- animals[animals$Genotype == g, , drop = FALSE]
  cm  <- correlate(sub, analysis_analytes, morph_here, subgroup = g,
                   fdr_scope = "by_row")
  cm_all[[g]] <- cm
  cm_key[[g]] <- correlate(sub, analysis_analytes, morph_here, subgroup = g,
                           fdr_scope = "key", restrict = key_cytokines[[g]])
  pal <- if (g == "WT") flicker_colors_wt else flicker_colors_fad
  sig <- cm[!is.na(cm$p_raw) & cm$p_raw < corr_p_threshold, ]
  save_scatters(sub, sig, pal, g, "cytokine_morphology")
}
write_stats(do.call(rbind, cm_all), "corr_cytokine_morphology_byrow")
write_stats(do.call(rbind, cm_key), "corr_cytokine_morphology_keyset")

# ----------------------------------------------------------------------------
# 3. Console summary: what passes p (raw) and what survives the heatmap-panel FDR.
# ----------------------------------------------------------------------------
report <- function(tbl, name) {
  if (is.null(tbl)) return(invisible())
  praw <- sum(tbl$p_raw < corr_p_threshold, na.rm = TRUE)
  pfdr <- sum(tbl$p_adj < corr_p_threshold, na.rm = TRUE)
  cat(sprintf("[07] %-30s pairs p<%.2f raw: %3d | survive panel-FDR: %3d\n",
              name, corr_p_threshold, praw, pfdr))
}
report(cyt_amy,                  "cytokine~amyloid (panel FDR)")
report(cyt_amy_key,              "cytokine~amyloid (key set)")
report(do.call(rbind, cm_all),   "cytokine~morphology (panel)")
report(do.call(rbind, cm_key),   "cytokine~morphology (key set)")

surv <- cyt_amy[!is.na(cyt_amy$p_adj) & cyt_amy$p_adj < alpha, ]
cat("[07] cytokine~amyloid surviving panel FDR:\n")
if (nrow(surv)) print(surv[, c("target","analyte","rho","p_raw","p_adj")], row.names = FALSE) else
  cat("   (none)\n")
cat("[07] correlations done ->", corr_dir, "\n")

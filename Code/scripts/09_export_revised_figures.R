# ============================================================================
# 09_export_revised_figures.R
# Corrected figures rendered to MATCH the manuscript's existing style so they
# drop into the composites coherently:
#   - PLS-DA: ellipse scores plot (Wood-lab scores.plot.ellipse), "Scores on
#     LV1/LV2", big flicker-colored points. The invalid LV-score p-value violins
#     are dropped; the Garza-style permutation result goes in the subtitle.
#   - Loadings: bars with LOOCV error bars (mean +/- SD over 1000 iterations),
#     the Wood-lab format, with the top-3-per-LV picks highlighted.
#   - Bars (amyloid, key cytokines): flicker-colored bar + SEM + jittered dots,
#     all three pairwise CORRECTED p's printed above brackets.
#   - Fig 4C: rho heatmap grid (cytokine x Abeta), red-white-blue, asterisks by
#     the heatmap-panel FDR -- the manuscript's format, not scatterplots.
# Colors follow the manuscript (WT yellow/blue/pink, 5xFAD orange/navy/purple).
# Output -> results/revised_figures_png/.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(ggpubr); library(heatmap3); library(gplots); library(ropls)
})
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")
source("lib/WoodLabFunctions.R")
source("lib/plsda_helpers.R")

pdir <- file.path(results_dir, "processed")
animals <- readRDS(file.path(pdir, "animal_table.rds"))
adj     <- readRDS(file.path(pdir, "animal_table_sexadj.rds"))
analysis_analytes <- readRDS(file.path(pdir, "analysis_analytes.rds"))
key_cytokines     <- readRDS(file.path(pdir, "key_cytokines.rds"))
fits              <- readRDS(file.path(pdir, "plsda_fits.rds"))
key_tbl <- read.csv(file.path(stats_dir, "key_cytokines.csv"), check.names = FALSE)
outdir  <- file.path(results_dir, "revised_figures_png")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
DPI <- 220
zscore <- function(x) { s <- sd(x, na.rm = TRUE); if (s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }
stars  <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))

# ---------------------------------------------------------------------------
# manuscript-style theme + palettes
# ---------------------------------------------------------------------------
pal_of <- function(g) if (g == "WT") flicker_colors_wt else flicker_colors_fad
theme_ms <- function(base = 13) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none")
}

# bar + SEM + dots, corrected pairwise p's above brackets (matches Fig 4B / Supp 5)
ms_bar <- function(d, var, g, ylab = "a.u.") {
  pal <- pal_of(g)
  ph  <- compare_flicker(d, var, subgroup = g, variable = var)$posthoc
  # base ggplot with .data[[var]] so hyphenated names (IL-1a, IP-10) are safe
  p <- ggplot(d, aes(x = Flicker, y = .data[[var]], fill = Flicker)) +
    stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.85) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(aes(shape = Sex), width = 0.12, size = 1.7, stroke = 0.4) +
    scale_fill_manual(values = pal) + scale_shape_manual(values = sex_shapes) +
    labs(x = NULL, y = ylab, title = var) + theme_ms()
  if (!is.null(ph) && nrow(ph)) {
    rng <- max(d[[var]], na.rm = TRUE) - min(0, min(d[[var]], na.rm = TRUE))
    ph$y.position <- max(d[[var]], na.rm = TRUE) + rng * (0.08 + 0.13 * (seq_len(nrow(ph)) - 1))
    ph$lab <- paste0(fmt(ph$p_adj), stars(ph$p_adj))
    p <- p + stat_pvalue_manual(ph, label = "lab", xmin = "group1", xmax = "group2",
                                tip.length = 0.01, size = 3.4, bracket.size = 0.4)
  }
  p
}

# ---------------------------------------------------------------------------
# 1. PLS-DA scores + loadings (Fig 2C/D, Supp 4B/C/D)
#    The model, its axis orientation and its LOOCV loadings all come from stage
#    04 (results/processed/plsda_fits.rds), so nothing is refit here and the
#    figures cannot drift from the exported statistics.
# ---------------------------------------------------------------------------
perm <- read.csv(file.path(stats_dir, "plsda_permutation.csv"), check.names = FALSE)
for (g in genotypes) {
  f <- fits[[g]]; if (is.null(f)) next
  rot <- f$rot; y <- f$y
  pr  <- perm[perm$subgroup == g, ]

  p <- scores.plot.ellipse(rot$T1, rot$T2, as.character(y), "Flicker",
                           analysis = "PLSDA", plot.title = "PLSDA") +
    scale_color_manual(values = pal_of(g)) +
    labs(subtitle = sprintf(
      "permutation (1000x): LOOCV accuracy = %s vs chance %s, p = %s%s",
      fmt(pr$loocv_accuracy), fmt(pr$chance), fmt(pr$p_permute),
      if (isTRUE(pr$p_permute >= alpha)) " (n.s.)" else "")) +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 10))
  fn <- if (g == "WT") "Fig2C_plsda_scores_WT.png" else "Supp4B_plsda_scores_5XFAD.png"
  ggsave(file.path(outdir, fn), p, width = 6, height = 5.2, dpi = DPI)

  # "Signals in LV1/LV2" from the SAME model, with LOOCV error bars. Top-3 picks
  # for that LV are drawn in a darker fill so the selection is visible.
  lbl <- if (g == "WT") "Fig2D" else "Supp4B"
  for (lv in 1:2) {
    L    <- if (lv == 1) rot$P1 else rot$P2
    SD   <- f$loocv$sd[, lv]
    keys <- key_tbl$analyte[key_tbl$genotype == g & key_tbl$LV == paste0("LV", lv)]
    pl <- loadings_chart_sd_gg(L, SD, sprintf("LV%d", lv), x_size = 9,
                               highlight = keys)
    ggsave(file.path(outdir, sprintf("%s_loadings_%s_LV%d.png", lbl, g, lv)),
           pl, width = 7, height = 3.8, dpi = DPI)
  }
}

# ---------------------------------------------------------------------------
# 2. Amyloid bars (Fig 4B) -- flicker-colored, corrected KW/Dunn p's
# ---------------------------------------------------------------------------
fad <- animals[animals$Genotype == "5XFAD", , drop = FALSE]
ab <- ggarrange(ms_bar(fad, "Abeta 40", "5XFAD", "Aβ 1-40 (a.u.)"),
                ms_bar(fad, "Abeta 42", "5XFAD", "Aβ 1-42 (a.u.)"),
                ncol = 2, common.legend = TRUE, legend = "right")
ab <- annotate_figure(ab, top = text_grob("Pathology Markers", face = "bold", size = 15))
ggsave(file.path(outdir, "Fig4B_amyloid_bars.png"), ab, width = 8, height = 4.3, dpi = DPI)

# ---------------------------------------------------------------------------
# 3. Fig 4C -- Spearman rho heatmap grid (cytokine x Abeta).
#    [REVISED SPEC] Rows are the 5xFAD key cytokines PLUS any retained analyte
#    reaching raw p < threshold on either marker -- "the p values are the driver
#    for what is included in the figure". The old hand-picked row list included
#    analytes outside the assay's linear range and is gone.
#    Each cell prints rho (in biology a rho of ~.46 is a reportable relationship,
#    so the coefficient belongs on the figure) with an asterisk where the pair
#    survives the heatmap-panel FDR.
# ---------------------------------------------------------------------------
ca <- read.csv(file.path(stats_dir, "corr_cytokine_amyloid_byrow.csv"), check.names = FALSE)
fig4c_rows <- unique(c(key_cytokines[["5XFAD"]],
                       ca$analyte[!is.na(ca$p_raw) & ca$p_raw < corr_p_threshold]))
fig4c_rows <- cytokines_all[cytokines_all %in% fig4c_rows]
cc <- ca[ca$analyte %in% fig4c_rows & ca$target %in% c("Abeta 40", "Abeta 42"), ]
cc$target  <- factor(cc$target, levels = c("Abeta 40", "Abeta 42"),
                     labels = c("Aβ 1-40", "Aβ 1-42"))
cc$analyte <- factor(cc$analyte, levels = rev(fig4c_rows))
cc$lab <- paste0(formatC(cc$rho, format = "f", digits = 2), stars(cc$p_adj))
n_surv <- sum(!is.na(cc$p_adj) & cc$p_adj < alpha)
p4c <- ggplot(cc, aes(target, analyte, fill = rho)) +
  geom_tile(color = "black", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation\nCoefficient (ρ)") +
  labs(x = NULL, y = NULL, title = "Spearman Correlation\nCytokines vs Pathology Markers",
       caption = sprintf(
         "Rows: 5xFAD key cytokines + any analyte with p < %.2f\n* q < %.2f, BH across the %d retained analytes per marker%s",
         corr_p_threshold, alpha, length(unique(ca$analyte)),
         if (n_surv == 0) " — no pair survives" else "")) +
  # no coord_fixed: match the original's WIDE cells / ~1.39:1 panel aspect
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        axis.text = element_text(color = "black", face = "bold"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 8, color = "grey30"),
        plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(outdir, "Fig4C_cytokine_amyloid_heatmap.png"), p4c,
       width = 7.0, height = 5.0, dpi = DPI)

# ---------------------------------------------------------------------------
# 4. Key-cytokine bar panels (Supp 5) -- top 3 on LV1 + top 3 on LV2 per
#    genotype, corrected Tukey / Dunn p's.
# ---------------------------------------------------------------------------
for (g in genotypes) {
  d    <- animals[animals$Genotype == g, , drop = FALSE]
  keys <- intersect(key_cytokines[[g]], analysis_analytes)
  if (!length(keys)) next
  plots <- lapply(keys, function(a) ms_bar(d, a, g, paste0(a, " (a.u.)")))
  nc  <- min(3, length(plots)); nr <- ceiling(length(plots) / nc)
  pan <- ggarrange(plotlist = plots, ncol = nc, nrow = nr,
                   common.legend = TRUE, legend = "top")
  ggsave(file.path(outdir, sprintf("Fig_key_cytokines_%s.png", g)), pan,
         width = 3.7 * nc, height = 3.5 * nr, dpi = DPI)
}

# ---------------------------------------------------------------------------
# 5. Sex-batched heatmaps (heatmap3 -- already manuscript-consistent)
# ---------------------------------------------------------------------------
for (g in genotypes) {
  df <- animals[animals$Genotype == g, , drop = FALSE]
  matZ <- apply(as.matrix(df[, analysis_analytes]), 2, zscore)
  grp  <- interaction(df$Sex, df$Flicker, sep = " ", lex.order = TRUE)
  gm   <- t(t(apply(matZ, 2, function(c) tapply(c, grp, mean, na.rm = TRUE))))
  hc   <- hclust(dist(t(gm)), "average")
  side <- c(M = "#111111", F = "#FFFFFF")[substr(rownames(gm), 1, 1)]
  br   <- c(-200, seq(-1.5, 1.5, 0.01), 200)
  png(file.path(outdir, sprintf("Fig_sexbatched_heatmap_%s.png", g)),
      width = 11, height = 5, units = "in", res = DPI)
  heatmap3(gm, Rowv = NA, Colv = as.dendrogram(hc), scale = "none",
           col = colorpanel(length(br) - 1, "blue", "white", "red"), breaks = br,
           RowSideColors = side, cexRow = 1, cexCol = 0.9, margins = c(7, 8),
           main = sprintf("%s cytokines by Sex x Flicker (z-scored)", g),
           legendfun = function() showLegend(legend = c("Male", "Female"), col = c("#111111", "#BBBBBB")))
  dev.off()
}

# stale figures from the previous (myeloid-panel) round would otherwise sit next
# to the new ones and get pasted into the manuscript by mistake
for (f in file.path(outdir, sprintf("Fig_myeloid_panel_%s.png", genotypes)))
  if (file.exists(f)) file.remove(f)

cat("[09] revised figures (manuscript style) ->", outdir, "\n")
cat("    ", paste(list.files(outdir), collapse = "\n     "), "\n")

# ============================================================================
# 10_figure5_correlations.R
# Figure 5 (5xFAD): pathology / cytokines vs glial morphology, rebuilt in the
# ORIGINAL Wood-lab scatter style (AshleyDataAnalysis_SBversion_Final.R L910-947):
#   shape-21 points filled by flicker with a black border, a black longdash lm
#   line WITH the grey CI ribbon (se=TRUE), Spearman rho + p annotated on top,
#   5xFAD palette (#4E5491 20Hz / #863775 40Hz / #B76945 No Flicker).
# [CORRECTION] Changes from the original:
#   - males/females are shown (circle = M, triangle = F), initialComments.md L26;
#   - the statistic is Spearman (it already was), the line stays a linear/Pearson
#     trend for display only, and p-values drive the asterisks -- no rho cutoff;
#   - [REVISED SPEC] the cytokines shown are the per-genotype KEY CYTOKINES from
#     the PLS-DA loadings (top 3 on LV1 + top 3 on LV2, stage 04b). The original
#     panels hard-coded Eotaxin and IL-15, both outside the assay's linear range;
#   - each cytokine panel also carries q, the BH-adjusted p across the retained
#     heatmap panel for that morphology measure (stage 07).
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(ggpubr) })
if (!exists("combined_xlsx")) source("config.R")

animals <- readRDS(file.path(results_dir, "processed", "animal_table.rds"))
key_cytokines <- readRDS(file.path(results_dir, "processed", "key_cytokines.rds"))
cm <- read.csv(file.path(stats_dir, "corr_cytokine_morphology_byrow.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
fad <- animals[animals$Genotype == "5XFAD", , drop = FALSE]
wt  <- animals[animals$Genotype == "WT", , drop = FALSE]

# BH-adjusted p for one cytokine x morphology pair, from the stage-07 table
q_of <- function(g, analyte, target) {
  i <- cm$subgroup == g & cm$analyte == analyte & cm$target == target
  if (!any(i)) NA_real_ else cm$p_adj[which(i)[1]]
}
outdir <- file.path(results_dir, "revised_figures_png")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
DPI <- 220
star <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))

# one scatter, matched to the original style; sex shown by point shape.
# ribbon = TRUE reproduces Fig 5 (se ribbon), FALSE reproduces Supp 7 (line only).
ms_scatter <- function(df, xcol, ycol, xlab, ylab, pal, ribbon = TRUE, q = NA_real_) {
  d <- data.frame(x = df[[xcol]], y = df[[ycol]],
                  Flicker = factor(df$Flicker, levels = flicker_levels),
                  Sex = factor(df$Sex, levels = sexes))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  ct  <- suppressWarnings(stats::cor.test(d$x, d$y, method = "spearman"))
  rho <- unname(ct$estimate); p <- ct$p.value
  # p drives whether a pair is shown; the asterisks follow the CORRECTED value
  # (q) wherever a correction applies, so the stars can't outrun the FDR.
  sub <- sprintf("ρ = %.2f, p = %s%s %s", rho,
                 formatC(p, format = "fg", digits = 3),
                 if (is.na(q)) "" else sprintf(", q = %s", formatC(q, format = "fg", digits = 3)),
                 star(if (is.na(q)) p else q))
  ggplot(d, aes(x = x, y = y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = ribbon, color = "black",
                linetype = "longdash", linewidth = 0.7) +
    geom_point(aes(fill = Flicker, shape = Sex), size = 3.6, stroke = 0.8,
               color = "black") +
    scale_shape_manual(values = c(M = 21, F = 24), name = "Sex") +
    scale_fill_manual(values = pal, name = "Treatment") +
    labs(x = xlab, y = ylab, subtitle = sub) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
           shape = guide_legend(override.aes = list(fill = "grey40", size = 4))) +
    theme(panel.background = element_rect(fill = "white"),
          text = element_text(size = 13),
          axis.text = element_text(color = "black", size = 11),
          axis.line = element_line(color = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
          panel.grid = element_blank(),
          plot.subtitle = element_text(hjust = 0, size = 11),
          plot.margin = margin(8, 12, 6, 8))
}

# a labeled row of panels sharing one legend + a bold row title
row_fig <- function(df, specs, title, ncol, file, w, h, pal, ribbon = TRUE) {
  plots <- lapply(specs, function(s)
    ms_scatter(df, s$x, s$y, s$xlab, s$ylab, pal, ribbon,
               q = if (is.null(s$q)) NA_real_ else s$q))
  g <- ggarrange(plotlist = plots, ncol = ncol, nrow = 1,
                 common.legend = TRUE, legend = "right")
  g <- annotate_figure(g, top = text_grob(title, face = "bold", size = 15))
  ggsave(file.path(outdir, file), g, width = w, height = h, dpi = DPI)
  cat(sprintf("[10] %-42s -> %s\n", title, file))
}

AST_TILE <- "AST mean_Distance to Nearest Neighbor"
MG_TILE  <- "MG mean_Distance to Nearest Neighbor"
AST_N    <- "AST mean_Number of Cells"
MG_N     <- "MG mean_Number of Cells"

# ============================ Figure 5 (5xFAD, with CI ribbon) ==============
row_fig(fad, list(
  list(x = AST_TILE, y = "Abeta 40", xlab = "Astrocyte Average\nDistance to Nearest Neighbor", ylab = "Aβ 40 (a.u.)"),
  list(x = AST_TILE, y = "Abeta 42", xlab = "Astrocyte Average\nDistance to Nearest Neighbor", ylab = "Aβ 42 (a.u.)")),
  "Pathology vs Avg Astrocyte Tiling", 2, "Fig5A_amyloid_vs_astro_tiling.png", 8, 4, flicker_colors_fad)
row_fig(fad, list(
  list(x = MG_TILE, y = "Abeta 40", xlab = "Microglia Average\nDistance to Nearest Neighbor", ylab = "Aβ 40 (a.u.)"),
  list(x = MG_TILE, y = "Abeta 42", xlab = "Microglia Average\nDistance to Nearest Neighbor", ylab = "Aβ 42 (a.u.)")),
  "Pathology vs Avg Microglia Tiling", 2, "Fig5B_amyloid_vs_mg_tiling.png", 8, 4, flicker_colors_fad)
# [REVISED SPEC] panels C/D now show the 5xFAD key cytokines, not the original's
# hard-coded RANTES/Eotaxin/IL-15/IL-10 (Eotaxin and IL-15 are out of range).
keys_fad <- key_cytokines[["5XFAD"]]
nk <- length(keys_fad)
row_fig(fad, lapply(keys_fad, function(cy)
  list(x = AST_N, y = cy, xlab = "Average Number of Astrocyte", ylab = paste0(cy, " (a.u.)"),
       q = q_of("5XFAD", cy, AST_N))),
  "Key Cytokines vs Avg Number of Astrocytes", nk, "Fig5C_cytokine_vs_astro_number.png", 3.5 * nk, 4, flicker_colors_fad)
row_fig(fad, lapply(keys_fad, function(cy)
  list(x = MG_N, y = cy, xlab = "Average Number of Microglia", ylab = paste0(cy, " (a.u.)"),
       q = q_of("5XFAD", cy, MG_N))),
  "Key Cytokines vs Avg Number of Microglia", nk, "Fig5D_cytokine_vs_mg_number.png", 3.5 * nk, 4, flicker_colors_fad)

# ============================ Supp Fig 7 (WT, no ribbon, Spearman) ==========
# [CORRECTION] original used Pearson ("R="); now Spearman (ρ) + males/females.
keys_wt <- key_cytokines[["WT"]]
nw <- length(keys_wt)
supp <- function(target, xlab, title, file)
  row_fig(wt, lapply(keys_wt, function(cy)
    list(x = target, y = cy, xlab = xlab, ylab = cy, q = q_of("WT", cy, target))),
    title, nw, file, 3.5 * nw, 4, flicker_colors_wt, ribbon = FALSE)
supp(AST_TILE, "AST mean_Distance to Nearest Neighbor",
     "Key Cytokines vs Astrocyte Tiling",           "Supp7A_cyt_vs_astro_tiling_WT.png")
supp(MG_TILE,  "MG mean_Distance to Nearest Neighbor",
     "Key Cytokines vs Microglia Tiling",           "Supp7B_cyt_vs_mg_tiling_WT.png")
supp(AST_N,    "AST mean_Number of Cells",
     "Key Cytokines vs Avg Number of Astrocytes",   "Supp7C_cyt_vs_astro_number_WT.png")
supp(MG_N,     "MG mean_Number of Cells",
     "Key Cytokines vs Avg Number of Microglia",    "Supp7D_cyt_vs_mg_number_WT.png")

cat("[10] Figure 5 + Supp 7 rebuilt in original style ->", outdir, "\n")

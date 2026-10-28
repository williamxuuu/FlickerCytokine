# ============================================================================
# 04_panels.R  —  the published Figure 2 panel language, as functions.
#
# Figure 2 (WT) and Supplementary Figure 4 (5xFAD) are the SAME FIGURE drawn on
# two genotypes: cytokine heatmap, PLS-DA scores, LV score violins, LV1/LV2
# loadings. They are defined here exactly once; stages 10 and 13 supply only a
# genotype, so the two figures cannot drift into different visual languages.
#
# Everything here READS THE CACHED PLS-DA from stage 02. No panel refits a
# model, so a figure cannot disagree with the statistics reported beside it.
#
# MATCHED from the published panels (manuscript/ManuscriptMostRecent.docx,
# image2.png = Figure 2, image10.png = Supplementary Figure 4):
#   heatmap   rows = animals grouped No Flicker / 40Hz / 20Hz (that order), group
#             name on the LEFT; Flicker + Sex colour strips; analyte names at the
#             BOTTOM, rotated; no dendrogram drawn; z-score bar -1.5..1.5 labelled
#             "z-scored"; Flicker legend as coloured dots, M/F as filled/open
#             squares; title "Inflammatory Cytokines and Chemokine Profile".
#   scores    "PLSDA", solid black axes through the origin, group ellipses,
#             filled circles, axes labelled "Scores on LV1/LV2 (R2X%)".
#   violins   one panel per LV, stacked; violin filled by flicker with a black
#             outline; black mean point with SD whiskers inside; boxed panel, no
#             grid; 45-degree x labels; no legend.
#   loadings  "Signals in LV1/LV2", vertical light-blue bars with black outline,
#             ascending, loadings normalised so |max| = 1, names at the bottom,
#             the high/low flicker group named at each end.
#
# THE VIOLINS CARRY NO STATISTICS. The published panels print p-values and
# significance brackets on them (Fig 2C: 0.0031, 0.06, 0.00013, 0.00068). Those
# are tests run on a latent variable that was FITTED to separate the groups, so
# they are circular -- the Wood lab asked for them to go (initialComments.md:8,
# "cannot do stats on linked variable"), so the violins stay as a picture only.
# Model credibility is reported instead as LOOCV accuracy
# against a label-permutation null, in the scores subtitle, in words.
# Do not add brackets, asterisks or p-values to these panels.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

MS_ORDER  <- c("NoFlicker", "40Hz", "20Hz")     # heatmap rows, top-to-bottom, as published
MS_XORDER <- c("NoFlicker", "20Hz", "40Hz")     # violin x-axis, as published
MS_LABEL  <- c(NoFlicker = "No\nFlicker", "40Hz" = "40Hz", "20Hz" = "20Hz")
MS_PRETTY <- c(NoFlicker = "No Flicker", "20Hz" = "20Hz", "40Hz" = "40Hz")

zsc <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s
}

# Annotation-column ramps, sampled from the published Supp 4A colour bars.
AB_RAMPS <- list(
  "Abeta 40" = c("#F9F3DC", "#F7D9A1", "#F7C76D", "#F3AF38", "#EB9D1F"),
  "Abeta 42" = c("#E3E0C5", "#CDC494", "#B8AA66", "#A29037", "#8A7A2E"))

# rows ordered the published way, so every panel for one genotype agrees
ms_rows <- function(g) {
  d <- ANIMALS[ANIMALS$Genotype == g, , drop = FALSE]
  d$Flicker <- factor(as.character(d$Flicker), levels = MS_ORDER)
  d[order(d$Flicker, d$Sex), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Heatmap (Fig 2B / Supp 4A).
#
# `annot` names continuous per-animal columns to draw as extra annotation
# columns between the group strips and the analytes -- Supp 4A carries Abeta 42
# and Abeta 40 there, as the published panel does. Figure 2B passes none: WT has
# no amyloid, and the published Figure 2 has no such columns.
#
# `col_order` picks the analyte (column) order: "cluster" is the average-linkage
# order the published panel used; "lv1" sorts the columns by the LV1 loadings,
# ascending, so they line up left-to-right with the "Signals in LV1" panel
# (collaborator request for Figure 2B).
# ---------------------------------------------------------------------------
ms_heatmap <- function(g, figure, stem, annot = character(), width = 12.2, height = 5.4,
                       col_order = c("cluster", "lv1")) {
  col_order <- match.arg(col_order)
  cols <- pal(g)
  d <- ms_rows(g)
  Z <- apply(as.matrix(d[, PANEL]), 2, zsc)
  ord <- if (col_order == "lv1") {
    match(ms_lv_order(readRDS(file.path(STATS, sprintf("plsda_%s.rds", g))), 1), PANEL)
  } else hclust(dist(t(Z), method = "euclidean"), method = "average")$order
  d$row <- factor(seq_len(nrow(d)), levels = rev(seq_len(nrow(d))))

  long <- data.frame(
    row     = rep(d$row, times = length(PANEL)),
    Flicker = rep(d$Flicker, times = length(PANEL)),
    analyte = factor(rep(PANEL[ord], each = nrow(d)), levels = PANEL[ord]),
    z       = as.vector(Z[, ord]))

  hm <- ggplot(long, aes(analyte, row, fill = z)) +
    geom_tile(colour = "black", linewidth = 0.18) +
    facet_grid(Flicker ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(Flicker = MS_LABEL)) +
    scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1.5, 1.5),
                         oob = scales::squish, name = "z-scored",
                         breaks = c(-1.5, -0.75, 0, 0.75, 1.5)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title = "Inflammatory Cytokines and Chemokine Profile", x = NULL, y = NULL) +
    theme_ms(11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, colour = INK),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line = element_blank(), panel.grid = element_blank(),
          panel.spacing.y = unit(1.6, "pt"),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"),
          plot.title = element_text(hjust = 0.5),
          legend.key.height = unit(26, "pt"), legend.key.width = unit(9, "pt"))

  # --- Flicker / Sex strips, plus any annotation columns -------------------
  # Annotation fills are resolved to hex here and drawn through
  # scale_fill_identity behind a new_scale_fill, for two reasons: a fill vector
  # passed outside aes() is not split correctly when the plot is facetted, and
  # an `aes(x = 2 + i)` inside a loop resolves `i` at render time, so every
  # column would land on the last one. Both of those were live bugs here.
  strip <- ggplot(d, aes(x = 1, y = row)) +
    geom_tile(aes(fill = Flicker), width = 0.9) +
    geom_tile(aes(x = 2, fill = Sex), width = 0.9, colour = "black", linewidth = 0.18) +
    scale_fill_manual(values = c(cols, M = "black", F = "white"), guide = "none")

  ab_keys <- list()
  if (length(annot)) {
    ann <- do.call(rbind, lapply(seq_along(annot), function(i) {
      v  <- d[[annot[i]]]
      rp <- AB_RAMPS[[annot[i]]]
      if (is.null(rp)) rp <- c("#F2F2F0", "#6E6E6E")
      fillv <- scales::gradient_n_pal(rp)(scales::rescale(v))
      fillv[is.na(fillv)] <- "#FFFFFF"
      ab_keys[[annot[i]]] <<- list(ramp = rp, lo = min(v, na.rm = TRUE),
                                   hi = max(v, na.rm = TRUE), lab = annot[i])
      data.frame(x = 2 + i, row = d$row, Flicker = d$Flicker, col = fillv)
    }))
    strip <- strip + ggnewscale::new_scale_fill() +
      geom_tile(data = ann, aes(x = x, y = row, fill = col),
                width = 0.9, colour = "black", linewidth = 0.18, inherit.aes = FALSE) +
      scale_fill_identity(guide = "none")
  }

  strip <- strip +
    facet_grid(Flicker ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    theme_void() + theme(strip.text = element_blank(),
                         panel.spacing.y = unit(1.6, "pt"),
                         plot.margin = margin(0, 1, 0, 0))

  # --- right-hand key: Flicker dots, M/F squares, annotation colour bars ---
  klab <- MS_PRETTY
  kf <- data.frame(y = 3:1, lab = unname(klab[MS_ORDER]), col = unname(cols[MS_ORDER]))
  ks <- data.frame(y = c(-0.6, -1.4), lab = c("M", "F"), fill = c("black", "white"))
  leg <- ggplot() +
    geom_point(data = kf, aes(0, y), colour = kf$col, size = 4) +
    geom_text(data = kf, aes(0.28, y, label = lab), hjust = 0, size = 3.2, colour = INK) +
    annotate("text", x = 0, y = 4.0, label = "Flicker", hjust = 0, fontface = "bold",
             size = 3.5, colour = INK) +
    geom_tile(data = ks, aes(0, y), fill = ks$fill, colour = "black",
              width = 0.34, height = 0.55, linewidth = 0.3) +
    geom_text(data = ks, aes(0.28, y, label = lab), hjust = 0, size = 3.2, colour = INK)

  # Each annotation bar is stacked tiles with literal fills, so it needs no
  # scale of its own and renders the same on ragg and base pdf.
  y_bot <- -2.6
  for (nm in names(ab_keys)) {
    k <- ab_keys[[nm]]
    nseg <- 48; h <- 2.6
    seg <- data.frame(y = y_bot - h + (seq_len(nseg) - 0.5) * h / nseg,
                      col = scales::gradient_n_pal(k$ramp)(seq(0, 1, length.out = nseg)))
    leg <- leg +
      geom_tile(data = seg, aes(0, y), fill = seg$col, width = 0.34,
                height = h / nseg, inherit.aes = FALSE) +
      annotate("rect", xmin = -0.17, xmax = 0.17, ymin = y_bot - h, ymax = y_bot,
               fill = NA, colour = "black", linewidth = 0.3) +
      annotate("text", x = 0.28, y = y_bot, hjust = 0, vjust = 1, size = 2.7,
               colour = INK, label = signif(k$hi, 2)) +
      annotate("text", x = 0.28, y = y_bot - h, hjust = 0, vjust = 0, size = 2.7,
               colour = INK, label = signif(k$lo, 2)) +
      annotate("text", x = 0, y = y_bot + 0.30, hjust = 0, vjust = 0, size = 2.9,
               colour = INK, label = ab_lab(k$lab, ""))
    y_bot <- y_bot - h - 1.5
  }

  y_lo <- if (length(ab_keys)) y_bot else -2.0
  leg <- leg + scale_x_continuous(limits = c(-0.3, if (length(ab_keys)) 2.6 else 2.1)) +
    scale_y_continuous(limits = c(y_lo, 4.4)) + theme_void()

  save_fig((strip | hm | leg) +
             plot_layout(widths = c(0.035 + 0.018 * length(annot), 1,
                                    0.085 + 0.045 * length(ab_keys))),
           figure, stem, width, height)
}

# ---------------------------------------------------------------------------
# PLS-DA scores, violins and loadings (Fig 2C/2D, Supp 4B/4C/4D).
# All three read the cached fit from stage 02.
# ---------------------------------------------------------------------------

# --- scores ----------------------------------------------------------------
ms_scores_plot <- function(g, cc) {
  cols <- pal(g)
  r2 <- round(cc$r2x * 100, 1)
  # MS_XORDER, not MS_ORDER: the heatmap's row order exists to group the rows,
  # but a legend should read No Flicker -> 20Hz -> 40Hz.
  # Sex is point shape, as in every other scatter (03_theme.R): the model is fit
  # on the sex-z-scored table, and the reader should see the sexes interleave.
  sc <- data.frame(LV1 = cc$scores[, "LV1"], LV2 = cc$scores[, "LV2"],
                   Flicker = factor(cc$flicker, levels = MS_XORDER),
                   Sex     = factor(cc$sex, levels = SEXES))
  ggplot(sc, aes(LV1, LV2)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    stat_ellipse(aes(colour = Flicker), type = "norm", level = 0.95,
                 linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(fill = Flicker, shape = Sex), colour = "black",
               size = 3.4, stroke = 0.45) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_fill_manual(values = cols, name = NULL, labels = MS_PRETTY) +
    scale_shape_manual(values = SEX_SHAPE, name = "Sex") +
    guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 3.4)),
           shape = guide_legend(order = 2,
                     override.aes = list(fill = "grey40", colour = "black", size = 3.4))) +
    labs(title = "PLSDA",
         subtitle = sprintf("%d analytes - LOOCV %.0f%% vs %.0f%% chance - permutation p = %.2f",
                            length(cc$analytes), 100 * cc$obs, 100/3, cc$pperm),
         x = sprintf("Scores on LV1 (%.1f%%)", r2[1]),
         y = sprintf("Scores on LV2 (%.1f%%)", r2[2])) +
    theme_ms(11) + theme(plot.title = element_text(hjust = 0.5),
                         # name the children: theme_ms sets panel.grid.major, so
                         # clearing the parent alone leaves the grid showing
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black", linewidth = 0.4))
}

# --- LV score violins, NO STATISTICS (see header) --------------------------
# Matched to the published panel: violin filled by flicker with a black
# outline, a black mean point with SD whiskers inside, boxed panel, no grid,
# 45-degree x labels, no legend. What is deliberately absent is the bracket,
# the p-value and the asterisk the published panel carries.
ms_violin_plot <- function(g, cc, k) {
  cols <- pal(g)
  r2 <- round(cc$r2x * 100, 1)
  v  <- data.frame(y = cc$scores[, k], Flicker = factor(cc$flicker, levels = MS_XORDER))
  mu <- do.call(rbind, lapply(split(v$y, v$Flicker), function(z)
          data.frame(m = mean(z), s = stats::sd(z))))
  mu$Flicker <- factor(rownames(mu), levels = MS_XORDER)

  ggplot(v, aes(Flicker, y)) +
    geom_violin(aes(fill = Flicker), colour = "black", linewidth = 0.4,
                width = 0.85, trim = FALSE, scale = "width") +
    geom_errorbar(data = mu, aes(x = Flicker, ymin = m - s, ymax = m + s),
                  width = 0.10, linewidth = 0.45, colour = "black", inherit.aes = FALSE) +
    geom_point(data = mu, aes(x = Flicker, y = m), inherit.aes = FALSE,
               size = 1.9, colour = "black") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_discrete(labels = MS_PRETTY[MS_XORDER]) +
    labs(x = NULL, y = sprintf("Scores on LV%d (%.1f%%)", k, r2[k])) +
    theme_ms(10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
          axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, colour = INK),
          axis.title.y = element_text(size = 8.5, colour = INK))
}

# --- loadings --------------------------------------------------------------
# signed, |max|-normalised LV k loadings; shared with the heatmap's "lv1" order
ms_lv_loadings <- function(cc, k) {
  v <- cc$loadings[, k] * cc$sign[k]; v / max(abs(v))
}
ms_lv_order <- function(cc, k) { v <- ms_lv_loadings(cc, k); names(v)[order(v)] }

ms_loading_plot <- function(g, cc, k) {
  v <- ms_lv_loadings(cc, k)
  mu <- tapply(cc$scores[, k], cc$flicker, mean)
  ends <- c(hi = names(which.max(mu)), lo = names(which.min(mu)))
  df <- data.frame(analyte = factor(names(v), ms_lv_order(cc, k)), loading = v)
  ggplot(df, aes(analyte, loading)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_col(fill = "#BDD7EE", colour = "black", linewidth = 0.3, width = 0.72) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = c(0, 0)) +
    annotate("text", x = 0.6, y = 0.85, hjust = 0, size = 3.1, colour = INK,
             label = MS_PRETTY[[ends[["hi"]]]]) +
    annotate("text", x = 0.6, y = -0.85, hjust = 0, size = 3.1, colour = INK,
             label = MS_PRETTY[[ends[["lo"]]]]) +
    labs(title = sprintf("Signals in LV%d", k), x = NULL, y = NULL) +
    theme_ms(11) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, colour = INK),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line.y = element_line(colour = "black", linewidth = 0.4))
}

# ---------------------------------------------------------------------------
# Driver: writes the scores, the violins (alone and stacked), the published
# scores+violins composite, and both loadings panels.
# `stems` names the output files so Fig 2 and Supp 4 keep their own numbering.
# ---------------------------------------------------------------------------
ms_plsda <- function(g, figure, stems) {
  cc <- readRDS(file.path(STATS, sprintf("plsda_%s.rds", g)))

  sco <- ms_scores_plot(g, cc)
  save_fig(sco, figure, stems$scores, 5.4, 5.0)

  v1 <- ms_violin_plot(g, cc, 1); v2 <- ms_violin_plot(g, cc, 2)
  vio <- (v1 / v2) + plot_layout(heights = c(1, 1))
  save_fig(vio, figure, stems$violins, 2.9, 5.0)

  # published layout: scores on the left, the two violins stacked to the right
  save_fig((sco | vio) + plot_layout(widths = c(1.85, 1)),
           figure, stems$composite, 8.3, 5.0)

  for (k in 1:2) save_fig(ms_loading_plot(g, cc, k), figure, sprintf(stems$loadings, k), 6.2, 4.4)
}

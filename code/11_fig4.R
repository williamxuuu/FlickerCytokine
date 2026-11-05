# ============================================================================
# 11_fig4.R  —  Figure 4 B / C, drawn to match the published Figure 4
# (manuscript/ManuscriptMostRecent.docx, word/media/image4.png).
#
# Panel 4A is the BioRender paradigm schematic and is not generated here.
#
#   4B  amyloid bars, Abeta 1-40 and Abeta 1-42, with every pairwise p-value
#       bracketed as published and the omnibus test stated beneath each title.
#   4C  cytokine ~ amyloid Spearman rho heatmap. ROWS are chosen, not listed:
#       the 5xFAD key cytokines from the stage-02 PLS-DA, plus any analyte
#       reaching permutation p < 0.05 on either marker. ORDER is presentation
#       only -- strongest to weakest on Abeta 1-40, the way the published panel
#       reads.
#
# Asterisks follow BH q across the 28-analyte panel, not the raw p the
# published panel starred. Westfall-Young is the stricter FWER cross-check and
# lives in stats/corr_cytokine_amyloid.csv alongside n_distinct and prop_tied;
# read the star with that table, not instead of it.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "03_theme.R"))
source(file.path(HERE, "07_figure_helpers.R"))
source(file.path(HERE, "lib", "stats_helpers.R"))

fad <- ANIMALS[ANIMALS$Genotype == "5XFAD", , drop = FALSE]

# ============================ Fig 4B — amyloid bars =========================
r  <- compare_flicker_many(fad, AMYLOID_TARGETS, subgroup = "5XFAD")
ps <- lapply(AMYLOID_TARGETS, function(v)
  bar_panel(fad, v, "5XFAD", r$omnibus[r$omnibus$variable == v, ], r$posthoc,
            title_expr = ab_lab(v, "")))

fig4b <- wrap_plots(ps, ncol = 2, guides = "collect") +
  plot_annotation(title = "Pathology Markers",
    theme = theme_ms(12) +
      theme(text = element_text(family = MS_FONT),
            plot.title = element_text(family = MS_FONT, face = "bold", size = 15,
                                      hjust = 0, colour = INK, margin = margin(b = 4))))
fig4b <- fig4b & theme(legend.position = "right",
                       legend.text = element_text(family = MS_FONT, size = 9, colour = INK))
save_fig(fig4b, "Fig4", "Fig4B_amyloid_bars", 8.6, 5.6)

# ============================ Fig 4C — rho heatmap ==========================
CA <- read.csv(file.path(STATS, "corr_cytokine_amyloid.csv"), check.names = FALSE)

rows  <- unique(c(key_cytokines("5XFAD"), CA$analyte[CA$p_perm < 0.05]))
ord40 <- CA[CA$target == "Abeta 40" & CA$analyte %in% rows, ]
rows  <- ord40$analyte[order(ord40$rho, decreasing = TRUE)]

h <- CA[CA$analyte %in% rows, ]
h$analyte <- factor(h$analyte, levels = rev(rows))
h$target  <- factor(h$target, levels = AMYLOID_TARGETS)

fig4c <- ggplot(h, aes(target, analyte, fill = rho)) +
  geom_tile(colour = "black", linewidth = 0.7) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), vapply(p_adj, stars, "")),
                colour = abs(rho) > 0.62),
            size = 3.3, show.legend = FALSE, family = MS_FONT) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradientn(colours = MS_DIVERGING, limits = c(-1, 1),
                       breaks = c(-1, 0, 1),
                       name = expression(paste("Correlation Coefficient (", rho, ")")),
                       guide = guide_colourbar(
                         title.position = "right", frame.colour = "black",
                         frame.linewidth = 0.4, ticks = FALSE,
                         barheight = unit(112, "pt"), barwidth = unit(9, "pt"),
                         title.theme = element_text(family = MS_FONT, angle = -90,
                                                    size = 9, colour = INK, hjust = 0.5))) +
  scale_x_discrete(position = "bottom",
                   labels = c(ab_lab("Abeta 40", ""), ab_lab("Abeta 42", ""))) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Spearman Correlation\nCytokines vs Pathology Markers", x = NULL, y = NULL,
       caption = "* q < 0.05, ** q < 0.01, *** q < 0.001\nBenjamini-Hochberg across the 28-analyte panel") +
  theme_fig4(11) +
  theme(panel.border = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 1, face = "plain", size = 10, colour = "black"),
        axis.text.x = element_text(face = "plain", size = 11, colour = "black",
                                   margin = margin(t = 3)),
        plot.title = element_text(family = MS_FONT, face = "bold", size = 13,
                                  hjust = 0.5, lineheight = 1.15, margin = margin(b = 8)),
        legend.title = element_text(family = MS_FONT),
        plot.caption = element_text(family = MS_FONT, size = 7.5, colour = INK2,
                                    hjust = 0.5, margin = margin(t = 8)),
        plot.margin = margin(6, 6, 6, 6))

save_fig(fig4c, "Fig4", "Fig4C_cytokine_amyloid_heatmap",
         6.0, 0.34 * nlevels(h$analyte) + 2.2)

cat("[11] Figure 4 done; 4C rows:", paste(rows, collapse = ", "), "\n")

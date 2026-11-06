# ============================================================================
# 03_theme.R  —  shared figure styling and output. Source after 00_setup.R.
#
# Sex is shown by point shape (circle = male, triangle = female) in every
# per-animal plot, so readers can see that the sexes are interleaved after the
# within-sex z-score rather than taking the adjustment on trust.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

SEX_SHAPE <- c(M = 21, F = 24)          # filled circle / filled triangle
INK       <- "#1A1A1A"                  # primary text
INK2      <- "#4A4A4A"                  # secondary text
INK3      <- "#8A8A8A"                  # muted / grid
SURFACE   <- "#FFFFFF"

# Diverging ramp for correlation (rho) heatmaps.
DIVERGING <- c("#4053A4", "#8E9BC9", "#FFFFFF", "#F09B92", "#ED2724")

theme_ms <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      text             = element_text(colour = INK),
      plot.title       = element_text(face = "bold", size = base + 2, colour = INK,
                                      margin = margin(b = 2)),
      plot.subtitle    = element_text(size = base - 1.5, colour = INK2,
                                      margin = margin(b = 8)),
      plot.caption     = element_text(size = base - 2.5, colour = INK3, hjust = 0),
      axis.title       = element_text(size = base - 0.5, colour = INK2),
      axis.text        = element_text(size = base - 1.5, colour = INK2),
      axis.line        = element_line(colour = INK3, linewidth = 0.3),
      axis.ticks       = element_line(colour = INK3, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#EDEDEA", linewidth = 0.3),
      legend.title     = element_text(size = base - 1.5, colour = INK2),
      legend.text      = element_text(size = base - 1.5, colour = INK2),
      legend.key.size  = unit(10, "pt"),
      strip.text       = element_text(face = "bold", size = base - 0.5, colour = INK),
      plot.margin      = margin(10, 12, 8, 10)
    )
}

# format p / q to 3 significant figures, with a floor for permutation p-values
fmt_p <- function(p, B = 50000) {
  if (is.na(p)) return("NA")
  if (p <= 1 / (B + 1)) sprintf("< %.0e", 1 / (B + 1)) else signif(p, 3)
}

# Abeta axis labels via plotmath -- renders identically on both devices.
ab_lab <- function(v, unit = " (z-score)")
  parse(text = sprintf('paste("A", beta, " %s%s")', sub("Abeta ", "1-", v), unit))

# Every measure downstream is a within-sex z-score (00_setup.R), so axes say so
# rather than "(a.u.)", which would imply raw MFI.
zlab <- function(v) paste0(v, " (z-score)")

# OUTPUT DEVICES. Every figure is written three ways:
#   PNG  300 dpi raster (ragg)
#   PDF  vector (base pdf device)
#   SVG  vector with live, editable text (svglite; fix_text_size = FALSE omits
#        the textLength attribute so text is not stretched when an editor
#        substitutes a font). Sizes are in points, so the SVG opens at the
#        physical size below.
# The base pdf() device cannot render Greek letters or arrows, so all figure
# text is ASCII or built with plotmath; the three formats are glyph-identical.
#
# FIG_SCALE scales every canvas while font, point and line sizes stay fixed in
# points, i.e. it sets text size relative to the panel. Stages pass nominal
# canvas sizes; the written size is nominal x FIG_SCALE.
FIG_SCALE <- 0.60

write_fig <- function(plot, f, width, height) {
  width <- width * FIG_SCALE; height <- height * FIG_SCALE
  ggsave(paste0(f, ".png"), plot, width = width, height = height, dpi = DPI,
         device = ragg::agg_png, bg = SURFACE)
  grDevices::pdf(paste0(f, ".pdf"), width = width, height = height, bg = SURFACE)
  print(plot); invisible(grDevices::dev.off())
  svglite::svglite(paste0(f, ".svg"), width = width, height = height, bg = SURFACE,
                   fix_text_size = FALSE)
  print(plot); invisible(grDevices::dev.off())
}

save_fig <- function(plot, figure, stem, width, height) {
  write_fig(plot, fig_path(figure, stem), width, height)
  cat(sprintf("     -> %s/%s.{png,pdf,svg}\n", figure, stem))
}

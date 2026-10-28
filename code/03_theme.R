# ============================================================================
# 03_theme.R  —  shared figure language. Source after 00_setup.R.
#
# The flicker palette lives in 00_setup.R and is the PUBLISHED manuscript's own
# set, so regenerated panels drop in next to untouched ones without looking
# foreign. See the note there on its colourblind trade-off, and why identity in
# this pipeline is never carried by colour alone.
#
# SEX IS POINT SHAPE in every scatter. This is a deliberate departure from the
# published panels, which use circles throughout. It is kept because sex is the
# largest nuisance effect in this dataset -- 13 of 56 genotype-analyte pairs
# show a sex difference at p < 0.05 before correction -- and a reader should be
# able to see the sexes interleaved rather than take the z-score on trust.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

SEX_SHAPE <- c(M = 21, F = 24)          # filled circle / filled triangle
INK       <- "#1A1A1A"                  # primary text
INK2      <- "#4A4A4A"                  # secondary text
INK3      <- "#8A8A8A"                  # muted / grid
SURFACE   <- "#FFFFFF"

# Diverging ramp for rho heatmaps, sampled from the published Figure 4 colour
# bar so Fig 4C keeps its own look.
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

# NOTE ON DEVICES. This R build reports capabilities("cairo") == TRUE but
# cairo_pdf() fails with "failed to load cairo DLL", so it silently writes
# nothing. We therefore use the base pdf() device -- which cannot render Greek
# letters or arrows (rho, beta, up-arrow all fail; em/en dashes are silently
# swapped for hyphens). Every figure string is kept ASCII-safe, or built with
# plotmath, so the PDF is glyph-identical to the PNG.
save_fig <- function(plot, figure, stem, width, height) {
  f <- fig_path(figure, stem)
  ggsave(paste0(f, ".png"), plot, width = width, height = height, dpi = DPI,
         device = ragg::agg_png, bg = SURFACE)
  grDevices::pdf(paste0(f, ".pdf"), width = width, height = height, bg = SURFACE)
  print(plot); invisible(grDevices::dev.off())
  cat(sprintf("     -> %s/%s.{png,pdf}\n", figure, stem))
}

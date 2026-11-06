# ============================================================================
# 07_figure_helpers.R  —  shared bar and scatter panels (Figures 4-5, Supp 5-7).
#
# One definition per panel type, used by every figure stage:
#   bars     group mean +/- SEM with one point per animal; Tukey / Dunn
#            adjusted pairwise p-values as brackets; the omnibus test (ANOVA
#            or Kruskal-Wallis) and its p printed beneath the title.
#   scatter  one point per animal with a dashed least-squares line (display
#            only); Spearman rho, permutation p and BH q printed on the panel.
#
# Conventions:
#   - Sex is point shape (03_theme.R).
#   - Pairwise contrasts are starred only when the omnibus test is
#     significant (p < 0.05); otherwise the adjusted p is printed unstarred.
#   - Correlation significance markers are based on BH q, not raw p.
#     Westfall-Young p-values are reported in the stats tables.
#   - Analytes measured at 6 or fewer distinct raw values are labelled
#     "[k levels]" on scatter panels, as a warning of coarse resolution.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# Helvetica is a base-14 PDF font AND a real system font, so the base pdf()
# device and ragg agree glyph for glyph. See the device note in 03_theme.R.
MS_FONT <- "Helvetica"

# diverging ramp for the Figure 4C rho heatmap, at -1, -0.5, 0, +0.5, +1
MS_DIVERGING <- c("#4053A4", "#A38EC3", "#FFFFFF", "#F79F86", "#ED2724")

FLICKER_MS <- c(NoFlicker = "No\nFlicker", "20Hz" = "20Hz", "40Hz" = "40Hz")

# Significance stars from adjusted p / q, defined once for every panel.
stars <- function(q) if (is.na(q)) "" else if (q < 0.001) " ***" else
                     if (q < 0.01) " **" else if (q < 0.05) " *" else ""

# pale bar fill: the group colour blended toward white
tint <- function(col, w = 0.42) {
  r <- grDevices::col2rgb(col)
  grDevices::rgb(t(255 - (255 - r) * w), maxColorValue = 255)
}

# panel theme: boxed panel, no grid, bold tick labels
theme_fig4 <- function(base = 11) {
  theme_ms(base) +
    theme(
      text             = element_text(family = MS_FONT, colour = INK),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major = element_blank(),      # must name the child: theme_ms set it
      panel.grid.minor = element_blank(),
      axis.line        = element_blank(),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(colour = "black", face = "bold", size = base - 1),
      plot.title       = element_text(family = MS_FONT, face = "plain", size = base,
                                      hjust = 0.5, colour = INK, margin = margin(b = 1)),
      plot.subtitle    = element_text(family = MS_FONT, size = base - 3, colour = INK2,
                                      hjust = 0.5, margin = margin(b = 4))
    )
}

# ---------------------------------------------------------------------------
# BAR PANEL (Fig 4B, Supp 5, Supp 6)
#
# `title_expr` lets Figure 4B pass a plotmath Abeta label where the cytokine
# panels pass a plain analyte name.
# `all_brackets` = TRUE draws every pairwise comparison with its p-value;
# FALSE draws only those with adjusted p < 0.05.
# ---------------------------------------------------------------------------
bar_panel <- function(d, v, g, om, ph, title_expr = v, ylab = NULL,
                      all_brackets = TRUE) {
  cols <- pal(g)
  dd <- data.frame(y = d[[v]], Flicker = factor(d$Flicker, FLICKER_LEVELS),
                   Sex = factor(d$Sex, SEXES))
  dd <- dd[complete.cases(dd), ]
  mu <- aggregate(y ~ Flicker, dd, function(z) c(m = mean(z), se = sd(z)/sqrt(length(z))))
  mu <- data.frame(Flicker = mu$Flicker, m = mu$y[, "m"], se = mu$y[, "se"])

  top  <- max(c(dd$y, mu$m + mu$se), na.rm = TRUE)
  base <- min(0, min(dd$y, na.rm = TRUE))
  span <- top - base; if (span <= 0) span <- 1
  data_top <- top          # keep the data range; `top` grows to cover the brackets

  p <- ggplot(dd, aes(Flicker, y)) +
    geom_col(data = mu, aes(y = m, fill = Flicker), width = 0.66,
             colour = "black", linewidth = 0.45) +
    geom_errorbar(data = mu, aes(y = m, ymin = m - se, ymax = m + se),
                  width = 0.17, linewidth = 0.5, colour = "black") +
    # points carry the full group colour; only the bar is tinted.
    # fill is passed outside aes() so it does not enter the bar's fill scale.
    geom_point(aes(shape = Sex), fill = cols[as.character(dd$Flicker)],
               colour = "black", size = 2.4, stroke = 0.45,
               position = position_jitter(width = 0.11, height = 0, seed = 1)) +
    scale_fill_manual(values = tint(cols), guide = "none") +
    scale_shape_manual(values = SEX_SHAPE, name = NULL) +
    scale_x_discrete(labels = FLICKER_MS) +
    labs(x = NULL, y = ylab, title = title_expr,
         subtitle = if (is.null(om) || !nrow(om)) NULL
                    else sprintf("%s  p = %s", om$test, signif(om$p, 3)))

  if (!is.null(ph)) {
    s <- ph[ph$variable == v, , drop = FALSE]
    if (!all_brackets) s <- s[s$p_adj < 0.05, , drop = FALSE]
    if (nrow(s)) {
      x1 <- match(s$group1, FLICKER_LEVELS); x2 <- match(s$group2, FLICKER_LEVELS)
      # adjacent up, span down; spaced for FIG_SCALE so the labels clear each other
      lvl <- ifelse(abs(x2 - x1) == 1, 0.33, 0.10)
      yb  <- top + span * lvl
      tick <- span * 0.018
      # Stars only behind a significant omnibus test: with omnibus p >= 0.05
      # the adjusted contrast p is printed but not starred.
      gate <- !is.null(om) && nrow(om) && isTRUE(om$p[1] < 0.05)
      lab <- ifelse(gate & s$p_adj < 0.05,
                    sprintf("%s *", signif(s$p_adj, 2)), sprintf("%s", signif(s$p_adj, 2)))
      br <- data.frame(x = x1, xend = x2, y = yb)
      p <- p +
        geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = 0.4, colour = INK3) +
        geom_segment(data = br, aes(x = x, xend = x, y = y, yend = y - tick),
                     inherit.aes = FALSE, linewidth = 0.4, colour = INK3) +
        geom_segment(data = br, aes(x = xend, xend = xend, y = y, yend = y - tick),
                     inherit.aes = FALSE, linewidth = 0.4, colour = INK3) +
        geom_text(data = cbind(br, lab = lab),
                  aes(x = (x + xend)/2, y = y, label = lab), inherit.aes = FALSE,
                  vjust = -0.45, size = 3.1, colour = INK, family = MS_FONT)
      top <- max(yb)
    }
  }

  # breaks are chosen on the data range, not the bracket headroom, so the
  # axis keeps a consistent tick density.
  p + coord_cartesian(ylim = c(base, top + span * 0.30), clip = "off") +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4)(c(base, data_top)),
                       expand = expansion(mult = c(0.01, 0))) +
    theme_fig4(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     face = "bold", size = 9.5, colour = "black"),
          plot.margin = margin(4, 8, 2, 6))
}

# ---------------------------------------------------------------------------
# SCATTER PANEL (Fig 5, Supp 7)
#
# `k` is the analyte's number of distinct raw values. At 6 or fewer, Spearman's
# rho is computed on a handful of levels and can be inflated by quantisation,
# so the panel is labelled "[k levels]".
# ---------------------------------------------------------------------------
scat <- function(d, xv, yv, xlab, ylab, g, rho, p, q = NA, ribbon = TRUE,
                 ann_pos = "above", k = NA) {
  cols <- pal(g)
  # long axis titles are wrapped so neighbouring panels cannot collide
  if (is.character(xlab)) xlab <- paste(strwrap(gsub("\n", " ", xlab), 28), collapse = "\n")
  kn <- if (!is.na(k) && k <= 6) sprintf("  [%d levels]", k) else ""
  # two lines, so neighbouring panels' statistics cannot run into each other
  ann <- sprintf("rho = %.2f  p = %s%s%s", rho, fmt_p(p),
                 if (is.na(q)) "" else sprintf("\nq = %s%s", signif(q, 3), stars(q)), kn)
  ggplot(d, aes(.data[[xv]], .data[[yv]])) +
    geom_smooth(method = "lm", formula = y ~ x, se = ribbon, colour = "black",
                fill = "grey80", linewidth = 0.5, linetype = "dashed", alpha = 0.55) +
    geom_point(aes(fill = Flicker, shape = Sex), colour = "black",
               size = 3.1, stroke = 0.45) +
    scale_fill_manual(values = cols, name = "Treatment",
                      labels = c(NoFlicker = "No Flicker", "20Hz" = "20Hz", "40Hz" = "40Hz")) +
    scale_shape_manual(values = SEX_SHAPE, name = "Sex") +
    # the "inside" annotation sits top-left, so buy it clear space rather than
    # letting it land on the data
    scale_y_continuous(expand = expansion(
      mult = c(0.05, if (ann_pos == "inside") 0.38 else 0.05))) +
    guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 3.4)),
           shape = guide_legend(order = 2,
                     override.aes = list(fill = "grey40", colour = "black", size = 3.4))) +
    labs(x = xlab, y = ylab,
         subtitle = if (ann_pos == "above") ann else NULL,
         title    = if (ann_pos == "inside") stars(q) else NULL) +
    (if (ann_pos == "inside")
       annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.15,
                label = sprintf("rho = %.2g\np = %s%s", rho, fmt_p(p), sub("^  ", "\n", kn)),
                size = 2.8, fontface = "italic", colour = INK)
     else NULL) +
    theme_ms(10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
          axis.line = element_blank(),
          plot.title = element_text(size = 13, colour = INK, hjust = 0.5,
                                    margin = margin(b = 0)),
          plot.subtitle = element_text(size = 8.4, colour = INK, hjust = 0.5,
                                       margin = margin(b = 3)),
          plot.margin = margin(4, 6, 4, 4))
}

# ---------------------------------------------------------------------------
# Composition helpers
# ---------------------------------------------------------------------------
# Subtitles are wrapped to the final (FIG_SCALE'd) width so they never run off
# the canvas; ~15 characters per inch at the 8.4 pt subtitle size.
wrap_sub <- function(x, w) paste(strwrap(x, width = floor(15 * w * FIG_SCALE)), collapse = "\n")

row_fig <- function(panels, title, subtitle, figure, stem, w, h, nrow = 1) {
  subtitle <- wrap_sub(subtitle, w)
  p <- patchwork::wrap_plots(panels, nrow = nrow, guides = "collect") +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
      theme = theme_ms() + theme(plot.title = element_text(size = 15),
                                 plot.subtitle = element_text(size = 8.4, colour = INK2),
                                 plot.margin = margin(8, 10, 4, 8)))
  p <- p & theme(legend.position = "right")
  save_fig(p, figure, stem, w, h)
}

grid_fig <- function(panels, title, subtitle, figure, stem, ncol,
                     panel_w = 3.7, panel_h = 4.2) {
  nr <- ceiling(length(panels) / ncol)
  subtitle <- wrap_sub(subtitle, panel_w * ncol + 1.2)
  p <- patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
      theme = theme_ms(12) +
        theme(text = element_text(family = MS_FONT),
              plot.title = element_text(family = MS_FONT, face = "bold", size = 15,
                                        hjust = 0, colour = INK, margin = margin(b = 4)),
              plot.subtitle = element_text(family = MS_FONT, size = 8.4, colour = INK2,
                                           hjust = 0, margin = margin(b = 6))))
  p <- p & theme(legend.position = "right",
                 legend.text = element_text(family = MS_FONT, size = 9, colour = INK))
  save_fig(p, figure, stem, panel_w * ncol + 1.2, panel_h * nr + 0.9)
}

# Look up one row of a correlation table. Returns NULL when the pair is absent
# rather than a zero-row data.frame, so callers can test it directly.
lookup <- function(tbl, g, analyte, target) {
  r <- tbl[tbl$subgroup == g & tbl$analyte == analyte & tbl$target == target, ]
  if (!nrow(r)) NULL else r[1, ]
}

# The PLS-DA key cytokines for a genotype, from the stage-02 cache.
key_cytokines <- function(g) readRDS(file.path(STATS, sprintf("plsda_%s.rds", g)))$key

# ============================================================================
# 06_amyloid.R
# Aβ40 / Aβ42 (and plaque load) across flicker in 5xFAD (Fig 4B).
# [CORRECTION] The original drew a bare global ANOVA p (no post-hoc) on the
# amyloid bars. Here it's Anderson-Darling -> ANOVA+Tukey or KW+Dunn's, pooled
# and split by sex, with corrected pairwise p on the plot and in the stats table.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(ggpubr) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")

processed_dir <- file.path(results_dir, "processed")
animals <- readRDS(file.path(processed_dir, "animal_table.rds"))

# amyloid is only measured in 5xFAD
fad <- animals[animals$Genotype == "5XFAD", , drop = FALSE]
amy <- intersect(amyloid_vars, colnames(fad))

# ----------------------------------------------------------------------------
# 1. Stats: pooled sex + by sex.
# ----------------------------------------------------------------------------
om <- list(); ph <- list()
for (st in c(list(list(s = NA)), lapply(sexes, function(s) list(s = s)))) {
  d <- fad; lab <- "5XFAD"
  if (!is.na(st$s)) { d <- d[d$Sex == st$s, , drop = FALSE]; lab <- paste("5XFAD", st$s) }
  r <- compare_flicker_many(d, amy, subgroup = lab)
  om[[lab]] <- r$omnibus; ph[[lab]] <- r$posthoc
}
write_stats(do.call(rbind, om), "amyloid_omnibus")
write_stats(do.call(rbind, ph[!vapply(ph, is.null, logical(1))]), "amyloid_posthoc")

# ----------------------------------------------------------------------------
# 2. Bar/violin plots with corrected pairwise p.
# ----------------------------------------------------------------------------
plot_amyloid <- function(d, v, lab) {
  r <- compare_flicker(d, v, subgroup = lab, variable = v)
  g <- ggplot(d, aes(x = Flicker, y = .data[[v]])) +
    geom_bar(aes(fill = Flicker), stat = "summary", fun = mean, alpha = 0.6,
             show.legend = FALSE) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(aes(shape = Sex), width = 0.12, size = 1.8) +
    scale_fill_manual(values = flicker_colors_fad) +
    scale_shape_manual(values = sex_shapes) +
    labs(x = NULL, y = v, title = v,
         subtitle = sprintf("5XFAD: %s (p=%s)", r$omnibus$test, fmt(r$omnibus$p))) +
    theme_bw(base_size = 11)
  if (!is.null(r$posthoc)) {
    sig <- r$posthoc[!is.na(r$posthoc$p_adj) & r$posthoc$p_adj < alpha, ]
    if (nrow(sig)) {
      sig$y.position <- max(d[[v]], na.rm = TRUE) * (1 + 0.08 * seq_len(nrow(sig)))
      sig$p.lab <- paste0("p=", fmt(sig$p_adj))
      g <- g + stat_pvalue_manual(sig, label = "p.lab", xmin = "group1",
                                  xmax = "group2", tip.length = 0.01, size = 3)
    }
  }
  g
}

plots <- lapply(amy, function(v) plot_amyloid(fad, v, "5XFAD"))
ggsave(file.path(amy_dir, "amyloid_bars_5XFAD.pdf"),
       ggpubr::ggarrange(plotlist = plots, ncol = 3,
                         nrow = ceiling(length(plots) / 3), common.legend = TRUE),
       width = 12, height = 3.4 * ceiling(length(plots) / 3))
cat("[06] amyloid:", length(amy), "measures ->", basename(amy_dir), "\n")

# ============================================================================
# 11_supp_myeloid_exploratory.R
# EXPLORATORY supplement: the four myeloid analytes the linear-range prune
# removes (G-CSF, GM-CSF, M-CSF, LIF).
#
# Why this figure exists. The published version of this manuscript built its
# myeloid story on these four analytes. They are absent from the Wood lab's own
# in-linear-range panel, and in these samples they sit on the assay floor: each
# takes 5-9 distinct values across all 43 animals, spanning 2.5-5 MFI above
# background. A correlation or group test on a variable with six levels is
# dominated by ties, so they are excluded from the primary analysis.
#
# They are shown here rather than deleted, so that (a) the myeloid thread stays
# visible and (b) the detection problem is stated by us rather than discovered by
# a reviewer. Everything in this figure is EXPLORATORY and labelled as such.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(ggpubr) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")

pdir    <- file.path(results_dir, "processed")
animals <- readRDS(file.path(pdir, "animal_table.rds"))
lr      <- read.csv(file.path(stats_dir, "analyte_linear_range.csv"), check.names = FALSE)
outdir  <- file.path(results_dir, "revised_figures_png")
DPI     <- 220

# the myeloid analytes the prune removed, in the order the manuscript names them
excluded_myeloid <- intersect(c("GM-CSF", "G-CSF", "LIX", "LIF", "M-CSF"),
                              lr$analyte[!lr$used])
# a retained analyte to anchor the detection panel visually
anchor <- "IL-2"

pal_of <- function(g) if (g == "WT") flicker_colors_wt else flicker_colors_fad
stars  <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***",
                      ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
theme_ms <- function(base = 12) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

# ---------------------------------------------------------------------------
# 1. Group comparisons, same corrected tests as the primary analysis
# ---------------------------------------------------------------------------
rows <- list()
for (g in genotypes) {
  d <- animals[animals$Genotype == g, , drop = FALSE]
  r <- compare_flicker_many(d, excluded_myeloid, subgroup = g)
  rows[[paste0(g, "_om")]] <- r$omnibus
  rows[[paste0(g, "_ph")]] <- r$posthoc
}
om <- do.call(rbind, rows[grepl("_om$", names(rows))])
ph <- do.call(rbind, rows[grepl("_ph$", names(rows))])
write_stats(om, "myeloid_exploratory_omnibus")
write_stats(ph, "myeloid_exploratory_posthoc")

ms_bar <- function(d, var, g) {
  p <- compare_flicker(d, var, subgroup = g, variable = var)$posthoc
  gg <- ggplot(d, aes(x = Flicker, y = .data[[var]], fill = Flicker)) +
    stat_summary(fun = mean, geom = "col", color = "black", width = 0.7, alpha = 0.85) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(aes(shape = Sex), width = 0.12, size = 1.6, stroke = 0.4) +
    scale_fill_manual(values = pal_of(g)) + scale_shape_manual(values = sex_shapes) +
    labs(x = NULL, y = paste0(var, " (a.u.)"), title = paste(g, var)) + theme_ms()
  if (!is.null(p) && nrow(p)) {
    rng <- max(d[[var]], na.rm = TRUE) - min(0, min(d[[var]], na.rm = TRUE))
    p$y.position <- max(d[[var]], na.rm = TRUE) + rng * (0.10 + 0.16 * (seq_len(nrow(p)) - 1))
    p$lab <- paste0(fmt(p$p_adj), stars(p$p_adj))
    gg <- gg + stat_pvalue_manual(p, label = "lab", xmin = "group1", xmax = "group2",
                                  tip.length = 0.01, size = 3.1, bracket.size = 0.35)
  }
  gg
}

plots <- unlist(lapply(genotypes, function(g) {
  d <- animals[animals$Genotype == g, , drop = FALSE]
  lapply(excluded_myeloid, function(a) ms_bar(d, a, g))
}), recursive = FALSE)
nc  <- length(excluded_myeloid)
pan <- ggarrange(plotlist = plots, ncol = nc, nrow = 2, common.legend = TRUE, legend = "top")
pan <- annotate_figure(pan,
  top = text_grob("EXPLORATORY — myeloid analytes outside the assay's linear range",
                  face = "bold", size = 14),
  bottom = text_grob(paste(
    "These analytes are excluded from the primary analysis. Each takes only 5-9 distinct values across all 43 animals",
    "(see detection panel), so group tests and correlations on them are dominated by ties. Shown for continuity with the",
    "previously reported myeloid results; p-values are Tukey / Dunn-BH, uncorrected across analytes.", sep = "\n"),
    size = 8.5, color = "grey25"))
ggsave(file.path(outdir, "Supp9_myeloid_exploratory_bars.png"), pan,
       width = 3.4 * nc, height = 8.2, dpi = DPI)

# ---------------------------------------------------------------------------
# 2. Detection panel — why they are excluded, in one picture
# ---------------------------------------------------------------------------
show <- c(excluded_myeloid, anchor)
dd <- animals %>%
  select(Genotype, all_of(show)) %>%
  pivot_longer(-Genotype, names_to = "Analyte", values_to = "MFI") %>%
  mutate(Analyte = factor(Analyte, levels = show),
         panel = ifelse(Analyte == anchor,
                        "retained (in linear range)", "excluded (floor)"))
lab <- lr %>% filter(analyte %in% show) %>%
  transmute(Analyte = factor(analyte, levels = show),
            txt = sprintf("%d distinct\nvalues", n_distinct))

p_det <- ggplot(dd, aes(x = Analyte, y = MFI, color = panel)) +
  geom_jitter(width = 0.18, height = 0, size = 1.5, alpha = 0.75) +
  geom_text(data = lab, aes(x = Analyte, y = Inf, label = txt), inherit.aes = FALSE,
            vjust = 1.2, size = 3, color = "grey25") +
  scale_color_manual(values = c("excluded (floor)" = "#B2182B",
                                "retained (in linear range)" = "#2166AC"), name = NULL) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 50, 200, 800),
                     expand = expansion(mult = c(0.04, 0.22))) +
  labs(x = NULL, y = "Background-subtracted MFI (log scale)",
       title = "Detection range of the excluded myeloid analytes",
       subtitle = sprintf("one point = one animal (n = %d); the reader quantises in 0.5 MFI steps",
                          nrow(animals))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
        axis.text = element_text(color = "black"),
        legend.position = "top")
ggsave(file.path(outdir, "Supp9_myeloid_exploratory_detection.png"), p_det,
       width = 7.5, height = 5, dpi = DPI)

cat("[11] exploratory myeloid supplement:", paste(excluded_myeloid, collapse = ", "), "\n")
cat("[11] -> Supp9_myeloid_exploratory_{bars,detection}.png\n")

# ---------------------------------------------------------------------------
# 3. The analytes that DO respond to flicker (replaces the old Supp Fig 6,
#    which was built almost entirely from analytes now outside the linear range:
#    Eotaxin, IL-5, IL-6, IL-7 and IL-15 are all excluded).
#    Within the retained panel only two contrasts clear 0.05, and both are here.
# ---------------------------------------------------------------------------
wt <- animals[animals$Genotype == "WT", , drop = FALSE]
p_ifny <- ms_bar(wt, "IFNy", "WT") + labs(title = "WT IFN-γ (sexes pooled)")

lixF <- wt[wt$Sex == "F", , drop = FALSE]
p_lix <- ms_bar(lixF, "LIX", "WT") + labs(title = "WT LIX (females only)")

pan6 <- ggarrange(p_ifny, p_lix, ncol = 2, common.legend = TRUE, legend = "top")
pan6 <- annotate_figure(pan6,
  top = text_grob("Flicker-responsive analytes within the assay's linear range",
                  face = "bold", size = 13),
  bottom = text_grob(paste(
    "The only two contrasts reaching p < 0.05 among the 16 retained analytes.",
    "IFN-γ: No Flicker vs 20Hz and vs 40Hz, both p = 0.044 (Kruskal-Wallis, Dunn-BH).",
    "LIX: 20Hz vs 40Hz in females, p = 0.048 (ANOVA, Tukey). Sexes pooled unless stated.", sep = "\n"),
    size = 8.5, color = "grey25"))
ggsave(file.path(outdir, "Supp6_flicker_responsive_WT.png"), pan6,
       width = 8.5, height = 5.4, dpi = DPI)
cat("[11] -> Supp6_flicker_responsive_WT.png\n")

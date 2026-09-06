# ============================================================================
# 05_cytokine_group_stats.R
# Per-analyte flicker comparisons within genotype (Supp 5-6), the core of the
# cytokine corrections.
# [CORRECTION] The original used raw pairwise t-tests / a bare global ANOVA with
# no post-hoc, and drew uncorrected p on the plots. Here each analyte goes
# through: Anderson-Darling -> ANOVA + Tukey (normal) or Kruskal-Wallis + Dunn's
# (non-normal), and only the *corrected* pairwise p is ever plotted. All tests,
# statistics, N and correction method land in results/stats/ so the Supplementary
# Statistics table can be regenerated instead of hand-copied.
# Uses RAW (not sex-adjusted) values; sex is handled by splitting, per the
# manuscript, and pooled-sex results are reported alongside.
#
# [REVISED SPEC] Two changes:
#   - the figure panel is now the per-genotype KEY CYTOKINES from the PLS-DA
#     loadings (stage 04b), not the a priori myeloid panel;
#   - Tukey / Dunn IS the multiple-comparison correction. No second FDR layer
#     across analytes is used for inference. p_omnibus_fdr stays in the CSV as a
#     diagnostic column only, and nothing in the manuscript should cite it.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(ggpubr) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")

processed_dir     <- file.path(results_dir, "processed")
animals           <- readRDS(file.path(processed_dir, "animal_table.rds"))
analysis_analytes <- readRDS(file.path(processed_dir, "analysis_analytes.rds"))
key_cytokines     <- readRDS(file.path(processed_dir, "key_cytokines.rds"))

# ----------------------------------------------------------------------------
# 1. Stats: pooled-sex within genotype, and split by sex.
# ----------------------------------------------------------------------------
all_om <- list(); all_ph <- list()
strata <- c(list(list(g = "WT", s = NA), list(g = "5XFAD", s = NA)),
            unlist(lapply(genotypes, function(g)
              lapply(sexes, function(s) list(g = g, s = s))), recursive = FALSE))

for (st in strata) {
  d <- animals[animals$Genotype == st$g, , drop = FALSE]
  lab <- st$g
  if (!is.na(st$s)) { d <- d[d$Sex == st$s, , drop = FALSE]; lab <- paste(st$g, st$s) }
  r <- compare_flicker_many(d, analysis_analytes, subgroup = lab)
  all_om[[lab]] <- r$omnibus; all_ph[[lab]] <- r$posthoc
}
omnibus <- do.call(rbind, all_om)
posthoc <- do.call(rbind, all_ph[!vapply(all_ph, is.null, logical(1))])

# DIAGNOSTIC ONLY: BH-FDR across analytes within each stratum on the omnibus p.
# Reported so the extra multiplicity layer is visible, but per the revised spec
# the corrected p used for inference is the Tukey / Dunn pairwise value.
omnibus$p_omnibus_fdr <- NA_real_
for (lab in unique(omnibus$subgroup)) {
  i <- omnibus$subgroup == lab & !is.na(omnibus$p)
  omnibus$p_omnibus_fdr[i] <- p.adjust(omnibus$p[i], method = "BH")
}
write_stats(omnibus, "cytokine_omnibus")
write_stats(posthoc, "cytokine_posthoc")

# ----------------------------------------------------------------------------
# 2. Plots: violin + box across flicker, per genotype (pooled sex), annotated
#    with the CORRECTED pairwise p (Tukey/Dunn). The key cytokines (top 3 on
#    LV1 + top 3 on LV2, per genotype) get their own figure.
# ----------------------------------------------------------------------------
plot_analyte <- function(d, analyte, pal, lab) {
  ph <- compare_flicker(d, analyte, subgroup = lab, variable = analyte)
  g  <- ggplot(d, aes(x = Flicker, y = .data[[analyte]])) +
    geom_violin(aes(fill = Flicker), alpha = 0.35, scale = "width",
                trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(shape = Sex), width = 0.12, size = 1.6) +
    scale_fill_manual(values = pal) + scale_shape_manual(values = sex_shapes) +
    labs(x = NULL, y = analyte,
         title = analyte,
         subtitle = sprintf("%s: %s (p=%s)", lab,
                            ph$omnibus$test, fmt(ph$omnibus$p))) +
    theme_bw(base_size = 11)
  # annotate only the corrected pairwise p's that are significant
  if (!is.null(ph$posthoc)) {
    sig <- ph$posthoc[!is.na(ph$posthoc$p_adj) & ph$posthoc$p_adj < alpha, ]
    if (nrow(sig)) {
      sig$y.position <- max(d[[analyte]], na.rm = TRUE) *
        (1 + 0.08 * seq_len(nrow(sig)))
      sig$p.lab <- paste0("p=", fmt(sig$p_adj))
      g <- g + stat_pvalue_manual(sig, label = "p.lab",
                                  xmin = "group1", xmax = "group2",
                                  tip.length = 0.01, size = 3)
    }
  }
  g
}

for (g in genotypes) {
  d    <- animals[animals$Genotype == g, , drop = FALSE]
  pal  <- if (g == "WT") flicker_colors_wt else flicker_colors_fad
  keys <- intersect(key_cytokines[[g]], analysis_analytes)
  plots <- lapply(keys, function(a) plot_analyte(d, a, pal, g))
  ggsave(file.path(cyt_dir, sprintf("key_cytokine_panel_%s.pdf", g)),
         ggpubr::ggarrange(plotlist = plots, ncol = 3,
                           nrow = ceiling(length(plots) / 3), common.legend = TRUE),
         width = 12, height = 3.2 * ceiling(length(plots) / 3))
  cat(sprintf("[05] %-6s key-cytokine panel: %s\n", g, paste(keys, collapse = ", ")))
  # drop the previous round's myeloid panel so it can't be mistaken for current
  old <- file.path(cyt_dir, sprintf("myeloid_panel_%s.pdf", g))
  if (file.exists(old)) file.remove(old)
}

# console summary: the corrected pairwise results for the key cytokines, which
# are the numbers the manuscript reports.
cat("[05] key-cytokine pairwise contrasts (Tukey / Dunn-BH), p <", alpha, ":\n")
key_ph <- posthoc[posthoc$subgroup %in% genotypes &
                  mapply(function(sg, v) v %in% key_cytokines[[sg]],
                         posthoc$subgroup, posthoc$variable) &
                  !is.na(posthoc$p_adj) & posthoc$p_adj < alpha,
                  c("subgroup", "variable", "group1", "group2", "p_adj", "posthoc")]
if (nrow(key_ph)) print(key_ph, row.names = FALSE) else cat("   (none)\n")
cat("[05] cytokine group stats done ->", stats_dir, "\n")

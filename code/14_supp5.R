# ============================================================================
# 14_supp5.R  —  Supplementary Figure 5: key cytokine group comparisons,
# A = WT, B = 5xFAD. Styled to match the published supplementary bar panels.
#
# ANALYTES ARE THE PLS-DA KEY CYTOKINES for that genotype, read from the
# stage-02 cache -- top 3 by |loading| on LV1 and on LV2. This figure replaces
# the published myeloid panel entirely: G-CSF / GM-CSF / M-CSF / LIF were
# presented there as a curated biological group, which is not a selection rule
# the analysis can defend. Whichever analytes carry the separation is.
#
# Statistics are stage 06's: Anderson-Darling, then ANOVA + Tukey or
# Kruskal-Wallis + Dunn(BH), on the sex-z-scored pooled table. Every pairwise
# comparison is bracketed with its adjusted p, significant or not, as the
# published panels do.
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

for (g in GENOTYPES) {
  d  <- ANIMALS[ANIMALS$Genotype == g, , drop = FALSE]
  ks <- key_cytokines(g)
  rr <- compare_flicker_many(d, ks, subgroup = g)
  ps <- lapply(ks, function(v)
    bar_panel(d, v, g, rr$omnibus[rr$omnibus$variable == v, ], rr$posthoc,
              title_expr = v, ylab = NULL))
  sig <- rr$posthoc[rr$posthoc$p_adj < 0.05, , drop = FALSE]
  grid_fig(ps, sprintf("%s Key Cytokines", g),
           sprintf("PLS-DA top 3 on LV1 and LV2; sex z-scored; AD -> ANOVA+Tukey or KW+Dunn(BH). %s",
                   if (nrow(sig)) sprintf("%d contrast(s) reach 0.05", nrow(sig))
                   else "Nothing reaches 0.05."),
           "Supp5", sprintf("Supp5_key_cytokines_%s", g), ncol = 3)
  cat(sprintf("[14] %-6s keys: %s\n", g, paste(ks, collapse = ", ")))
}
cat("[14] Supplementary Figure 5 done\n")

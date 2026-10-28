# ============================================================================
# 16_supp7.R  —  Supplementary Figure 7: WT cytokine ~ morphology correlations.
#
# Supplementary Figure 7 IS Figure 5C/D on WT, so it uses the same rule: the
# WT PLS-DA KEY CYTOKINES from the stage-02 cache, against each of the four
# morphology targets. The published panel named its analytes directly; this one
# derives them, so the selection follows the analysis.
#
# Published style: no CI ribbon, dashed black trend line, statistics italic
# inside the panel at top-left, significance star centred above, black panel
# box, no grid. Sex is point shape, as everywhere here.
#
# The subtitle states the best q on each target so "key" is never mistaken for
# "significant" -- the two are different claims and this figure often shows
# analytes that carry separation but no correlation.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "03_theme.R"))
source(file.path(HERE, "07_figure_helpers.R"))

CM <- read.csv(file.path(STATS, "corr_cytokine_morphology.csv"), check.names = FALSE)
wt <- ANIMALS[ANIMALS$Genotype == "WT", , drop = FALSE]
KEYW <- key_cytokines("WT")
cat("[16] WT key cytokines:", paste(KEYW, collapse = ", "), "\n")

for (cfg in list(
  list(v = "AST mean_Distance to Nearest Neighbor", ttl = "Key Cytokines vs Astrocyte Tiling",
       lab = "Astrocyte Avg Distance to Nearest Neighbor", stem = "Supp7A_cyt_vs_astro_tiling_WT"),
  list(v = "MG mean_Distance to Nearest Neighbor",  ttl = "Key Cytokines vs Microglia Tiling",
       lab = "Microglia Avg Distance to Nearest Neighbor", stem = "Supp7B_cyt_vs_mg_tiling_WT"),
  list(v = "AST mean_Number of Cells", ttl = "Key Cytokines vs Avg Number of Astrocytes",
       lab = "Average Number of Astrocytes", stem = "Supp7C_cyt_vs_astro_number_WT"),
  list(v = "MG mean_Number of Cells",  ttl = "Key Cytokines vs Avg Number of Microglia",
       lab = "Average Number of Microglia", stem = "Supp7D_cyt_vs_mg_number_WT"))) {

  ps <- lapply(KEYW, function(cy) {
    r <- lookup(CM, "WT", cy, cfg$v)
    scat(wt, cfg$v, cy, cfg$lab, zlab(cy), "WT", r$rho, r$p_perm, r$p_adj,
         ribbon = FALSE, ann_pos = "inside", k = r$n_distinct)
  })
  qs <- CM$p_adj[CM$subgroup == "WT" & CM$target == cfg$v & CM$analyte %in% KEYW]
  best <- min(qs, na.rm = TRUE)
  sub <- sprintf(paste("PLS-DA key cytokines (top 3 on LV1 and LV2); BH q across the",
                       "%d-analyte family, best q here = %.2f%s"),
                 length(PANEL), best,
                 if (best >= 0.05) " - nothing reaches q < 0.05, so no panel is starred" else "")
  grid_fig(ps, cfg$ttl, sub, "Supp7", cfg$stem, ncol = 3, panel_w = 3.1, panel_h = 3.3)
  cat(sprintf("[16] %-34s best q = %.3f\n", cfg$stem, best))
}
cat("[16] Supplementary Figure 7 done\n")

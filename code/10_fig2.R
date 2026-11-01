# ============================================================================
# 10_fig2.R  —  Figure 2 B / C / D, WT.
#
# The drawing lives in 04_panels.R so Figure 2 and Supplementary Figure 4
# cannot drift apart -- they are the same figure on two genotypes. Read that
# file's header for what was matched from the published panel.
#
# Figure 2 passes no heatmap annotation columns: WT has no amyloid, and the
# published Figure 2 has none. Supp 4A does -- see 13_supp4.R.
#
# Outputs, all from the cached stage-02 fit:
#   Fig2B_heatmap_WT            28-analyte cytokine heatmap
#   Fig2C_plsda_scores_WT       scores plot alone
#   Fig2C_plsda_violins_WT      LV1 over LV2 violins, NO statistics
#   Fig2C_plsda_composite_WT    scores + violins, the published panel C layout
#   Fig2D_loadings_WT_LV1/LV2   signals in each latent variable
#   Fig2E_LV1_scores_by_group_WT  LV1 score per animal by group, NO statistics
#                               (our version of Garza et al. Fig 1E)
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "03_theme.R"))
source(file.path(HERE, "04_panels.R"))

ms_heatmap("WT", "Fig2", "Fig2B_heatmap_WT", col_order = "lv1")   # columns follow Fig 2D LV1
ms_plsda("WT", "Fig2", list(
  scores    = "Fig2C_plsda_scores_WT",
  violins   = "Fig2C_plsda_violins_WT",
  composite = "Fig2C_plsda_composite_WT",
  loadings  = "Fig2D_loadings_WT_LV%d",
  lv1_groups = "Fig2E_LV1_scores_by_group_WT"))

cat("[10] Figure 2 done\n")

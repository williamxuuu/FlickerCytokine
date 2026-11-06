# ============================================================================
# 13_supp4.R  —  Supplementary Figure 4 A-E: 5xFAD cytokine profile and PLS-DA.
#
# The 5xFAD counterpart of Figure 2, drawn by the same functions in
# 04_panels.R.
#
#   4A  28-analyte heatmap, columns in LV1 order, with Abeta 1-42 and
#       Abeta 1-40 annotation columns
#   4B  LV1 vs LV2 scores, with LV1 / LV2 score distributions
#   4C / 4D  loadings on LV1 / LV2
#   4E  LV1 score per animal, mean +/- SEM
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "03_theme.R"))
source(file.path(HERE, "04_panels.R"))

ms_heatmap("5XFAD", "Supp4", "Supp4A_heatmap_5XFAD", annot = c("Abeta 42", "Abeta 40"),
           col_order = "lv1", height = 8.6)   # taller: legend carries two Abeta colour bars
ms_plsda("5XFAD", "Supp4", list(
  scores    = "Supp4B_plsda_scores_5XFAD",
  violins   = "Supp4B_plsda_violins_5XFAD",
  composite = "Supp4B_plsda_composite_5XFAD",
  loadings  = "Supp4CD_loadings_5XFAD_LV%d",
  lv1_groups = "Supp4E_LV1_scores_by_group_5XFAD"))

cat("[13] Supplementary Figure 4 done\n")

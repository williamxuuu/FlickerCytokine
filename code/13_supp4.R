# ============================================================================
# 13_supp4.R  —  Supplementary Figure 4 A / B / C / D, 5xFAD.
#
# Supplementary Figure 4 IS Figure 2 drawn on 5xFAD: the same cytokine heatmap,
# the same PLS-DA scores, the same LV score violins, the same LV1/LV2 loadings.
# Every panel comes from the shared 04_panels.R, so the two figures cannot
# diverge in style again.
#
#   4A  cytokine heatmap, 28 analytes, columns in LV1-loading order (matching
#       4C), with the Abeta 1-42 and Abeta 1-40
#       annotation columns the published 4A carries, on their own ramps
#       sampled from that panel's colour bars.
#   4B  PLS-DA scores (+ the score violins, no statistics).
#   4C / 4D  Signals in LV1 / LV2.
#
# WORTH KNOWING. The published caption is "Flicker reduces expression of IL5 in
# 5xFAD mice", but IL-5 is one of the four analytes the >80%-above-background
# gate drops, so it is not in this panel and cannot be. That is a claim to
# settle in the text; the figure cannot resolve it. See 00_setup.R.
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
           col_order = "lv1")   # columns follow Supp 4C LV1
ms_plsda("5XFAD", "Supp4", list(
  scores    = "Supp4B_plsda_scores_5XFAD",
  violins   = "Supp4B_plsda_violins_5XFAD",
  composite = "Supp4B_plsda_composite_5XFAD",
  loadings  = "Supp4CD_loadings_5XFAD_LV%d"))

cat("[13] Supplementary Figure 4 done\n")

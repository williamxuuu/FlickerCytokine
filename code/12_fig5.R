# ============================================================================
# 12_fig5.R  —  Figure 5 A / B / D + supplement (5xFAD), matched to the published
# Figure 5 (manuscript/ManuscriptMostRecent.docx, word/media/image5.png).
#
#   5A / 5B  amyloid burden vs glial tiling. These are NAMED comparisons -- the
#            figure's whole claim is about Abeta against astrocyte and
#            microglial tiling -- so the targets are fixed by the question, not
#            selected. That is not a hardcoded analyte list.
#   5D       cytokines vs microglial number, chosen by multiple-testing-
#            corrected significance (see the selection rule below). 5C is
#            vacated; the astrocyte-number and key-cytokine panels are in
#            figures/Fig5_supp/.
#
# ANALYTE SELECTION IS NOT HARDCODED. 5D reads its analytes from the stage-05
# statistics, and the supplement's key-cytokine panels from the stage-02
# PLS-DA, so if the panel, the z-score or the model changes, these follow.
#
# Statistics come from stats/corr_cytokine_morphology.csv (stage 05):
# permutation Spearman on the sex-z-scored table, BH q across the 28-analyte
# family for that target. Trend lines are linear for display only, which is the
# Wood lab's convention (initialComments.md:15).
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "03_theme.R"))
source(file.path(HERE, "07_figure_helpers.R"))

set.seed(42); B <- 50000
CM  <- read.csv(file.path(STATS, "corr_cytokine_morphology.csv"), check.names = FALSE)
fad <- ANIMALS[ANIMALS$Genotype == "5XFAD", , drop = FALSE]

AST_TILE <- "AST mean_Distance to Nearest Neighbor"
MG_TILE  <- "MG mean_Distance to Nearest Neighbor"
AST_N    <- "AST mean_Number of Cells"
MG_N     <- "MG mean_Number of Cells"

# Permutation p for a pair that is not in the cytokine x target tables -- used
# only by 5A/5B, which correlate amyloid against morphology. Free permutation,
# matching stage 05: sex is already removed by the z-score.
perm_p <- function(d, xv, yv) {
  ok <- !is.na(d[[xv]]) & !is.na(d[[yv]]); x <- d[[xv]][ok]; y <- d[[yv]][ok]
  n <- sum(ok)
  r  <- suppressWarnings(cor(x, y, method = "spearman"))
  rx <- rank(x); ry <- rank(y)
  nul <- replicate(B, cor(rx, ry[sample.int(n)]))
  list(rho = r, p = (1 + sum(abs(nul) >= abs(r) - 1e-12)) / (B + 1), n = n)
}

# ============================ Fig 5A / 5B ===================================
for (cfg in list(list(v = AST_TILE, lab = "Astrocyte Average\nDistance to Nearest Neighbor",
                      ttl = "Pathology vs Avg Astrocyte Tiling",
                      stem = "Fig5A_amyloid_vs_astro_tiling"),
                 list(v = MG_TILE,  lab = "Microglia Average\nDistance to Nearest Neighbor",
                      ttl = "Pathology vs Avg Microglia Tiling",
                      stem = "Fig5B_amyloid_vs_mg_tiling"))) {
  ps <- lapply(AMYLOID_TARGETS, function(ab) {
    s <- perm_p(fad, cfg$v, ab)
    scat(fad, cfg$v, ab, cfg$lab, ab_lab(ab), "5XFAD", s$rho, s$p)
  })
  row_fig(ps, cfg$ttl, "sex z-scored; permutation Spearman, B = 50,000",
          "Fig5", cfg$stem, 9.8, 5.6)
}

# ============================ Fig 5D + supplement ===========================
# SELECTION RULE (2026-11, PI decision): Fig 5D shows the cytokines whose
# correlation with microglial number survives BOTH Benjamini-Hochberg q < 0.05
# across all 28 analytes AND Westfall-Young FWER p < 0.05 (stage 05). All 28
# were tested and both corrections count all 28 tests, so selecting on them is
# not cherry-picking -- the caption must state the rule. Pairs passing BH only
# go to the supplement. The PLS-DA key-cytokine panels that used to be 5C/5D
# (all null at q < 0.05) move to the supplement as well, for transparency.
#
# 5C is vacated: no cytokine passes BH against astrocyte number after the sex
# z-score, so that panel now lives in the supplement. Re-letter in Illustrator.
unlink(Sys.glob(paste0(fig_path("Fig5", "Fig5C_cytokine_vs_astro_number"), ".*")))

fam <- CM[CM$subgroup == "5XFAD" & CM$target == MG_N, ]
fam <- fam[order(fam$p_perm), ]
BOTH    <- fam$analyte[fam$p_adj < 0.05 & fam$p_wy <  0.05]
BH_ONLY <- fam$analyte[fam$p_adj < 0.05 & fam$p_wy >= 0.05]
cat("[12] MG number, BH + WY :", paste(BOTH,    collapse = ", "), "\n")
cat("[12] MG number, BH only :", paste(BH_ONLY, collapse = ", "), "\n")

cyt_grid <- function(cys, v, lab, ttl, sub, figure, stem, ncol) {
  if (!length(cys)) { cat("[12]", stem, "- no analytes pass, not drawn\n"); return(invisible()) }
  ps <- lapply(cys, function(cy) {
    r <- lookup(CM, "5XFAD", cy, v)
    scat(fad, v, cy, lab, zlab(cy), "5XFAD", r$rho, r$p_perm, r$p_adj, k = r$n_distinct)
  })
  grid_fig(ps, ttl, sub, figure, stem, ncol = min(ncol, length(ps)),
           panel_w = 4.0, panel_h = 4.4)
}

cyt_grid(BOTH, MG_N, "Average Number of Microglia", "Cytokines vs Avg Number of Microglia",
         sprintf(paste("Cytokines passing BH q < 0.05 across all %d analytes and",
                       "Westfall-Young p < 0.05"), length(PANEL)),
         "Fig5", "Fig5D_cytokine_vs_mg_number", ncol = 2)

cyt_grid(BH_ONLY, MG_N, "Average Number of Microglia", "Cytokines vs Avg Number of Microglia",
         sprintf("Cytokines passing BH q < 0.05 across all %d analytes, not Westfall-Young",
                 length(PANEL)),
         "Fig5_supp", "Fig5S1_cytokine_vs_mg_number_BH_only", ncol = 2)

# --- PLS-DA key cytokines, moved from 5C / 5D --------------------------------
KEYF <- key_cytokines("5XFAD")
cat("[12] 5xFAD key cytokines:", paste(KEYF, collapse = ", "), "\n")
for (cfg in list(list(v = AST_N, lab = "Average Number of Astrocytes",
                      ttl = "Key Cytokines vs Avg Number of Astrocytes",
                      stem = "Fig5S2_key_cytokines_vs_astro_number"),
                 list(v = MG_N,  lab = "Average Number of Microglia",
                      ttl = "Key Cytokines vs Avg Number of Microglia",
                      stem = "Fig5S3_key_cytokines_vs_mg_number"))) {
  best <- min(CM$p_adj[CM$subgroup == "5XFAD" & CM$target == cfg$v &
                       CM$analyte %in% KEYF], na.rm = TRUE)
  cyt_grid(KEYF, cfg$v, cfg$lab, cfg$ttl,
           sprintf(paste("PLS-DA key cytokines (top 3 on LV1 and LV2);",
                         "BH q across the %d-analyte family, best q here = %.3g"),
                   length(PANEL), best),
           "Fig5_supp", cfg$stem, ncol = 3)
}

cat("[12] Figure 5 done\n")

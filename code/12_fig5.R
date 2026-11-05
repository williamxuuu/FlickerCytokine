# ============================================================================
# 12_fig5.R  —  Figure 5 A / B / C / D (5xFAD), matched to the published
# Figure 5 (manuscript/ManuscriptMostRecent.docx, word/media/image5.png).
#
#   5A / 5B  amyloid burden vs glial tiling. These are NAMED comparisons -- the
#            figure's whole claim is about Abeta against astrocyte and
#            microglial tiling -- so the targets are fixed by the question, not
#            selected. That is not a hardcoded analyte list.
#   5C / 5D  KEY CYTOKINES vs astrocyte and microglial density.
#
# ANALYTE SELECTION IS NOT HARDCODED. Both panels show the 5xFAD key cytokines
# published by stage 02 -- top 3 by |loading| on LV1 and on LV2 -- so if the
# panel, the z-score or the model changes, these figures change with them.
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

# ============================ Fig 5C / 5D ===================================
KEYF <- key_cytokines("5XFAD")
cat("[12] 5xFAD key cytokines:", paste(KEYF, collapse = ", "), "\n")

for (cfg in list(list(v = AST_N, lab = "Average Number of Astrocytes",
                      ttl = "Key Cytokines vs Avg Number of Astrocytes",
                      stem = "Fig5C_cytokine_vs_astro_number"),
                 list(v = MG_N,  lab = "Average Number of Microglia",
                      ttl = "Key Cytokines vs Avg Number of Microglia",
                      stem = "Fig5D_cytokine_vs_mg_number"))) {
  ps <- lapply(KEYF, function(cy) {
    r <- lookup(CM, "5XFAD", cy, cfg$v)
    scat(fad, cfg$v, cy, cfg$lab, zlab(cy), "5XFAD",
         r$rho, r$p_perm, r$p_adj, k = r$n_distinct)
  })
  best <- min(CM$p_adj[CM$subgroup == "5XFAD" & CM$target == cfg$v &
                       CM$analyte %in% KEYF], na.rm = TRUE)
  sub <- sprintf(paste("PLS-DA key cytokines (top 3 on LV1 and LV2);",
                       "BH q across the %d-analyte family, best q here = %.3g"),
                 length(PANEL), best)
  grid_fig(ps, cfg$ttl, sub, "Fig5", cfg$stem, ncol = 3, panel_w = 4.0, panel_h = 4.4)
}

cat("[12] Figure 5 done\n")

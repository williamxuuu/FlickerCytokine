# ============================================================================
# 12_fig5.R  —  Figure 5 and supplementary panels: amyloid, cytokines and glial
# organisation (5xFAD).
#
#   5A / 5B  Abeta 1-40 and Abeta 1-42 vs astrocyte (5A) and microglia (5B)
#            tiling (mean distance to nearest neighbour). Permutation Spearman,
#            B = 50,000, on the sex-z-scored table.
#   5C       cytokines vs microglial number, selected by corrected significance
#            (rule below).
#   Fig5_supp/
#     S1     cytokines vs microglial number passing BH only (not Westfall-Young)
#     S2/S3  PLS-DA key cytokines vs astrocyte (S2) and microglial (S3) number
#
# SELECTION RULE FOR 5C. All 28 analytes were correlated with microglial number
# (stage 05). Panel 5C shows those passing BOTH Benjamini-Hochberg q < 0.05 and
# Westfall-Young p < 0.05, each computed across all 28 tests; those passing BH
# only are shown in S1. Because both corrections account for every analyte
# tested, selecting on them introduces no additional multiplicity.
#
# Analytes are read from the stage-02 and stage-05 outputs, never listed by
# name. Trend lines are least-squares fits for display only; the reported
# statistic is Spearman's rho.
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

# Permutation p for amyloid ~ morphology pairs (5A/5B), which are not in the
# stage-05 cytokine tables. Unrestricted permutation, as in stage 05.
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

# ============================ Fig 5C + supplement ===========================
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
         "Fig5", "Fig5C_cytokine_vs_mg_number", ncol = 2)

cyt_grid(BH_ONLY, MG_N, "Average Number of Microglia", "Cytokines vs Avg Number of Microglia",
         sprintf("Cytokines passing BH q < 0.05 across all %d analytes, not Westfall-Young",
                 length(PANEL)),
         "Fig5_supp", "Fig5S1_cytokine_vs_mg_number_BH_only", ncol = 2)

# --- PLS-DA key cytokines vs glial number (S2 / S3) -------------------------
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

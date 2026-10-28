# ============================================================================
# 15_supp6.R  —  Supplementary Figure 6: cytokine profiling, A = WT, B = 5xFAD.
# Matched to the published panel (word/media/image12.png): one FIXED set of six
# analytes shown for BOTH genotypes, 2 rows x 3 columns, all three pairwise
# comparisons bracketed, analyte name as the panel title, no y-axis title.
#
# ---------------------------------------------------------------------------
# DELIBERATE EXCEPTION TO "PLS-DA KEYS EVERYWHERE" -- OPEN QUESTION.
#
# Every other downstream figure selects the PLS-DA key cytokines. This one
# cannot: Supplementary Figure 5 already shows exactly those analytes, per
# genotype, as bar panels with the same statistics, so the same rule here would
# make Supp 6 a pixel-level duplicate of Supp 5.
#
# Supp 6 therefore keeps the published figure's job: ONE analyte set shown
# across BOTH genotypes, so WT and 5xFAD can be read against each other. The set
# is chosen by flicker responsiveness -- the six analytes with the smallest
# omnibus p from stage 06, ranked on their best result across the two genotypes.
# Still data-driven and still cascading (it moves when the panel or the z-score
# moves); it simply ranks on the group test rather than on the LV loadings.
#
# Unresolved: whether to keep Supp 6 as-is, cut it, or point it at the key
# cytokines and cut Supp 5 instead. Changing the rule is a one-line edit below.
# ---------------------------------------------------------------------------
#
# Departure from the published panel: sex is point shape. Published uses circles.
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

N_SUPP6 <- 6

# --- the analyte set: most flicker-responsive, across both genotypes ---------
om <- read.csv(file.path(STATS, "group_omnibus_cytokines.csv"), check.names = FALSE)
best <- tapply(om$p, om$variable, function(p) min(p, na.rm = TRUE))
SUPP6 <- names(sort(best))[seq_len(N_SUPP6)]
SUPP6 <- SUPP6[SUPP6 %in% PANEL]

cat("[15] Supp 6 analytes (most flicker-responsive):", paste(SUPP6, collapse = ", "), "\n")
for (a in SUPP6) cat(sprintf("       %-12s best omnibus p = %.4f\n", a, best[[a]]))

for (g in GENOTYPES) {
  d  <- ANIMALS[ANIMALS$Genotype == g, , drop = FALSE]
  rr <- compare_flicker_many(d, SUPP6, subgroup = g)
  ps <- lapply(SUPP6, function(v)
    bar_panel(d, v, g, rr$omnibus[rr$omnibus$variable == v, ], rr$posthoc,
              title_expr = v, ylab = NULL))
  grid_fig(ps, sprintf("%s Cytokine Profiling", g),
           sprintf("the %d most flicker-responsive analytes of %d, ranked on the omnibus test across both genotypes; sex z-scored",
                   length(SUPP6), length(PANEL)),
           "Supp6", sprintf("Supp6%s_cytokine_profiling_%s",
                            if (g == "WT") "A" else "B", g), ncol = 3)
  sig <- rr$posthoc[rr$posthoc$p_adj < 0.05, , drop = FALSE]
  cat(sprintf("[15] %-6s significant post-hoc: %s\n", g,
      if (nrow(sig)) paste(sprintf("%s %s/%s q=%.3g", sig$variable, sig$group1,
                                   sig$group2, sig$p_adj), collapse = "; ") else "none"))
}
cat("[15] Supplementary Figure 6 done\n")

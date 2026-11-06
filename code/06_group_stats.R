# ============================================================================
# 06_group_stats.R  —  flicker group comparisons, within genotype.
#
# For every panel analyte (and, in 5xFAD, Abeta 40 / Abeta 42), on the
# sex-z-scored table with sexes pooled (WT n = 20, 5xFAD n = 23):
#
#   Anderson-Darling test on model residuals (lib/stats_helpers.R)
#     normal      -> one-way ANOVA, Tukey HSD post hoc
#     non-normal  -> Kruskal-Wallis, Dunn post hoc (Benjamini-Hochberg)
#
# The post hoc procedure controls error across the three pairwise flicker
# contrasts within an analyte. No additional correction is applied across
# analytes; omnibus p-values should be read with the 28-analyte panel size in
# mind.
#
# Outputs: stats/group_omnibus_{cytokines,amyloid}.csv,
#          stats/group_posthoc_{cytokines,amyloid}.csv.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "lib", "stats_helpers.R"))

oms <- list(); phs <- list()
for (g in GENOTYPES) {
  d <- ANIMALS[ANIMALS$Genotype == g, , drop = FALSE]
  r <- compare_flicker_many(d, PANEL, subgroup = g)
  oms[[g]] <- r$omnibus; phs[[g]] <- r$posthoc
}
omnibus <- do.call(rbind, oms)
posthoc <- do.call(rbind, phs)

write.csv(omnibus, file.path(STATS, "group_omnibus_cytokines.csv"), row.names = FALSE)
write.csv(posthoc, file.path(STATS, "group_posthoc_cytokines.csv"), row.names = FALSE)

# --- amyloid (5xFAD only) ----------------------------------------------------
fad <- ANIMALS[ANIMALS$Genotype == "5XFAD", , drop = FALSE]
ra  <- compare_flicker_many(fad, AMYLOID_TARGETS, subgroup = "5XFAD")
write.csv(ra$omnibus, file.path(STATS, "group_omnibus_amyloid.csv"), row.names = FALSE)
write.csv(ra$posthoc, file.path(STATS, "group_posthoc_amyloid.csv"), row.names = FALSE)

# --- summary -----------------------------------------------------------------
cat("\n=== FLICKER EFFECTS, 28-analyte panel, sex z-scored ===\n")
for (g in GENOTYPES) {
  om <- omnibus[omnibus$subgroup == g, ]
  ph <- posthoc[posthoc$subgroup == g & posthoc$p_adj < 0.05, , drop = FALSE]
  cat(sprintf("\n%s: omnibus p < 0.05 on %d of %d analytes\n",
              g, sum(om$p < 0.05, na.rm = TRUE), nrow(om)))
  if (nrow(ph)) {
    for (i in seq_len(nrow(ph)))
      cat(sprintf("   %-12s %-10s vs %-10s  %s p_adj = %.4f\n", ph$variable[i],
                  ph$group1[i], ph$group2[i], ph$posthoc[i], ph$p_adj[i]))
  } else cat("   no post-hoc contrast reaches 0.05\n")
}

cat("\n=== AMYLOID (5xFAD) ===\n")
for (i in seq_len(nrow(ra$omnibus)))
  cat(sprintf("   %-10s %s p = %.4f\n", ra$omnibus$variable[i],
              ra$omnibus$test[i], ra$omnibus$p[i]))
if (!is.null(ra$posthoc)) {
  sg <- ra$posthoc[ra$posthoc$p_adj < 0.05, , drop = FALSE]
  cat(if (nrow(sg)) "" else "   no amyloid post-hoc contrast reaches 0.05\n")
  for (i in seq_len(nrow(sg)))
    cat(sprintf("   %-10s %s vs %s  p_adj = %.4f\n", sg$variable[i],
                sg$group1[i], sg$group2[i], sg$p_adj[i]))
}
cat("\n[06] ->", STATS, "\n")

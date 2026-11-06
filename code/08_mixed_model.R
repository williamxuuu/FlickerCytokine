# ============================================================================
# 08_mixed_model.R  —  sensitivity analysis: sex modelled rather than removed.
#
# Stage 06 adjusts for sex by z-scoring within sex and pooling. This stage
# addresses the same question on the RAW values with a linear mixed model:
#
#   value ~ Flicker + Sex + (1 | litter)
#
# (lme4 / lmerTest; falls back to lm(value ~ Flicker + Sex) when the litter
# variance is singular or there are fewer than 3 litters -- the `model` column
# records which was used). The Flicker x Sex interaction, which the z-score
# approach cannot assess, comes from a separate Flicker * Sex fit; pairwise
# contrasts are sex-averaged estimated marginal means (emmeans, Tukey). Models
# are fitted on raw MFI and on log1p(MFI), since MFI is right-skewed;
# `p_resid_normal` reports a Shapiro-Wilk test on the residuals.
#
# Results are written alongside the stage-06 p-values for comparison; no figure
# depends on this stage.
#
# Outputs: stats/mixed_model_omnibus.csv, stats/mixed_model_posthoc.csv.
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
  d <- ANIMALS_RAW[ANIMALS_RAW$Genotype == g, , drop = FALSE]
  vars <- c(PANEL, if (g == "5XFAD") AMYLOID_TARGETS)
  for (tr in c("none", "log1p")) for (v in vars) {
    r <- compare_flicker_mixed(d, v, subgroup = g, transform = tr)
    oms[[paste(g, tr, v)]] <- r$omnibus; phs[[paste(g, tr, v)]] <- r$posthoc
  }
}
om <- do.call(rbind, oms); ph <- do.call(rbind, phs)
rownames(om) <- NULL; rownames(ph) <- NULL

# stage-06 omnibus p (z-score approach), for direct comparison
z <- rbind(read.csv(file.path(STATS, "group_omnibus_cytokines.csv"), check.names = FALSE),
           read.csv(file.path(STATS, "group_omnibus_amyloid.csv"), check.names = FALSE))
om$p_zscore_route <- z$p[match(paste(om$subgroup, om$variable),
                               paste(z$subgroup, z$variable))]

write.csv(om, file.path(STATS, "mixed_model_omnibus.csv"), row.names = FALSE)
write.csv(ph, file.path(STATS, "mixed_model_posthoc.csv"), row.names = FALSE)

cat("\n[08] model used (raw scale):\n")
print(table(om$subgroup[om$scale == "none"], om$model[om$scale == "none"]))
for (tr in c("none", "log1p")) {
  o <- om[om$scale == tr, ]
  cat(sprintf("\n[08] scale = %s\n", tr))
  for (g in GENOTYPES) {
    oc <- o[o$subgroup == g & o$variable %in% PANEL, ]
    cat(sprintf("  %-6s flicker p<0.05: %2d/%d (z-score route %d) | sex p<0.05: %2d | FlickerxSex p<0.05: %2d | resid non-normal: %2d\n",
                g, sum(oc$p_flicker < 0.05, na.rm = TRUE), nrow(oc),
                sum(oc$p_zscore_route < 0.05, na.rm = TRUE),
                sum(oc$p_sex < 0.05, na.rm = TRUE), sum(oc$p_flickerXsex < 0.05, na.rm = TRUE),
                sum(oc$p_resid_normal < 0.05, na.rm = TRUE)))
  }
}
cat("[08] -> stats/mixed_model_omnibus.csv, stats/mixed_model_posthoc.csv\n")

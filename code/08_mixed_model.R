# ============================================================================
# 08_mixed_model.R  —  sex-controlled mixed model, a CROSS-CHECK on stage 06.
#
# Stage 06 handles sex by z-scoring within sex, then pooling. This stage asks
# the same question the other way, as the Wood lab suggested ("Multiple effects
# model ... to control for sex", initialComments.md): keep the RAW values and
# put Sex in the model, with litter as a random intercept because littermates
# are not independent. See compare_flicker_mixed() in lib/stats_helpers.R.
#
# It runs on the RAW table. On the z-scored table Sex is zero by construction,
# which is why stage 06 dropped this model; here that is the point.
#
# Two things only this model can report:
#   - Flicker x Sex: does the flicker effect itself differ between the sexes?
#     The z-score cannot see this; pooling averages over it.
#   - litter: whether litter carries variance at all (the `model` column says
#     which fit was used).
#
# Fitted on raw MFI and, as a sensitivity, on log1p(MFI) -- MFI is skewed and
# a linear model assumes roughly normal residuals (p_resid_normal reports it).
#
# NOTHING DOWNSTREAM READS THIS. It is reported beside stage 06, not instead of
# it; no figure changes.
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

# beside it, the stage-06 omnibus p on the z-scored table, for direct comparison
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

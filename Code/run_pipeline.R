# ============================================================================
# run_pipeline.R — top-level driver for the corrected cytokine/morphology
# pipeline. Usage:  Rscript run_pipeline.R
# Working directory must be this folder (cytokine_analysis/).
# ============================================================================

t_start <- Sys.time()
source("config.R")

steps <- c(
  "scripts/01_load_and_prep.R",          # read + background-subtract + factors
  "scripts/02_linear_range_and_batch.R", # prune to linear panel + sex-adjust
  "scripts/03_cytokine_heatmaps.R",      # heatmap3 per genotype
  "scripts/04_plsda.R",                  # PLS-DA + Garza-style permutation, LOOCV loadings
  "scripts/04b_key_cytokines.R",         # top 3 on LV1 + top 3 on LV2, per genotype
  "scripts/05_cytokine_group_stats.R",   # AD -> ANOVA+Tukey / KW+Dunn per analyte
  "scripts/05b_cytokine_mixed_effects.R",# sex-controlled mixed/multiple-effects
  "scripts/06_amyloid.R",                # Abeta bars, corrected
  "scripts/07_correlations.R",           # Spearman, p-driven, row/myeloid FDR
  "scripts/08_morphology_stats.R",       # per-cell violins (re-derivation; see hdr)
  "scripts/09_export_revised_figures.R", # corrected figures as PNG for the manuscript
  "scripts/10_figure5_correlations.R",   # Figure 5 correlations in the original style
  "scripts/11_supp_myeloid_exploratory.R" # exploratory myeloid supplement + detection panel
)

for (s in steps) {
  cat("\n==========================================================\n")
  cat("[run] sourcing:", s, "\n")
  cat("==========================================================\n")
  t0 <- Sys.time()
  tryCatch(
    source(s, echo = FALSE),
    error = function(e) {
      cat("[run] ERROR in", s, ":\n", conditionMessage(e), "\n"); stop(e)
    })
  cat(sprintf("[run] %s done in %.1fs\n", s,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
cat(sprintf("\n[run] pipeline finished in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

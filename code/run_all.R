# ============================================================================
# run_all.R  —  runs the pipeline in dependency order.
#
#   Rscript run_all.R            everything
#   Rscript run_all.R stats      setup + QC + PLS-DA + statistics, no figures
#   Rscript run_all.R figures    figures only (statistics must already exist)
#
# Each stage runs in a separate Rscript process (each sources 00_setup.R), so
# a failure stops the run at that stage and no state leaks between stages.
#
# Stage 02 (PLS-DA permutation + LOOCV) takes ~4 min per genotype on first run
# and is cached thereafter (see 02_plsda.R). Stage 05 takes ~1 min.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()

STATS_STAGES <- c(
  "01_sex_zscore_qc.R",     # QC of the within-sex z-score
  "02_plsda.R",             # PLS-DA, permutation test, key cytokines
  "05_correlations.R",      # cytokine ~ amyloid / morphology
  "06_group_stats.R",       # flicker group comparisons
  "08_mixed_model.R",       # sensitivity: mixed model on raw values
  "09_corr_selection.R")    # correlations under two selection rules

FIGURE_STAGES <- c(
  "10_fig2.R",              # Figure 2: WT heatmap and PLS-DA
  "11_fig4.R",              # Figure 4: amyloid, cytokine ~ amyloid
  "12_fig5.R",              # Figure 5: amyloid / cytokines ~ glia
  "13_supp4.R",             # Supp 4: 5xFAD heatmap and PLS-DA
  "14_supp5.R",             # Supp 5: key cytokine group comparisons
  "15_supp6.R",             # Supp 6: most responsive cytokines
  "16_supp7.R",             # Supp 7: WT cytokines ~ glia
  "17_data_folders.R")      # results filed by topic (runs last)

what <- commandArgs(trailingOnly = TRUE)
stages <- if (!length(what)) {
  c(STATS_STAGES, FIGURE_STAGES)
} else if (identical(what[1], "stats")) {
  STATS_STAGES
} else if (identical(what[1], "figures")) {
  FIGURE_STAGES
} else {
  stop("usage: Rscript run_all.R [stats|figures]")
}

rscript <- file.path(R.home("bin"), "Rscript")
t0 <- Sys.time()
for (s in stages) {
  cat("\n", strrep("=", 72), "\n", s, "\n", strrep("=", 72), "\n", sep = "")
  st <- system2(rscript, file.path(HERE, s))
  if (st != 0) stop("stage failed: ", s, " (exit ", st, ")")
}
cat(sprintf("\n[run_all] %d stages in %.1f min\n",
            length(stages), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

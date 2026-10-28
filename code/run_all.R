# ============================================================================
# run_all.R  —  the whole cascade, in dependency order.
#
#   Rscript run_all.R            everything
#   Rscript run_all.R stats      setup + QC + PLS-DA + statistics, no figures
#   Rscript run_all.R figures    figures only (statistics must already exist)
#
# Each stage is a separate Rscript process, so a failure stops the run at the
# stage that failed and cannot leave a half-updated global environment behind
# for the next one.
#
# Stage 02 (PLS-DA) is the slow one, ~4 min per genotype for the permutation
# and LOOCV. It caches on a fingerprint of its actual input, so re-runs are
# instant unless the panel, the z-score or the data changed. Stage 05 is ~1 min
# for 50,000 permutations per family.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()

STATS_STAGES <- c(
  "01_sex_zscore_qc.R",     # is the z-score safe here? balance + before/after
  "02_plsda.R",             # THE HUB: model -> key cytokines
  "05_correlations.R",      # cytokine ~ amyloid / morphology
  "06_group_stats.R",       # flicker group comparisons
  "08_mixed_model.R")       # cross-check: raw values, Sex + (1|litter) in the model

FIGURE_STAGES <- c(
  "10_fig2.R",              # WT heatmap, PLS-DA scores + violins, loadings
  "11_fig4.R",              # amyloid bars, cytokine ~ amyloid heatmap
  "12_fig5.R",              # amyloid ~ morphology, key cytokines ~ morphology
  "13_supp4.R",             # 5xFAD version of Figure 2
  "14_supp5.R",             # key cytokine group panels
  "15_supp6.R",             # cytokine profiling, both genotypes
  "16_supp7.R",             # WT key cytokines ~ morphology
  "17_data_folders.R")      # re-file results by topic into data/ (runs last)

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

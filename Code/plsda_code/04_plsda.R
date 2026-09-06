# ============================================================================
# 04_plsda.R
# Supervised PLS-DA (ropls) of cytokines by flicker, per genotype (Fig 2/Supp 4).
# [CORRECTION] Three changes from the original:
#   1. "cannot do stats on a linked variable" -- the original ran t-test/ANOVA on
#      the LV1/LV2 SCORES. That is removed. Group significance comes from the
#      per-analyte tests (stage 05); PLS-DA is a weighting, not a group test.
#   2. "Permutation analysis as in Kristie's studies" -- validation now follows
#      Garza et al. 2020 (J Neurosci 40:1211) / Sci Adv 2023 eadf5672: LOOCV
#      classification accuracy, compared against the accuracies of models fit to
#      randomly reassigned groups (p_permute). ropls' own R2Y/Q2 permutation is
#      kept as a secondary metric.
#   3. Loadings are now exported (results/stats/plsda_loadings.csv) as mean +/- SD
#      over 1000 LOOCV iterations, so "the top 3 analytes on LV1" is reproducible
#      rather than read off a bar chart. Stage 04b does that selection.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(ropls) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/WoodLabFunctions.R")   # scores.plot, loadings.chart2
source("lib/plsda_helpers.R")

processed_dir     <- file.path(results_dir, "processed")
adj               <- readRDS(file.path(processed_dir, "animal_table_sexadj.rds"))
analysis_analytes <- readRDS(file.path(processed_dir, "analysis_analytes.rds"))
set.seed(rng_seed)

perm_rows <- list(); load_rows <- list(); fits <- list()

for (g in genotypes) {
  df <- adj[adj$Genotype == g, , drop = FALSE]
  X  <- plsda_matrix(df, analysis_analytes)
  y  <- factor(df$Flicker, levels = flicker_levels)
  pal <- if (g == "WT") flicker_colors_wt else flicker_colors_fad
  cat(sprintf("[04] %-6s PLS-DA on %d animals x %d analytes\n", g, nrow(X), ncol(X)))

  # --- model + secondary ropls permutation (R2Y / Q2) -----------------------
  mod <- tryCatch(plsda_fit(X, y, permI = plsda_ropls_permI),
                  error = function(e) { cat("[04]", g, "opls error:", conditionMessage(e), "\n"); NULL })
  if (is.null(mod)) next
  rot <- plsda_orient(mod, label = g)

  # --- LOOCV loadings (mean +/- SD over 1000 iterations) --------------------
  cat(sprintf("[04] %-6s LOOCV loadings, %d iterations ...\n", g, plsda_loocv_runs))
  lo <- plsda_loocv_loadings(X, y, runs = plsda_loocv_runs, ref = mod)
  # apply the same axis orientation to the LOOCV means
  lo_mean <- sweep(lo$mean, 2, rot$sign, "*")

  load_rows[[g]] <- data.frame(
    genotype = g, analyte = colnames(X),
    LV1 = rot$P1, LV2 = rot$P2,
    LV1_loocv_mean = lo_mean[, "LV1"], LV1_loocv_sd = lo$sd[, "LV1"],
    LV2_loocv_mean = lo_mean[, "LV2"], LV2_loocv_sd = lo$sd[, "LV2"],
    loocv_runs = lo$runs, row.names = NULL, stringsAsFactors = FALSE)

  # --- permutation test on the model (the "Kristie" statistic) --------------
  cat(sprintf("[04] %-6s label permutation, %d iterations ...\n", g, plsda_n_perm))
  pm  <- plsda_permutation(X, y, n_perm = plsda_n_perm)
  sdf <- mod@summaryDF
  perm_rows[[g]] <- data.frame(
    subgroup = g,
    R2X = sdf$`R2X(cum)`, R2Y = sdf$`R2Y(cum)`, Q2 = sdf$`Q2(cum)`,
    loocv_accuracy = pm$accuracy, chance = 1 / nlevels(y),
    perm_null_mean = pm$null_mean, perm_null_sd = pm$null_sd,
    p_permute = pm$p_permute, n_perm_loocv = pm$n_perm,
    p_R2Y_ropls = sdf$pR2Y, p_Q2_ropls = sdf$pQ2, n_perm_ropls = plsda_ropls_permI,
    n_animals = nrow(X), n_analytes = ncol(X),
    note = "p_permute = LOOCV accuracy vs randomly reassigned groups; NO stats on LV scores",
    stringsAsFactors = FALSE)

  # --- figures (visualisation only) -----------------------------------------
  p_sc <- scores.plot(rot$T1, rot$T2, as.character(y), "Flicker",
                      analysis = "PLSDA", plot.title = paste0(g, " PLS-DA")) +
    ggplot2::scale_color_manual(values = pal)
  ggplot2::ggsave(file.path(cyt_dir, sprintf("plsda_scores_%s.pdf", g)),
                  p_sc, width = 6, height = 5)

  for (lv in 1:2) {
    L  <- if (lv == 1) rot$P1 else rot$P2
    SD <- lo$sd[, lv]
    ggplot2::ggsave(
      file.path(cyt_dir, sprintf("plsda_loadings_%s_LV%d.pdf", g, lv)),
      loadings_chart_sd_gg(L, SD, sprintf("%s LV%d", g, lv)),
      width = 8, height = 5)
  }

  fits[[g]] <- list(X = X, y = y, mod = mod, rot = rot, loocv = lo, perm = pm)
}

perm_tbl <- do.call(rbind, perm_rows)
load_tbl <- do.call(rbind, load_rows)
if (!is.null(perm_tbl)) {
  write.csv(perm_tbl, file.path(stats_dir, "plsda_permutation.csv"), row.names = FALSE)
  cat("[04] model validity (Garza-style permutation):\n")
  print(perm_tbl[, c("subgroup", "loocv_accuracy", "chance", "perm_null_mean",
                     "p_permute", "R2Y", "Q2", "p_R2Y_ropls")], row.names = FALSE)
}
if (!is.null(load_tbl))
  write.csv(load_tbl, file.path(stats_dir, "plsda_loadings.csv"), row.names = FALSE)
saveRDS(fits, file.path(processed_dir, "plsda_fits.rds"))
cat("[04] PLS-DA done ->", cyt_dir, "\n")

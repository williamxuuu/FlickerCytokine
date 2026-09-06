# ============================================================================
# lib/plsda_helpers.R
# PLS-DA support shared by stages 04, 04b and 09, so the model, its orientation
# and its loadings are identical everywhere they appear.
#
# [CORRECTION] "Permutation analysis as in Kristie's studies" — the validation
# scheme of Garza et al. 2020 (J Neurosci 40:1211) and Sci Adv 2023 eadf5672:
#   - orthogonal rotation in the LV1-LV2 plane;
#   - leave-one-out cross-validation, loadings reported as mean +/- SD across
#     iterations (that is what the error bars on a Wood-lab loadings chart are);
#   - animals randomly reassigned to groups and the PLS-DA refit many times;
#     p_permute is the fraction of permuted models that match or beat the true
#     model, and a model is only called valid at p_permute < 0.05.
# No statistic is ever computed on the LV scores themselves ("cannot do stats on
# a linked variable") — group differences come from the per-analyte tests.
# ============================================================================

suppressPackageStartupMessages({ library(ropls); library(matrixStats) })

# z-score and drop constant columns; the matrix the model actually sees.
plsda_matrix <- function(df, analytes) {
  X <- as.matrix(df[, analytes, drop = FALSE])
  X <- apply(X, 2, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s
  })
  X[, apply(X, 2, function(z) sd(z) > 0), drop = FALSE]
}

plsda_fit <- function(X, y, permI = 0) {
  suppressWarnings(opls(X, y, predI = plsda_n_lv, permI = permI,
                        scaleC = "standard", fig.pdfC = "none", info.txtC = "none"))
}

# --- deterministic axis orientation -----------------------------------------
# The sign of a latent variable is arbitrary: flipping an LV flips its scores and
# its loadings together and changes nothing about the model. Fix it so repeated
# runs and the different stages draw the same picture: on each LV, the analyte
# with the largest absolute loading is made positive. `plsda_flip` in config.R
# can override this per genotype if a published orientation must be matched.
plsda_orient <- function(mod, label = "") {
  L <- mod@loadingMN
  flip <- vapply(1:2, function(k) L[which.max(abs(L[, k])), k] < 0, logical(1))
  ov <- if (exists("plsda_flip") && !is.null(plsda_flip[[label]]))
          plsda_flip[[label]] else NULL
  if (!is.null(ov)) flip <- as.logical(ov)
  sgn <- ifelse(flip, -1, 1)
  list(T1 = mod@scoreMN[, 1] * sgn[1], T2 = mod@scoreMN[, 2] * sgn[2],
       P1 = L[, 1] * sgn[1],           P2 = L[, 2] * sgn[2],
       sign = sgn)
}

# --- LOOCV loadings: mean +/- SD --------------------------------------------
# Garza et al.: "the mean and SDs after iteratively excluding single samples (a
# leave-one-out cross-validation), one at a time, and recalculating the PLSDA
# 1000 times". Each iteration's loading vector is sign-aligned to the full model
# before averaging, because the LV sign is arbitrary per fit.
plsda_loocv_loadings <- function(X, y, runs = plsda_loocv_runs, ref = NULL) {
  full <- if (is.null(ref)) plsda_fit(X, y) else ref
  Lref <- full@loadingMN
  P <- array(NA_real_, dim = c(runs, ncol(X), 2),
             dimnames = list(NULL, colnames(X), c("LV1", "LV2")))
  n <- nrow(X)
  for (i in seq_len(runs)) {
    out <- sample.int(n, 1)
    Xi  <- X[-out, , drop = FALSE]
    ok  <- colSds(Xi) > 1e-10
    mi  <- tryCatch(plsda_fit(Xi[, ok, drop = FALSE], y[-out]), error = function(e) NULL)
    if (is.null(mi)) next
    for (k in 1:2) {
      v <- rep(0, ncol(X)); v[ok] <- mi@loadingMN[, k]
      P[i, , k] <- v * sign(sum(v * Lref[, k]))     # align to the full model
    }
  }
  list(mean = cbind(LV1 = colMeans(P[, , 1], na.rm = TRUE),
                    LV2 = colMeans(P[, , 2], na.rm = TRUE)),
       sd   = cbind(LV1 = colSds(P[, , 1], na.rm = TRUE),
                    LV2 = colSds(P[, , 2], na.rm = TRUE)),
       runs = runs)
}

# --- LOOCV classification accuracy ------------------------------------------
# Systematic leave-one-out (every animal held out exactly once), so the statistic
# is deterministic and the observed model and every permuted model are scored the
# same way.
plsda_loocv_accuracy <- function(X, y) {
  n <- nrow(X)
  hit <- vapply(seq_len(n), function(i) {
    Xi <- X[-i, , drop = FALSE]
    ok <- colSds(Xi) > 1e-10
    mi <- tryCatch(plsda_fit(Xi[, ok, drop = FALSE], y[-i]), error = function(e) NULL)
    if (is.null(mi)) return(NA)
    pr <- tryCatch(as.character(predict(mi, newdata = X[i, ok, drop = FALSE])),
                   error = function(e) NA_character_)
    if (is.na(pr)) NA else pr == as.character(y[i])
  }, logical(1))
  mean(hit, na.rm = TRUE)
}

# --- permutation test on the model ------------------------------------------
# Animals are randomly reassigned to flicker groups and the whole PLS-DA
# (including its LOOCV scoring) is rerun. p_permute = P(permuted accuracy >=
# observed accuracy), the fraction of shuffled models that do as well or better.
plsda_permutation <- function(X, y, n_perm = plsda_n_perm) {
  obs  <- plsda_loocv_accuracy(X, y)
  null <- numeric(n_perm)
  for (i in seq_len(n_perm)) null[i] <- plsda_loocv_accuracy(X, sample(y))
  list(accuracy   = obs,
       null_mean  = mean(null, na.rm = TRUE),
       null_sd    = sd(null, na.rm = TRUE),
       # +1 in numerator and denominator: a permutation p can never be exactly 0
       p_permute  = (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1),
       n_perm     = n_perm,
       null       = null)
}

# --- loadings chart with LOOCV error bars ------------------------------------
# ggplot twin of the Wood lab's loadings.chart.sd (lib/WoodLabFunctions.R:208),
# styled like loadings.chart2 so it drops into the manuscript figures unchanged.
# Bars are the full-model loading scaled by the largest absolute loading; error
# bars are the LOOCV SD on the same scale.
# `highlight` names the analytes to draw in the darker fill (the top-N picks for
# that LV), so the selection stage 04b made is visible on the chart itself.
loadings_chart_sd_gg <- function(loading, sdev, component_str, x_size = 10,
                                 fill = "lightblue", highlight = character(0),
                                 highlight_fill = "#2C7FB8") {
  suppressPackageStartupMessages(library(ggplot2))
  s   <- max(abs(loading))
  ord <- order(loading)
  d <- data.frame(Analyte  = factor(names(loading)[ord], levels = names(loading)[ord]),
                  Loading  = loading[ord] / s,
                  SD       = if (is.null(sdev)) 0 else sdev[ord] / s)
  d$Pick <- ifelse(as.character(d$Analyte) %in% highlight, "top", "other")
  ggplot(d, aes(x = Analyte, y = Loading, fill = Pick)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(top = highlight_fill, other = fill), guide = "none") +
    geom_errorbar(aes(ymin = Loading - SD, ymax = Loading + SD),
                  width = 0.25, linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(x = NULL, y = NULL, title = paste0("Signals in ", component_str)) +
    theme_classic(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = x_size, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(color = "black"))
}

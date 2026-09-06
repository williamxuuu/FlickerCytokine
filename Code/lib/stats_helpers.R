# ============================================================================
# lib/stats_helpers.R
# The corrected statistics, factored out so every stage uses the same rules.
# These functions ARE the Wood-lab corrections from ../initialComments.md:
#   - group comparisons: Anderson-Darling -> ANOVA+Tukey or Kruskal-Wallis+Dunn
#   - correlations: Spearman for the stat, Pearson line only, p-value driven,
#     FDR by row (or restricted to the myeloid panel)
#   - scatterplots: males/females labeled
# Every function returns TIDY tables (one row per test) so results/stats/ can be
# assembled into the manuscript's Supplementary Statistics table directly.
# ============================================================================

suppressPackageStartupMessages({
  library(nortest)   # Anderson-Darling
  library(rstatix)   # dunn_test, tukey_hsd
  library(ggplot2)
})

# --- normality ---------------------------------------------------------------
# [CORRECTION] manuscript pre-registers Anderson-Darling to choose parametric vs
# not. Tested on within-group-centered residuals (the ANOVA/KW assumption).
# Small groups (n<8) fall back to Shapiro; n<3 -> treat as non-normal.
ad_is_normal <- function(values, groups, alpha = normality_alpha) {
  ok <- !is.na(values) & !is.na(groups)
  values <- values[ok]; groups <- as.factor(groups[ok])
  if (length(values) < 3 || length(unique(values)) < 3) return(FALSE)
  resid <- values - ave(values, groups, FUN = function(z) median(z))
  p <- tryCatch({
    if (length(resid) >= 8) nortest::ad.test(resid)$p.value
    else shapiro.test(resid)$p.value
  }, error = function(e) 0)
  isTRUE(p >= alpha)
}

# --- one variable across the 3 flicker groups --------------------------------
# Returns list(omnibus = <1-row tibble>, posthoc = <tibble of pairwise>).
# Tukey p is already family-wise adjusted; Dunn's is adjusted with
# posthoc_padjust (config). No raw pairwise p is ever reported as final.
compare_flicker <- function(df, value_col, subgroup = "", variable = value_col,
                            group_col = "Flicker") {
  d <- data.frame(y = df[[value_col]], g = factor(df[[group_col]],
                                                  levels = flicker_levels))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d <- droplevels(d)
  ns <- table(d$g)
  base <- data.frame(subgroup = subgroup, variable = variable,
                     n_total = nrow(d),
                     n_groups = paste(names(ns), ns, sep = "=", collapse = ","),
                     stringsAsFactors = FALSE)

  # not enough data to test
  if (nlevels(d$g) < 2 || any(ns < 2) || nrow(d) < 4) {
    om <- cbind(base, test = "none", statistic = NA, df = NA,
                p = NA, method = "insufficient-n")
    return(list(omnibus = om, posthoc = NULL))
  }

  normal <- ad_is_normal(d$y, d$g)
  if (normal) {
    fit <- stats::aov(y ~ g, data = d)
    s   <- summary(fit)[[1]]
    om  <- cbind(base, test = "ANOVA",
                 statistic = s[["F value"]][1], df = s[["Df"]][1],
                 p = s[["Pr(>F)"]][1], method = "AD-normal -> ANOVA+Tukey")
    ph  <- tryCatch(rstatix::tukey_hsd(d, y ~ g), error = function(e) NULL)
    if (!is.null(ph)) {
      ph <- data.frame(subgroup = subgroup, variable = variable,
                       group1 = ph$group1, group2 = ph$group2,
                       estimate = ph$estimate, p_adj = ph$p.adj,
                       posthoc = "TukeyHSD", stringsAsFactors = FALSE)
    }
  } else {
    kw  <- stats::kruskal.test(y ~ g, data = d)
    om  <- cbind(base, test = "Kruskal-Wallis",
                 statistic = unname(kw$statistic), df = unname(kw$parameter),
                 p = kw$p.value, method = "AD-nonnormal -> KW+Dunn")
    ph  <- tryCatch(
      rstatix::dunn_test(d, y ~ g, p.adjust.method = posthoc_padjust),
      error = function(e) NULL)
    if (!is.null(ph)) {
      ph <- data.frame(subgroup = subgroup, variable = variable,
                       group1 = ph$group1, group2 = ph$group2,
                       estimate = ph$statistic, p_adj = ph$p.adj,
                       posthoc = paste0("Dunn(", posthoc_padjust, ")"),
                       stringsAsFactors = FALSE)
    }
  }
  list(omnibus = om, posthoc = ph)
}

# Run compare_flicker over many variables; bind into two tidy tables.
compare_flicker_many <- function(df, value_cols, subgroup = "") {
  oms <- list(); phs <- list()
  for (v in value_cols) {
    if (!v %in% colnames(df)) next
    r <- compare_flicker(df, v, subgroup = subgroup, variable = v)
    oms[[v]] <- r$omnibus; phs[[v]] <- r$posthoc
  }
  list(omnibus = do.call(rbind, oms),
       posthoc = if (length(phs)) do.call(rbind, phs) else NULL)
}

# --- sex-controlled mixed / multiple-effects model ---------------------------
# [CORRECTION, initialComments.md line 14] "Multiple effects model for bar plots
# ... to control for sex". This tests the flicker effect while controlling for
# sex, and reports whether the flicker effect itself differs by sex.
#
# Model choice (kept honest for the small N):
#   primary  : linear mixed model  y ~ Flicker + Sex + (1 | litter)
#              litter is the biologically-standard clustering variable here
#              (5-6 litters, sizes 2-9, deliberately split across flicker), so a
#              random intercept accounts for littermate correlation while Sex is
#              a fixed covariate. Flicker tested by Type-III ANOVA with
#              Satterthwaite df (lmerTest); flicker pairwise by estimated
#              marginal means, Tukey-adjusted, averaged over sex (emmeans).
#   fallback : if the litter variance is ~0 (singular fit -- common with only
#              5-6 groups) or lmer fails, drop the random term and use the
#              fixed-effects model y ~ Flicker + Sex (Type-II ANOVA, which is
#              contrast-invariant). The model actually used is reported per
#              analyte so nothing is hidden.
# Cohort is NOT used as a random effect: it has 2 badly-imbalanced levels here,
# which would make the random term degenerate.
compare_flicker_mixed <- function(df, value_col, subgroup = "",
                                  variable = value_col, group = "litter") {
  suppressPackageStartupMessages({
    library(lmerTest); library(emmeans); library(car)
  })
  d <- data.frame(y = df[[value_col]],
                  Flicker = factor(df$Flicker, levels = flicker_levels),
                  Sex = factor(df$Sex, levels = sexes),
                  grp = factor(df[[group]]))
  d <- droplevels(d[stats::complete.cases(d), , drop = FALSE])
  base <- data.frame(subgroup = subgroup, variable = variable,
                     n_total = nrow(d), n_litter = nlevels(d$grp),
                     stringsAsFactors = FALSE)
  na_out <- function(model) list(
    omnibus = cbind(base, model = model, F_flicker = NA, df1 = NA, df2 = NA,
                    p_flicker = NA, p_sex = NA, p_flickerXsex = NA),
    posthoc = NULL)

  # need both sexes and all three flicker groups populated
  if (nlevels(d$Flicker) < 2 || nlevels(d$Sex) < 2 ||
      nrow(d) < 6 || any(table(d$Flicker) < 2))
    return(na_out("insufficient-n"))

  # --- fit primary mixed model, detect singular litter variance -------------
  ctrl <- lme4::lmerControl(check.conv.singular = "ignore")
  m <- tryCatch(lmerTest::lmer(y ~ Flicker + Sex + (1 | grp), data = d,
                               control = ctrl, REML = TRUE),
                error = function(e) NULL)
  use_lmer <- !is.null(m) && nlevels(d$grp) >= 3 && !lme4::isSingular(m)

  if (use_lmer) {
    model <- "lmer: Flicker+Sex+(1|litter)"
    an  <- anova(m)                                  # Type III, Satterthwaite
    Ff  <- an["Flicker", "F value"]; pf <- an["Flicker", "Pr(>F)"]
    d1  <- an["Flicker", "NumDF"];   d2 <- an["Flicker", "DenDF"]
    ps  <- an["Sex", "Pr(>F)"]
  } else {
    m <- stats::lm(y ~ Flicker + Sex, data = d)      # fixed-effects fallback
    model <- if (nlevels(d$grp) < 3) "lm: Flicker+Sex (few litters)"
             else "lm: Flicker+Sex (litter var~0)"
    an  <- car::Anova(m, type = 2)                   # contrast-invariant
    Ff  <- an["Flicker", "F value"]; pf <- an["Flicker", "Pr(>F)"]
    d1  <- an["Flicker", "Df"];      d2 <- m$df.residual
    ps  <- an["Sex", "Pr(>F)"]
  }

  # --- flicker effect differ by sex? (interaction diagnostic) ---------------
  pfs <- tryCatch({
    if (use_lmer) {
      mi <- lmerTest::lmer(y ~ Flicker * Sex + (1 | grp), data = d, control = ctrl)
      anova(mi)["Flicker:Sex", "Pr(>F)"]
    } else {
      car::Anova(stats::lm(y ~ Flicker * Sex, data = d), type = 2)["Flicker:Sex", "Pr(>F)"]
    }
  }, error = function(e) NA_real_)

  # --- flicker pairwise, sex-averaged marginal means, Tukey-adjusted --------
  ph <- tryCatch({
    em <- emmeans::emmeans(m, ~ Flicker)
    pr <- as.data.frame(emmeans::contrast(em, method = "pairwise", adjust = "tukey"))
    parts <- do.call(rbind, strsplit(as.character(pr$contrast), " - ", fixed = TRUE))
    data.frame(subgroup = subgroup, variable = variable,
               group1 = trimws(parts[, 1]), group2 = trimws(parts[, 2]),
               estimate = pr$estimate, p_adj = pr$p.value,
               posthoc = "emmeans-Tukey (sex-averaged)", stringsAsFactors = FALSE)
  }, error = function(e) NULL)

  list(omnibus = cbind(base, model = model, F_flicker = Ff, df1 = d1, df2 = d2,
                       p_flicker = pf, p_sex = ps, p_flickerXsex = pfs),
       posthoc = ph)
}

compare_flicker_mixed_many <- function(df, value_cols, subgroup = "") {
  oms <- list(); phs <- list()
  for (v in value_cols) {
    if (!v %in% colnames(df)) next
    r <- compare_flicker_mixed(df, v, subgroup = subgroup, variable = v)
    oms[[v]] <- r$omnibus; phs[[v]] <- r$posthoc
  }
  list(omnibus = do.call(rbind, oms),
       posthoc = if (length(phs)) do.call(rbind, phs) else NULL)
}

# --- correlations ------------------------------------------------------------
# [CORRECTION] Spearman rho + p is the statistic; the p-value drives inclusion
# (no |rho| cutoff). [REVISED SPEC] "FDR correct based on pre-selected cytokines
# based on the heatmap - not all 32 possible cytokines": the family is the
# retained heatmap panel, corrected BY ROW (per target). "myeloid" and "key"
# restrict the family further and exist for sensitivity tables only. FDR is never
# applied blanket across the whole matrix.
corr_one <- function(x, y, method = corr_stat_method) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(c(rho = NA, p = NA, n = length(x)))
  ct <- suppressWarnings(stats::cor.test(x, y, method = method))
  c(rho = unname(ct$estimate), p = unname(ct$p.value), n = length(x))
}

# analytes: character vec of columns; targets: character vec of columns.
# fdr_scope: "by_row" (the retained heatmap panel, per target), "key" (only the
# PLS-DA-selected key cytokines, passed in via `restrict`), or "myeloid".
correlate <- function(df, analytes, targets, subgroup = "",
                      method = corr_stat_method,
                      fdr_scope = corr_fdr_scope, restrict = NULL) {
  if (identical(fdr_scope, "myeloid"))
    analytes <- intersect(analytes, cytokines_myeloid)
  if (identical(fdr_scope, "key") && !is.null(restrict))
    analytes <- intersect(analytes, restrict)
  rows <- list()
  for (tg in targets) {
    for (an in analytes) {
      if (!all(c(tg, an) %in% colnames(df))) next
      r <- corr_one(df[[an]], df[[tg]], method = method)
      rows[[paste(tg, an)]] <- data.frame(
        subgroup = subgroup, target = tg, analyte = an,
        rho = r["rho"], p_raw = r["p"], n = r["n"],
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(out)
  # FDR by row = within each target, across the analytes tested
  out$p_adj <- NA_real_
  for (tg in unique(out$target)) {
    i <- out$target == tg & !is.na(out$p_raw)
    out$p_adj[i] <- stats::p.adjust(out$p_raw[i], method = "BH")
  }
  out$fdr_scope <- fdr_scope
  rownames(out) <- NULL
  out
}

# --- corrected scatter (Spearman stat, Pearson line, M/F labeled) ------------
# [CORRECTION] statistic annotated is Spearman (rho + p); the line is a Pearson
# linear trend for display only; sex is shown by point shape + color. `q` is the
# BH-adjusted p within the analyte's target row, annotated when supplied so the
# corrected value travels with the figure.
scatter_sex <- function(df, xcol, ycol, subgroup = "",
                        colors = flicker_colors_wt, q = NA_real_) {
  d <- data.frame(x = df[[xcol]], y = df[[ycol]],
                  Sex = factor(df[["Sex"]], levels = sexes))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  st <- corr_one(d$x, d$y, method = "spearman")
  lab <- sprintf("Spearman rho = %s, p = %s%s", fmt(st["rho"]), fmt(st["p"]),
                 if (is.na(q)) "" else sprintf(", q = %s", fmt(q)))
  ggplot(d, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30",
                linetype = "longdash", linewidth = 0.8) +  # Pearson line (display)
    geom_point(aes(shape = Sex, color = Sex), size = 3) +
    scale_shape_manual(values = sex_shapes) +
    scale_color_manual(values = c(M = "#1f78b4", F = "#e31a1c")) +
    labs(x = xcol, y = ycol,
         title = paste(subgroup, ycol, "vs", xcol),
         subtitle = lab) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9))
}

# --- tidy stats sink ---------------------------------------------------------
# Append a tidy stats table to results/stats/<name>.csv (one row per test).
write_stats <- function(tbl, name) {
  if (is.null(tbl) || !nrow(tbl)) return(invisible(NULL))
  path <- file.path(stats_dir, paste0(name, ".csv"))
  utils::write.csv(tbl, path, row.names = FALSE)
  cat(sprintf("[stats] wrote %-32s (%d rows)\n", basename(path), nrow(tbl)))
}

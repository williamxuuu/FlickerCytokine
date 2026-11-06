# ============================================================================
# lib/stats_helpers.R
# Flicker group comparisons, shared so every stage applies the same test.
# Sourced by stages 06, 08, 11, 14 and 15.
#
# Decision rule: Anderson-Darling test for normality, then one-way ANOVA with
# Tukey HSD when normal, or Kruskal-Wallis with Dunn's test (Benjamini-Hochberg
# adjusted) when not. The post hoc procedure controls error across the three
# pairwise flicker contrasts; no further correction is applied across
# variables.
#
# All functions return tidy tables (one row per test).
#
# Globals expected from 00_setup.R: flicker_levels, normality_alpha,
# posthoc_padjust.
# ============================================================================

suppressPackageStartupMessages({
  library(nortest)   # Anderson-Darling
  library(rstatix)   # dunn_test, tukey_hsd
})

# --- normality ---------------------------------------------------------------
# Anderson-Darling on within-group-centred residuals (the ANOVA assumption).
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
# Returns list(omnibus = <1-row df>, posthoc = <df of pairwise>).
# Tukey p is already family-wise adjusted; Dunn's is adjusted with
# posthoc_padjust. Only adjusted pairwise p-values are reported.
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

# --- sex-adjusted mixed model (stage 08 sensitivity analysis) ----------------
# Per variable, within genotype, on the RAW table (on the z-scored table the
# sex effect is zero by construction):
#     y ~ Flicker + Sex + (1 | litter)
# falling back to lm(y ~ Flicker + Sex) when the litter variance is singular or
# there are < 3 litters. Reports the sex-controlled flicker effect, the sex main
# effect, the Flicker x Sex interaction (from a separate y ~ Flicker * Sex fit),
# a Shapiro test on the residuals, and sex-averaged emmeans-Tukey contrasts.
# `transform` = "log1p" refits on log(1 + y), for skewed MFI.
compare_flicker_mixed <- function(df, value_col, subgroup = "", variable = value_col,
                                  group = "litter", transform = c("none", "log1p")) {
  transform <- match.arg(transform)
  suppressPackageStartupMessages({ library(lmerTest); library(emmeans); library(car) })
  y <- df[[value_col]]
  if (transform == "log1p") y <- log1p(pmax(y, 0))
  d <- data.frame(y = y, Flicker = factor(df$Flicker, levels = flicker_levels),
                  Sex = factor(df$Sex, levels = c("M", "F")), grp = factor(df[[group]]))
  d <- droplevels(d[stats::complete.cases(d), , drop = FALSE])
  base <- data.frame(subgroup = subgroup, variable = variable, scale = transform,
                     n_total = nrow(d), n_litter = nlevels(d$grp), stringsAsFactors = FALSE)
  if (nlevels(d$Flicker) < 3 || nlevels(d$Sex) < 2 || any(table(d$Flicker) < 2) ||
      stats::sd(d$y) == 0)
    return(list(omnibus = cbind(base, model = "insufficient", F_flicker = NA, df1 = NA,
                                df2 = NA, p_flicker = NA, p_sex = NA, p_flickerXsex = NA,
                                p_resid_normal = NA), posthoc = NULL))

  ctrl <- lme4::lmerControl(check.conv.singular = "ignore")
  m <- tryCatch(lmerTest::lmer(y ~ Flicker + Sex + (1 | grp), data = d, control = ctrl),
                error = function(e) NULL)
  use_lmer <- !is.null(m) && nlevels(d$grp) >= 3 && !lme4::isSingular(m)

  if (use_lmer) {
    model <- "lmer: Flicker+Sex+(1|litter)"
    an <- anova(m)                                   # Satterthwaite
    Ff <- an["Flicker", "F value"]; pf <- an["Flicker", "Pr(>F)"]
    d1 <- an["Flicker", "NumDF"];   d2 <- an["Flicker", "DenDF"]
    ps <- an["Sex", "Pr(>F)"]
    pfs <- tryCatch(anova(lmerTest::lmer(y ~ Flicker * Sex + (1 | grp), data = d,
                                         control = ctrl))["Flicker:Sex", "Pr(>F)"],
                    error = function(e) NA_real_)
  } else {
    m <- stats::lm(y ~ Flicker + Sex, data = d)
    model <- if (nlevels(d$grp) < 3) "lm: Flicker+Sex (few litters)"
             else "lm: Flicker+Sex (litter var ~ 0)"
    an <- car::Anova(m, type = 2)
    Ff <- an["Flicker", "F value"]; pf <- an["Flicker", "Pr(>F)"]
    d1 <- an["Flicker", "Df"];      d2 <- m$df.residual
    ps <- an["Sex", "Pr(>F)"]
    pfs <- tryCatch(car::Anova(stats::lm(y ~ Flicker * Sex, data = d),
                               type = 2)["Flicker:Sex", "Pr(>F)"],
                    error = function(e) NA_real_)
  }
  pn <- tryCatch(stats::shapiro.test(stats::residuals(m))$p.value, error = function(e) NA_real_)

  ph <- tryCatch({
    pr <- as.data.frame(emmeans::contrast(emmeans::emmeans(m, ~ Flicker),
                                          method = "pairwise", adjust = "tukey"))
    parts <- do.call(rbind, strsplit(as.character(pr$contrast), " - ", fixed = TRUE))
    data.frame(subgroup = subgroup, variable = variable, scale = transform,
               group1 = trimws(parts[, 1]), group2 = trimws(parts[, 2]),
               estimate = pr$estimate, p_adj = pr$p.value,
               posthoc = "emmeans-Tukey (sex-averaged)", stringsAsFactors = FALSE)
  }, error = function(e) NULL)

  list(omnibus = cbind(base, model = model, F_flicker = Ff, df1 = d1, df2 = d2,
                       p_flicker = pf, p_sex = ps, p_flickerXsex = pfs,
                       p_resid_normal = pn),
       posthoc = ph)
}

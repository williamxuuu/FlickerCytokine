# ============================================================================
# 02_linear_range_and_batch.R
# Revised-spec order:
#   (a) batch-adjust cytokines for sex within genotype, so pooled-sex figures
#       aren't driven by a sex shift;
#   (b) remove analytes on the linear-range criteria — "suspicious ones that
#       show no change".
# (a) is done on all 32 analytes and (b) is a column subset, so the two steps
# commute: sex adjustment and z-scoring are both per-analyte, and pruning first
# or last gives identical numbers. The spec's order is used here so the script
# reads the way the analysis is described.
#
# [REVISED SPEC] The prune is now a real removal (config: keep_myeloid_regardless
# is FALSE). Two independent lines of evidence are re-derived and written to
# results/stats/analyte_linear_range.csv so the call is documented, not asserted:
#   1. loading-curve fit — MFI vs amount of lysate loaded, from the BCA linear-
#      range plate (this is the data original/LR Code.R plotted but never scored);
#   2. floor / "shows no change" screen on the 43 animals — an analyte sitting on
#      the assay floor has few distinct values (quantised readings), a narrow
#      range, and/or many animals clipped to zero.
# The decision itself stays anchored on the Wood lab's curated in-linear-range
# panel (config: cytokines_linear), which matches their own "Only Linear
# Cytokines" folder; the no-change screen can remove further analytes from it.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(readxl) })
if (!exists("combined_xlsx")) source("config.R")

processed_dir <- file.path(results_dir, "processed")
animals <- readRDS(file.path(processed_dir, "animal_table.rds"))

# ----------------------------------------------------------------------------
# 1. Sex adjustment WITHIN genotype: remove the additive sex shift per analyte
#    while preserving the flicker effect. (Equivalent to a one-factor
#    removeBatchEffect; done by hand so it's transparent and dependency-free.)
#    Use these values only for pooled-sex figures (heatmaps, PLSDA). The
#    by-sex comparisons and all correlations keep the raw values.
# ----------------------------------------------------------------------------
adj <- animals
for (a in cytokines_all) {
  for (g in genotypes) {
    idx <- adj$Genotype == g
    if (!any(idx)) next
    v  <- adj[[a]][idx]
    sx <- adj$Sex[idx]
    sex_mean   <- ave(v, sx, FUN = function(z) mean(z, na.rm = TRUE))
    grand_mean <- mean(v, na.rm = TRUE)
    adj[[a]][idx] <- v - sex_mean + grand_mean
  }
}
cat("[02] sex-adjusted", length(cytokines_all), "analytes within genotype\n")

# ----------------------------------------------------------------------------
# 2a. Loading-curve evidence: does the signal rise with the amount of lysate?
#     Sheet layout: header on row 2, 8 loading levels (Sample = uL loaded),
#     then two Background0 wells. Background-subtract, clip at 0, fit
#     MFI ~ loading per analyte per genotype.
# ----------------------------------------------------------------------------
lr_fit_sheet <- function(sheet) {
  d <- suppressMessages(read_excel(lr_xlsx, sheet = sheet, range = "A2:AH12"))
  d <- d[!is.na(d$Sample), , drop = FALSE]
  bg    <- d[d$Sample == "Background0", , drop = FALSE]
  curve <- d[d$Sample != "Background0", , drop = FALSE]
  load  <- suppressWarnings(as.numeric(curve$Sample))
  keep  <- !is.na(load); curve <- curve[keep, , drop = FALSE]; load <- load[keep]

  raw_names <- setdiff(colnames(d), c("Location", "Sample"))
  out <- lapply(raw_names, function(nm) {
    y_bg <- mean(suppressWarnings(as.numeric(bg[[nm]])), na.rm = TRUE)
    y    <- pmax(suppressWarnings(as.numeric(curve[[nm]])) - y_bg, 0)
    if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0)
      return(data.frame(raw = nm, slope = 0, r2 = 0, span = 0))
    f <- stats::lm(y ~ load)
    data.frame(raw   = nm,
               slope = unname(coef(f)[2]),
               r2    = summary(f)$r.squared,
               span  = max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
  })
  out <- do.call(rbind, out)
  # translate loading-curve spelling -> combined-file spelling
  out$analyte <- ifelse(out$raw %in% names(lr_name_map),
                        lr_name_map[out$raw], out$raw)
  out
}

lr <- lapply(names(lr_sheets), function(g) {
  f <- lr_fit_sheet(lr_sheets[[g]])
  setNames(f[, c("analyte", "slope", "r2")],
           c("analyte", paste0(tolower(g), "_slope"), paste0(tolower(g), "_r2")))
})
lr_tab <- Reduce(function(a, b) merge(a, b, by = "analyte", all = TRUE), lr)
missing_lr <- setdiff(cytokines_all, lr_tab$analyte)
if (length(missing_lr))
  cat("[02] WARNING: no loading-curve column for:", paste(missing_lr, collapse = ", "), "\n")

# ----------------------------------------------------------------------------
# 2b. Floor / "shows no change" screen on the animals themselves.
# ----------------------------------------------------------------------------
screen <- do.call(rbind, lapply(cytokines_all, function(a) {
  v <- animals[[a]]; v <- v[!is.na(v)]
  data.frame(analyte    = a,
             median     = median(v),
             min        = min(v),
             max        = max(v),
             cv         = if (mean(v) > 0) sd(v) / mean(v) else NA_real_,
             pct_zero   = mean(v <= 0),
             n_distinct = length(unique(v)))
}))
screen$no_change_flag <- with(screen,
  n_distinct <= nochange_max_distinct |
  (max - min) <= nochange_max_range |
  pct_zero   >= nochange_max_pctzero)

# ----------------------------------------------------------------------------
# 3. The decision.
# ----------------------------------------------------------------------------
tab <- merge(screen, lr_tab, by = "analyte", all.x = TRUE)
tab$curated_linear <- tab$analyte %in% cytokines_linear
tab$is_myeloid     <- tab$analyte %in% cytokines_myeloid
tab$lr_pass <- with(tab, pmax(wt_r2, `5xfad_r2`, na.rm = TRUE) >= lr_min_r2 &
  ifelse(wt_r2 >= `5xfad_r2`, wt_slope, `5xfad_slope`) > lr_min_slope)

tab$used <- tab$curated_linear
if (isTRUE(apply_nochange_screen)) tab$used <- tab$used & !tab$no_change_flag
if (isTRUE(keep_myeloid_regardless)) tab$used <- tab$used | tab$is_myeloid

tab$reason <- with(tab, ifelse(
  used, "in curated linear range; passes no-change screen",
  ifelse(!curated_linear & no_change_flag, "outside linear range; floor / no change",
  ifelse(!curated_linear,                  "outside linear range",
                                           "in curated linear range but flagged: floor / no change"))))

analysis_analytes <- cytokines_all[cytokines_all %in% tab$analyte[tab$used]]
dropped           <- setdiff(cytokines_all, analysis_analytes)

tab <- tab[match(cytokines_all, tab$analyte),
           c("analyte", "wt_slope", "wt_r2", "5xfad_slope", "5xfad_r2", "lr_pass",
             "median", "min", "max", "cv", "pct_zero", "n_distinct", "no_change_flag",
             "curated_linear", "is_myeloid", "used", "reason")]
write.csv(tab, file.path(stats_dir, "analyte_linear_range.csv"), row.names = FALSE)

cat("[02] analytes used:", length(analysis_analytes), "of", length(cytokines_all), "\n")
cat("[02]   kept   :", paste(analysis_analytes, collapse = ", "), "\n")
cat("[02]   dropped:", paste(dropped, collapse = ", "), "\n")

# Where the two derived criteria disagree with the curated call — printed so the
# prune can be argued (or overridden) rather than taken on faith.
d1 <- tab$analyte[tab$curated_linear & tab$no_change_flag]
if (length(d1))
  cat("[02]   curated-linear but flagged as floor/no-change (removed):",
      paste(d1, collapse = ", "), "\n")
d2 <- tab$analyte[!is.na(tab$lr_pass) & tab$lr_pass != tab$curated_linear]
if (length(d2))
  cat("[02]   loading-curve fit disagrees with the curated call for:",
      paste(d2, collapse = ", "),
      "\n         (8-point curve, noisy — reported as evidence, not used to decide)\n")

# ----------------------------------------------------------------------------
# 4. Persist for downstream stages.
# ----------------------------------------------------------------------------
saveRDS(analysis_analytes, file.path(processed_dir, "analysis_analytes.rds"))
saveRDS(adj, file.path(processed_dir, "animal_table_sexadj.rds"))
write.csv(adj, file.path(processed_dir, "animal_table_sexadj.csv"), row.names = FALSE)
cat("[02] wrote analysis_analytes.rds + animal_table_sexadj.{rds,csv}\n")

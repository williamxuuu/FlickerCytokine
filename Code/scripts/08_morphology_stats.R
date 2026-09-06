# ============================================================================
# 08_morphology_stats.R
# Per-cell glial morphology across flicker (astrocyte + microglia), within
# genotype and by sex (Fig 1/3, Supp 2/3).
#
# NOTE: the manuscript's morphology stats were computed in MATLAB. This stage
# RE-DERIVES them in R from the per-cell "outliersremoved" exports so everything
# lives in one place and uses the same AD -> ANOVA+Tukey / KW+Dunn's rule.
# Coverage is partial: the per-cell files here are tiling + process/filament
# metrics + soma volume. The headline cell-VOLUME, DENSITY and CONVEX-HULL
# violins came from a MATLAB export not included per-cell, so those are not
# reproduced here -- cross-check the canonical MATLAB numbers for those.
# ============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
if (!exists("combined_xlsx")) source("config.R")
source("lib/stats_helpers.R")

# ----------------------------------------------------------------------------
# 1. Discover per-cell files and load each metric (both sexes) into one frame.
# ----------------------------------------------------------------------------
files <- c(
  list.files(morph_dir_astro, pattern = "outliersremoved.*\\.csv$", full.names = TRUE),
  list.files(morph_dir_mg,    pattern = "^HPC-CA1_.*outliersremoved.*\\.csv$",
             full.names = TRUE))

metric_of <- function(f) {
  b <- basename(f)
  b <- sub("^HPC[-_]?(CA1_)?", "", b)
  b <- sub("(?i)(fe)?males?_outliersremoved.*$", "", b, perl = TRUE)
  trimws(b)
}
celltype_of <- function(f) if (grepl("mCytokine", f)) "Microglia" else "Astrocyte"

load_metric <- function(paths) {
  d <- do.call(rbind, lapply(paths, function(p) {
    x <- read.csv(p, check.names = FALSE)
    data.frame(value = x[[6]],                       # metric is always col 6
               Flicker = x$Flicker, Genotype = x$Genotype,
               Sex = substr(x$Sex, 1, 1),            # "Female"/"F" -> "F"
               censor = if ("OutlierCensor" %in% names(x)) x$OutlierCensor else 0,
               stringsAsFactors = FALSE)
  }))
  d <- d[is.na(d$censor) | d$censor == 0, ]          # honor outlier flags
  d$Flicker  <- factor(d$Flicker,  levels = flicker_levels)
  d$Genotype <- factor(ifelse(d$Genotype %in% c("5XFAD","5xFAD"), "5XFAD", "WT"),
                       levels = genotypes)
  d$Sex <- factor(d$Sex, levels = sexes)
  d[stats::complete.cases(d[, c("value","Flicker","Genotype","Sex")]), ]
}

# group files by (celltype, metric)
keys <- paste(vapply(files, celltype_of, ""), vapply(files, metric_of, ""), sep = " :: ")
groups <- split(files, keys)
cat("[08] per-cell metrics found:", length(groups), "\n")

# ----------------------------------------------------------------------------
# 2. Corrected group stats per metric: pooled sex within genotype, and by sex.
# ----------------------------------------------------------------------------
om <- list()
for (k in names(groups)) {
  parts <- strsplit(k, " :: ", fixed = TRUE)[[1]]
  ct <- parts[1]; metric <- parts[2]
  d  <- load_metric(groups[[k]])
  strata <- c(list(list(g = "WT", s = NA), list(g = "5XFAD", s = NA)),
              unlist(lapply(genotypes, function(g) lapply(sexes,
                     function(s) list(g = g, s = s))), recursive = FALSE))
  for (st in strata) {
    dd <- d[d$Genotype == st$g, , drop = FALSE]
    lab <- paste(ct, st$g)
    if (!is.na(st$s)) { dd <- dd[dd$Sex == st$s, ]; lab <- paste(ct, st$g, st$s) }
    if (nrow(dd) < 4) next
    r <- compare_flicker(dd, "value", subgroup = lab, variable = metric)
    om[[paste(k, lab)]] <- r$omnibus
    # posthoc rows carry the pairwise corrected p
    if (!is.null(r$posthoc)) om[[paste(k, lab, "ph")]] <- NULL  # (kept in plots)
  }
}
morph_stats <- do.call(rbind, om)
write_stats(morph_stats, "morphology_percell_omnibus")

# ----------------------------------------------------------------------------
# 3. Violin plots (per-cell), the astrocyte tiling one being the key Fig 1/3
#    comparison. One PDF per (celltype, metric), genotype-faceted, corrected p.
# ----------------------------------------------------------------------------
for (k in names(groups)) {
  parts <- strsplit(k, " :: ", fixed = TRUE)[[1]]; ct <- parts[1]; metric <- parts[2]
  d <- load_metric(groups[[k]])
  subt <- sapply(genotypes, function(g) {
    r <- compare_flicker(d[d$Genotype == g, ], "value", variable = metric)
    sprintf("%s: %s p=%s", g, r$omnibus$test, fmt(r$omnibus$p))
  })
  p <- ggplot(d, aes(x = Flicker, y = value, fill = Flicker)) +
    geom_violin(alpha = 0.4, scale = "width", trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5) +
    facet_wrap(~ Genotype) +
    scale_fill_manual(values = flicker_colors_wt) +
    labs(x = NULL, y = metric, title = paste(ct, "-", metric),
         subtitle = paste(subt, collapse = "   |   ")) +
    theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8))
  fn <- file.path(morph_dir, sprintf("%s_%s.pdf", ct, make.names(metric)))
  ggsave(fn, p, width = 8, height = 4)
}
cat("[08] morphology per-cell stats + violins done ->", morph_dir, "\n")
cat("[08] (re-derivation of MATLAB morphology; volume/density/convex-hull not in per-cell data)\n")

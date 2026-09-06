# ============================================================================
# 03_cytokine_heatmaps.R
# heatmap3 of z-scored cytokines per genotype (Fig 2 WT, Fig 3 5xFAD), rows =
# animals ordered by flicker (Flicker + Sex side bars), columns = analytes
# clustered by euclidean/average. 5xFAD also carries Abeta40/42 side bars.
# Faithful to the original; the only change is the analyte set is now the
# linear-range panel and values are sex-adjusted (see 02). No p-values here --
# those were never on the heatmap; the "remove p-values" correction is about the
# per-analyte violins in stage 05.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(heatmap3); library(gplots)
})
if (!exists("combined_xlsx")) source("config.R")

processed_dir     <- file.path(results_dir, "processed")
adj               <- readRDS(file.path(processed_dir, "animal_table_sexadj.rds"))
analysis_analytes <- readRDS(file.path(processed_dir, "analysis_analytes.rds"))

zscore   <- function(x) { s <- sd(x, na.rm = TRUE); if (s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }
grad_col <- function(v, lo = "white", hi = "#6A51A3") {
  if (all(is.na(v))) return(rep("grey80", length(v)))
  r <- rank(v, na.last = "keep"); r <- (r - min(r, na.rm = TRUE)) /
    max(1, diff(range(r, na.rm = TRUE)))
  colorRampPalette(c(lo, hi))(100)[pmax(1, ceiling(r * 100))]
}

cyt_heatmap <- function(df, colors, title, outfile, add_amyloid = FALSE) {
  df  <- df[order(df$Flicker, df$Sex), , drop = FALSE]
  mat <- as.matrix(df[, analysis_analytes, drop = FALSE])
  matZ <- apply(mat, 2, zscore)                 # z-score within analyte (column)
  rownames(matZ) <- df$BlindName

  # column (analyte) dendrogram, euclidean/average -- matches original
  hc <- hclust(dist(t(matZ), method = "euclidean"), method = "average")

  side <- cbind(Flicker = colors[as.character(df$Flicker)],
                Sex     = c(M = "#333333", F = "#BBBBBB")[as.character(df$Sex)])
  if (add_amyloid && all(c("Abeta 40", "Abeta 42") %in% colnames(df))) {
    side <- cbind(side,
                  `Ab40` = grad_col(df[["Abeta 40"]]),
                  `Ab42` = grad_col(df[["Abeta 42"]]))
  }

  breaks <- c(-200, seq(-1.5, 1.5, 0.01), 200)
  cols   <- colorpanel(length(breaks) - 1, "blue", "white", "red2")
  pdf(outfile, width = 11, height = 8)
  heatmap3(matZ, Rowv = NA, Colv = as.dendrogram(hc), scale = "none",
           col = cols, breaks = breaks, RowSideColors = side,
           cexRow = 0.6, cexCol = 0.9, margins = c(7, 6), main = title,
           legendfun = function() showLegend(legend = names(colors),
                                             col = unname(colors)))
  dev.off()
  cat(sprintf("[03] %-28s %d animals x %d analytes -> %s\n",
              title, nrow(matZ), ncol(matZ), basename(outfile)))
}

for (g in genotypes) {
  sub <- adj[adj$Genotype == g, , drop = FALSE]
  pal <- if (g == "WT") flicker_colors_wt else flicker_colors_fad
  cyt_heatmap(sub, pal,
              title   = sprintf("%s cytokine profile (z-scored, sex-adjusted)", g),
              outfile = file.path(cyt_dir, sprintf("heatmap_%s.pdf", g)),
              add_amyloid = (g == "5XFAD"))
}

# ----------------------------------------------------------------------------
# [CORRECTION, initialComments.md line 11] "Are there sex-specific differences?
# Can show batched heat map if needed." A sex-batched heatmap: RAW (not
# sex-adjusted) values z-scored per analyte, then averaged within each of the 6
# Sex x Flicker groups, so male vs female patterns are directly visible.
# ----------------------------------------------------------------------------
raw <- readRDS(file.path(processed_dir, "animal_table.rds"))

sex_batched_heatmap <- function(df, title, outfile) {
  matZ <- apply(as.matrix(df[, analysis_analytes, drop = FALSE]), 2, zscore)
  grp  <- interaction(df$Sex, df$Flicker, sep = " ", lex.order = TRUE)
  gm   <- t(apply(matZ, 2, function(col) tapply(col, grp, mean, na.rm = TRUE)))
  gm   <- t(gm)                                   # rows = Sex x Flicker groups
  hc   <- hclust(dist(t(gm), method = "euclidean"), method = "average")
  side <- c(M = "#333333", F = "#BBBBBB")[substr(rownames(gm), 1, 1)]
  breaks <- c(-200, seq(-1.5, 1.5, 0.01), 200)
  cols   <- colorpanel(length(breaks) - 1, "blue", "white", "red2")
  pdf(outfile, width = 11, height = 5)
  heatmap3(gm, Rowv = NA, Colv = as.dendrogram(hc), scale = "none",
           col = cols, breaks = breaks, RowSideColors = side,
           cexRow = 1, cexCol = 0.9, margins = c(7, 8), main = title,
           legendfun = function() showLegend(legend = c("Male", "Female"),
                                             col = c("#333333", "#BBBBBB")))
  dev.off()
  cat(sprintf("[03] sex-batched %-6s -> %s\n", title, basename(outfile)))
}

for (g in genotypes) {
  sub <- raw[raw$Genotype == g, , drop = FALSE]
  sex_batched_heatmap(sub, g,
    file.path(cyt_dir, sprintf("heatmap_sexbatched_%s.pdf", g)))
}
cat("[03] cytokine heatmaps done ->", cyt_dir, "\n")

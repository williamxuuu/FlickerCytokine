# ============================================================================
# 00_setup.R  —  paths, analyte panel, and the shared analysis table.
# Sourced by every other stage.
#
# PIPELINE STRUCTURE. Each stage reads only the outputs of earlier stages, so a
# change upstream propagates to every downstream table and figure:
#
#   00 setup      raw MFI --> detectability filter --> analyte PANEL
#                 raw MFI --> within-sex z-score   --> ANIMALS (analysis table)
#   02 plsda      PANEL + ANIMALS --> LV1/LV2 loadings --> key cytokines
#   05/06 stats   PANEL + ANIMALS --> correlations, group comparisons
#   10-16 figures key cytokines + statistics --> figure panels
#
# No script contains a hardcoded list of selected cytokines; every analyte
# selection is derived from the data by an explicit rule.
#
# ---------------------------------------------------------------------------
# 1. ANALYTE PANEL (detectability filter)
#
#   An analyte is retained if MORE THAN 80% of animals have background-
#   subtracted MFI > 0. Background is the mean of the blank wells, with
#   negative values set to 0 (prep/01_load_and_prep.R). In this study the
#   filter retains 28 of 32 analytes (IL-5, IL-12(p40), MCP-1 and TNFa fall
#   below it).
#
#   The filter is computed on raw values, not on z-scores: "above background"
#   is a property of the assay measurement, which standardisation removes.
#
#   The plate was run without standard curves, so values are relative MFI
#   rather than concentrations, and no LOD/LLOQ is defined.
#
# ---------------------------------------------------------------------------
# 2. SEX NORMALISATION -- z-score within sex, then pool.
#
#   With n = 3-4 animals per sex x flicker cell, the sexes cannot be analysed
#   separately with useful power. Pooling raw values instead would confound
#   flicker with sex, because several measures (e.g. IL-10, RANTES, LIX, VEGF,
#   IFNy, astrocyte cell number) differ strongly between sexes.
#
#   Each measure is therefore standardised within sex, within genotype:
#
#       z = (v - mean(v | sex, genotype)) / sd(v | sex, genotype)
#
#   This removes the sex main effect by construction and lets males and females
#   be pooled (WT n = 20, 5xFAD n = 23). Flicker effects are preserved provided
#   flicker groups are balanced across sex; stage 01 verifies that balance and
#   stops the pipeline if it fails.
#
#   The same transform is applied to cytokines, amyloid and morphology, so that
#   correlations compare like with like: correlating a z-scored cytokine with a
#   raw, sex-dimorphic morphology measure would reintroduce the sex confound.
#
#   Downstream values are unitless z-scores; figure axes are labelled
#   accordingly.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
CODE    <- .here()
ROOT    <- normalizePath(file.path(CODE, ".."))          # repository root
DATA    <- file.path(ROOT, "data", "processed")
STATS   <- file.path(ROOT, "stats")
FIG_DIR <- file.path(ROOT, "figures")
DPI     <- 300
for (d in c(STATS, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

GENOTYPES      <- c("WT", "5XFAD")
FLICKER_LEVELS <- c("NoFlicker", "20Hz", "40Hz")
SEXES          <- c("M", "F")

# --- settings lib/stats_helpers.R expects as globals --------------------------
flicker_levels   <- FLICKER_LEVELS
normality_alpha  <- 0.05        # Anderson-Darling: p >= this => parametric
posthoc_padjust  <- "BH"        # Dunn's pairwise adjustment

# --- palette ----------------------------------------------------------------
# Flicker-group colours used across all manuscript figures. 20Hz and 40Hz are
# close under protanopia, so group identity is never carried by colour alone:
# groups are also labelled on the axis, and sex is shown by point shape.
COLORS_WT  <- c(NoFlicker = "#E8A941", "20Hz" = "#5C81BD", "40Hz" = "#C85EAD")
COLORS_FAD <- c(NoFlicker = "#B76945", "20Hz" = "#4E5491", "40Hz" = "#863775")
pal <- function(g) if (g == "WT") COLORS_WT else COLORS_FAD

CYTOKINES_ALL <- c(
  "G-CSF","Eotaxin","GM-CSF","IFNy","IL-1a","IL-1B","IL-2","IL-4","IL-3","IL-5",
  "IL-6","IL-7","IL-9","IL-10","IL-12(p40)","IL-12(p70)","LIF","IL-13","LIX",
  "IL-15","IL-17","IP-10","KC","MCP-1","MIP-1a","MIP-1b","M-CSF","MIP-2","MIG",
  "RANTES","VEGF","TNFa")

AMYLOID_TARGETS <- c("Abeta 40", "Abeta 42")
MORPH_TARGETS   <- c(
  "AST mean_Distance to Nearest Neighbor", "AST median_Distance to Nearest Neighbor",
  "mean_Astrocyte Cell Volume",            "median_Astrocyte Cell Volume",
  "AST mean_Number of Cells",              "AST median_Number of Cells",
  "MG mean_Distance to Nearest Neighbor",  "MG median_Distance to Nearest Neighbor",
  "mean_Microglia Cell Volume",            "median_Microglia Cell Volume",
  "MG mean_Number of Cells",               "MG median_Number of Cells")

# ---------------------------------------------------------------------------
# Load raw, gate the panel on it.
# ---------------------------------------------------------------------------
ANIMALS_RAW <- read.csv(file.path(DATA, "animal_table.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
ANIMALS_RAW$Genotype <- factor(ANIMALS_RAW$Genotype, levels = GENOTYPES)
ANIMALS_RAW$Flicker  <- factor(ANIMALS_RAW$Flicker,  levels = FLICKER_LEVELS)
ANIMALS_RAW$Sex      <- factor(ANIMALS_RAW$Sex,      levels = SEXES)

DETECT_MIN <- 0.80
frac_above_bg <- sapply(CYTOKINES_ALL, function(a) {
  v <- ANIMALS_RAW[[a]]; v <- v[!is.na(v)]; sum(v > 0) / length(v)
})
PANEL   <- CYTOKINES_ALL[frac_above_bg > DETECT_MIN]
DROPPED <- setdiff(CYTOKINES_ALL, PANEL)

# ---------------------------------------------------------------------------
# The analysis table: z-scored within sex, within genotype.
# Rebuilt on every run (it is deterministic and fast) so no stage can read a
# stale copy; a CSV is written alongside for inspection.
# ---------------------------------------------------------------------------
Z_COLS <- c(CYTOKINES_ALL, AMYLOID_TARGETS, MORPH_TARGETS)

sex_zscore <- function(df, cols, min_n = 3) {
  out <- df
  for (v in cols) {
    if (!v %in% names(df)) next
    for (g in GENOTYPES) for (s in SEXES) {
      i <- which(df$Genotype == g & df$Sex == s)
      if (!length(i)) next
      x  <- df[[v]][i]
      ok <- !is.na(x)
      # Too few observations, or no spread, to standardise on: centre only, so
      # the cell still contributes at mean 0 instead of becoming NaN.
      if (sum(ok) < min_n) { out[[v]][i] <- x - mean(x, na.rm = TRUE); next }
      s_d <- stats::sd(x, na.rm = TRUE)
      out[[v]][i] <- if (is.na(s_d) || s_d == 0) x - mean(x, na.rm = TRUE)
                     else (x - mean(x, na.rm = TRUE)) / s_d
    }
  }
  out
}

ANIMALS <- sex_zscore(ANIMALS_RAW, Z_COLS)
write.csv(ANIMALS, file.path(DATA, "animal_table_sex_zscored.csv"), row.names = FALSE)

cat(sprintf("[00] panel  : %d of %d analytes (>%.0f%% above background)\n",
            length(PANEL), length(CYTOKINES_ALL), 100 * DETECT_MIN))
cat("[00] dropped:", paste(DROPPED, collapse = ", "), "\n")
cat(sprintf("[00] table  : sex z-scored within genotype, %d animals, %d measures\n",
            nrow(ANIMALS), sum(Z_COLS %in% names(ANIMALS))))
cat("[00] figures ->", FIG_DIR, "at", DPI, "dpi\n")

# --- shared figure device ----------------------------------------------------
# ragg embeds the resolution in the PNG, so it opens at 300 dpi in other
# software rather than the 72 dpi default. Falls back to grDevices::png.
open_png <- function(file, width, height) {
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(file, width = width, height = height, units = "in", res = DPI)
  else
    grDevices::png(file, width = width, height = height, units = "in", res = DPI)
}

# Panels are written to figures/<Figure>/<stem>.
fig_path <- function(figure, stem) {
  d <- file.path(FIG_DIR, figure)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, stem)
}

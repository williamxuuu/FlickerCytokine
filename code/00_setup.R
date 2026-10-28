# ============================================================================
# 00_setup.R  —  paths, the analyte panel, and THE one analysis table.
# Source this first; every other stage depends on it.
#
# THE CASCADE. This analysis is strictly downstream-flowing. Each stage's output
# is the next stage's input, and nothing is hardcoded from a previous round:
#
#   00 setup      raw MFI --> detectability gate --> 28-analyte PANEL
#                 raw MFI --> sex z-score --> ANIMALS          (the one table)
#   02 plsda      PANEL + ANIMALS --> LV1/LV2 loadings --> KEY[[genotype]]
#   05/06 stats   PANEL + ANIMALS --> correlations, group comparisons
#   10-16 figures KEY + stats --> every panel
#
# Change the gate here and all 28-analyte claims move together. Change the
# PLS-DA and every downstream figure's analyte selection moves with it. No
# script contains a literal list of cytokine names.
#
# ---------------------------------------------------------------------------
# 1. ANALYTE PANEL (detectability gate)
#
#   An analyte is retained if MORE THAN 80% of animals have background-
#   subtracted MFI > 0. Background = mean of the two Background0 wells,
#   negatives clipped to 0 (upstream, prep/01_load_and_prep.R).
#
#   Retains 28 of 32; drops IL-5, IL-12(p40), MCP-1, TNFa.
#
#   The gate is computed on RAW values, never on the z-scored table: "above
#   background" is a statement about the assay, and z-scoring destroys it.
#
#   NOTE: this plate was run WITHOUT standards, so there is no LOD/LLOQ. Values
#   are MFI, not pg/mL, and "background" is two blank wells. State in Methods.
#
# ---------------------------------------------------------------------------
# 2. SEX HANDLING (PI's call, 2026-09) -- z-score within sex, then pool.
#
#   Per sex cell this study has n = 3-4 animals, which is not enough to analyse
#   the sexes separately. Pooling them raw is also wrong: IL-10, RANTES, LIX,
#   VEGF, IFNy and astrocyte cell number are all strongly sex-dimorphic, so a
#   pooled test would be reading a sex difference as a flicker effect.
#
#   So each measure is z-scored WITHIN sex, WITHIN genotype:
#
#       z = (v - mean(v | sex, genotype)) / sd(v | sex, genotype)
#
#   which sets every sex to mean 0, SD 1 and lets M and F be pooled into one
#   group of n = 20 (WT) / 23 (5xFAD). The sex main effect is cancelled by
#   construction; the flicker effect is preserved, because flicker is close to
#   balanced across sex within every genotype (WT 3-4 F and 3 M per arm; 5xFAD
#   4 F and 3-4 M per arm), so removing a sex mean does not remove a flicker
#   mean. Stage 01 asserts that balance and stops if it ever breaks.
#
#   This is applied to CYTOKINES, AMYLOID AND MORPHOLOGY alike -- a correlation
#   between a z-scored cytokine and a raw morphology measure would put the
#   confound straight back in through the target.
#
#   Consequence to know: values downstream are unitless z-scores, not MFI, so
#   figure axes read "(z-score)". That is the price of pooling the sexes.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
CODE    <- .here()
ROOT    <- normalizePath(file.path(CODE, ".."))          # cytokine/
PROJECT <- normalizePath(file.path(ROOT, ".."))          # AshleyPrichard/
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

# --- palette: the PUBLISHED manuscript hexes ---------------------------------
# New panels must be indistinguishable from the untouched ones already in the
# manuscript, so these are the figures' own colours, not a substitute set.
# KNOWN TRADE-OFF: 20Hz (#5C81BD) vs 40Hz (#C85EAD) separate by only dE 3.8
# under protanopia, and NoFlicker (#E8A941) sits at 2.0:1 contrast on white.
# Identity is therefore never carried by colour alone anywhere in this pipeline
# -- sex is point SHAPE and groups are direct-labelled or faceted.
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
# THE one analysis table: z-scored within sex, within genotype.
# Deterministic and cheap, so it is rebuilt here rather than read from disk --
# no stage can silently run on a stale copy. A CSV is written for the record.
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
# ragg writes the pHYs chunk, so the PNG actually reports 300 dpi to Word /
# Illustrator / Box instead of defaulting to 72. Falls back to grDevices::png.
open_png <- function(file, width, height) {
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(file, width = width, height = height, units = "in", res = DPI)
  else
    grDevices::png(file, width = width, height = height, units = "in", res = DPI)
}

# Every panel lands in figures/<Figure>/ so subfigures stay together.
fig_path <- function(figure, stem) {
  d <- file.path(FIG_DIR, figure)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, stem)
}

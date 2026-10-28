# ============================================================================
# config.R — settings for the upstream loader (01_load_and_prep.R).
#
# This file configures ONE step: reading the raw plate workbook and writing
# data/processed/animal_table.csv. Everything downstream of that — the analyte
# panel, the sex z-score, the PLS-DA, the statistics, the figures — is
# configured in code/00_setup.R, not here.
#
# Column names below are spelled exactly as they appear in the workbook.
# ============================================================================

# --- Paths -------------------------------------------------------------------
# Self-locating, so the loader runs from any working directory.
.prep_here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
PREP_DIR     <- .prep_here()                                    # code/prep/
project_root <- normalizePath(file.path(PREP_DIR, "..", ".."))  # package root
data_dir     <- file.path(project_root, "data", "raw")
results_dir  <- file.path(project_root, "data")    # -> data/processed/
stats_dir    <- file.path(project_root, "stats")
dir.create(file.path(results_dir, "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

# The one combined per-animal table the analysis runs off: metadata + 32
# cytokines + amyloid + glial-morphology summaries, one row per mouse.
# 45 rows in range = 43 animals + 2 Background0 wells.
combined_xlsx  <- file.path(data_dir, "Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx")
combined_sheet <- "Simplified"
combined_range <- "A2:BP47"      # header on row 2

# Fail loudly rather than resolving somewhere wrong. .prep_here() falls back to
# getwd() when this file is source()d interactively instead of run through
# Rscript, which would otherwise point data_dir at the wrong place silently.
if (!file.exists(combined_xlsx))
  stop("raw workbook not found at:\n  ", combined_xlsx,
       "\nExpected <package>/data/raw/. If you source()d this file by hand, run\n",
       "  Rscript code/prep/01_load_and_prep.R\nfrom the package instead.")

# --- Experimental design -----------------------------------------------------
# NoFlicker is the reference level; WT before 5XFAD.
genotypes      <- c("WT", "5XFAD")
flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
sexes          <- c("M", "F")

# --- Analytes ----------------------------------------------------------------
# All 32 Luminex analytes. The loader background-subtracts every one of them;
# the 28-analyte panel is selected downstream in 00_setup.R, not here.
cytokines_all <- c(
  "G-CSF", "Eotaxin", "GM-CSF", "IFNy", "IL-1a", "IL-1B", "IL-2", "IL-4",
  "IL-3", "IL-5", "IL-6", "IL-7", "IL-9", "IL-10", "IL-12(p40)", "IL-12(p70)",
  "LIF", "IL-13", "LIX", "IL-15", "IL-17", "IP-10", "KC", "MCP-1", "MIP-1a",
  "MIP-1b", "M-CSF", "MIP-2", "MIG", "RANTES", "VEGF", "TNFa"
)

# --- Amyloid + morphology columns, carried through unchanged -----------------
amyloid_vars <- c("Abeta 40", "Abeta 42", "Number of Plaques",
                  "Total Plaque Volume", "Ratio 42 to 40")

# Per-animal morphology summaries used for the cytokine<->morphology
# correlations: tiling (distance to nearest neighbor) and density (number of
# cells), each as mean and median.
morph_vars <- c(
  "AST mean_Distance to Nearest Neighbor", "AST median_Distance to Nearest Neighbor",
  "mean_Astrocyte Cell Volume",            "median_Astrocyte Cell Volume",
  "AST mean_Number of Cells",              "AST median_Number of Cells",
  "MG mean_Distance to Nearest Neighbor",  "MG median_Distance to Nearest Neighbor",
  "mean_Microglia Cell Volume",            "median_Microglia Cell Volume",
  "MG mean_Number of Cells",               "MG median_Number of Cells"
)

cat("[config] data dir:", data_dir, "\n")
cat("[config] analytes:", length(cytokines_all), "\n")

# ============================================================================
# config.R — central configuration for the flicker cytokine/morphology pipeline
# Edit values here before running run_pipeline.R. No need to touch the scripts.
#
# This is a corrected re-implementation of the Wood-lab analysis
# (original: SARA Temp/AshleyDataAnalysis_SBversion_Final.R, S. Bitarafan 2025).
# Corrections applied per the Wood-lab feedback in ../initialComments.md and the
# REVISED spec (2026-08) — every change is tagged  # [CORRECTION] ...  where it
# happens, so it's easy to diff against the original intent.
#
# REVISED SPEC (2026-08), in order:
#   1. heatmaps, batch-corrected for sex
#   2. remove cytokines by linear-range criteria (suspicious / show no change)
#   3. PLS-DA, no stats on the latent variables
#   4. permutation analysis as in Garza et al. (the "Kristie" papers)
#   5. FDR over the pre-selected heatmap panel, not all 32 analytes
#   6. key cytokines = top 3 on LV1 + top 3 on LV2
#   7. group comparisons on the key cytokines: ANOVA+Tukey / KW+Dunn / LMM
#   8. Spearman correlations, p-value driven, males/females marked
# ============================================================================

# --- Paths --------------------------------------------------------------------
project_root <- normalizePath(".", mustWork = TRUE)
# PACKAGE LAYOUT: the folders mirror the lab working copy ("SARA Temp/"), so each
# raw input sits where that folder keeps it rather than in one flat directory.
# The files are byte-identical to the originals; only these paths changed.
data_dir     <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
apf_dir      <- file.path(data_dir, "Ashley Prichard Flicker Paper 2")

# The one combined per-animal table the original analysis runs off of:
# metadata + 32 cytokines + amyloid + glial-morphology summaries, one row/mouse.
combined_xlsx  <- file.path(data_dir, "Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx")
combined_sheet <- "Simplified"
combined_range <- "A2:BP47"      # header on row 2, 45 animals (matches original)

# Linear-range (loading-curve) source: 8-point lysate loading series per genotype
# plus two Background0 wells, 32 analytes. This is the data original/LR Code.R
# plotted but never scored.
lr_xlsx    <- file.path(apf_dir, "BCA", "BCA_Analyzed_PrichardSamples_11425_V2.xlsx")
lr_sheets  <- c(WT = "WT", `5XFAD` = "FAD")

# Per-cell morphology (one row = one cell) for the violin-level comparisons.
# Two vintages exist; the 08-May-2025 astrocyte set is the newer one.
morph_dir_astro <- file.path(apf_dir, "Morphology",
                             "pca_morphology_data_2Dspace",
                             "pca_morphology_data_2Dspace")
morph_dir_mg    <- file.path(apf_dir, "mCytokine")   # HPC-CA1_*.csv, as in the lab folder

results_dir <- file.path(project_root, "results")
cyt_dir     <- file.path(results_dir, "cytokine")     # heatmaps, PLSDA, per-analyte
amy_dir     <- file.path(results_dir, "amyloid")      # Abeta bars
corr_dir    <- file.path(results_dir, "correlations") # Spearman cyt~amyloid, cyt~morphology
morph_dir   <- file.path(results_dir, "morphology")   # per-cell violins
stats_dir   <- file.path(results_dir, "stats")        # tidy stats tables (one row per test)
for (d in c(results_dir, cyt_dir, amy_dir, corr_dir, morph_dir, stats_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Experimental design ------------------------------------------------------
# 12 groups: Genotype (2) x Flicker (3) x Sex (2). Comparisons are made WITHIN
# genotype (and, where morphology diverges, split by sex) as the Wood lab asked.
genotypes      <- c("WT", "5XFAD")
flicker_levels <- c("NoFlicker", "20Hz", "40Hz")   # NoFlicker is the reference
flicker_ref    <- "NoFlicker"
sexes          <- c("M", "F")

# --- Analytes -----------------------------------------------------------------
# All 32 Luminex analytes, spelled exactly as the combined-file columns.
cytokines_all <- c(
  "G-CSF", "Eotaxin", "GM-CSF", "IFNy", "IL-1a", "IL-1B", "IL-2", "IL-4",
  "IL-3", "IL-5", "IL-6", "IL-7", "IL-9", "IL-10", "IL-12(p40)", "IL-12(p70)",
  "LIF", "IL-13", "LIX", "IL-15", "IL-17", "IP-10", "KC", "MCP-1", "MIP-1a",
  "MIP-1b", "M-CSF", "MIP-2", "MIG", "RANTES", "VEGF", "TNFa"
)

# Analyte names differ between the loading-curve workbook and the combined file.
# Map loading-curve spelling -> combined-file spelling. Anything not listed is
# assumed identical.
lr_name_map <- c(
  "IFN-g"    = "IFNy",
  "IL-1b"    = "IL-1B",
  "IL-12p40" = "IL-12(p40)",
  "IL-12p70" = "IL-12(p70)",
  "TNF-a"    = "TNFa"
)

# [CORRECTION] "Prune out targeted based on linear range criteria" and (revised
# spec) "remove cytokines based on the linear range criteria — suspicious ones
# that show no change". The Wood lab's own curated in-linear-range panel; it
# matches their folder SARA Temp/.../DotBox plots/Zscored/All/Only Linear
# Cytokines/ exactly, so it is the authoritative call. Stage 02 re-derives the
# supporting evidence (loading-curve fit + a floor/no-change screen) and reports
# it alongside, so the prune is documented rather than asserted.
cytokines_linear <- c(
  "IFNy", "IL-10", "IL-13", "IL-17", "IL-1B", "IL-1a", "IL-2", "IL-4", "IL-9",
  "IP-10", "KC", "LIX", "MIG", "MIP-1a", "MIP-1b", "RANTES", "VEGF"
)

# [REVISED SPEC] Strict prune. TRUE would re-add the myeloid panel even when it
# falls outside linear range (the previous round's compromise); the revised spec
# selects key cytokines from the PLS-DA loadings instead, so the myeloid panel no
# longer needs protecting and the prune can be applied as written.
keep_myeloid_regardless <- FALSE

# Myeloid-lineage analytes the manuscript used to focus on. Kept as a LABEL only
# (to annotate which selected cytokines happen to be myeloid) — no longer drives
# analyte selection or the FDR family.
cytokines_myeloid <- c("G-CSF", "GM-CSF", "M-CSF", "LIX", "LIF", "MIG")

# --- Linear-range criteria (stage 02) -----------------------------------------
# Loading-curve criterion: signal should rise with the amount of lysate loaded.
lr_min_r2       <- 0.50   # R^2 of MFI ~ loading, in the better-fitting genotype
lr_min_slope    <- 0      # slope must be positive
# Floor / "shows no change" screen on the 43 animals (background-subtracted MFI).
# An analyte sitting on the assay floor has a tiny median, few distinct values
# (quantised readings) and a large share of animals clipped to zero.
apply_nochange_screen <- TRUE  # set FALSE to prune on the curated list alone
nochange_max_distinct <- 10    # distinct values across the 43 animals
nochange_max_range    <- 5     # max - min MFI above background
nochange_max_pctzero  <- 0.20  # fraction of animals clipped to zero

# --- Amyloid + morphology columns (names as they appear in the combined file) -
amyloid_vars <- c("Abeta 40", "Abeta 42", "Number of Plaques",
                  "Total Plaque Volume", "Ratio 42 to 40")

# Morphology summaries used for the cytokine<->morphology correlations. The
# manuscript leans on tiling (distance to nearest neighbor) and density (number
# of cells); mean and median are both available per animal.
morph_vars <- c(
  "AST mean_Distance to Nearest Neighbor", "AST median_Distance to Nearest Neighbor",
  "mean_Astrocyte Cell Volume",            "median_Astrocyte Cell Volume",
  "AST mean_Number of Cells",              "AST median_Number of Cells",
  "MG mean_Distance to Nearest Neighbor",  "MG median_Distance to Nearest Neighbor",
  "mean_Microglia Cell Volume",            "median_Microglia Cell Volume",
  "MG mean_Number of Cells",               "MG median_Number of Cells"
)

# --- PLS-DA (stage 04) --------------------------------------------------------
# [CORRECTION] "Permutation analysis as in Kristie's studies" — Garza et al. 2020
# J Neurosci 40:1211 and Sci Adv 2023 eadf5672:
#   - orthogonal rotation in the LV1-LV2 plane
#   - leave-one-out cross-validation, loadings reported as mean +/- SD
#   - animals randomly reassigned to groups and the PLS-DA refit; p_permute is
#     the fraction of permuted models matching or beating the true model.
plsda_n_lv         <- 2      # latent variables
plsda_loocv_runs   <- 1000   # LOOCV iterations for loadings mean +/- SD
plsda_n_perm       <- 1000   # label permutations for p_permute
plsda_ropls_permI  <- 500    # ropls' own R2Y/Q2 permutation (secondary metric)

# [REVISED SPEC] Key cytokines = top N by |loading| on LV1 and top N on LV2,
# selected per genotype from that genotype's own model.
key_n_per_lv <- 3

# --- Statistics ---------------------------------------------------------------
alpha            <- 0.05      # significance threshold
sig_digits       <- 3         # [CORRECTION] Singer note: report ~3 sig figs everywhere
normality_alpha  <- 0.05      # Anderson-Darling: p >= this => treat as normal (parametric)

# [CORRECTION] Per-analyte flicker comparisons (N ~ 4-8/group, within genotype):
#   normal      -> one-way ANOVA + Tukey HSD
#   non-normal  -> Kruskal-Wallis + Dunn's, corrected for multiple comparisons
# (comments.md: "ANOVA and multiple comparisons correction (Tukeys) and if the
#  data is not parametric use KW and Dunns Post Hoc Corrections")
# [REVISED SPEC] Tukey/Dunn IS the correction — no second FDR layer across
# analytes is used for inference. The p_omnibus_fdr column stays in the CSVs as a
# diagnostic only.
posthoc_padjust  <- "BH"      # Dunn's / pairwise p-adjust method (BH = FDR)

# [CORRECTION] Correlations are Spearman for the STAT, Pearson line only as a
# visual trend (comments.md). The p-value drives inclusion, not a rho cutoff.
corr_stat_method <- "spearman"   # rho + p
corr_line_method <- "lm"         # Pearson/linear trend line for display only
corr_p_threshold <- 0.05         # p-driven inclusion (no |rho| cutoff)
# [REVISED SPEC] "FDR correct based on pre-selected cytokines based on the
# heatmap - not all 32 possible cytokines": BH across the retained heatmap panel,
# separately per target. "myeloid" and "key" remain available for sensitivity
# tables but are not the primary output.
corr_fdr_scope   <- "by_row"

# --- Colors -------------------------------------------------------------------
# Flicker palette to match the manuscript figure legends
# (NoFlicker warm, 20Hz blue, 40Hz pink/purple). Tune hex to the final figures.
flicker_colors_wt   <- c(NoFlicker = "#E8A941", "20Hz" = "#5C81BD", "40Hz" = "#C85EAD")
flicker_colors_fad  <- c(NoFlicker = "#B76945", "20Hz" = "#4E5491", "40Hz" = "#863775")
sex_shapes          <- c(M = 16, F = 17)   # [CORRECTION] label M/F on scatterplots

# --- Misc ---------------------------------------------------------------------
rng_seed <- 42
set.seed(rng_seed)

# Round-and-format helper so every table/figure obeys sig_digits consistently.
fmt <- function(x, d = sig_digits) formatC(signif(x, d), format = "g", digits = d)

cat("[config] project root:", project_root, "\n")
cat("[config] data dir    :", data_dir, "\n")
cat("[config] analytes: all =", length(cytokines_all),
    "| curated linear =", length(cytokines_linear),
    "| strict prune =", !keep_myeloid_regardless, "\n")

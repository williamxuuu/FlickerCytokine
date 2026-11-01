# ============================================================================
# 17_data_folders.R  —  the 28-analyte results, filed by topic.
#
# figures/ is organised by manuscript figure and stats/ is one flat folder.
# This stage re-files the same results by TOPIC, in the layout the round-2
# (16-analyte) results/ folders used, so anyone working from that layout finds
# the current numbers where they expect them:
#
#   data/cytokine/      per-animal values (raw MFI + sex z-scored), the panel
#                       gate, flicker group stats, PLS-DA loadings / scores /
#                       key cytokines; plots/ = heatmaps, PLS-DA, key panels
#   data/amyloid/       per-animal Abeta + plaque values, flicker group stats;
#                       plots/ = amyloid bars
#   data/morphology/    per-animal glial summaries (raw + sex z-scored)
#   data/correlations/  cytokine ~ amyloid / morphology tables, the BH q < 0.05
#                       pairs, one scatter per such pair, and plots/ = the
#                       correlation figures
#
# NOTHING IS COMPUTED HERE THAT A STAGE BEFORE IT DID NOT ALREADY DECIDE. Every
# number is either the analysis table from 00_setup.R, a stats/ CSV, or the
# stage-02 PLS-DA cache, and every plot is a copy of a figures/ panel -- except
# the per-pair scatters, which are drawn with the same scat() as Fig 5 from the
# stage-05 statistics. So this stage must run LAST, after the figure stages.
#
# No morphology group statistics or plots: the per-cell exports those need are
# not in the source data (README "Known limits"). The per-animal summaries here
# are the correlation targets, nothing more.
#
# Each folder is cleared before it is written, so a pair that stops passing BH
# on a re-run disappears instead of lingering from the last run.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- .here()
source(file.path(HERE, "00_setup.R"))
source(file.path(HERE, "03_theme.R"))
source(file.path(HERE, "04_panels.R"))           # ms_lv_loadings(), as Fig 2D draws it
source(file.path(HERE, "07_figure_helpers.R"))

OUT  <- file.path(ROOT, "data")
META <- c("Sample Number", "Subject", "BlindName", "Genotype", "Sex", "Flicker",
          "Cohort", "litter")
AMYLOID_ALL <- c(AMYLOID_TARGETS, "Number of Plaques", "Total Plaque Volume", "Ratio 42 to 40")

topic_dir <- function(topic) {
  d <- file.path(OUT, topic)
  unlink(d, recursive = TRUE)
  dir.create(file.path(d, "plots"), recursive = TRUE, showWarnings = FALSE)
  d
}
put_csv <- function(df, dir, name) {
  write.csv(df, file.path(dir, name), row.names = FALSE)
  cat(sprintf("     -> %s/%s\n", basename(dir), name))
}
put_stats <- function(dir, files) for (f in files)
  stopifnot(file.copy(file.path(STATS, f), file.path(dir, f), overwrite = TRUE))
# copy a figures/ panel (png + pdf + svg) into <topic>/plots/
put_plots <- function(dir, figure, stems) for (s in stems) for (ext in c(".png", ".pdf", ".svg"))
  stopifnot(file.copy(paste0(fig_path(figure, s), ext),
                      file.path(dir, "plots", paste0(s, ext)), overwrite = TRUE))
values <- function(tbl, cols, genos = GENOTYPES) {
  d <- tbl[tbl$Genotype %in% genos, c(META, cols), drop = FALSE]
  d[order(d$Genotype, d$Flicker, d$Sex), , drop = FALSE]
}

# ================================ cytokine ==================================
d <- topic_dir("cytokine")
put_csv(values(ANIMALS_RAW, PANEL), d, "cytokine_values_raw_MFI.csv")
put_csv(values(ANIMALS,     PANEL), d, "cytokine_values_sex_zscored.csv")
put_csv(data.frame(analyte = CYTOKINES_ALL,
                   frac_above_background = unname(frac_above_bg),
                   in_panel = CYTOKINES_ALL %in% PANEL),
        d, "cytokine_panel_gate.csv")
put_stats(d, c("group_omnibus_cytokines.csv", "group_posthoc_cytokines.csv",
               "plsda_key_cytokines.csv"))

pl <- lapply(GENOTYPES, function(g) readRDS(file.path(STATS, sprintf("plsda_%s.rds", g))))
names(pl) <- GENOTYPES
put_csv(do.call(rbind, lapply(GENOTYPES, function(g) {
  cc <- pl[[g]]
  data.frame(genotype = g, analyte = cc$analytes,
             LV1 = cc$loadings[, 1] * cc$sign[1], LV2 = cc$loadings[, 2] * cc$sign[2],
             LV1_normalised = ms_lv_loadings(cc, 1)[cc$analytes],
             LV2_normalised = ms_lv_loadings(cc, 2)[cc$analytes],
             key_LV1 = cc$analytes %in% cc$key_by_lv[[1]],
             key_LV2 = cc$analytes %in% cc$key_by_lv[[2]], row.names = NULL)
})), d, "plsda_loadings.csv")
put_csv(do.call(rbind, lapply(GENOTYPES, function(g) {
  cc <- pl[[g]]; i <- as.integer(rownames(cc$scores))
  data.frame(genotype = g, Subject = ANIMALS$Subject[i], Sex = cc$sex,
             Flicker = cc$flicker, LV1 = cc$scores[, 1], LV2 = cc$scores[, 2],
             row.names = NULL)
})), d, "plsda_scores.csv")
put_csv(do.call(rbind, lapply(GENOTYPES, function(g) {
  cc <- pl[[g]]
  data.frame(genotype = g, n_animals = nrow(cc$scores), n_analytes = length(cc$analytes),
             loocv_accuracy = cc$obs, chance = 1 / 3, permutation_p = cc$pperm,
             R2X_LV1 = cc$r2x[1], R2X_LV2 = cc$r2x[2])
})), d, "plsda_model_summary.csv")

put_plots(d, "Fig2", c("Fig2B_heatmap_WT", "Fig2C_plsda_composite_WT",
                       "Fig2D_loadings_WT_LV1", "Fig2D_loadings_WT_LV2",
                       "Fig2E_LV1_scores_by_group_WT"))
put_plots(d, "Supp4", c("Supp4A_heatmap_5XFAD", "Supp4B_plsda_composite_5XFAD",
                        "Supp4CD_loadings_5XFAD_LV1", "Supp4CD_loadings_5XFAD_LV2",
                        "Supp4E_LV1_scores_by_group_5XFAD"))
put_plots(d, "Supp5", c("Supp5_key_cytokines_WT", "Supp5_key_cytokines_5XFAD"))
put_plots(d, "Supp6", c("Supp6A_cytokine_profiling_WT", "Supp6B_cytokine_profiling_5XFAD"))

# ================================ amyloid ===================================
d <- topic_dir("amyloid")
put_csv(values(ANIMALS_RAW, AMYLOID_ALL,     "5XFAD"), d, "amyloid_values_raw.csv")
put_csv(values(ANIMALS,     AMYLOID_TARGETS, "5XFAD"), d, "amyloid_values_sex_zscored.csv")
put_stats(d, c("group_omnibus_amyloid.csv", "group_posthoc_amyloid.csv"))
put_plots(d, "Fig4", "Fig4B_amyloid_bars")

# =============================== morphology =================================
d <- topic_dir("morphology")
unlink(file.path(d, "plots"), recursive = TRUE)      # nothing reproducible to plot
put_csv(values(ANIMALS_RAW, MORPH_TARGETS), d, "morphology_values_raw.csv")
put_csv(values(ANIMALS,     MORPH_TARGETS), d, "morphology_values_sex_zscored.csv")

# ============================== correlations ================================
d <- topic_dir("correlations")
put_stats(d, c("corr_cytokine_amyloid.csv", "corr_cytokine_morphology.csv",
               "corr_morphology_by_selection.csv"))
CR <- list(cytokine_amyloid    = read.csv(file.path(STATS, "corr_cytokine_amyloid.csv"),
                                          check.names = FALSE),
           cytokine_morphology = read.csv(file.path(STATS, "corr_cytokine_morphology.csv"),
                                          check.names = FALSE))
sig <- do.call(rbind, lapply(names(CR), function(f)
  cbind(family = f, CR[[f]][CR[[f]]$p_adj < 0.05, , drop = FALSE])))
put_csv(sig[order(sig$p_adj, sig$p_wy), ], d, "corr_significant_BH.csv")

# one scatter per BH-significant pair, in the round-2 per-pair layout
for (i in seq_len(nrow(sig))) {
  r  <- sig[i, ]
  sd <- file.path(d, r$family); dir.create(sd, showWarnings = FALSE)
  xl <- if (r$target %in% AMYLOID_TARGETS) ab_lab(r$target) else zlab(r$target)
  p  <- scat(ANIMALS[ANIMALS$Genotype == r$subgroup, , drop = FALSE], r$target, r$analyte,
             xl, zlab(r$analyte), r$subgroup, r$rho, r$p_perm, r$p_adj, k = r$n_distinct) +
        labs(title = sprintf("%s: %s vs %s", r$subgroup, r$analyte, r$target)) +
        theme(plot.title = element_text(size = 10, face = "bold"))
  f <- file.path(sd, gsub("[^A-Za-z0-9.()-]+", "_",
                          sprintf("%s_%s_vs_%s", r$subgroup, r$analyte, r$target)))
  write_fig(p, f, 5.2, 4.2)
}
cat(sprintf("     -> correlations/: %d per-pair scatters (BH q < 0.05)\n", nrow(sig)))

put_plots(d, "Fig4", "Fig4C_cytokine_amyloid_heatmap")
put_plots(d, "Fig5", c("Fig5A_amyloid_vs_astro_tiling", "Fig5B_amyloid_vs_mg_tiling",
                       "Fig5C_cytokine_vs_astro_number", "Fig5D_cytokine_vs_mg_number"))
put_plots(d, "Supp7", c("Supp7A_cyt_vs_astro_tiling_WT", "Supp7B_cyt_vs_mg_tiling_WT",
                        "Supp7C_cyt_vs_astro_number_WT", "Supp7D_cyt_vs_mg_number_WT"))

cat("[17] data/{cytokine,amyloid,morphology,correlations} done\n")

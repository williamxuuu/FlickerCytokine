# ============================================================================
# 01_load_and_prep.R
# Read the combined per-animal table, background-subtract the cytokines, set up
# factors, and write one tidy animal-level table the rest of the pipeline uses.
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
})

# [PACKAGING FIX 2026-09-26] Source the config by absolute path. Previously
# `source("config.R")`, which (a) required the working directory to be this
# folder and (b) named a file that did not exist -- the config shipped as
# `config_round2.R`, so this stage could not run at all. Renamed to config.R
# and located relative to this script. Loader logic below is unchanged.
.prep_here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
if (!exists("combined_xlsx")) source(file.path(.prep_here(), "config.R"))

processed_dir <- file.path(results_dir, "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 1. Read the combined "Simplified" sheet (metadata + cytokines + amyloid +
#    morphology summaries; last 2 rows are Background0 wells).
# ----------------------------------------------------------------------------
raw <- read_excel(combined_xlsx, sheet = combined_sheet, range = combined_range)
raw <- as.data.frame(raw, check.names = FALSE)
cat("[01] read", nrow(raw), "rows x", ncol(raw), "cols\n")

# Background wells = the rows with no genotype (Sample == "Background0").
is_bg   <- is.na(raw$Genotype) | raw$Sample %in% c("Background0", "background0")
bg      <- raw[is_bg, , drop = FALSE]
animals <- raw[!is_bg, , drop = FALSE]
cat("[01] background wells:", nrow(bg), "| animals:", nrow(animals), "\n")

# ----------------------------------------------------------------------------
# 2. Background-subtract the cytokines, clip negatives to 0.
#    (manuscript: "Cytokine data was subtracted from background values, and any
#     negative values were set to 0".) Amyloid/morphology are left as-is.
# ----------------------------------------------------------------------------
stopifnot(all(cytokines_all %in% colnames(animals)))
bg_means <- vapply(cytokines_all, function(a) mean(bg[[a]], na.rm = TRUE),
                   numeric(1))
for (a in cytokines_all) {
  animals[[a]] <- pmax(animals[[a]] - bg_means[a], 0)
}
cat("[01] subtracted background (mean over", nrow(bg), "wells); clipped <0 to 0\n")

# ----------------------------------------------------------------------------
# 3. Factors. NoFlicker is the reference flicker level; keep WT before 5XFAD.
# ----------------------------------------------------------------------------
animals$Genotype <- factor(animals$Genotype, levels = genotypes)
animals$Flicker  <- factor(animals$Flicker,  levels = flicker_levels)
animals$Sex      <- factor(animals$Sex,      levels = sexes)
if ("Cohort" %in% colnames(animals)) animals$Cohort <- factor(animals$Cohort)

# convenient combined labels the original script leaned on
animals$Condition   <- paste(animals$Genotype, animals$Sex, animals$Flicker)
animals$GenotypeSex <- paste(animals$Genotype, animals$Sex)

# ----------------------------------------------------------------------------
# 4. Sanity: group sizes, and note them for the stats tables (small-N design).
# ----------------------------------------------------------------------------
grp <- animals %>%
  count(Genotype, Flicker, Sex, name = "n") %>%
  arrange(Genotype, Flicker, Sex)
cat("[01] group sizes (Genotype x Flicker x Sex):\n")
print(grp, row.names = FALSE)
cat("[01] per genotype x flicker (sexes pooled):\n")
print(animals %>% count(Genotype, Flicker, name = "n"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Write the tidy animal table + a fast-reload RDS.
# ----------------------------------------------------------------------------
keep_cols <- c("Sample Number", "Subject", "BlindName", "Genotype", "Sex",
               "Flicker", "Cohort", "litter", "Condition", "GenotypeSex",
               cytokines_all, amyloid_vars, morph_vars)
keep_cols <- intersect(keep_cols, colnames(animals))
animals   <- animals[, keep_cols, drop = FALSE]

write.csv(animals, file.path(processed_dir, "animal_table.csv"), row.names = FALSE)
saveRDS(animals,  file.path(processed_dir, "animal_table.rds"))
write.csv(grp, file.path(stats_dir, "group_sizes.csv"), row.names = FALSE)
cat("[01] wrote processed/animal_table.{csv,rds} (", nrow(animals), "animals)\n")

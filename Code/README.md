# Flicker Cytokine / Morphology Pipeline (corrected)

A cleaned, staged re-implementation of the Wood-lab hippocampus cytokine +
amyloid + glial-morphology analysis behind the flicker manuscript
(*"Noninvasive audiovisual neurostimulation differentially alters glial
responses at 20Hz and 40Hz"*). It reproduces the original figures and then
applies the statistical corrections the Wood lab asked for (see
[../initialComments.md](../initialComments.md)).

The original analysis is one 3,652-line script
(`original/AshleyDataAnalysis_SBversion_Final.R`, S. Bitarafan, 2025). It's kept
verbatim in [original/](original/) for reference; nothing here overwrites it.
This pipeline breaks that script into config-driven stages and marks every
statistical change with a `# [CORRECTION]` tag pointing back to the comment.

## Design

12 groups: **Genotype** (WT, 5XFAD) × **Flicker** (NoFlicker, 20Hz, 40Hz) ×
**Sex** (M, F). All comparisons are made **within genotype**, and split by sex
where the morphology diverges. The cytokine cohort is small — n≈6–8 per flicker
group within a genotype (sex pooled), n≈3–4 when split by sex — which is exactly
why the original blanket FDR across 32 analytes was too strict.

## Inputs

- `../SARA Temp/Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx`
  (sheet *Simplified*) — the combined per-animal table: metadata + 32 cytokines
  + amyloid (Aβ40/42, plaques) + astro/microglia morphology summaries.
- `../SARA Temp/.../HPC*_*outliersremoved*.csv` — per-cell morphology (one row
  per cell) for the violin-level comparisons.
- Original Wood-lab scripts + helpers in [original/](original/); the plotting /
  PLSDA helpers are sourced from [lib/WoodLabFunctions.R](lib/WoodLabFunctions.R).

## Corrections applied (from initialComments.md)

| Wood-lab note | What the pipeline does |
|---|---|
| Prune analytes by linear range | **strict removal** to the in-range panel; the loading curve is re-fit and a floor/no-change screen applied, both written to `results/stats/analyte_linear_range.csv` (stage 02) |
| Focus analytes | **top 3 on LV1 + top 3 on LV2, per genotype**, from the PLS-DA loadings (stage 04b). Replaces the a priori myeloid panel |
| No stats on PLSDA latent variables | PLSDA kept for visualization + a **Garza-style permutation test** (LOOCV accuracy vs randomly reassigned groups); per-analyte bar/violin carry the stats |
| Remove p-values on the violins-by-heatmap | those annotations stripped |
| Batch-adjust for sex (after linear range) | sex adjustment within genotype |
| 3 flicker groups, N≈4–8/group | Anderson-Darling → **ANOVA + Tukey** (normal) or **Kruskal-Wallis + Dunn's** (non-normal) |
| Correlations: Spearman stat, Pearson line | Spearman ρ + p for the stat, linear trend line for display only |
| p-value drives inclusion (not a rho cutoff) | `corr_p_threshold`, no |ρ| gate |
| FDR "across everything" is too much | FDR **by row over the retained heatmap panel** (`corr_fdr_scope = "by_row"`); key-set-restricted tables written as a sensitivity analysis |
| Label males/females on scatterplots | sex encoded by point shape on every scatter |
| ~3 significant figures throughout | `fmt()` / `sig_digits` applied to tables + figure text |

## Configure

Edit [config.R](config.R): analyte lists, the linear-range criteria
(`apply_nochange_screen`, `nochange_max_*`, `lr_min_r2`), the PLS-DA validation
settings (`plsda_loocv_runs`, `plsda_n_perm`, `key_n_per_lv`), colors, and every
statistical knob (`normality_alpha`, `posthoc_padjust`, `corr_stat_method`,
`corr_fdr_scope`, `sig_digits`, …).

The two knobs most likely to need revisiting:

- `apply_nochange_screen` — `TRUE` removes analytes that sit on the assay floor
  (few distinct values, narrow range, many readings at zero) on top of the
  curated linear-range panel. `FALSE` prunes on the curated panel alone.
- `key_n_per_lv` — how many top-weighted analytes to take from each latent
  variable. The revised spec asks for 3.

## Run

```bash
cd "/Users/williamxu/Documents/Research/BALES/AshleyPrichard/cytokine_analysis"
Rscript run_pipeline.R
```

Or a single stage interactively:

```r
setwd("/Users/williamxu/Documents/Research/BALES/AshleyPrichard/cytokine_analysis")
source("config.R"); source("lib/WoodLabFunctions.R")
source("scripts/07_correlations.R")
```

## Outputs

```
results/
├── cytokine/      heatmaps, PLSDA scores/loadings (LOOCV error bars), key-cytokine panels
├── amyloid/       Abeta40/42 bars + tests
├── correlations/  cytokine~amyloid + cytokine~morphology (Spearman, M/F-labeled)
├── morphology/    per-cell violins (volume, tiling, density, convex hull)
└── stats/         one tidy CSV per analysis — test, statistic, raw p, adj p, N, method
```

Every reported statistic also lands in a `results/stats/*.csv` row (test used,
post-hoc, correction method, N, comparisons) so the manuscript's Supplementary
Statistics table can be regenerated instead of hand-copied.

## Sanity checks after running

1. `results/stats/` tables are non-empty and N per group matches the design above.
2. Corrected cytokine p-values are less extreme than the original raw ones but
   not all wiped out (the point of narrowing FDR scope).
3. Spot-check one correlation (e.g. MIG ~ Aβ40) against the original
   `../SARA Temp/SpearmanCorrelation_*_SummaryTable_*.csv`.
4. `results/stats/key_cytokines.csv` has `key_n_per_lv` rows per LV per genotype,
   and every analyte in it appears in `analyte_linear_range.csv` with `used=TRUE`.
5. `results/stats/plsda_permutation.csv`: `p_permute` is the model's validity.
   Where it is **not** below 0.05, the loadings ranking is descriptive only —
   the figures and the key-cytokine table say so, and the manuscript must too.

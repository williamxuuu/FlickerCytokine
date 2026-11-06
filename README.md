# FlickerCytokine

Analysis code for the hippocampal cytokine, amyloid and glia–cytokine
correlation analyses in:

> **Noninvasive audiovisual neurostimulation differentially alters glial
> responses at 20Hz and 40Hz.** (Manuscript in preparation.)

Wild-type (WT) and 5xFAD mice (6 months old, both sexes) received chronic
audiovisual flicker: No Flicker, 20 Hz or 40 Hz, 1 h/day for 7 days. The
right hippocampus was assayed with a 32-plex Luminex cytokine/chemokine panel
and, in 5xFAD mice, an Aβ1-40 / Aβ1-42 immunoassay. Glial morphology (astrocyte
and microglia tiling, volume and number) was quantified separately. Per-animal
summaries of those measures are used here as correlation targets.

The pipeline takes the raw per-animal workbook and produces every statistic
and every cytokine-related figure panel in the manuscript. It is
written so that the same analysis can be run on a new dataset with the same
structure.

**This repository contains code only.** The study data are not distributed
here. To run the pipeline, supply a workbook in the format described in
[`data/raw/README.md`](data/raw/README.md). `.gitignore` excludes `data/`,
`stats/` and `figures/`, so generated outputs are never committed.

---

## Contents

1. [Requirements](#requirements)
2. [Quick start](#quick-start)
3. [Analysis](#analysis)
4. [Pipeline stages](#pipeline-stages)
5. [Outputs](#outputs)
6. [Reproducibility](#reproducibility)
7. [Running on your own data](#running-on-your-own-data)
8. [Repository layout](#repository-layout)

---

## Requirements

R ≥ 4.5 (developed and tested on R 4.5.2, macOS). Package versions used:

| Package | Version | Used for |
|---|---|---|
| ropls (Bioconductor) | 1.42.0 | PLS-DA |
| rstatix | 0.7.3 | Tukey HSD, Dunn's test |
| nortest | 1.0.4 | Anderson-Darling normality test |
| matrixStats | 1.5.0 | resampling summaries |
| lme4, lmerTest | 2.0.1, 3.2.1 | mixed-model sensitivity analysis |
| emmeans, car | 2.0.2, 3.1.5 | mixed-model contrasts and tests |
| readxl | 1.4.5 | reading the raw workbook |
| tidyverse | 2.0.0 | data handling |
| ggplot2, patchwork, ggnewscale, scales | 4.0.2, 1.3.2, 0.5.2, 1.4.0 | figures |
| ragg, svglite | 1.5.1, 2.2.2 | PNG and SVG output |

```r
install.packages(c("tidyverse", "readxl", "matrixStats", "ggplot2", "patchwork",
                   "ggnewscale", "scales", "ragg", "svglite", "rstatix", "nortest",
                   "lmerTest", "lme4", "emmeans", "car", "BiocManager"))
BiocManager::install("ropls")
```

---

## Quick start

```bash
# 1. Place the raw workbook in data/raw/ and build the per-animal table
Rscript code/prep/01_load_and_prep.R

# 2. Run the analysis
Rscript code/run_all.R            # statistics and figures
Rscript code/run_all.R stats      # statistics only
Rscript code/run_all.R figures    # figures only (requires existing statistics)
```

A first run takes about 15 minutes on a laptop, most of it in the PLS-DA
permutation test. The results are cached, and later runs take about 3 minutes.

---

## Analysis

All analyses are performed **within genotype**, comparing the three flicker
groups.

### 1. Preprocessing

Analyte values are background-subtracted median fluorescence intensities
(MFI). Background is the mean of the blank wells on the plate, and values below
background are set to 0. The plate was run without standard curves, so values
are relative rather than absolute concentrations.
*(`code/prep/01_load_and_prep.R`)*

### 2. Analyte detectability filter

An analyte is analysed only if **more than 80% of animals** have values above
background. The filter is applied to raw values, before any transformation. In
this study it retains 28 of 32 analytes. IL-5, IL-12(p40), MCP-1 and TNF-α are
excluded. *(`code/00_setup.R`)*

### 3. Sex normalisation

Each sex × flicker cell has only 3–4 animals, too few to analyse the sexes
separately. Several measures also differ strongly between sexes, for example
IL-10, RANTES and astrocyte number. Pooling raw values would therefore
confound flicker with sex. Instead, every measure (cytokines, amyloid and
morphology) is **z-scored within sex within genotype**, and the sexes are then
pooled:

```
z = (x − mean(x | sex, genotype)) / sd(x | sex, genotype)
```

This removes sex differences in baseline by construction. It preserves flicker
effects only if sex is balanced across the flicker groups.
`01_sex_zscore_qc.R` checks this:

- **Balance.** It stops the pipeline if any group is more unbalanced than
  80/20.
- **Sex effect removed.** It reports the sex effect for each analyte before
  and after normalisation.
- **Flicker effect retained.** It reports the flicker effect for each analyte
  before and after normalisation.

All downstream values are unitless z-scores.

### 4. PLS-DA

A two-component partial least squares discriminant analysis (ropls) is
fitted with flicker group as the outcome. Each analyte is scaled to unit
variance.

- **Model validation.** Leave-one-out cross-validated (LOOCV) classification
  accuracy is compared with a null distribution from 1,000 random permutations
  of the group labels:
  p = (#{permuted accuracy ≥ observed} + 1) / 1,001.
- **Loading uncertainty.** This is the SD of each loading across 1,000 refits,
  each leaving out one randomly chosen animal. Loadings are aligned in sign
  before the SD is taken.
- **Key cytokines.** These are the three analytes with the largest absolute
  loading on LV1 and the three on LV2. They are used to select analytes for
  Supplementary Figures 5 and 7.
- **No tests on latent-variable scores.** The latent variables are constructed
  to separate the groups, so testing group differences on their scores would be
  circular. Score plots are descriptive. Evidence for the model as a whole
  comes from the permutation test.

*(`code/02_plsda.R`)*

### 5. Group comparisons

For each analyte, and for Aβ1-40 and Aβ1-42 in 5xFAD, the test depends on an
Anderson-Darling test of the residuals:

| Residuals | Omnibus test | Post hoc |
|---|---|---|
| normal (p ≥ 0.05) | one-way ANOVA | Tukey HSD |
| non-normal | Kruskal-Wallis | Dunn's test, Benjamini-Hochberg adjusted |

The post hoc procedure controls error across the three pairwise contrasts. No
further correction is applied across analytes. In figures, a pairwise
contrast is starred only when its omnibus test is significant.
*(`code/06_group_stats.R`, `code/lib/stats_helpers.R`)*

### 6. Correlations

Spearman's ρ is computed between each analyte and each target (Aβ1-40,
Aβ1-42 and per-animal glial measures), within genotype, on the sex-normalised
values.

- **P values** come from 50,000 permutations of the target. This makes the null
  distribution exact for the tie structure of the data, which matters at
  n ≈ 20 with coarsely measured analytes.
- **Multiple testing.** A family is one genotype × one target, which means 28
  tests. The same permutations are applied to every analyte in a family. Two
  corrections are reported:
  - `p_adj`: Benjamini-Hochberg false discovery rate (primary).
  - `p_wy`: Westfall-Young single-step maxT family-wise error rate (stricter).
    It is computed on studentised ρ, so analytes with many ties do not
    dominate the maximum.
- **Selection for display.** Figure 5C shows the cytokines whose correlation
  with microglial number passes both corrections. Pairs passing
  Benjamini-Hochberg only are shown in the supplement. Both corrections count
  all 28 tests, so selecting on them adds no hidden multiplicity.
- **Measurement resolution.** Every correlation row reports `n_distinct` and
  `prop_tied`, computed on raw MFI. On scatter panels, analytes measured at 6
  or fewer distinct values are labelled `[k levels]`, because quantisation can
  inflate ρ at that resolution.

Regression lines on scatter plots are least-squares fits shown for display
only.
*(`code/05_correlations.R`)*

### 7. Sensitivity analyses

| Analysis | Question | Stage |
|---|---|---|
| Sex-stratified permutation (`p_perm_sexstrat`) | Do correlation p-values depend on permuting across sexes? | 05 |
| Linear mixed model on raw values, `y ~ Flicker + Sex + (1 \| litter)`, raw and log1p scale | Does modelling sex, rather than removing it, change the flicker results? Is there a Flicker × Sex interaction? | 08 |
| Correlations under two selection rules (PLS-DA key cytokines vs most flicker-responsive analytes) | Does the choice of displayed analytes change the conclusions? | 09 |

---

## Pipeline stages

Each stage runs as a separate R process and reads only the outputs of earlier
stages. No script contains a hardcoded list of selected cytokines. Every
analyte selection is derived from the data by the rules above, so a change to
the input or to any upstream setting propagates through all tables and
figures.

```
raw workbook ─► prep/01_load_and_prep.R ─► animal table
                                              │
                         00_setup.R: detectability filter, within-sex z-score
                                              │
            ┌─────────────────┬───────────────┼────────────────┬──────────────┐
            ▼                 ▼               ▼                ▼              ▼
   01 z-score QC      02 PLS-DA ──►    05 correlations   06 group tests   08 mixed model
                      key cytokines          │                │
                            │                ▼                │
                            │         09 selection rules      │
                            └────────────┬────────────────────┘
                                         ▼
                              10–16 figure stages ─► 17 results by topic
```

| Stage | Purpose | Main outputs (`stats/`) |
|---|---|---|
| `prep/01_load_and_prep.R` | read workbook, background subtraction | `data/processed/animal_table.csv`, `group_sizes.csv` |
| `00_setup.R` | paths, detectability filter, within-sex z-score (sourced by every stage) | `data/processed/animal_table_sex_zscored.csv` |
| `01_sex_zscore_qc.R` | balance check; sex and flicker effects before vs after | `qc_sex_balance.csv`, `qc_sex_zscore_effect.csv` |
| `02_plsda.R` | PLS-DA, permutation test, loading SDs, key cytokines | `plsda_<genotype>.rds`, `plsda_key_cytokines.csv` |
| `05_correlations.R` | permutation Spearman, BH and Westfall-Young | `corr_cytokine_amyloid.csv`, `corr_cytokine_morphology.csv` |
| `06_group_stats.R` | flicker group comparisons | `group_omnibus_*.csv`, `group_posthoc_*.csv` |
| `08_mixed_model.R` | mixed-model sensitivity analysis | `mixed_model_omnibus.csv`, `mixed_model_posthoc.csv` |
| `09_corr_selection.R` | correlations under two selection rules | `corr_morphology_by_selection.csv` |
| `10`–`16` | figure panels (see below) | `figures/` |
| `17_data_folders.R` | files all results by topic under `data/` | `data/{cytokine,amyloid,morphology,correlations}/` |

---

## Outputs

### Figures

Every panel is written as a 300 dpi PNG, a vector PDF and an SVG with editable
text. Panels are written to `figures/<Figure>/`.

| Figure | Panels | Stage |
|---|---|---|
| Figure 2 (WT) | B heatmap (columns ordered by LV1 loading) · C PLS-DA scores and score distributions · D LV1/LV2 loadings · E LV1 scores by group | `10_fig2.R` |
| Figure 4 (5xFAD) | B Aβ1-40 / Aβ1-42 by flicker group · C cytokine × amyloid Spearman ρ heatmap | `11_fig4.R` |
| Figure 5 (5xFAD) | A/B amyloid vs astrocyte / microglia tiling · C cytokines vs microglial number (BH and Westfall-Young significant) | `12_fig5.R` |
| Figure 5, supplementary panels | S1 cytokines vs microglial number (BH only) · S2/S3 PLS-DA key cytokines vs astrocyte / microglial number | `12_fig5.R` |
| Supplementary Figure 4 (5xFAD) | A–E, as Figure 2, with amyloid annotation on the heatmap | `13_supp4.R` |
| Supplementary Figure 5 | key-cytokine group comparisons, WT and 5xFAD | `14_supp5.R` |
| Supplementary Figure 6 | six most flicker-responsive cytokines, shown for both genotypes | `15_supp6.R` |
| Supplementary Figure 7 (WT) | key cytokines vs astrocyte / microglia tiling and number | `16_supp7.R` |

Experimental schematics (Figure 2A, Figure 4A) and the glial morphology
figures are not generated by this code. Glial morphology group statistics come
from separate cell-level analyses.

### Statistics

All statistics are written as tidy CSV files in `stats/`, one row per test,
with the columns needed to check each result. Correlation tables carry `rho`,
`n`, `p_perm`, `p_adj`, `p_wy`, `p_perm_sexstrat`, `n_distinct` and
`prop_tied`. Group tables carry the test used, the statistic, degrees of
freedom, group sizes and adjusted pairwise p-values.

---

## Reproducibility

- **Seeds.** Every stage with a random component sets `set.seed(42)`, so
  identical input gives identical output.
- **Monte Carlo error.** Permutation p-values carry sampling error. With
  B = 50,000 (correlations), a p-value near 0.04 is stable to about ±0.002.
  With B = 1,000 (PLS-DA), the resolution is 0.001.
- **Caching.** The PLS-DA permutation and LOOCV results are cached in
  `stats/plsda_<genotype>.rds`. The cache is keyed on a fingerprint of the
  input matrix, so any change to the data, analyte panel or transform triggers
  recomputation. Delete the `.rds` files to force it.
- **Determinism checks.** Rebuilding everything from the raw workbook in a
  clean directory, with no cache, reproduces every statistics table and every
  PNG and SVG figure byte-for-byte. PDFs differ only in their embedded
  creation timestamp.

---

## Running on your own data

The pipeline assumes one row per animal, with genotype, sex and flicker
group, and the analytes as columns. To adapt it:

1. **Workbook location and layout.** Set `combined_xlsx`, `combined_sheet`
   and `combined_range` in `code/prep/config.R`.
2. **Column names.** Set `cytokines_all`, `amyloid_vars` and `morph_vars` in
   `code/prep/config.R`. The analyte list is repeated as `CYTOKINES_ALL` in
   `code/00_setup.R`, and the two must match.
3. **Design.** `GENOTYPES`, `FLICKER_LEVELS` and `SEXES` are set in
   `code/00_setup.R` (and in `config.R` for the loader).
4. **Required columns.** `Genotype`, `Sex`, `Flicker`, and `Sample` (blank
   wells are identified by a missing genotype or `Sample == "Background0"`).
   `litter` is used only by the mixed-model sensitivity analysis.
5. **Analysis settings**, all named constants:

   | Setting | Value here | Location |
   |---|---|---|
   | `DETECT_MIN` | 0.80 | `00_setup.R` |
   | `N_PERM` | 1,000 | `02_plsda.R` |
   | `N_KEY` | 3 | `02_plsda.R` |
   | `B` | 50,000 | `05_correlations.R` |
   | `normality_alpha` | 0.05 | `00_setup.R` |
   | `posthoc_padjust` | `"BH"` | `00_setup.R` |

The figure stages use this study's figure numbering and morphology column
names. The statistics stages (00–09) are general.

---

## Repository layout

```
FlickerCytokine/
├── README.md
├── code/
│   ├── run_all.R              runs all stages in order
│   ├── 00_setup.R             paths, detectability filter, within-sex z-score
│   ├── 01_sex_zscore_qc.R     z-score quality control
│   ├── 02_plsda.R             PLS-DA, permutation test, key cytokines
│   ├── 03_theme.R             figure styling; PNG / PDF / SVG output
│   ├── 04_panels.R            PLS-DA panels (Figure 2, Supplementary Figure 4)
│   ├── 05_correlations.R      permutation Spearman, BH, Westfall-Young
│   ├── 06_group_stats.R       flicker group comparisons
│   ├── 07_figure_helpers.R    shared bar and scatter panels
│   ├── 08_mixed_model.R       mixed-model sensitivity analysis
│   ├── 09_corr_selection.R    correlations under two selection rules
│   ├── 10_fig2.R … 16_supp7.R one script per figure
│   ├── 17_data_folders.R      results filed by topic
│   ├── lib/stats_helpers.R    group-comparison functions
│   └── prep/                  raw workbook → per-animal table
└── data/raw/README.md         expected input format

Generated on run (not tracked):
├── data/processed/            per-animal tables (raw and sex-normalised)
├── data/{cytokine,amyloid,morphology,correlations}/   results by topic
├── stats/                     all statistics as CSV
└── figures/                   all figure panels (PNG, PDF, SVG)
```

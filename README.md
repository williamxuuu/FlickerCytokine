# 28-Cytokine Analysis — code

Chronic audiovisual flicker (No Flicker / 20 Hz / 40 Hz, 1 h/day × 7 d) in
6-month 5xFAD and WT mice, both sexes. This package holds the **complete
cytokine analysis**: the 28-analyte panel, the within-sex z-score, the PLS-DA,
the group statistics, the correlations, and every manuscript figure panel that
comes out of them.

Packaged 2026-09-26 from `AshleyPrichard/cytokine/`. Everything here was
re-run and checked before packaging — see **[VERIFICATION.md](VERIFICATION.md)**
for what was tested and what was found.

> **This repository is code only.** The raw plate workbook, the per-animal
> tables, the statistics and the figure files are **not** included — the data
> are not public. Put the workbook in `data/raw/` (format in
> [data/raw/README.md](data/raw/README.md)) and the pipeline regenerates every
> table and figure. `.gitignore` keeps `data/`, `stats/` and `figures/` out of
> the repository so they cannot be committed by accident.

---

## Run it

```bash
# 1. once: place the plate workbook in data/raw/, then build the analysis table
Rscript code/prep/01_load_and_prep.R     # data/raw/*.xlsx -> data/processed/animal_table.csv

# 2. the analysis
cd code
Rscript run_all.R            # everything      (~13 min the first time, PLS-DA permutation)
Rscript run_all.R stats      # statistics only
Rscript run_all.R figures    # figures only (statistics must already exist)
```

**Requirements.** R 4.5.2 and sixteen packages — every one of them actually
loaded by the code, nothing carried along (lmerTest, lme4, emmeans and car are
used only by the stage-08 mixed-model cross-check):

```r
install.packages(c("tidyverse", "readxl", "matrixStats", "ggplot2", "patchwork",
                   "ggnewscale", "scales", "ragg", "svglite", "rstatix", "nortest",
                   "lmerTest", "lme4", "emmeans", "car"))
BiocManager::install("ropls")          # Bioconductor
```

The first run computes the PLS-DA permutation and LOOCV and caches them in
`stats/plsda_*.rds`; later runs reuse the cache (~1.5 min) unless the input
changed. The seed is fixed at 42, so the same input reproduces the same values.

---

## What's in here

```
28CytokineAnalysisCode/
├── README.md              this file
├── VERIFICATION.md        what was checked before packaging, and what was found
├── code/
│   ├── run_all.R          runs the cascade in dependency order
│   ├── 00_setup.R         paths, palette, the 28-analyte gate, the sex z-score
│   ├── 01_sex_zscore_qc.R is the z-score safe here? balance + before/after
│   ├── 02_plsda.R         THE HUB — model, permutation, LOOCV, key cytokines
│   ├── 03_theme.R         shared figure language
│   ├── 04_panels.R        Fig 2 / Supp 4 panels
│   ├── 05_correlations.R  cytokine ~ amyloid / morphology
│   ├── 06_group_stats.R   flicker group comparisons
│   ├── 07_figure_helpers.R Fig 4 / Fig 5 panel language
│   ├── 08_mixed_model.R   cross-check: sex + litter in the model, raw values
│   ├── 09_corr_selection.R cytokine ~ morphology: key vs most flicker-responsive
│   ├── 10_fig2.R … 16_supp7.R   one stage per manuscript figure
│   ├── 17_data_folders.R  re-files the results by topic into data/ (runs last)
│   ├── lib/stats_helpers.R      the flicker group comparison
│   └── prep/              raw .xlsx -> data/processed/animal_table.csv
│
│   ── generated locally, not in this repository ──
├── data/
│   ├── raw/               the source plate workbook (README only here)
│   ├── processed/         per-animal tables (raw + sex z-scored)
│   ├── cytokine/          28-analyte values, panel gate, group stats, PLS-DA + plots
│   ├── amyloid/           Abeta / plaque values, group stats + bars
│   ├── morphology/        per-animal glial summaries (no group stats; see Known limits)
│   └── correlations/      rho tables, the BH q < 0.05 pairs, one scatter per pair + plots
├── stats/                 every statistic as a tidy CSV
└── figures/               manuscript-ready panels, 300 dpi PNG, vector PDF, editable SVG
```

**Every file in `code/` runs.** There is no dead code, no unused helper and no
superseded stage in this package: `run_all.R` executes stages 00–17, stages 06, 08,
11, 14 and 15 source `lib/stats_helpers.R`, and `prep/` builds the analysis
table from the raw workbook. Nothing else is carried along.

Round-1 and round-2 material — the Wood lab's original function library, the
earlier PLS-DA helpers, the linear-range analyte screen and their outputs —
lives in the project folder (`AshleyPrichard/cytokine/` and
`AshleyPrichard/archive/`) and is deliberately not duplicated here.

---

## The cascade

The analysis flows strictly downstream. Each stage's output is the next stage's
input, and **no script contains a literal list of cytokine names.**

```
  raw MFI ──┬─► detectability gate ──► PANEL (28 analytes)      00_setup.R
            └─► z-score within sex ──► ANIMALS (the one table)  00_setup.R
                          │
                          ▼
              PLS-DA per genotype ──► KEY CYTOKINES             02_plsda.R
                          │            (top 3 |loading| LV1 + LV2)
                          │
       ┌──────────────────┼──────────────────┐
       ▼                  ▼                  ▼
  correlations      group stats         every figure           05, 06, 10-16
   (05)              (06)                selects on KEY
```

Change the gate and all 28-analyte claims move together. Change the PLS-DA and
every downstream figure's analyte selection moves with it.

---

## The two decisions everything rests on

### 1. The panel: 28 analytes, >80% above background

An analyte is kept if **more than 80% of animals** have background-subtracted
MFI > 0. Keeps 28; drops **IL-5, IL-12(p40), MCP-1, TNFa**.

One line in Methods, literature precedent (Olink-style "% above LOD"), applied
to all 32 without exception. It replaces the Wood-lab curated list plus the
ad-hoc "no-change" screen.

**There is no LOD or LLOQ.** The plate was run without standards — 43 unknowns
and two Background0 wells, `AnalysisType = None`, blank units. Values are MFI,
not pg/mL, and "background" is two blank wells. This must be stated in Methods;
a real LOD would need a re-run with standards.

> **Margin worth knowing.** IL-5 is the closest call: 34 of 43 animals above
> background = 0.791, against a 0.80 threshold. One more animal (35/43 = 0.814)
> and it would be in the panel. This matters because Supplementary Figure 4's
> published caption is specifically about IL-5 — see "What the manuscript still
> says" below.

### 2. Sex: z-score within sex, then pool (PI's call)

Per sex cell there are only **n = 3–4** animals. The sexes cannot be analysed
separately, and pooling them raw would read sex differences as flicker effects.
So every measure — cytokines, amyloid **and** morphology — is z-scored within
sex within genotype, then M and F are pooled into one group of n = 20 (WT) /
23 (5xFAD).

`01_sex_zscore_qc.R` checks this and refuses to continue if it stops holding:

- **Balance** — the z-score is only valid if sex is balanced across flicker
  arms, otherwise subtracting a sex mean removes flicker signal too. Worst arm
  here is **57/43** (hard limit 80/20). Safe.
- **What it removed** — analytes with a sex effect at p < 0.05: **13 → 0**.
  The strongest were 5xFAD IL-10 (p = 7.5e-05), RANTES (1.3e-04), WT IL-10
  (0.0022).
- **What it kept** — flicker p-values before vs after correlate at
  **Spearman 0.915**; 7 → 5 analytes at p < 0.05.

Downstream values are unitless z-scores, so **figure axes read "(z-score)" and
bars can go below zero**. That is the visible price of pooling the sexes.

---

## Results

### PLS-DA — descriptive only, both genotypes

| | LOOCV | chance | permutation p |
|---|---|---|---|
| WT | 0.400 | 0.333 | **0.202** |
| 5xFAD | 0.348 | 0.333 | **0.320** |

Neither model beats chance. The scores and loadings are shown as a picture and
a descriptive weighting; **no statistics are run on the latent variables** and
the score violins carry no p-values, brackets or asterisks (the Wood lab asked
for this — `initialComments.md:8`, "cannot do stats on linked variable").

The z-score did improve WT LOOCV from 0.250 to 0.400, but not enough.

### Key cytokines (top 3 |loading| per latent variable)

| | LV1 | LV2 |
|---|---|---|
| **WT** | IL-15, IL-6, LIF | IL-2, IL-17, LIX |
| **5xFAD** | KC, IP-10, IL-17 | LIF, M-CSF, IL-12(p70) |

These drive analyte selection in Fig 4C, Fig 5C/D, Supp 5 and Supp 7.

### Flicker group effects — WT responds, 5xFAD does not

Omnibus p < 0.05 on **5 of 28** analytes in WT, **0 of 28** in 5xFAD.

| Analyte | Contrast | Adjusted p |
|---|---|---|
| IL-15 | NoFlicker vs 20Hz | 0.0034 |
| IL-6 | NoFlicker vs 20Hz | 0.0086 |
| M-CSF | NoFlicker vs 20Hz | 0.0211 |
| Eotaxin | NoFlicker vs 20Hz | 0.0325 |
| IL-15 | NoFlicker vs 40Hz | 0.0363 |
| IFNy | NoFlicker vs 20Hz, vs 40Hz | 0.0381 |

Two of these (IL-15, IL-6) are also WT LV1 key cytokines — the PLS-DA and the
group test agree on where the WT signal is.

*Caveat:* 5xFAD MIG shows a Tukey contrast at p = 0.043 with a non-significant
omnibus (p = 0.051). Reported for completeness; it should not be claimed.

### Amyloid — Abeta 40 is significant again

Kruskal-Wallis **p = 0.039**; 40Hz vs No Flicker **p = 0.048**, 40Hz vs 20Hz
**p = 0.048**. Abeta 42 is null (ANOVA p = 0.534).

Round 2 had Abeta 40 at p = 0.060 (n.s.). The sex z-score recovered it. The
manuscript's "40Hz reduces Abeta 40" claim is back in play.

### Correlations — the microglia result is the solid one

Nine pairs survive BH across the 28-analyte family; **five also survive
Westfall–Young** (FWER, the strict cross-check):

| Pair | rho | q (BH) | p (WY) |
|---|---|---|---|
| GM-CSF ~ Abeta 42 | 0.65 | 0.030 | **0.024** |
| IL-3 ~ microglia number | 0.65 | 0.012 | **0.026** |
| IL-15 ~ microglia number | 0.63 | 0.012 | **0.036** |
| Eotaxin ~ microglia number | 0.62 | 0.012 | **0.040** |
| IL-10 ~ microglia number | 0.62 | 0.012 | **0.041** |

All 5xFAD. IL-7, MIG, MIP-1b and MIP-2 also clear BH but fail WY.

**The astrocyte correlations are gone.** Round 2 reported astrocyte density ~
IL-10 (rho = −0.61) and ~ RANTES (−0.56). Both vanish under the sex z-score —
they were the sex confound, exactly as suspected. That disappearance is
evidence the correction is doing its job, and it is a result worth stating.

GM-CSF ~ Abeta 42 now passes both corrections (it passed BH but failed WY at
p = 0.25 before the z-score). Read it with its resolution columns: **5 distinct
values, 33% ties**.

### Correlations by selection rule (stage 09, 2026-11)

Stage 05 correlates all 28 analytes with every target, so choosing a different
analyte set changes no rho or p, only the multiplicity family. Stage 09 puts the
PLS-DA key cytokines next to the six most flicker-responsive analytes (smallest
omnibus p) on the four figure targets. The responsive set is not stronger: WT
has 0 nominal hits (key: 2), 5xFAD has 4 (key: 3), and one pair (5xFAD MIG ~
microglia number) clears BH across the 28 but fails Westfall-Young (p = 0.18).
The pairs that pass both corrections (IL-3, IL-15, Eotaxin, IL-10 ~ microglia
number) are in neither set. See `stats/corr_morphology_by_selection.csv`.

---

## What the manuscript still says (not yet reconciled — PI's call)

`manuscript/7.21.26ManuscriptWilliam.docx` was written against round 2:

- Methods says **16 analytes**; the panel is **28**.
- G-CSF, GM-CSF, M-CSF and LIF are written up as *excluded* for linear range.
  All four are back in, and LIF and M-CSF are 5xFAD key cytokines.
- **Supp Fig 9** exists only to justify excluding those four. Its premise is
  void; it should be cut.
- Abeta 40 is written as n.s. (p = 0.060). It is now p = 0.039.
- Supp 4's caption is "Flicker reduces expression of IL5 in 5xFAD mice". IL-5
  is one of the four analytes the gate drops, so it is not in the panel and
  cannot be. Either the caption changes or the gate has to be argued around.

---

## Known limits

- **Morphology group statistics (Fig 1/3/6, Supp 2/3/8) are not reproducible**
  from the per-cell files in the original `data/SARA Temp/` — those figures were
  made in MATLAB from an export not in that folder. They are **not part of this
  package**, and the correlation targets inherit that blocker. Needs canonical
  source data.
- **Permutation p-values carry Monte Carlo noise.** Correlations use B = 50,000
  (a p near 0.04 is stable to about ±0.002); the PLS-DA permutation uses
  B = 1,000, which is coarser. The seed is fixed (42), so a given run reproduces.
- **Coarse analytes.** Every correlation row carries `n_distinct` and
  `prop_tied`, measured on raw MFI. An analyte read at 3–5 levels can throw a
  large rho from quantisation alone; that is what the BH/WY gap detects.
  Report `p_adj` and `p_wy` together.
- **Supp 6 does not use the key cytokines** — it would be a pixel duplicate of
  Supp 5. It ranks by flicker responsiveness instead. See the header of
  `15_supp6.R`; this is flagged for a decision.

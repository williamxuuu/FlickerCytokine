# PLS-DA

Partial least squares discriminant analysis of the 16 retained cytokines against
flicker group (NoFlicker / 20Hz / 40Hz), fit separately per genotype.

- `code/` — `04_plsda.R` (fit + LOOCV + permutation), `04b_key_cytokines.R`
  (picks the focus cytokines from the loadings), `plsda_helpers.R`.
- `stats/plsda_loadings.csv` — LV1/LV2 loadings with LOOCV standard deviations.
- `stats/key_cytokines.csv` — top 3 per latent variable per genotype.
- `stats/plsda_fits.rds` — the saved models. Every downstream figure reads this,
  so figures cannot drift from the reported statistics.
- `figures/` — scores and loadings plots (PDF vector + manuscript PNG panels).

**Key cytokines (top 3 on LV1, top 3 on LV2):**
- WT: IP-10, IL-10, IL-1B / VEGF, MIP-1a, KC
- 5xFAD: KC, IL-1a, IL-17 / MIP-1a, IL-10, IP-10

**Important caveat:** both models fail permutation testing (see
`../Permutation Analysis/`). The loadings are used descriptively — to choose which
cytokines to display — and no inferential claim is made from the LV scores. The
earlier LV-score p-values (0.0031 etc.) were invalid and have been removed.

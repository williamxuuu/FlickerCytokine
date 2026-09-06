# Corrections traceability — every change tied back to initialComments.md

The point of this file: confirm that the edits map to what the Wood lab actually
asked for in [../initialComments.md](../initialComments.md), and flag the few
that go beyond it so they're a conscious choice, not a silent addition.

## Every requested change → the comment that asked for it

| initialComments.md | What it asked | Where it's implemented |
|---|---|---|
| **L1** "figures … reporting uncorrected P values … FDR … too strict" | the overall motivation | everything below |
| **L3** "the code re run (can be altered to be Matlab or Python if easier)" | re-run allowed in another language | whole pipeline re-run in R; Methods note "in R (cytokine…)" |
| **L6** "Prune out targeted based on linear range criteria" | drop out-of-range analytes | stage 02 (17-analyte linear panel); Methods; Supp 6 flagged (5/6 out of range) |
| **L7** "ones that are significant vs pulling out myeloid related ones" | focus analytes | myeloid panel (`cytokines_myeloid`); Supp 5 / Fig 4C rows |
| **L8** "PLSDA – cannot do stats on linked variable … remove p values on the violin plots by heat map" | drop LV-score stats | **Fig 2C & Supp 4B: the LV1/LV2 score-violin p-value panels removed**; stage 04 |
| **L9** "Permutation analysis is statistically sound following PLSDA" | validate PLS-DA by permutation | stage 04 permutation (R²Y p, Q²); shown in Fig 2C / Supp 4B captions |
| **L10** "batch adjustment for sex after the linear range" | adjust for sex | stage 02 (sex-adjust within genotype) |
| **L11** "Are there sex-specific differences? Can show batched heat map" | sex differences + batched heatmap | sex-batched heatmap (stage 03); mixed model reports sex effect |
| **L12–13** "3 flicker groups … ANOVA + Tukey … if not parametric KW + Dunn's" | correct the 3-group tests | stage 05; **Supp 5 myeloid p's**, **Fig 4B amyloid p's** |
| **L14** "Alternatively: Multiple effects model … to control for sex?" | sex-controlled model | stage 05b `Flicker + Sex + (1\|litter)` |
| **L15, L25** "statistical analysis is Spearman, the correlation line is Pearson" | Spearman stat, Pearson line | stage 07; Fig 4C, Fig 5; Supp 7 "Pearson"→"Spearman" |
| **L16–17, L24** "Use the p values" / "p value is the driver" | p drives inclusion, not a ρ cutoff | stage 07 (`corr_p_threshold`, no ρ gate) |
| **L18–21** "FDR across everything … too much … by row (32) … or Myeloid" | narrow the FDR | stage 07 (`by_row` / `myeloid`); Fig 4C asterisks |
| **L22–25** "Spearman … analyte and a morphology parameter" | cytokine↔morphology correlations | stage 07; Fig 5 |
| **L26** "Labeled Males and Females on scatter plots" | label M/F on scatters | Fig 5 scatters (M/F by shape) |

## Changes that go BEYOND the letter of initialComments.md — your call

These follow from the requests but weren't spelled out. Keep or revert:

1. **Amyloid (Aβ40/42) bars corrected (Fig 4B).** initialComments.md talks about
   the *cytokine* panel; amyloid isn't named. I applied the same KW+Dunn's
   correction (L12–13) because it's the identical 3-flicker-group comparison and
   the same "uncorrected p-values" concern (L1). Consequence: Aβ40 40Hz-vs-No
   Flicker goes 0.024 → 0.060 (n.s.). If you'd rather leave amyloid on its
   original test, say so.

2. **Males/females added to the BAR plots** (Supp 5 myeloid, Fig 4B amyloid).
   L26 asks for M/F only on *scatter plots*. I added M/F shapes to the bars too,
   for consistency. Easy to drop if you want the bars plain.

3. **Wording softened to match the new stats** ("significantly reduced" →
   "reduced"; "greater expression" → "trend toward"; "is more effective" → "may
   be more effective"). Not a separate request — it's required so the prose
   doesn't contradict the now-non-significant numbers (downstream of L1).

4. **"~3 significant figures"** on the reported values. This came from a
   co-author comment inside the manuscript (A. Singer), **not** from
   initialComments.md. Harmless but noting the source.

5. **"Addison-Darling" → "Anderson-Darling"** typo fix. Not requested; just the
   correct name of the test.


## Revised spec (2026-08) — the second round

Ashley's revised note re-orders the analysis and changes how the focus analytes
are chosen. Each line below is quoted from that note.

| Revised-spec line | What it asked | Where it's implemented |
|---|---|---|
| "revising the heatmaps to batch correct for females and males" | sex batch correction, applied before the heatmaps | stage 02 (now adjusts all 32 analytes, then prunes — the two steps commute) |
| "remove cytokines based on the linear range criteria (suspicious ones that show no change)" | a real removal, not a flag | stage 02: `keep_myeloid_regardless = FALSE`; the loading curve is re-fit and a floor/no-change screen added; 16 analytes retained |
| "then PLSDA without stats (cannot do stats on a linked variable)" | unchanged from round 1 | stage 04 — no test on LV scores |
| "then permutation analysis as in Kristie's studies" | Garza et al. 2020 (J Neurosci 40:1211) / Sci Adv 2023 eadf5672 | stage 04: LOOCV loadings mean ± SD over 1000 iterations, and 1000 label permutations scored by LOOCV classification accuracy (`p_permute`) |
| "FDR correct based on pre-selected cytokines based on the heatmap - not all 32" | narrow the FDR family to the heatmap panel | stage 07: `corr_fdr_scope = "by_row"` over the 16 retained analytes, per target |
| "Identify top weighted cytokines based on LV1 (top 3) and on LV2 (top 3)" | data-driven focus set, replacing the myeloid panel | **stage 04b** (new): top 3 by \|LOOCV mean loading\| per LV, per genotype |
| "comparing between groups for the key cytokines by using an ANOVA (tukey correction) or KW (dunn correction) OR LMM" | unchanged tests, new analyte set | stages 05 / 05b, now driven by the key set |
| "Use p values to determine significance for Spearman correlations" / "in biology .46 is okay … but the p values are the driver" | no ρ gate | stage 07 (`corr_p_threshold`, no ρ cutoff) |
| "males and females need to be marked on the scatterplots (different symbols)" | unchanged from round 1 | stage 10 (M = circle, F = triangle) |

### Decisions taken on this round (William + Ashley)

1. **Linear range: re-derive, then prune strictly.** The Wood lab's curated
   17-analyte panel is kept as the anchor (it matches their own "Only Linear
   Cytokines" folder), with two derived evidence columns written alongside it.
   The floor/no-change screen removes one further analyte (**MIP-1b**: median 2.0,
   8 distinct values across 43 animals, range 0–4). Final panel = **16 analytes**.
2. **Correlation FDR family = the heatmap panel, per target** (not the key set,
   not all 32). A key-set-restricted table is written as a sensitivity analysis.
3. **Key cytokines are chosen per genotype**, each from its own model.
4. **Tukey / Dunn is the only multiplicity correction** applied to the group
   comparisons. `p_omnibus_fdr` remains in the CSVs as a diagnostic column and is
   not cited anywhere.

### Beyond the letter of the revised spec — flagging these

- **The no-change screen is mine, not the Wood lab's.** The spec says "suspicious
  ones that show no change" without giving a rule; I implemented it as
  `n_distinct <= 10 OR range <= 5 OR >= 20% of animals at zero` (thresholds in
  `config.R`). Set `apply_nochange_screen <- FALSE` to prune on the curated list
  alone, which puts MIP-1b back.
- **Selection is by |loading|, not signed loading.** An LV separates groups in
  both directions, so a strongly negative loading is as informative as a positive
  one. Selecting the top 3 *positive* loadings instead would change the 5xFAD LV2
  set (IP-10 has a negative loading there).
- **The stability column.** `key_cytokines.csv` reports whether each pick's LOOCV
  error bar clears zero. Not requested; it costs nothing and makes an unstable
  pick visible.

## Deliberately NOT changed (initialComments.md didn't ask, so untouched)

- **Morphology figures** — Fig 1, 3, 6, Supp 2, 3, 8. initialComments.md covers
  cytokines and cytokine↔morphology *correlations* only; the morphology group
  comparisons already use KW+Dunn's. Left as published.
- Stage 08 re-derives per-cell morphology stats for our own checking only; it is
  **not** applied to any manuscript figure.

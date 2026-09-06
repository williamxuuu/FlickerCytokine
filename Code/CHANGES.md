# Corrected cytokine / morphology re-analysis — what changed and why

This folder is a clean re-implementation of the Wood-lab flicker analysis
(`original/AshleyDataAnalysis_SBversion_Final.R`, S. Bitarafan 2025) with the
statistical corrections the Wood lab asked for in
[../initialComments.md](../initialComments.md). The original 3,652-line script
is preserved untouched in [original/](original/); nothing here overwrites it.

The short version: the original reported **uncorrected p-values** on the figures
(FDR existed only inside a few CSVs), ran **stats on PLS-DA latent variables**,
used **pairwise t-tests / bare ANOVA with no post-hoc**, and **never actually
applied the linear-range filter** it QC'd. Fixing those changes several of the
manuscript's headline numbers. The details are below.

## How it's organized

The original is ~6 copy-pasted blocks (WT / WT-M / WT-F / 5xFAD / 5xFAD-M …) of
the same recipe. I collapsed that into config-driven stages that loop over
subgroups. Every statistical change is tagged `# [CORRECTION]` at the point it
happens.

| Stage | File | What it does |
|---|---|---|
| config | `config.R` | analyte lists, linear-range panel, colors, every stat knob |
| helpers | `lib/stats_helpers.R` | the corrected stats (AD→ANOVA/Tukey vs KW/Dunn, Spearman+FDR, sex-labeled scatter) |
| helpers | `lib/plsda_helpers.R` | PLS-DA matrix, deterministic axis orientation, LOOCV loadings, permutation test |
| 01 | `scripts/01_load_and_prep.R` | read combined table, background-subtract cytokines, factors |
| 02 | `scripts/02_linear_range_and_batch.R` | sex-adjust within genotype, then prune to the linear panel (re-derived + no-change screen) |
| 03 | `scripts/03_cytokine_heatmaps.R` | heatmap3 per genotype (Fig 2/3) |
| 04 | `scripts/04_plsda.R` | PLS-DA + Garza-style permutation + LOOCV loadings, no stats on scores (Fig 2) |
| 04b | `scripts/04b_key_cytokines.R` | key cytokines = top 3 on LV1 + top 3 on LV2, per genotype |
| 05 | `scripts/05_cytokine_group_stats.R` | per-analyte flicker tests (Supp 5/6) |
| 05b | `scripts/05b_cytokine_mixed_effects.R` | sex-controlled mixed/multiple-effects model |
| 06 | `scripts/06_amyloid.R` | Aβ40/42 bars (Fig 4B) |
| 07 | `scripts/07_correlations.R` | cytokine~amyloid + cytokine~morphology (Fig 4/5) |
| 08 | `scripts/08_morphology_stats.R` | per-cell morphology (Fig 1/3 — see caveat) |

Every test also writes a tidy row to `results/stats/` (test, statistic, raw p,
adjusted p, N, correction method) so the Supplementary Statistics table can be
regenerated instead of hand-copied.

## What each correction changed

| Wood-lab note (initialComments.md) | Original behavior | Now |
|---|---|---|
| Prune analytes by linear range | never applied in code (visual QC only) | keep the 17-analyte linear panel; log kept/dropped (`02`) |
| Focus on significant + myeloid analytes | all 32 shown | myeloid panel (`cytokines_myeloid`) called out for the focused stats |
| No stats on PLS-DA latent variables | t-test/ANOVA on LV1/LV2 scores | removed; PLS-DA gets a **permutation** test instead (`04`) |
| Remove p-values on the violins-by-heatmap | raw p drawn | those annotations dropped; corrected p only on the dedicated per-analyte figures (`05`) |
| Batch-adjust for sex | sex handled by running separate subgroups | additive sex shift removed within genotype for pooled figures (`02`) |
| "Multiple/mixed-effects model to control for sex" | not done | `y ~ Flicker + Sex + (1\|litter)` (lmer), fixed-effects fallback on singular fit; emmeans Tukey; reports sex + Flicker×Sex (`05b`) |
| "Show a sex-batched heatmap" | not done | Sex×Flicker group-mean heatmap per genotype (`03`) |
| 3 flicker groups, N≈4–8 | pairwise t-test / bare ANOVA, no post-hoc | Anderson-Darling → ANOVA+Tukey or Kruskal-Wallis+Dunn's (`05`,`06`,`08`) |
| Correlations: Spearman stat, Pearson line | mixed Pearson/Spearman, raw p on plot | Spearman ρ+p for the stat, linear line for display only (`07`) |
| p drives inclusion, not a ρ cutoff | — | `corr_p_threshold`, no \|ρ\| gate (`07`) |
| FDR "across everything" is too much | FDR across the whole matrix (in CSV only) | FDR **by row** (per target) or **restricted to the myeloid panel** (`07`) |
| Label males/females on scatterplots | not labeled | sex = point shape + color on every scatter (`07`) |
| ~3 significant figures throughout | up to 15 digits | `fmt()` / `sig_digits = 3` everywhere |

## What the corrections do to the results

This is the part that matters for the manuscript. Numbers are from the current
run; they move a little with the exact linear-range / FDR choices in `config.R`.

**WT myeloid cytokines (Supp 5).** Proper KW+Dunn's / ANOVA+Tukey collapses most
of the WT myeloid story:

| Comparison | Manuscript (uncorrected) | Corrected |
|---|---|---|
| M-CSF, 20Hz vs NoFlicker | p = 0.008 | **p = 0.014 — holds** |
| G-CSF vs NoFlicker | p = 0.028 / 0.05 | p = 0.058 — lost |
| LIX, 20Hz vs 40Hz | p = 0.045 | p = 0.145 — lost |
| GM-CSF, 20Hz vs 40Hz | p = 0.047 | p = 0.130 — lost |
| any analyte, FDR across the panel | (not done) | nothing survives |

So the "20Hz raises myeloid cytokines in WT" claim rests almost entirely on
**M-CSF** once the pairwise tests are corrected. The sex-controlled mixed model
(`05b`) agrees: M-CSF is the only WT analyte with a flicker effect after
controlling for sex (F = 5.89, p = 0.012; 20Hz vs No Flicker p = 0.0095). That
model also flags genuine **sex main effects** (IL-10, RANTES in 5xFAD) and
**Flicker×Sex interactions** (IL-9, MIP-1a/b in WT — effects that only appear
when the sexes are split), which is why some comparisons are reported by sex.

**PLS-DA (Fig 2C).** The original t-tested the LV scores and reported group
p-values (e.g. NoFlicker vs 20Hz p = 0.0031). Those aren't a valid group test.
The permutation test on the model itself is **not significant** (WT: Q² = −0.53,
p(R²Y) = 0.37; 5xFAD: Q² = −0.09, p = 0.39). The PLS-DA separates the plot but
does not pass permutation — it should be presented as visualization only.

**Amyloid (Fig 4B).** Aβ40 40Hz vs NoFlicker goes from p = 0.024 to **p = 0.060**
under KW+Dunn's — it no longer clears 0.05.

**Cytokine ↔ amyloid (Fig 4C).** The Spearman ρ's reproduce the manuscript
(GM-CSF~Aβ42 ρ = 0.67 vs 0.62; LIX~Aβ42 ρ = 0.57 vs 0.54). What changes is
multiplicity:
- **By-row FDR** (32 analytes/target): only GM-CSF~Aβ42 and LIX~Aβ42 survive.
- **Myeloid-restricted FDR** (the Wood lab's suggested route): GM-CSF, LIX, MIG
  and G-CSF ~Aβ42 all survive (p_adj 0.003–0.029).

The myeloid~amyloid finding is defensible **if** the correlations are restricted
to the pre-selected myeloid panel before correcting — which is exactly what the
Wood lab recommended. Reported blanket across all 32, it mostly doesn't hold.

## Decisions that need Ashley's call

1. **Linear range vs myeloid focus.** The curated linear-range set (17 analytes)
   drops G-CSF, GM-CSF, M-CSF and LIF — four of the manuscript's headline myeloid
   cytokines. `config.R` currently re-adds the myeloid panel
   (`keep_myeloid_regardless = TRUE`) so the story is testable; set it `FALSE`
   for a strict linear-range prune. Needs a call on whether the myeloid analytes
   are reported despite being outside linear range.

2. **FDR scope for correlations — DECIDED: myeloid-restricted.** `config.R` is set
   to `corr_fdr_scope = "myeloid"`, following the Wood lab's explicit note that
   FDR "across everything… is probably too much" and to "limit to pre-selected /
   myeloid" analytes (initialComments.md L19/L21). This keeps GM-CSF/LIX/MIG/G-CSF
   ~Aβ42 significant. State the scope in Methods. (`by_row` remains a one-flag
   fallback if a reviewer wants it stricter.)

3. **Morphology source data (blocker).** Stage 08 re-derives the per-cell
   morphology stats, but the per-cell files in `SARA Temp/` **do not reproduce**
   the manuscript's morphology figures. Example: WT astrocyte tiling in the
   provided file has medians 84.7 / 87.8 / 93.3 µm across NoFlicker/20Hz/40Hz
   (essentially no effect, 40Hz slightly *more* dispersed), whereas Fig 1D reports
   H(2) = 102.9, p = 4.6e-23 in the opposite direction. Those figures were made in
   MATLAB from an export that isn't in this folder. The morphology group-stats
   can't be corrected until the canonical source data is provided — the cytokine /
   amyloid / correlation results above are unaffected by this.

## Manuscript revisions

Edits are on a **copy**, `../ManuscriptMostRecent_REVISED_stats.docx` (the
original is untouched). Three kinds of change, all tracked:

- **Statistics substituted in the prose** (31 edits): the changed p-values are
  now in the text (e.g. Aβ40 0.024 → 0.060; the WT myeloid p's; correlation p's
  swapped to FDR), and the wording softened to match (e.g. "significantly
  reduced" → "reduced"; "greater expression" → "trend toward"). The PLS-DA
  LV-score p-values in the Fig 2 legend are replaced by the permutation result.
  The Methods stats paragraph is rewritten to describe the corrected approach.
- **Margin comments** (7, author "W. Xu (re-analysis)"), each giving the old
  stat, the new stat, and the reasoning / which stage produced it.
- **Corrected figures embedded** as a labeled addendum at the end (8 PNGs,
  rendered to match the manuscript's own style — flicker-colored bar+SEM+dot
  plots, PLS-DA ellipse scores, and the Fig 4C ρ-heatmap grid: revised Fig 2C
  PLS-DA, Fig 4B amyloid bars, Fig 4C cytokine~amyloid ρ-heatmap, Supp 5 myeloid
  panels, and the two new sex-batched heatmaps). These are DATA panels:
  the published figures are composites that also carry microscopy / schematic
  panels we can't regenerate, so they are provided for the figure designer to
  swap into the composites rather than overwriting them. Source PNGs are in
  `results/revised_figures_png/`.

Figures were subsequently rebuilt/extended and re-styled to match the originals:
Fig 2 (PLS-DA scores + Signals-in-LV loadings, oriented like the original),
Fig 4B/4C, **Fig 5** (5xFAD correlations, Sara's scatter style + M/F shapes),
Supp 4 (PLS-DA + relabeled C/D loadings), Supp 5, **Supp 7** (WT correlations,
Spearman + M/F), and the sex-batched heatmaps. All in `results/revised_figures_png/`.

## Final multi-agent review pass (working file `7.21.26ManuscriptWilliam.docx`)

A 7-dimension agent swarm cross-checked the finalized draft (leftover/uncorrected
stats, internal consistency, stats-vs-analysis, figure/panel refs, scientific
overclaims, phrasing, and my comments), each finding adversarially verified.
**64 issues confirmed.** The headline: the *narrative* (abstract, highlights,
intro, discussion, some figure titles) still told the pre-correction story while
the corrected stats no longer support it.

**Fixed in the working file (unambiguous):**
- Section heading "amyloid I-40" → "1-42" (the surviving cytokine↔amyloid
  correlations are with Aβ1-42, not Aβ1-40).
- Fig 4C sentence: "LIF" → "MIG" in the FDR-surviving list; "correlated with
  cytokine modulation induced by flicker" → "…with Aβ1-42 burden".
- "significantly correlated… in 5xFAD mice (Supp Fig 7b)" → "…in WT mice" (Supp 7
  is the WT/healthy figure).
- Supp 6 IFN-γ raw p's 0.018/0.041 → corrected 0.044 (n.s.).
- Typos: "ytokine"×2, "cyotkine", "make of up", "include Alzheimer" →
  "including", "we show here… here".
- My 9 margin comments rewritten in William's voice (first-person, plain).

**Flagged (comments, not rewritten — narrative/scientific calls for the authors):**
abstract & highlights myeloid claim (only M-CSF survives; no significant
20Hz-vs-40Hz difference), Fig 2 title "distinct clustering" (PLS-DA fails
permutation), amyloid "may be more effective than 20Hz" (Aβ40 n.s.), Results
"more pronounced… GM-CSF/G-CSF/LIX/LIF" (all n.s.; LIF out of range), 5xFAD
"robust proinflammatory response" (no significant flicker effect in 5xFAD), and
the mismatched Supp Fig 4 title.

**Reported separately (not touched):** grammar errors and a broken statistic
(F(2)=292, "p = 5.4115"; "p = < .001") inside the **morphology** supplements —
those live in the MATLAB-derived morphology figures, which are out of this
re-analysis's scope and still need the canonical source data.

## Running it

```bash
cd "/Users/williamxu/Documents/Research/BALES/AshleyPrichard/cytokine_analysis"
Rscript run_pipeline.R
```

Outputs land in `results/` (`cytokine/`, `amyloid/`, `correlations/`,
`morphology/`, `stats/`). See `README.md` for the layout and sanity checks.

---

# Round 2 — the revised spec (2026-08)

Ashley's revised note changed two things about the analysis and re-ordered the
rest. Everything already implemented for round 1 (sex batch correction, no stats
on the LVs, ANOVA/Tukey vs KW/Dunn, the mixed model, Spearman with p-driven
inclusion, M/F symbols) was kept unchanged and re-pointed at the new analyte set.

## The two substantive changes

**1. The linear-range prune is now a real removal.** Round 1 kept the four
out-of-range myeloid analytes (`keep_myeloid_regardless = TRUE`) so the
manuscript's myeloid story stayed testable. The revised spec asks to "remove
cytokines based on the linear range criteria — suspicious ones that show no
change", and no longer asks for a myeloid focus, so the flag is now `FALSE`.

Stage 02 also re-derives the evidence instead of asserting the panel. Two
independent columns are written to `results/stats/analyte_linear_range.csv`:

- **Loading-curve fit.** The BCA linear-range plate is an 8-point lysate loading
  series (0–8 µL) per genotype plus two background wells. Background-subtract,
  clip at zero, fit MFI ~ loading. This is the data `original/LR Code.R` plotted
  but never scored — it computes R² and p and leaves both commented out. With
  only 8 points the WT curve is noisy, so it is reported as supporting evidence,
  not used to decide; the script prints every analyte where it disagrees with the
  curated call.
- **Floor / no-change screen** on the 43 animals: an analyte sitting on the assay
  floor has few distinct values, a narrow range, and/or many readings clipped to
  zero. Rule (thresholds in `config.R`): `n_distinct <= 10` OR `range <= 5` OR
  `>= 20%` at zero.

The decision stays anchored on the Wood lab's curated 17-analyte panel — it
matches their own `DotBox plots/Zscored/All/Only Linear Cytokines/` folder
exactly, so it is their call, not an inference. The no-change screen removes one
further analyte. **Final panel: 16 analytes.**

| Dropped this round | median | max | distinct values | curated linear |
|---|---|---|---|---|
| G-CSF | 1.0 | 2.5 | 5 | no |
| GM-CSF | 1.5 | 3.0 | 6 | no |
| M-CSF | 1.5 | 3.5 | 6 | no |
| LIF | 2.0 | 5.0 | 9 | no |
| MIP-1b | 2.0 | 4.0 | 8 | **yes** — removed by the no-change screen |

For scale, a retained analyte like IL-2 has a median of 324 and 42 distinct
values across the same 43 animals.

**2. Key cytokines now come from the PLS-DA loadings, not the myeloid panel.**
New stage `04b`: top 3 by |LOOCV mean loading| on LV1 and top 3 on LV2, per
genotype, from that genotype's own model. Selection is by magnitude, because an
LV separates groups in both directions.

| Genotype | LV1 top 3 | LV2 top 3 |
|---|---|---|
| WT | IP-10 (0.451), IL-10 (0.434), IL-1B (0.419) | VEGF (0.425), MIP-1a (0.373), KC (0.369) |
| 5xFAD | KC (0.317), IL-1a (0.316), IL-17 (0.306) | MIP-1a (0.693), IL-10 (0.374), IP-10 (−0.310) |

All 12 picks are stable in the sense that their LOOCV error bar clears zero.

## Permutation test — now the Garza ("Kristie") form

Round 1 used ropls' built-in permutation, which shuffles Y and compares R²Y/Q².
Garza et al. 2020 (J Neurosci 40:1211) and Sci Adv 2023 eadf5672 do something
different and more interpretable: leave-one-out cross-validation with the
loadings reported as mean ± SD over 1000 iterations, and 1000 refits on randomly
reassigned groups, with `p_permute` the fraction of shuffled models matching or
beating the true one.

Both are now reported. `p_permute` is primary; the ropls R²Y/Q² permutation stays
as a secondary column.

| | LOOCV accuracy | chance | permuted mean | **p_permute** | R²Y | Q² | ropls p(R²Y) |
|---|---|---|---|---|---|---|---|
| WT | 0.250 | 0.333 | 0.261 | **0.585** | 0.266 | −0.517 | 0.428 |
| 5xFAD | 0.304 | 0.333 | 0.244 | **0.376** | 0.237 | −0.202 | 0.440 |

**Neither model classifies better than chance, and neither beats its permuted
null.** Both PLS-DA models are visualization plus a descriptive weighting. The
loadings charts now carry LOOCV error bars, the top-3 picks are highlighted, and
`key_cytokines.csv` carries `p_permute` and a `model_valid` flag on every row so
the table cannot be read without it.

## What round 2 does to the results

Every change below is attributable to the five dropped analytes — the amyloid
group stats, which do not depend on the cytokine panel, are byte-identical.

**Cytokine group comparisons.** Nothing among the six key cytokines differs
between flicker groups, in either genotype, by either route (ANOVA+Tukey /
KW+Dunn, or the sex-controlled mixed model). Across all 16 retained analytes only
two contrasts clear 0.05, and neither is a key cytokine:

| Contrast | Round 1 | Round 2 |
|---|---|---|
| WT IFNy, NoFlicker vs 20Hz (Dunn-BH) | 0.044 | **0.044 — holds** |
| WT IFNy, NoFlicker vs 40Hz (Dunn-BH) | 0.044 | **0.044 — holds** |
| WT female LIX, 20Hz vs 40Hz (Tukey) | 0.048 | **0.048 — holds** |
| WT M-CSF, NoFlicker vs 20Hz | 0.014 | analyte dropped (out of linear range) |
| WT male M-CSF, NoFlicker vs 20Hz | 0.041 | analyte dropped |
| WT male MIP-1b, NoFlicker vs 40Hz | 0.037 | analyte dropped |

M-CSF was round 1's single surviving myeloid result. It is gone — not because
the test changed, but because the analyte sits on the assay floor (median 1.5,
6 distinct values across 43 animals).

**Cytokine ↔ amyloid (Fig 4C).** Under the heatmap-panel FDR family (16
analytes per marker), **nothing survives**:

| Pair | ρ | raw p | q (16-analyte family) |
|---|---|---|---|
| LIX ~ Aβ42 | 0.571 | 0.0044 | 0.069 |
| MIG ~ Aβ42 | 0.534 | 0.0087 | 0.069 |
| MIG ~ Aβ40 | 0.507 | 0.014 | 0.203 |
| IL-10 ~ Aβ40 | 0.465 | 0.025 | 0.203 |
| IP-10 ~ Aβ40 | 0.421 | 0.045 | 0.243 |

Round 1 had GM-CSF~Aβ42 (q = 0.010) and LIX~Aβ42 (q = 0.046) surviving by row,
and four myeloid pairs surviving under the myeloid-restricted family. GM-CSF was
the smallest p in that row and carried the family; with it dropped, LIX and MIG
land just above threshold. Both are still substantial correlations (ρ ≈ 0.53–0.57)
with raw p < 0.01 — they are reportable as trends, but not as FDR-significant.

**Cytokine ↔ morphology (Fig 5).** This is where the signal is: 15 pairs survive
the panel FDR, all but one in 5xFAD.

| Pair | ρ | q |
|---|---|---|
| MG mean number of cells ~ MIG | 0.712 | 0.0022 |
| MG mean number of cells ~ IL-10 | 0.688 | 0.0023 |
| MG median number of cells ~ MIG | 0.688 | 0.0039 |
| MG mean tiling ~ MIG | −0.701 | 0.0032 |
| MG mean tiling ~ IL-10 | −0.665 | 0.0043 |
| MG median number of cells ~ RANTES | 0.652 | 0.0039 |
| WT median astrocyte volume ~ IL-2 | −0.644 | 0.044 |
| (plus 8 more: RANTES, IP-10, IL-10 with astro/MG density and tiling) | | |

Round 1 had 18; the seven lost all involved dropped analytes (M-CSF, LIF, MIP-1b,
GM-CSF). Four are new, surfacing because the smaller family makes the correction
less punishing (RANTES and IL-10 with astrocyte density, IL-2 with astrocyte
volume in WT).

**Sex effects.** The mixed model still finds real sex main effects — IL-10 in both
genotypes (WT p = 0.0011, 5xFAD p = 5.8e-05), RANTES (p = 4.6e-05), LIX and VEGF
in 5xFAD, IFNy in 5xFAD, IL-1B in WT — and Flicker×Sex interactions for WT IL-9
(p = 0.041) and WT MIP-1a (p = 0.023). This is what the sex batch correction and
the sex-batched heatmaps are for, and it is arguably the most robust cytokine
result in the dataset.

## Where this leaves the manuscript

Stated plainly, so the narrative rewrite has one source:

- **No cytokine differs between flicker groups** among the PLS-DA-selected key
  cytokines, in either genotype. The only surviving group differences in the whole
  retained panel are WT IFNy (lower in No Flicker than either flicker condition)
  and a WT-female LIX 20Hz-vs-40Hz difference.
- **Neither PLS-DA model passes permutation.** The scores plots and loadings are
  descriptive.
- **No cytokine↔amyloid correlation survives FDR.** LIX and MIG ~ Aβ42 are the
  closest (q = 0.069).
- **Cytokine↔microglial-morphology correlations are strong and survive FDR** in
  5xFAD — MIG, IL-10 and RANTES track microglial density and tiling with |ρ| of
  0.65–0.71. This is the result the cytokine figures can carry.
- **Sex is a real source of variance** in several analytes.

The abstract, highlights, and discussion still describe round 1's story (and in
places the pre-correction story). They need rewriting against the list above.
Per Annabelle's instruction, that is drafted and flagged for the authors, not
auto-applied.

## Round-2 config decisions

| Knob | Setting | Why |
|---|---|---|
| `keep_myeloid_regardless` | `FALSE` | revised spec: the prune is a real removal |
| `apply_nochange_screen` | `TRUE` | "suspicious ones that show no change"; set `FALSE` to keep MIP-1b |
| `corr_fdr_scope` | `"by_row"` | "FDR based on pre-selected cytokines based on the heatmap" |
| `key_n_per_lv` | `3` | "top 3" on each of LV1 and LV2 |
| `plsda_loocv_runs` / `plsda_n_perm` | `1000` / `1000` | matches Garza et al. |

## Still open after round 2

1. **Morphology source data** — unchanged blocker from round 1. Stage 08 cannot
   reproduce the manuscript's morphology figures from the per-cell files in
   `SARA Temp/`; the canonical MATLAB export is not in the folder.
2. **Narrative rewrite** — see "Where this leaves the manuscript" above.
3. **Figure swap** — the corrected panels in `results/revised_figures_png/` need
   pasting into the Canva composites. `Fig_myeloid_panel_*.png` has been deleted;
   the replacement is `Fig_key_cytokines_*.png`.
4. **Whether Supp 6 survives at all** — most of its analytes are now formally out
   of the analysis.

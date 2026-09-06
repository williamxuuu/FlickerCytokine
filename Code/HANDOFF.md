# Handoff — Ashley Prichard Flicker Manuscript: cytokine/morphology re-analysis

**Purpose of this doc:** everything needed to continue this project in a fresh
chat, self-contained (don't rely on prior chat memory). Written for William + a
new AI assistant.

---

## 1. What this project is

Revising the statistics + figures for a manuscript before submission:
**"Noninvasive audiovisual neurostimulation differentially alters glial responses
at 20Hz and 40Hz"** (Prichard / Singer lab / Wood lab). Chronic flicker
(NoFlicker / 20Hz / 40Hz, 1h/day×7d) in 6-mo 5xFAD vs WT mice, both sexes.
Readouts: glial **morphology** (IMARIS), **cytokines** (Luminex 32-plex),
**amyloid** (Aβ40/42 ELISA), and correlations among them.

The Wood lab asked for statistical corrections to the **cytokine + correlation**
analyses (NOT the morphology). Those requests are in
`initialComments.md` (top level). We re-ran that analysis cleanly in R and are
now updating the manuscript + figures to match.

**Round 2 (Aug 2026):** Ashley sent a revised spec that (a) makes the
linear-range prune a real removal and (b) picks the focus cytokines from the
PLS-DA loadings (top 3 on LV1, top 3 on LV2) instead of the myeloid panel. That
is implemented; §4 and §5 below are the round-2 numbers. The full before/after
is in `cytokine_analysis/CHANGES.md` ("Round 2 — the revised spec").

**Key working directory:** `/Users/williamxu/Documents/Research/BALES/AshleyPrichard/`

---

## 2. Environment / gotchas (read before editing anything)

- **Python:** use `/usr/local/bin/python3` (has `python-docx` 1.2.0, `openpyxl`).
  The project `.venv-lr` is broken.
- **R:** `/usr/local/bin/Rscript` (R 4.5.2). All needed packages installed
  (ropls, heatmap3, rstatix, lmerTest, emmeans, nortest, FSA, ggpubr, tidyverse…).
- **The active manuscript file is `7.21.26ManuscriptWilliam.docx`** — this is
  William's working copy. It has been opened/saved in **Word, which re-split the
  runs**, so any programmatic edit must be **content-based** (search for the text),
  NOT by paragraph/run index. Helper pattern: concatenate a paragraph's run texts,
  find the substring, splice across runs. (See `scratchpad` scripts referenced in
  §7 if still present, else re-derive.)
- **Always back up before editing** (`shutil.copyfile`). Existing backups:
  `7.21.26ManuscriptWilliam_BACKUP.docx`, `_BACKUP2.docx`.
- **Comments** are real Word tracked-comments. python-docx `document.add_comment([run], text, author, initials)` adds them; editing existing comment text needs direct XML on the comments part (find `w:comment` by `w:author`).
- Do NOT touch `ManuscriptMostRecent.docx` — that's the untouched original.

---

## 3. File map

```
AshleyPrichard/
├── ManuscriptMostRecent.docx          # ORIGINAL, untouched
├── 7.21.26ManuscriptWilliam.docx      # ACTIVE working file (edit this)
├── 7.21.26ManuscriptWilliam_BACKUP*.docx
├── ManuscriptMostRecent_REVISED_stats.docx  # my earlier reference copy (superseded by William's file)
├── initialComments.md                 # the Wood-lab corrections (source of truth for "what to change")
├── HANDOFF.md                         # this file
├── SARA Temp/                         # original Wood-lab data + code
│   ├── AshleyDataAnalysis_SBversion_Final.R   # the original 3652-line analysis (S. Bitarafan)
│   ├── WoodLabFunctions.R                     # scores.plot.ellipse, loadings.chart2, scatter helpers
│   ├── Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx  # combined per-animal data (sheet "Simplified", A2:BP47)
│   └── ...                                     # per-cell morphology CSVs, BCA/linear-range, etc.
└── cytokine_analysis/                 # THE RE-ANALYSIS PIPELINE (our work)
    ├── config.R                       # all knobs: analyte lists, colors, thresholds, FDR scope
    ├── run_pipeline.R                 # runs stages 01–10
    ├── README.md
    ├── CHANGES.md                     # detailed change log + review-pass summary  <-- read this
    ├── CORRECTIONS_TRACEABILITY.md    # every change mapped to an initialComments.md line
    ├── lib/  (stats_helpers.R, WoodLabFunctions.R)
    ├── original/  (verbatim copies of the Wood-lab scripts)
    ├── scripts/  01_load … 10_figure5_correlations.R
    └── results/
        ├── stats/                     # tidy CSVs: one row per test (the authoritative corrected numbers)
        └── revised_figures_png/       # corrected figures to paste into the manuscript (Canva)
```

Run the whole pipeline: `cd cytokine_analysis && Rscript run_pipeline.R` (~2–3 min).

---

## 4. Config decisions already made (in `config.R`)

Round 2 settings — these are what the current numbers come from:

- `keep_myeloid_regardless = FALSE` — **strict linear-range prune.** G-CSF,
  GM-CSF, M-CSF and LIF are removed (all sit on the assay floor). The revised
  spec picks focus cytokines from the LVs, so the myeloid panel no longer needs
  protecting.
- `apply_nochange_screen = TRUE` — a floor/no-change screen on top of the curated
  panel (few distinct values, narrow range, many readings at zero). Removes
  MIP-1b as well. **Final panel: 16 analytes.** Set FALSE to put MIP-1b back.
- `corr_fdr_scope = "by_row"` — **DECIDED (round 2).** Correlation FDR is BH
  across the 16 retained heatmap analytes, per target ("FDR correct based on
  pre-selected cytokines based on the heatmap"). Key-set-restricted tables are
  written as a sensitivity analysis.
- `key_n_per_lv = 3` — top 3 per latent variable, per genotype (stage 04b).
- `plsda_loocv_runs = plsda_n_perm = 1000` — matches Garza et al.
- Group comparisons use **Tukey / Dunn only**; no second FDR layer across
  analytes. `p_omnibus_fdr` stays in the CSVs as a diagnostic and is not cited.

---

## 5. The corrected results (round 2 — what the manuscript must match)

- **Analyte panel: 16.** Dropped for linear range / floor: G-CSF, GM-CSF, M-CSF,
  LIF, MIP-1b, plus the 11 already outside the curated panel.
- **Key cytokines (top 3 LV1 + top 3 LV2, per genotype):**
  WT = IP-10, IL-10, IL-1B / VEGF, MIP-1a, KC.
  5xFAD = KC, IL-1a, IL-17 / MIP-1a, IL-10, IP-10.
- **Group comparisons: nothing significant among the key cytokines**, either
  genotype, either route (ANOVA+Tukey / KW+Dunn, or the sex-controlled mixed
  model). Across all 16 retained analytes only **WT IFNy** (No Flicker vs 20Hz
  and vs 40Hz, Dunn-BH p=0.044 both) and **WT-female LIX** (20Hz vs 40Hz,
  Tukey p=0.048) clear 0.05. M-CSF (round 1's only survivor) is gone because the
  analyte is out of linear range, not because the test changed.
- **PLS-DA: fails permutation in both genotypes.** WT LOOCV accuracy 0.250 vs
  chance 0.333, p_permute = 0.585; 5xFAD 0.304, p_permute = 0.376. Visualization
  and a descriptive weighting only. The old LV-score p-values (0.0031 etc.) were
  invalid and stay removed.
- **Amyloid: unchanged from round 1** (it doesn't depend on the cytokine panel).
  Aβ40 40Hz vs No Flicker = **p=0.060 (n.s.)** (was 0.024).
- **Cytokine↔amyloid (Fig 4C): nothing survives FDR** over the 16-analyte panel.
  Closest: LIX~Aβ42 (ρ=0.571, p=0.0044, q=0.069) and MIG~Aβ42 (ρ=0.534, p=0.0087,
  q=0.069). Round 1's GM-CSF~Aβ42 carried that row and is now excluded.
- **Cytokine↔morphology (Fig 5): 15 pairs survive**, nearly all 5xFAD. Microglia
  density ~ MIG (ρ=0.71, q=0.0022), IL-10 (0.69, q=0.0023), RANTES (0.65);
  microglia tiling ~ MIG (−0.70), IL-10 (−0.66); astrocyte density ~ IL-10
  (−0.61), RANTES (−0.56); WT astrocyte volume ~ IL-2 (−0.64). **This is the
  cytokine result that holds up.**
- **Sex effects are real:** IL-10 (WT p=0.0011, 5xFAD p=5.8e-05), RANTES
  (4.6e-05), LIX, VEGF, IFNy in 5xFAD, IL-1B in WT; Flicker×Sex interactions for
  WT IL-9 (0.041) and WT MIP-1a (0.023).
- **Morphology group stats (Fig 1/3/6):** still NOT reproducible from the
  per-cell files in `SARA Temp/` — those were made in MATLAB from an export not
  in the folder. Left untouched. Needs canonical source before any correction.

Full per-test numbers: `cytokine_analysis/results/stats/*.csv`.
Before/after against round 1: `cytokine_analysis/CHANGES.md` § "Round 2".

---

## 6. Figure swap guide (paste PNGs into the composites in Canva)

Corrected data panels are in `cytokine_analysis/results/revised_figures_png/`.
They are DATA panels only — the published figures are composites with microscopy /
schematic panels that can't be regenerated, so swap the data sub-panels in.

| Manuscript panel | Use file | Notes |
|---|---|---|
| Fig 2C (PLS-DA scores) | `Fig2C_plsda_scores_WT.png` | replaces the scores plot AND deletes the LV-score violin panels; subtitle carries LOOCV accuracy + p_permute |
| Fig 2D (Signals LV1/LV2) | `Fig2D_loadings_WT_LV1.png`, `_LV2.png` | 16-analyte model, LOOCV error bars, top-3 picks in darker fill |
| Fig 4B (amyloid bars) | `Fig4B_amyloid_bars.png` | Aβ40 now n.s. |
| Fig 4C (ρ-heatmap) | `Fig4C_cytokine_amyloid_heatmap.png` | rows = 5xFAD key cytokines; **no asterisks** — nothing survives FDR |
| Fig 5 (correlations) | `Fig5A/B/C/D_*.png` | Sara's scatter style + M/F shapes; C/D now show the key cytokines (was RANTES/Eotaxin/IL-15/IL-10 — two of those are out of range) |
| Supp 4 B (5xFAD PLS-DA) | `Supp4B_plsda_scores_5XFAD.png` | delete its LV-score violins |
| Supp 4 C/D (loadings) | `Supp4B_loadings_5XFAD_LV1.png`, `_LV2.png` | C = LV2, D = LV1 (relabel already in legend) |
| Supp 5 (WT cytokine bars) | `Fig_key_cytokines_WT.png` | **replaces `Fig_myeloid_panel_WT.png`, which has been deleted.** Nothing significant |
| Supp 7 (WT correlations) | `Supp7A/B/C/D_*.png` | Spearman + M/F, no CI ribbon, key cytokines |
| (new, optional) sex heatmaps | `Fig_sexbatched_heatmap_WT/5XFAD.png` | worth promoting — sex is a real effect here |

Regenerate any of these: `Rscript scripts/09_export_revised_figures.R` (and
`scripts/10_figure5_correlations.R` for Fig 5 / Supp 7). Both read the PLS-DA
model saved by stage 04 (`results/processed/plsda_fits.rds`), so they cannot
drift from the exported statistics.

---

## 7. Manuscript edit state

`7.21.26ManuscriptWilliam.docx` has been updated to the **round-2** numbers
(backup: `7.21.26ManuscriptWilliam_BACKUP4_pre_round2.docx`). Applied:
- **Methods**: analyte selection rewritten (16 analytes, floor screen, exploratory
  myeloid supplement); PLS-DA permutation described as Garza et al. (cites refs
  7/8, which were already in the list); key-cytokine definition added; FDR family
  stated as the 16 retained analytes per target; cytokine values noted as
  background-subtracted MFI, not standard-curve interpolated.
- **Abstract, highlights, Results, Discussion, Summary**: myeloid frequency claims
  replaced with what the data support (IFN-γ in WT; cytokines tracking glial
  morphology in 5xFAD).
- **Figure legends**: Fig 2, 4, 5 and Supp 4, 5, 6, 7 retitled and renumbered to
  the corrected stats. New **Supplementary Figure 9** (exploratory myeloid +
  detection-range panel) added after Supp Fig 8.
- **Supplementary Table 2**: every cytokine block regenerated (Fig 2, 4, 5;
  Supp 5, 6, 7) and a Supp Fig 9 block appended.
- **28 "W. Xu (re-analysis)" comments**: the 19 that cited round-1 numbers were
  rewritten in place; resolved ones say so.

An audit for stale round-1 strings (old permutation p's, "21 analytes",
"17 of 32", myeloid FDR p's, "distinct clustering") returns clean.

---

## 8. Multi-agent review pass — findings (done 2026-07-22)

A 7-dimension agent swarm cross-checked the finalized draft; 64 issues confirmed.
Details in `cytokine_analysis/CHANGES.md` ("Final multi-agent review pass").

**Already fixed** in the working file: heading "I-40"→"1-42"; Fig 4C "LIF"→"MIG"
and "flicker"→"Aβ1-42 burden"; "5xFAD"→"WT" for Supp 7b; Supp 6 IFN-γ 0.018/0.041
→0.044 n.s.; typos (ytokine, cyotkine, make of up, include→including, doubled "here").

**Flagged (comments, not rewritten — need author decisions):**
- Abstract + Highlights: myeloid claim ("significantly differed 20Hz vs 40Hz") —
  only M-CSF (vs No Flicker) survives.
- Fig 2 title "distinct clustering" (PLS-DA fails permutation).
- Results "40Hz may be more effective than 20Hz" (Aβ40 n.s.).
- Results WT "more pronounced GM-CSF/G-CSF/LIX/LIF" (all n.s.; LIF out of range).
- Results 5xFAD "robust proinflammatory response" (no sig flicker effect in 5xFAD).
- Supp Fig 4 title copied from Supp 6 (retitle to match heatmap/PLS-DA).

**Reported, NOT touched (out of scope):** grammar + a broken stat
(`F(2)=292, p=5.4115`; `p = < .001`) in the **morphology** supplements (MATLAB
figures, pending source data).

---

## 9. Open TODOs (prioritized)

1. **[BIGGEST] Rewrite the narrative to match the corrected stats.** The abstract,
   Annabelle's highlights, intro preview, and discussion still tell the old story
   (strong myeloid frequency effects; 40Hz reduced amyloid; PLS-DA separation).
   **The round-2 reality:** no cytokine differs between flicker groups among the
   key cytokines (only WT IFNy across the whole retained panel); the myeloid
   analytes are out of the analysis entirely; amyloid reduction n.s.; PLS-DA
   fails permutation; no cytokine↔amyloid correlation survives FDR. What *does*
   hold is cytokine↔microglial-morphology in 5xFAD (MIG/IL-10/RANTES, |ρ| 0.65–0.71)
   and sex differences in several analytes — the paper's cytokine story probably
   needs to be rebuilt around those. See the flagging comments for suggested
   wordings. (Note in the doc: Annabelle asked that AI *flag* not rewrite the
   science — draft, don't auto-apply.)
   **The 19 existing "W. Xu (re-analysis)" comments cite round-1 numbers and need
   updating too.**
2. **Swap the corrected figures into the composites** in Canva (§6), and update the
   significance markers on Fig 4B / Supp 5 to match.
3. **Supp 6 decision:** keep-as-exploratory (with the caveat comment) or cut —
   most of its analytes are now formally outside the analysis.
4. ~~**Fig 5 Eotaxin / IL-15**~~ — resolved in round 2: Fig 5 C/D and Supp 7 now
   show the per-genotype key cytokines, so neither out-of-range analyte appears.
5. **Supp 4 title:** retitle to match its content (heatmap + PLS-DA).
6. **Morphology figures (Fig 1/3/6, Supp 2/3/8):** need the canonical MATLAB
   source data to (a) verify the H-statistics and (b) fix the grammar/broken-stat
   errors flagged in the review. Out of scope until that data is provided.

---

## 10. How to pick this up in a new chat

Point the assistant at this file and `cytokine_analysis/CHANGES.md` +
`CORRECTIONS_TRACEABILITY.md`. Emphasize: edit `7.21.26ManuscriptWilliam.docx`
**by content** (Word re-split the runs), back up first, use `/usr/local/bin/python3`.
The authoritative corrected numbers are in `cytokine_analysis/results/stats/`.

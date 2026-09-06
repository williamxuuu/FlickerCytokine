# Figure swap guide — round 2

Every data panel below is in `results/revised_figures_png/`. The published figures
are Canva composites containing microscopy and schematic panels that cannot be
regenerated, so swap the **data sub-panel** in and leave the rest.

| Manuscript panel | File | What changed this round |
|---|---|---|
| Fig 2B (heatmap) | `results/cytokine/heatmap_WT.pdf` | 16 analytes, sex batch-corrected (was all 32, unadjusted) |
| Fig 2C (PLS-DA scores) | `Fig2C_plsda_scores_WT.png` | subtitle now carries LOOCV accuracy + p(permutation); LV-score violins stay deleted |
| Fig 2D (Signals LV1/LV2) | `Fig2D_loadings_WT_LV1.png`, `_LV2.png` | LOOCV error bars added; top-3 picks in darker fill |
| Fig 4B (amyloid bars) | `Fig4B_amyloid_bars.png` | unchanged from round 1 (Aβ40 n.s., p = 0.060) |
| Fig 4C (ρ heatmap) | `Fig4C_cytokine_amyloid_heatmap.png` | new rows, ρ printed in-cell, **all asterisks removed** — nothing survives FDR |
| Fig 5A/B | `Fig5A_amyloid_vs_astro_tiling.png`, `Fig5B_amyloid_vs_mg_tiling.png` | unchanged |
| Fig 5C/D | `Fig5C_cytokine_vs_astro_number.png`, `Fig5D_cytokine_vs_mg_number.png` | now the 5xFAD key cytokines; Eotaxin and IL-15 removed (out of range) |
| Supp 4B (5xFAD PLS-DA) | `Supp4B_plsda_scores_5XFAD.png` | permutation result in subtitle |
| Supp 4C/D (loadings) | `Supp4B_loadings_5XFAD_LV1.png`, `_LV2.png` | C = LV2, D = LV1 (as already relabelled); error bars added |
| Supp 5 (WT + 5xFAD) | `Fig_key_cytokines_WT.png`, `Fig_key_cytokines_5XFAD.png` | **replaces the myeloid panel entirely.** The old `Fig_myeloid_panel_*.png` files were deleted so they cannot be pasted by mistake |
| Supp 6 | `Supp6_flicker_responsive_WT.png` | **replaces the old 6-analyte panel**; five of those six analytes are out of range |
| Supp 7A–D | `Supp7A_*.png` … `Supp7D_*.png` | key cytokines; ρ, p and q annotated |
| **Supp 9 (new)** | `Supp9_myeloid_exploratory_bars.png`, `Supp9_myeloid_exploratory_detection.png` | the four excluded myeloid analytes plus the detection-range panel justifying their exclusion |
| (optional) sex heatmaps | `Fig_sexbatched_heatmap_WT.png`, `_5XFAD.png` | worth promoting — sex is one of the more robust effects in the data |

## Significance markers that must change on the composites

- **Fig 4C**: remove every asterisk. Under FDR across the 16 retained analytes, no
  cytokine–amyloid correlation survives (LIX and MIG ~ Aβ42 are closest at q = 0.069).
- **Old Supp 5 (myeloid)**: the published panel carries LIX 0.045\*, G-CSF 0.028\*,
  GM-CSF 0.047\* and M-CSF 0.008\*/0.049\*. That panel is replaced; if any of it is
  reused in Supp 9, only WT M-CSF 20Hz vs No Flicker (p = 0.014) keeps a marker.
- **Fig 5C/D**: asterisks now follow q, not raw p.

Regenerate everything with `Rscript run_pipeline.R` (~7 min), or a single stage with
`Rscript -e 'source("config.R"); source("scripts/09_export_revised_figures.R")'`.

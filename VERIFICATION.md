# Verification report

What was checked before this package was assembled, what passed, and what was
found. Performed 2026-09-26 against `AshleyPrichard/cytokine/` on R 4.5.2
(`/usr/local/bin/Rscript`), macOS 25.6.0.

**Headline:** the analysis reproduces end to end, byte for byte. Six problems
were found — all in the *packaging and documentation*, none in the statistics or
the figures. All six are resolved: three path fixes, one documentation rewrite,
and two passes removing code and data the analysis never touches. Three further
items are recorded as things to know rather than defects.

Everything removed was verified unreachable first, and the pipeline was re-run
afterwards to confirm it still produces identical output.

---

## 1. Reproducibility — the whole chain, from the raw plate file

Each link was re-run in an isolated copy and compared against what was shipped.

| Link | Method | Result |
|---|---|---|
| raw `.xlsx` → `animal_table.csv` | re-ran the loader standalone | **byte-identical** (43 animals × 59 cols, from 45 rows = 43 animals + 2 background wells) |
| `animal_table.csv` → all statistics | full `run_all.R` in a clean copy | **all 9 CSVs byte-identical** |
| statistics → all figures | same run | **all 26 PNGs byte-identical**; all 26 PDFs identical apart from the embedded creation timestamp |
| PLS-DA, cache deleted | recomputed both genotypes from scratch (~13 min) | **exactly reproduces the cached values** — see below |

The PLS-DA is the one stage with a random component (1,000 label permutations,
1,000 LOOCV resamples). With the cache deleted and `set.seed(42)`:

| | LOOCV | permutation p | key LV1 | key LV2 |
|---|---|---|---|---|
| WT | 0.400 | 0.202 | IL-15, IL-6, LIF | IL-2, IL-17, LIX |
| 5xFAD | 0.348 | 0.320 | KC, IP-10, IL-17 | LIF, M-CSF, IL-12(p70) |

Identical to the cached fit in every value. The cache is keyed on a fingerprint
of the actual input matrix, so it cannot silently go stale.

### Clean-room test — the one that matters for Box

The finished package was copied to an unrelated directory, **every derived file
deleted** (all statistics, all figures, the processed table and both PLS-DA
caches), leaving only the code and the raw `.xlsx`, and rebuilt from nothing:

```
Rscript code/prep/01_load_and_prep.R     # raw .xlsx -> animal_table.csv
cd code && Rscript run_all.R             # 11 stages, 7.8 min
```

Against the shipped package:

| | Result |
|---|---|
| `animal_table.csv` | identical |
| all 10 statistics CSVs | identical |
| all 26 figure PNGs | identical |
| PLS-DA, both genotypes | LOOCV and permutation p reproduce exactly |

So someone who downloads this folder and has the R packages listed in the README
can regenerate every number and every panel from the raw plate file, with no
path editing and no reliance on anything outside the folder.

## 2. Every number in the README was re-derived from the data

Not copied forward — recomputed and compared. All matched.

- **Panel gate.** 28 of 32 analytes above the >80% threshold; dropped IL-5,
  IL-12(p40), MCP-1, TNFa. ✓
- **Group sizes.** WT n = 20, 5xFAD n = 23; per-arm sex counts as documented. ✓
- **Sex balance.** Worst arm 57.1% one sex, under the 80% hard stop. ✓
- **Z-score effect.** Sex effects at p < 0.05: 13 → 0. Flicker p-values before
  vs after: Spearman 0.915. Flicker p < 0.05: 7 → 5. ✓
- **Strongest sex effects removed.** 5xFAD IL-10 7.5e-05, RANTES 1.3e-04,
  WT IL-10 0.0022. ✓
- **Group effects.** WT 5 of 28 omnibus p < 0.05, 5xFAD 0 of 28; all six
  post-hoc contrasts match to the quoted digit (IL-15 0.0034, IL-6 0.0086,
  M-CSF 0.0211, Eotaxin 0.0325, IL-15/40Hz 0.0363, IFNy 0.0381). ✓
- **The MIG caveat.** 5xFAD MIG Tukey p = 0.0432 with omnibus p = 0.0506 —
  the README's "should not be claimed" is correct. ✓
- **Amyloid.** Abeta 40 KW p = 0.0393, both 40 Hz contrasts p = 0.0481;
  Abeta 42 ANOVA p = 0.534. ✓
- **Correlations.** 9 pairs clear BH, 5 also clear Westfall–Young, with the
  quoted rho / q / WY values and the "5 distinct values, 33% ties" note on
  GM-CSF ~ Abeta 42. All 5xFAD. ✓

Figures were also opened and read: Fig 2C states 28 analytes, LOOCV 40% vs 33%,
permutation p = 0.20, and carries **no statistics on the violins**, as the Wood
lab required. Fig 4C shows GM-CSF ~ Abeta 1-42 at 0.65 with exactly one star,
consistent with q = 0.030.

---

## Findings

### F1 — the prep stage could not run at all *(fixed)*

`prep/01_load_and_prep.R` did `source("config.R")`. No file of that name
shipped; the config was `config_round2.R`. Running the stage failed immediately:

```
cannot open file 'config.R': No such file or directory
```

This is the step that turns the raw plate workbook into the analysis table, so
as shipped the analysis could not be rebuilt from raw data by anyone.

**Fixed.** The config is now `code/prep/config.R`, and the loader sources it by
a path resolved relative to its own location, so it runs from any working
directory.

### F2 — the prep config pointed at a folder layout that no longer exists *(fixed)*

```r
project_root <- normalizePath(".", mustWork = TRUE)
data_dir     <- normalizePath(file.path(project_root, "..", "SARA Temp"), mustWork = TRUE)
```

This resolved the raw data as a *sibling* of the working directory. After the
2026-09-23 reorganisation the raw data sits at `AshleyPrichard/data/SARA Temp`,
so the path never resolves and `mustWork = TRUE` turns that into a hard error.

**Fixed.** Paths are now self-locating from the script's own position, and the
raw workbook ships inside the package at `data/raw/`.

### F3 — prep wrote its output where nothing reads it *(fixed)*

`results_dir` was `<project_root>/results`, so the loader wrote
`results/processed/animal_table.csv`. But `00_setup.R` reads
`data/processed/animal_table.csv`. Even with F1 and F2 fixed, a rebuild would
have written a table the pipeline then ignored, leaving the stale one in place —
the quiet version of this bug, and the more dangerous one.

**Fixed.** The loader now writes to the package's own `data/processed/` and
`stats/`, which is where the pipeline reads from.

> **These three fixes change paths only.** No analysis value, threshold or
> formula was touched. Proof: after the fix, `01_load_and_prep.R` regenerates
> `animal_table.csv` **byte-identically** to the shipped one.

### F4 — the figures manifest was stale and contradicted the results *(rewritten)*

`figures/MANIFEST.md` described an earlier round. Every one of these claims is
contradicted by the shipped `stats/`:

| MANIFEST claimed | Actually |
|---|---|
| correlations are "sex-stratified permutation Spearman" | free permutation; sex-stratified is a sensitivity column only |
| GM-CSF ~ Abeta 42 at q = 0.017 | q = 0.030 |
| that cell fails Westfall–Young at p = 0.25 | it **passes**, p = 0.024 |
| Supp 7 analytes "picked per panel by significance" | they are the PLS-DA key cytokines |
| Supp 7 best q = 0.72 / 0.58 / 0.32 / 0.97 | 0.681 / 0.656 / 0.221 / 0.989 |
| Fig 5D "only LIF significant" | nothing is significant; LIF and M-CSF tie at q = 0.12 |
| 5xFAD keys are KC, IL-1a, IL-17, IL-12(p70), LIF, G-CSF | KC, IP-10, IL-17, LIF, M-CSF, IL-12(p70) |
| the palette was changed to pass a colourblind check | the published hexes were deliberately kept; the trade-off is documented in `00_setup.R` and handled by using shape, not colour, for identity |

It also pointed at a `finalCode/` directory that no longer exists.

**Rewritten** from the shipped statistics. Had this gone to Box unchanged, a
reader would have had two documents disagreeing about the results.

### F5 — round-2 artifacts sat unlabelled beside the current ones *(removed)*

`data/processed/` held five files that **nothing in the live pipeline reads**:
`animal_table_sexadj.csv` / `.rds`, `analysis_analytes.rds`,
`key_cytokines.rds` and `plsda_fits.rds`. All are round-2 outputs, produced only
by the superseded `prep/02_linear_range_and_batch.R`.

The sharp edge: round 2's **sex mean-shift** table (`animal_table_sexadj.csv`)
sat directly beside round 3's **within-sex z-score** table
(`animal_table_sex_zscored.csv`), with nothing marking which was current. Both
are plausible-looking per-animal tables of the same shape, and picking the wrong
one silently reverts the single decision this analysis rests on.

**Removed from the package** (see finding F6). All five remain in
`AshleyPrichard/cytokine/data/processed/`.

### F6 — the package carried code the analysis never runs *(pruned)*

A reachability audit found 1,156 lines of R that nothing in the pipeline
executes, plus five unused data files:

| Removed | Lines | Why |
|---|---|---|
| `code/lib/WoodLabFunctions.R` | 816 | never sourced by any stage |
| `code/prep/02_linear_range_and_batch.R` | 194 | round-2 analyte screen, replaced by the gate in `00_setup.R`; also needs a workbook not shipped |
| `code/lib/plsda_helpers.R` | 146 | never sourced — stage 02 implements its own fit, LOOCV and permutation |
| `data/processed/superseded_round2/` | 5 files | round-2 outputs, read by nothing (see F5) |

Two files were also trimmed rather than removed, because part of each is live:

- **`code/prep/config.R`: 221 → 74 lines.** Only 11 of its 43 settings were
  reachable; the rest configured the round-2 linear-range screen, the "Hayden"
  analyte rule, PLS-DA parameters and colour palettes that `00_setup.R` now
  owns. A guard was added so a wrong path fails loudly instead of resolving
  somewhere unintended.
- **`code/lib/stats_helpers.R`: 289 → 105 lines.** Six of its nine functions
  were unreachable (`compare_flicker_mixed`, `compare_flicker_mixed_many`,
  `corr_one`, `correlate`, `scatter_sex`, `write_stats`). These were also the
  ones referencing globals `00_setup.R` does not define — `sexes`,
  `corr_fdr_scope`, `cytokines_myeloid`, `sex_shapes`, `flicker_colors_wt`,
  `fmt` — so any call would have failed with "object not found". Removing them
  made two globals in `00_setup.R` dead in turn (`stats_dir`,
  `corr_stat_method`); both were dropped.

Comment noise was removed in the same pass: notes addressed to a person,
requests for a decision, and refactor history of the form "this used to be X".
Scientific rationale was kept in full — why the 80% gate, why the z-score, why
no statistics on the latent variables, what was matched from the published
panels, and every caveat. The headers still carry the reasoning, just not the
changelog.

Net effect on `code/`: **3,527 → 1,961 lines (44% smaller)**, of which 358 were
comment lines. Comment *density* rose from 27.8% to 31.7%, which sounds backwards
but is arithmetic: the largest thing deleted (`WoodLabFunctions.R`, 816 lines)
was only 12% comment, so removing it raised the average for everything left.

**Nothing was lost.** Every removed file remains in
`AshleyPrichard/cytokine/` and `AshleyPrichard/archive/`; the package is a
distribution copy, not the master.

**Re-verified after pruning:** all 10 statistics CSVs and all 26 figure PNGs
byte-identical, and a second clean-room rebuild from the raw workbook with no
caches reproduced everything.

---

## Things to know (no change made)

**T1 — IL-5 misses the panel by one animal.** 34 of 43 above background =
0.791, against the 0.80 gate; 35/43 = 0.814 would retain it. This is the
tightest call in the gate and it lands on a contested claim — Supplementary
Figure 4's published caption is specifically about IL-5. Worth stating in
Methods rather than leaving for a reviewer to find. The next-closest analytes
are on the other side of the line: LIX 0.837 and IFNy 0.860 are retained.

**T2 — the two permutation budgets are different.** Correlations use
B = 50,000, where a p near 0.04 is stable to roughly ±0.002 — so the four
Westfall–Young values between 0.024 and 0.041 sit safely under 0.05. The PLS-DA
permutation uses B = 1,000, which is much coarser; it does not matter here
because both p-values (0.202, 0.320) are nowhere near a threshold.

**T3 — morphology group statistics remain unreproducible**, as the original
README states. Those figures were produced in MATLAB from an export that is not
in the source folder. They are not part of this package. The per-animal
morphology *summaries* used for the correlations are in the analysis table and
do reproduce; the per-cell group comparisons do not.

---

## Reproducing this verification

```bash
# 1. rebuild the analysis table from the raw workbook
Rscript code/prep/01_load_and_prep.R
# 2. force a full PLS-DA recompute, then run everything
rm stats/plsda_*.rds
cd code && Rscript run_all.R
```

Expect: 11 stages, about 13 minutes with the PLS-DA rebuild (1.5 minutes
with the cache). Statistics and PNGs should come out byte-identical; PDFs will
differ only in their embedded creation timestamp.

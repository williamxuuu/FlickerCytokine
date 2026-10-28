# Raw data

## `Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx`

The one combined per-animal workbook the analysis runs off. Sheet
**`Simplified`**, range **`A2:BP47`** (header on row 2).

45 rows = **43 animals + 2 `Background0` wells**. Columns: metadata, the 32
Luminex analytes, amyloid (Abeta 40 / 42, plaque counts), and the per-animal
glial morphology summaries.

`code/prep/01_load_and_prep.R` reads it, subtracts the mean of the two
background wells from each analyte, clips negatives to 0, sets the factor
levels, and writes `data/processed/animal_table.csv`. Amyloid and morphology
are carried through unchanged.

> **No standards were run on this plate.** 43 unknowns and two blanks,
> `AnalysisType = None`, blank units. So there is **no LOD or LLOQ**, values are
> **MFI, not pg/mL**, and "background" means those two blank wells. The
> 28-analyte detectability gate is defined against that background, not against
> a real detection limit. This needs saying in Methods.

## What is deliberately not here

- **The BCA loading-curve workbook** (`BCA_Analyzed_PrichardSamples_11425_V2.xlsx`)
  — needed only by the superseded round-2 analyte screen, which is not part of
  this package.
- **Per-cell morphology exports** — used by the MATLAB morphology figures, which
  are not part of this analysis and are not reproducible from the source folder
  (see T3 in `../../VERIFICATION.md`).

Both live in the lab's original `AshleyPrichard/data/SARA Temp/` folder.

# Input data

The pipeline starts from a single Excel workbook in this folder. The study
data are not distributed with this repository.

## Expected layout

One sheet, with a header row and **one row per animal**, followed by the
blank (background) wells from the Luminex plate. The location is set in
`code/prep/config.R`:

| Setting | Value used in this study |
|---|---|
| `combined_xlsx` | `data/raw/Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx` |
| `combined_sheet` | `Simplified` |
| `combined_range` | `A2:BP47` (header on row 2; 43 animals + 2 blank wells) |

## Required columns

| Column | Content |
|---|---|
| `Genotype` | `WT` or `5XFAD`. Left empty for blank wells |
| `Sex` | `M` or `F` |
| `Flicker` | `NoFlicker`, `20Hz` or `40Hz` |
| `Sample` | sample identifier. Blank wells are `Background0` |
| `litter` | litter identifier, used only by the mixed-model sensitivity analysis |
| one column per analyte | Luminex median fluorescence intensity (MFI), named as in `cytokines_all` (`config.R`) |
| `Abeta 40`, `Abeta 42` | amyloid measurements (5xFAD), named as in `amyloid_vars` |
| glial morphology summaries | per-animal mean / median tiling, cell volume and cell number, named as in `morph_vars` |

Optional identifier columns (`Sample Number`, `Subject`, `BlindName`, `Cohort`)
are carried through if present.

## Processing

`code/prep/01_load_and_prep.R` performs these steps:

1. Identifies the blank wells: rows with no genotype, or `Sample == "Background0"`.
2. Subtracts the mean of the blank wells from every analyte, and sets values
   below background to 0.
3. Writes `data/processed/animal_table.{csv,rds}`.

Amyloid and morphology values are carried through unchanged.

The plate in this study was run without standard curves. Analyte values are
therefore relative MFI rather than concentrations, and the detectability
filter (`code/00_setup.R`) is defined relative to the blank wells.

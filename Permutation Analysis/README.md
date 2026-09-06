# Permutation analysis (PLS-DA model validity)

Tests whether the PLS-DA models classify flicker group better than chance,
following Garza et al.: leave-one-out cross-validated classification accuracy is
compared against a null built by randomly reassigning the group labels.

- 1000 LOOCV iterations, 1000 label permutations, per genotype.
- `04_plsda.R` — the fit and permutation code (`plsda_loocv_runs` and
  `plsda_n_perm` are set in `../Code/config.R`).
- `plsda_permutation.csv` — the full result.

## Result

| Genotype | LOOCV accuracy | Chance | Null mean | p_permute | R2Y | Q2 |
|---|---|---|---|---|---|---|
| WT (n=20) | 0.250 | 0.333 | 0.261 | 0.585 | 0.266 | -0.517 |
| 5xFAD (n=23) | 0.304 | 0.333 | 0.244 | 0.376 | 0.237 | -0.202 |

**Neither model passes.** Both classify at or below chance, both p_permute are far
from significance, and both Q2 values are negative (the model predicts worse than
the mean). The PLS-DA is therefore reported as visualization and descriptive
weighting only — it does not support a claim of group separation. Any statement
that the groups cluster distinctly is not supported by this test.

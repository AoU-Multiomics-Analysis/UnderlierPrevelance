# RNA dual outlier definitions

## Goal

Generate RNA underlier prevalence results for both residual expression-z-score
calls without the haplo filter and calls that require the haplo expression-drop
filter as well as the z-score cutoff.

## Output contract

Keep the haplo-only outputs:

- `underliers_haplo.tsv.gz`
- `rna_outlier_prevalence_per_gene_haplo.tsv`

For every configured z-score cutoff, emit two explicitly named definitions:

- `underliers_z_<cutoff>.tsv.gz` and
  `rna_outlier_prevalence_per_gene_z_<cutoff>.tsv`: rows meeting only the
  residual expression-z-score condition, `expression_zscore < cutoff`.
- `underliers_haplo_z_<cutoff>.tsv.gz` and
  `rna_outlier_prevalence_per_gene_haplo_z_<cutoff>.tsv`: rows meeting both the
  haplo condition and `expression_zscore < cutoff`.

The haplo condition remains the existing per-gene adjusted-log2-CPM rule:

`expression_logcpm < gene_mean_logcpm - logcpm_drop`

Both comparisons remain strict (`<`). Prevalence continues to count unique
gene/subject pairs and uses the same tested-subject denominator.

## Implementation

Separate call construction from the definition loop so the z-only path starts
from all rows, while the intersection path starts from haplo-positive rows.
Reuse the existing z-cutoff configuration and prevalence calculation. For each
cutoff, write one z-only pair and one haplo-plus-z pair.

Update the RNA WDL task's artifact-count validation to expect:

`1 + 2 * length(z_cutoffs)`

underlier artifacts and the same number of prevalence artifacts. Update the
runtime documentation and tests to describe and assert the new meanings and
filenames.

## Testing

- Add focused tests for z-only calls, haplo-plus-z intersection calls, strict
  cutoff behavior, and prevalence counts.
- Update the RNA smoke-output test for both output pairs.
- Run the available focused test suite and static WDL validation. If the local
  environment lacks the RNA container dependencies, report that limitation and
  still run dependency-free tests and direct function-level assertions.

## Scope

This change does not alter normalization, PCA, residualization, z-score
calculation, cutoff configuration, or the haplo definition itself. It changes
only the call definitions, output names, artifact-count validation, tests, and
documentation.

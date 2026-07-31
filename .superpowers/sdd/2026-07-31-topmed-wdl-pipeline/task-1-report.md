# Task 1 implementation report

## Scope

Implemented only the Task 1 synthetic fixtures and validation tests in the repository. The fixtures provide:

- A standard GCT v1.2 count matrix with six genes and four samples (`S1`–`S4).
- A sample-keyed covariate TSV with `Genotype_PC1` and `Genotype_PC2`.
- Two protein-coding GFF3 gene records.
- A ClinVar VCF and cohort VCF sharing three pathogenic alleles.
- Input-contract tests for GCT headers/dimensions, covariate columns and alignment, annotation content, and shared pathogenic alleles.
- q² output-schema and variant-filtering tests for the scripts planned in later tasks.

## Changed files

- `tests/fixtures/counts.gct`
- `tests/fixtures/genotype_covariates.tsv`
- `tests/fixtures/annotation.gff3`
- `tests/fixtures/clinvar.vcf`
- `tests/fixtures/cohort.vcf`
- `tests/test_inputs.py`
- `tests/test_q2_script.py`
- `.superpowers/sdd/2026-07-31-topmed-wdl-pipeline/task-1-report.md`

## Verification

- `pytest -q tests/test_inputs.py` — could not run: `zsh:1: command not found: pytest` (exit 127).
- `python3 -m pytest -q tests/test_inputs.py` — could not run: `No module named pytest` (exit 1).
- `pytest -q tests/test_q2_script.py` — could not run: `zsh:1: command not found: pytest` (exit 127).
- `python3 -m py_compile tests/test_inputs.py tests/test_q2_script.py` — passed (exit 0).
- `git diff --check` — passed (exit 0).
- Standalone fixture assertions equivalent to the input tests — passed.

## Commit

Created commit `fee27da` (`test: add synthetic pipeline fixtures`).

## Concerns

The q² tests are intentionally command-line integration tests for `scripts/compute_q2_incidence.py` and `scripts/filter_variants.sh`, which are not part of Task 1 and do not yet exist. They will remain unavailable until the later q² implementation task adds those scripts and the environment provides pytest/bcftools.

## Fix-round summary

- Added helper-level expected-failure contracts for malformed GCT headers/dimensions, duplicate and missing sample IDs, missing genotype-PC columns, and nonnumeric genotype-PC values.
- Added an index-pair contract covering genotype VCF/index arrays and ClinVar VCF/index inputs without requiring binary `.tbi` files or tabix.
- Made all fixture and script paths resolve from the test file/repository root, independent of pytest’s working directory.
- Changed the filtering test threshold from `0.01` to `0.25`; the three singleton heterozygotes have AF `0.125`, while the homozygous alternate site has AF `1.0`, so the intended retained count is deterministically three after `+fill-tags`.
- Verification after this fix round: `py_compile`, `git diff --check`, and direct fixture/helper assertions pass. Pytest remains unavailable (`pytest` command missing and Python reports `No module named pytest`).

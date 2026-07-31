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
- `tests/fixtures/cohort_2.vcf`
- `tests/fixtures/cohort_2.vcf.tbi` (placeholder path file)
- `tests/fixtures/vcf_inputs.tsv`
- `tests/fixtures/cohort.vcf.tbi` (placeholder path file)
- `tests/fixtures/clinvar.vcf.tbi` (placeholder path file)
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

Created commits `fee27da` (`test: add synthetic pipeline fixtures`), `a234a37` (`test: strengthen pipeline fixture contracts`), and `804567c` (`test: close Task 1 fixture contract gaps`).

## Concerns

The q² tests are intentionally command-line integration tests for `scripts/compute_q2_incidence.py` and `scripts/filter_variants.sh`, which are not part of Task 1 and do not yet exist. They will remain unavailable until the later q² implementation task adds those scripts and the environment provides pytest/bcftools.

## Fix-round summary

- Added helper-level expected-failure contracts for malformed GCT headers/dimensions, duplicate and missing sample IDs, missing genotype-PC columns, and nonnumeric genotype-PC values.
- Added an index-pair contract covering genotype VCF/index arrays and ClinVar VCF/index inputs without requiring binary `.tbi` files or tabix.
- Made all fixture and script paths resolve from the test file/repository root, independent of pytest’s working directory.
- Changed the filtering test threshold from `0.01` to `0.25`; the three singleton heterozygotes have AF `0.125`, while the homozygous alternate site has AF `1.0`, so the intended retained count is deterministically three after `+fill-tags`.
- Verification after this fix round: `py_compile`, `git diff --check`, and direct fixture/helper assertions pass. Pytest remains unavailable (`pytest` command missing and Python reports `No module named pytest`).

## Fix-round 2 summary

- `read_gct_contract()` now checks every data-row field count and rejects empty Name or Description fields, with truncated-row and missing-metadata tests.
- `validate_covariates()` now checks every row’s field count before sample-ID and numeric-value validation, with a truncated-row test.
- Added `vcf_inputs.tsv` plus placeholder `.tbi` path files. Tests verify equal genotype VCF/index array lengths, ClinVar pairing, and path existence only; binary tabix validity remains the later WDL/bcftools task’s responsibility.
- q² helpers now load and validate the paired VCF/index manifest before invoking the later-task command interfaces, while retaining robust repository-root path resolution.
- Fix-round 2 verification: `python3 -m py_compile tests/test_inputs.py tests/test_q2_script.py`, `git diff --check`, and direct helper/fixture assertions pass. Pytest remains unavailable.

## Fix-round 3 summary

- Added `cohort_2.vcf` and its placeholder `.tbi` path, and expanded `vcf_inputs.tsv` to two genotype VCF/index pairs plus one ClinVar pair.
- The manifest parser now returns separate genotype VCFs, genotype indexes, ClinVar VCF, and ClinVar index values. Tests independently check list lengths and reject a shortened genotype index list.
- `run_q2()` validates and consumes the complete paired lists, then passes `cohort*.vcf` so both cohort VCFs are included. The `.tbi` files remain path placeholders only; binary tabix validation belongs to the later WDL/bcftools task.
- Fix-round 3 verification: `python3 -m py_compile tests/test_inputs.py tests/test_q2_script.py`, `git diff --check`, and direct helper/fixture assertions pass. Pytest remains unavailable.

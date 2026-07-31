# Task 5 Report: Reusable WDL Workflows

## Status

Implemented the RNA underlier, q² incidence, and conditional dispatcher WDL 1.0 workflows, plus a representative two-branch input JSON. All three workflows pass `miniwdl check` without lint notices.

The approved boundary is preserved: the RNA workflow accepts only GCT counts, genotype covariates, GENCODE annotation, RNA parameters, and runtime script/image dependencies. It contains no VCF, PLINK, or GDS inputs or processing. The q² workflow accepts paired genotype VCF/index arrays and a ClinVar VCF/index pair, validates and stages them with deterministic names, filters the genotype VCFs, and invokes the Python aggregator with the literal deterministic glob `cohort/*.filtered.vcf.gz`.

## Changed files

- `workflows/rna_underlier.wdl`
  - Converts GCT v1.2 to the existing counts-TSV contract with `scripts/convert_gct.py`.
  - Runs `scripts/run_rna_underlier.R` with the requested genotype-PC count, optional phenotype-PC noise override, connectivity threshold, log-CPM drop, and thread count.
  - Checks required artifacts and exposes selected-PC metadata, joined expression values, both underlier tables, both prevalence tables, and the complete RNA artifact glob.
  - Declares Docker image, CPU, memory, disk, and retry runtime settings for both tasks.
- `workflows/q2_incidence.wdl`
  - Rejects empty genotype input or unequal VCF/index array lengths in the staging command.
  - Validates `.tbi` basename pairing for every genotype and ClinVar pair and stages all pairs under deterministic names.
  - Calls the existing bcftools filtering script and emits `all.filtered.vcf.gz` with its index.
  - Stages the filtered pair as `cohort/0000.filtered.vcf.gz`, calls the existing Python aggregator with `cohort/*.filtered.vcf.gz`, and supports an optional gene whitelist.
  - Declares Docker image, CPU, memory, disk, and retry runtime settings for all three tasks.
- `workflows/main.wdl`
  - Imports both reusable sub-workflows.
  - Uses `if (run_rna)` and `if (run_q2)` call blocks.
  - Keeps branch data and repository-script inputs optional at the dispatcher boundary, requires them with `select_first` only inside the enabled call, and exposes optional branch outputs.
- `tests/test_wdl_inputs.json`
  - Supplies both branches with the RNA smoke fixtures, paired genotype/ClinVar fixtures, filtering parameters, existing scripts, and the Dockerfile's documented local image tag.

The four scripts are explicit `File` workflow dependencies because `envs/Dockerfile` installs their runtime dependencies but does not copy the repository scripts into the image. This also satisfies miniwdl's local-file allowlist; private `File` literals were tested and rejected by miniwdl before this correction.

## Red/green evidence

- Initial `miniwdl check` calls for all three workflow paths failed with `No such file or directory`, establishing the required RED state.
- The first dispatcher check caught and led to correction of the WDL scalar type from `Bool` to `Boolean`.
- A runtime probe caught miniwdl rejecting private repository-script `File` declarations as not expressly supplied workflow inputs. The scripts were promoted to explicit branch dependencies, after which execution reached container setup.
- Final `miniwdl check` passes for all three workflows.

## Tests and checks

Passed:

- `for wdl in workflows/*.wdl; do miniwdl check "$wdl"; done`
- `python3 -m json.tool tests/test_wdl_inputs.json`
- `miniwdl run -i tests/test_wdl_inputs.json workflows/main.wdl -j` (input names, types, and all local paths resolved)
- Disabled-dispatcher execution: both branches skipped and all branch outputs returned `null`.
- Enabled-branch validation probes: RNA-only and q²-only calls with omitted required branch inputs both failed at the enabled conditional call's `select_first`, as intended.
- Repository-script dependency probe: the RNA branch accepted all explicit file dependencies and reached `ConvertGct` container setup.
- `python3 -B -m py_compile scripts/convert_gct.py scripts/compute_q2_incidence.py tests/test_inputs.py tests/test_q2_script.py`
- `bash -n scripts/filter_variants.sh`
- Static boundary assertions: no `VCF`, `PLINK`, or `GDS` reference in `rna_underlier.wdl`; equal-length validation and the deterministic q² glob are present in `q2_incidence.wdl`.
- `git diff --check`

Unavailable or environment-blocked:

- `pytest -q`: `pytest` executable is not installed.
- `python3 -m pytest -q`: the active Python reports `No module named pytest`.
- Full containerized WDL smoke execution: miniwdl cannot access the Docker daemon in this sandbox (`PermissionError: Operation not permitted`).

## Concerns

- The Task 4 report records that the final Docker image build and runtime probes were not completed. The WDLs use the documented local tag `underlier-prevalence:test`, but a successful image build is still required before an end-to-end workflow run.
- The checked-in `.tbi` fixtures are placeholders by design. They exercise WDL pairing and staging paths but are not binary-valid tabix indexes for a production smoke run.
- The pre-existing modification to `.superpowers/sdd/2026-07-31-topmed-wdl-pipeline/task-1-report.md` was left untouched and excluded from Task 5 staging.

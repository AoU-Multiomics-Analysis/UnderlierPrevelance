# TOPMed WDL Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement independently runnable RNA-underlier and VCF/q² WDL workflows plus a conditional main dispatcher using GCT counts and precomputed genotype covariates.

**Architecture:** Keep analytical logic in versioned R/Python scripts and use WDL for staging, validation, branching, and declared outputs. The RNA branch consumes a sample-keyed genotype-PC TSV and never computes genotype PCs; the q² branch consumes indexed VCFs and ClinVar. A shared main workflow conditionally calls either branch.

**Tech Stack:** WDL 1.0, miniwdl validation, R/Bioconductor (`edgeR`, `limma`, `corral`, `PCAtools`, `scran`, `scuttle`, `WGCNA`, `rtracklayer`, `vroom`, `tidyverse`, `magrittr`), bcftools, Python 3 with pandas/numpy, Docker.

## Global Constraints

- GCT input must be standard v1.2 with `Name` and `Description` columns.
- Genotype covariates must be TSV columns `sample_id`, `Genotype_PC1`, `Genotype_PC2`, and so on.
- RNA processing must not invoke PLINK, GDS conversion, LD pruning, or genotype PCA.
- Phenotype-PC count must be selected with `PCAtools::chooseGavishDonoho`.
- Default GD noise variance is `1` under the unit-variance working convention; an explicit override is supported.
- q² processing must accept VCF/index arrays and ClinVar VCF/index inputs.
- Every WDL under `workflows/` must pass `miniwdl check`.

---

### Task 1: Add test fixtures and validation helpers

**Files:**
- Create: `tests/fixtures/counts.gct`
- Create: `tests/fixtures/genotype_covariates.tsv`
- Create: `tests/fixtures/annotation.gff3`
- Create: `tests/fixtures/clinvar.vcf`
- Create: `tests/fixtures/cohort.vcf`
- Create: `tests/test_inputs.py`
- Create: `tests/test_q2_script.py`

**Interfaces:**
- Consumes: synthetic GCT, covariate, annotation, and VCF fixtures.
- Produces: failing tests defining parser errors, sample alignment, covariate-column validation, and q² output schema.

- [ ] **Step 1: Write failing fixture tests**

```python
import pytest


def test_gct_fixture_has_standard_header():
    lines = Path("tests/fixtures/counts.gct").read_text().splitlines()
    assert lines[0] == "#1.2"
    assert lines[2].split("\t")[:2] == ["Name", "Description"]


def test_covariates_use_sample_id_and_genotype_pc_columns():
    header = Path("tests/fixtures/genotype_covariates.tsv").read_text().splitlines()[0].split("\t")
    assert header[:2] == ["sample_id", "Genotype_PC1"]
```

- [ ] **Step 2: Run tests to verify the fixture contract is enforced**

Run: `pytest -q tests/test_inputs.py`

Expected: FAIL because the fixture files and parser/validation implementation do not yet exist.

- [ ] **Step 3: Add minimal valid fixtures and validation assertions**

Use four samples, six genes, two genotype PCs, two protein-coding GFF records, and a cohort/ClinVar VCF pair with at least three shared pathogenic alleles so q² has a nonempty result.

- [ ] **Step 4: Run tests to verify the fixture tests pass**

Run: `pytest -q tests/test_inputs.py`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests
git commit -m "test: add synthetic pipeline fixtures"
```

### Task 2: Implement RNA input conversion and analysis script

**Files:**
- Create: `scripts/convert_gct.py`
- Create: `scripts/run_rna_underlier.R`
- Modify: `tools/lint_r.R`
- Test: `tests/test_inputs.py`

**Interfaces:**
- `convert_gct.py --input counts.gct --output counts.tsv` converts standard GCT to a tab-delimited matrix with `gene_id` followed by sample columns.
- `run_rna_underlier.R --counts counts.tsv --genotype-covariates genotype_covariates.tsv --gencode annotation.gff3 --out-dir out --n-geno-pcs N --phenotype-pc-noise VALUE` emits the RNA outputs listed in the design.
- `convert_gct.py` rejects malformed headers, wrong dimensions, duplicate gene IDs, and nonnumeric count cells.

- [ ] **Step 1: Write failing converter and schema tests**

```python
import subprocess
from pathlib import Path
import pytest


def run_converter(input_path, output_path):
    return subprocess.run(
        ["python3", "scripts/convert_gct.py", "--input", str(input_path), "--output", str(output_path)],
        check=True,
    )


def test_convert_gct_writes_gene_id_and_sample_columns(tmp_path):
    output = tmp_path / "counts.tsv"
    run_converter("tests/fixtures/counts.gct", output)
    header = output.read_text().splitlines()[0].split("\t")
    assert header == ["gene_id", "S1", "S2", "S3", "S4"]


def test_converter_rejects_duplicate_gene_ids(tmp_path):
    bad = tmp_path / "bad.gct"
    bad.write_text(Path("tests/fixtures/counts.gct").read_text().replace("GENE2", "GENE1", 1))
    with pytest.raises(subprocess.CalledProcessError):
        run_converter(bad, tmp_path / "out.tsv")
```

- [ ] **Step 2: Run the focused tests and confirm the expected failures**

Run: `pytest -q tests/test_inputs.py -k 'convert or duplicate'`

Expected: FAIL because `scripts/convert_gct.py` is absent.

- [ ] **Step 3: Implement the converter**

Parse the two GCT metadata lines, verify declared dimensions against the data, require `Name` and `Description`, coerce all count fields to finite numeric values, and write only `gene_id` plus sample columns.

- [ ] **Step 4: Implement the RNA analysis**

Port the supplied R pipeline’s expression filtering, normalization, connectivity QC, log2-CPM adjustment, underlier definitions, and prevalence outputs. Replace genotype-PCA code with strict TSV validation and covariate selection. Run `PCAtools::pca`, use default noise variance `1` under the unit-variance working convention unless an explicit override is supplied, call `PCAtools::chooseGavishDonoho`, clamp only to the available rank after rejecting zero-PC results, and write selected-PC metadata.

- [ ] **Step 5: Run the RNA unit/smoke test**

Run: `python3 scripts/convert_gct.py --input tests/fixtures/counts.gct --output /tmp/topmed-counts.tsv && Rscript scripts/run_rna_underlier.R --counts /tmp/topmed-counts.tsv --genotype-covariates tests/fixtures/genotype_covariates.tsv --gencode tests/fixtures/annotation.gff3 --out-dir /tmp/topmed-rna-test --n-geno-pcs 2`

Expected: exit 0 and output `selected_phenotype_pcs.tsv`, `expr_z_join.tsv.gz`, both underlier tables, and both prevalence tables.

- [ ] **Step 6: Commit**

```bash
git add scripts tools tests
git commit -m "feat: add GCT conversion and RNA underlier analysis"
```

### Task 3: Implement q² VCF filtering and aggregation

**Files:**
- Create: `scripts/filter_variants.sh`
- Create: `scripts/compute_q2_incidence.py`
- Test: `tests/test_q2_script.py`

**Interfaces:**
- `filter_variants.sh --input-dir DIR --output-dir DIR --af-max FLOAT --missing-max FLOAT --qual-min FLOAT --threads INT` writes filtered/indexed VCFs.
- `compute_q2_incidence.py --clinvar FILE --vcf-glob GLOB --out FILE [--genes FILE]` writes the documented q² table.

- [ ] **Step 1: Write failing q² schema and filtering tests**

```python
import subprocess
from pathlib import Path


def run_q2(output):
    subprocess.run(
        ["python3", "scripts/compute_q2_incidence.py", "--clinvar", "tests/fixtures/clinvar.vcf", "--vcf-glob", "tests/fixtures/cohort.vcf", "--out", str(output)],
        check=True,
    )
    return output


def run_filtering(output_dir):
    subprocess.run(
        ["bash", "scripts/filter_variants.sh", "--input-dir", "tests/fixtures", "--output-dir", str(output_dir), "--af-max", "0.01", "--missing-max", "0.1", "--qual-min", "125", "--threads", "1"],
        check=True,
    )
    return output_dir / "all.filtered.vcf.gz"


def read_variant_count(path):
    result = subprocess.run(["bcftools", "view", "-H", str(path)], check=True, capture_output=True, text=True)
    return len(result.stdout.splitlines())


def test_q2_output_contains_required_columns(tmp_path):
    output = tmp_path / "q2.tsv"
    run_q2(output)
    columns = output.read_text().splitlines()[0].split("\t")
    assert columns == ["gene", "n_plp_alleles_clinvar", "n_plp_alleles_observed", "q", "incidence", "carrier_freq"]


def test_filtering_removes_nonbiallelic_and_high_af_sites(tmp_path):
    filtered = run_filtering(tmp_path)
    assert read_variant_count(filtered) == 3
```

- [ ] **Step 2: Run tests to verify they fail for missing scripts**

Run: `pytest -q tests/test_q2_script.py`

Expected: FAIL because the new command-line scripts do not exist.

- [ ] **Step 3: Implement bcftools filtering**

Stage each input VCF and index, concatenate or process deterministically, run `bcftools +fill-tags` for AC/AN/AF/F_MISSING, apply biallelic/AF/missingness/QUAL filters, and create tabix indexes for every emitted VCF.

- [ ] **Step 4: Port q² aggregation**

Retain the supplied ClinVar P/LP small-variant matching, gene assignment, maximum duplicated-row AF behavior, q aggregation, `incidence = q * q`, and carrier-frequency calculation. Validate that the VCF glob matches at least one file.

- [ ] **Step 5: Run q² tests and the full focused test file**

Run: `pytest -q tests/test_q2_script.py`

Expected: PASS with a nonempty q² table and expected filtered-site count.

- [ ] **Step 6: Commit**

```bash
git add scripts tests
git commit -m "feat: add VCF filtering and q2 incidence analysis"
```

### Task 4: Add container definitions

**Files:**
- Modify: `envs/Dockerfile`
- Create: `envs/requirements-q2.txt`
- Create: `envs/bioconda-environment.yml`

**Interfaces:**
- Produces reproducible image tags used by RNA and q² WDL tasks.

- [ ] **Step 1: Add build smoke-test checks**

```bash
docker build -f envs/Dockerfile -t underlier-prevalence:test .
docker run --rm underlier-prevalence:test Rscript -e 'library(PCAtools); stopifnot(is.function(chooseGavishDonoho))'
docker run --rm underlier-prevalence:test bcftools --version
```

- [ ] **Step 2: Run the checks before implementation**

Run the commands above.

Expected: FAIL because the Dockerfile is currently empty and the dependencies are not declared.

- [ ] **Step 3: Implement the image definitions**

Use explicit R/Bioconductor package installation for the RNA dependencies and explicit Python/bcftools installation for q². Keep the WDL task image references in one WDL constant or input so deployments can override them.

- [ ] **Step 4: Build and verify the image**

Run the commands above.

Expected: PASS with both runtime probes succeeding.

- [ ] **Step 5: Commit**

```bash
git add envs
git commit -m "build: package RNA and q2 runtime dependencies"
```

### Task 5: Add the reusable WDL workflows

**Files:**
- Create: `workflows/rna_underlier.wdl`
- Create: `workflows/q2_incidence.wdl`
- Create: `workflows/main.wdl`
- Create: `tests/test_wdl_inputs.json`

**Interfaces:**
- `rna_underlier` consumes `counts_gct`, `genotype_covariates_tsv`, `gencode_gff`, `n_genotype_pcs`, `phenotype_pc_noise`, `conn_z`, `logcpm_drop`, and `threads`.
- `q2_incidence` consumes `genotype_vcfs`, `genotype_vcf_indexes`, `clinvar_vcf`, `clinvar_vcf_index`, filtering parameters, and optional `gene_whitelist`.
- `main` consumes branch flags and branch-specific optional inputs and exposes optional branch outputs.

- [ ] **Step 1: Write WDL validation tests and input JSON**

```bash
miniwdl check workflows/rna_underlier.wdl
miniwdl check workflows/q2_incidence.wdl
miniwdl check workflows/main.wdl
```

Expected: FAIL because the workflow files do not exist.

- [ ] **Step 2: Implement RNA WDL tasks**

Add conversion, analysis, and output-glob tasks. Declare runtime Docker image, CPU, memory, disk, and retry settings. Stage the converter output into the analysis task and emit selected-PC metadata plus all prevalence artifacts.

- [ ] **Step 3: Implement q² WDL tasks**

Add VCF staging, filtering, and q² aggregation tasks. Validate equal VCF/index array lengths in the command section and pass a deterministic glob to the Python script.

- [ ] **Step 4: Implement the main conditional dispatcher**

Use `if (run_rna)` and `if (run_q2)` call blocks. Require branch inputs inside the enabled call’s task command and expose optional outputs at the top level.

- [ ] **Step 5: Run miniwdl validation**

Run: `for wdl in workflows/*.wdl; do miniwdl check "$wdl"; done`

Expected: PASS for all three workflow files.

- [ ] **Step 6: Commit**

```bash
git add workflows tests/test_wdl_inputs.json
git commit -m "feat: add RNA, q2, and main WDL workflows"
```

### Task 6: Add CI and user documentation

**Files:**
- Modify: `.github/workflows/wdl-validation.yml`
- Modify: `README.md`
- Create: `docs/running.md`

**Interfaces:**
- CI validates every WDL and runs lightweight synthetic tests.
- Documentation explains inputs, outputs, the analysis interpretation, and a complete local smoke-test command.

- [ ] **Step 1: Add documentation/CI checks**

```bash
pytest -q
for wdl in workflows/*.wdl; do miniwdl check "$wdl"; done
git diff --check
```

- [ ] **Step 2: Run the checks before modifying CI/docs**

Expected: pytest and WDL checks fail or are incomplete because the workflow and test suite are not yet present.

- [ ] **Step 3: Implement CI and documentation**

Document the GCT v1.2 format, genotype-covariate TSV schema, exact RNA/q² output meanings, branch flags, GD noise override, and the distinction between haplo and strict calls. Extend CI to install the lightweight test dependencies and run `pytest` plus the existing `miniwdl check` loop.

- [ ] **Step 4: Run complete verification**

Run: `pytest -q && for wdl in workflows/*.wdl; do miniwdl check "$wdl"; done && git diff --check`

Expected: PASS with zero test failures, all WDLs valid, and no whitespace errors.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/wdl-validation.yml README.md docs/running.md
git commit -m "docs: document and validate TOPMed WDL workflows"
```

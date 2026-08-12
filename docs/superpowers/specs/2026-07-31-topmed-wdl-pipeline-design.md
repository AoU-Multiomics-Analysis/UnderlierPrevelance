# TOPMed Underlier and q² WDL Pipeline Design

## Goal

Convert the supplied TOPMed pipeline into reproducible WDL workflows in this repository, preserving both deliverables while changing the RNA workflow to accept precomputed genotype covariates and select phenotype PCs with the Gavish–Donoho method from `PCAtools`.

## Scope

The implementation contains two independently runnable workflows and one optional dispatcher:

1. `rna_underlier.wdl` produces RNA underlier calls and per-gene prevalence from a standard GCT count matrix.
2. `q2_incidence.wdl` filters genotype VCFs and computes ClinVar P/LP recessive-incidence (`q²`) tables.
3. `main.wdl` conditionally runs either or both workflows using `run_rna` and `run_q2` flags.

The RNA workflow does not compute genotype PCs, convert genotype VCFs to PLINK/GDS, perform LD pruning, or inspect the VCF. Sample identifiers are assumed to match between the GCT and the precomputed genotype-covariate TSV and are validated at runtime.

## Architecture and data flow

### RNA underlier workflow

Inputs:

- Standard GCT v1.2 count file: `#1.2`, dimensions, `Name`, `Description`, then one column per sample.
- Genotype-covariate TSV with `sample_id` followed by numeric `Genotype_PC1`, `Genotype_PC2`, … columns.
- GENCODE GFF/GFF3 annotation containing protein-coding gene records and gene symbols.
- Runtime parameters for connectivity QC, log2-CPM drop, genotype-PC count, GD noise estimation, and threads.

Processing:

1. Parse GCT into a gene-by-sample count matrix.
2. Validate unique sample IDs and exact sample alignment with the genotype-covariate table.
3. Restrict to protein-coding GENCODE genes.
4. Apply the supplied expression filtering and Freeman–Tukey normalization.
5. Remove connectivity outlier samples using the configured connectivity Z threshold.
6. Re-filter and re-normalize the QC-passing samples.
7. Run `PCAtools::pca` on the normalized expression matrix.
8. Select the number of phenotype PCs with `PCAtools::chooseGavishDonoho(x, var.explained = pca$sdev^2, noise = ...)`.
9. Regress out the selected phenotype PCs plus the requested genotype-PC columns.
10. Compute standardized residual expression z-scores, join them to PC-adjusted log2-CPM values, and emit haplo and strict underlier definitions.

The default GD noise variance is `1`, corresponding to the unit-variance working convention for the expression measurements. An explicit `phenotype_pc_noise` override is supported for analyses that have a validated noise estimate.

Outputs:

- Selected phenotype-PC count and GD diagnostic metadata.
- `expr_z_join.tsv.gz`.
- `underliers_haplo.tsv.gz` and the strict z-cutoff underlier table.
- `rna_outlier_prevalence_per_gene_haplo.tsv` and the strict prevalence TSV.
- Optional RDS intermediates for restartability and inspection.

### q² incidence workflow

Inputs:

- Array of genotype VCF/VCF.GZ files and an equally sized array of tabix indexes.
- ClinVar VCF/VCF.GZ and its tabix index.
- Variant filtering parameters: AF maximum, missingness maximum, QUAL minimum, and thread count.
- Optional gene-whitelist file.

Processing:

1. Stage the VCFs and indexes in a deterministic directory.
2. Use bcftools to retain biallelic variants, recompute AC/AN/AF/F_MISSING, and apply AF, missingness, and QUAL thresholds.
3. Index filtered VCF output.
4. Run the supplied q² aggregation logic against ClinVar pathogenic/likely-pathogenic small variants.

Outputs:

- Filtered/indexed VCFs.
- `q2_incidence.tsv` with gene, ClinVar allele counts, observed allele counts, q, incidence, and carrier frequency.

### Main dispatcher

`main.wdl` exposes `run_rna` and `run_q2`. Branch-specific inputs are optional at the top level and are required by WDL runtime validation when their branch is enabled. Branch outputs are optional in the main workflow and concrete in the sub-workflows.

## Runtime packaging

Use focused container environments:

- An R/Bioconductor image containing `edgeR`, `limma`, `corral`, `PCAtools`, `scran`, `scuttle`, `WGCNA`, `rtracklayer`, `vroom`, `tidyverse`, and `magrittr` for RNA tasks.
- A bcftools/Python image containing bcftools with `+fill-tags`, Python, pandas, and numpy for q² tasks.

The repository Dockerfile and workflow tasks will pin explicit base-image/package versions where practical. Analytical scripts remain in the repository and are invoked by WDL tasks, so command-line behavior and output names are version controlled.

## Validation and error handling

Tasks fail with actionable messages for:

- Invalid GCT headers, dimensions, duplicate gene IDs, or nonnumeric counts.
- Duplicate, missing, or nonmatching RNA/genotype sample IDs.
- Missing or nonnumeric required genotype-PC columns.
- Invalid requested genotype-PC count or GD result outside the available observation rank.
- VCF/index array length mismatches or no staged VCF inputs.
- Missing ClinVar annotations or no eligible ClinVar records.
- Missing expected output files after a task completes.

The workflows do not silently regenerate genotype PCs and do not silently drop unmatched samples.

## Testing strategy

Repository tests will include:

- Tiny synthetic GCT parsing and dimension-validation fixtures.
- Sample-alignment and genotype-covariate schema tests.
- A small PCA fixture that exercises the `PCAtools::chooseGavishDonoho` call and records the selected count.
- Synthetic VCF/ClinVar fixtures for filtering and q² aggregation.
- `miniwdl check` over every file under `workflows/`.

CI will run WDL syntax checks and lightweight tests. A documented local/container smoke test will run both complete branches on the synthetic fixtures; production-sized RNA and VCF inputs remain execution-backend dependent.

## Acceptance criteria

- Both sub-workflows validate with `miniwdl check`.
- The RNA workflow has no PLINK, GDS, genotype-PCA, or LD-pruning step.
- The RNA workflow consumes the agreed `sample_id`/`Genotype_PC*` TSV contract.
- Phenotype-PC count is selected by `PCAtools::chooseGavishDonoho` and emitted as an output.
- GCT input is accepted directly and converted inside the workflow.
- VCF input is accepted directly by the q² workflow, with required indexes and explicit filtering.
- The dispatcher can run either branch independently or both together.
- Synthetic tests cover parser, alignment, PC-selection wiring, q² filtering, and WDL validation.

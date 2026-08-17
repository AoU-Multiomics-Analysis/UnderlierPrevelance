# TOPMed RNA underlier prevalence and q² incidence workflows

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/UnderlierPrevelance/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->


This repository packages two reproducible WDL 1.0 analyses and a conditional
dispatcher:

- RNA underlier prevalence from a GCT expression matrix and precomputed
  genotype covariates.
- Gene-level recessive q² incidence from filtered cohort VCFs and ClinVar
  pathogenic/likely-pathogenic (P/LP) alleles.

`workflows/main.wdl` runs either branch or both. It is a dispatcher only; the
branch workflows can also be run directly.

## Quick start

For production, the WDL defaults to the image published at
`ghcr.io/aou-multiomics-analysis/underlierprevelance:main`. To run locally,
build a test image and override `docker_image` in the input JSON:

```bash
docker build -f envs/Dockerfile -t underlier-prevalence:test .
miniwdl run workflows/main.wdl -i tests/test_wdl_inputs.json
```

The complete input contracts, branch-specific commands, output schemas, and
runtime limitations are in [docs/running.md](docs/running.md).

## Inputs

RNA requires a standard GCT v1.2 count matrix, a protein-coding GENCODE
GFF/GFF3, and a precomputed genotype-covariate TSV. The workflow does not
derive genotype PCs, run PLINK, convert to GDS, or perform LD pruning.

q² requires paired cohort VCF/VCF.GZ and `.tbi` index arrays plus a paired
ClinVar VCF/VCF.GZ and index. The cohort array must contain chromosome or
region shards with identical sample IDs in identical order; independent
cohorts cannot be concatenated through this interface. It filters variants
before matching ClinVar P/LP small variants.

## Outputs

The RNA branch writes residual-expression z-scores, haplo-only calls, z-score-only
calls, haplo-plus-z-score intersection calls, their per-gene prevalence tables,
and metadata documenting Gavish–Donoho phenotype-PC selection. The q² branch
writes the filtered/indexed VCF and a per-gene table with `q`, `incidence = q²`,
and the supplied
`carrier_freq = 1 - (1 - q)²` semantics. See [docs/running.md](docs/running.md)
for precise definitions.

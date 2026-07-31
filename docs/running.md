# Running the TOPMed workflows

The repository provides two WDL 1.0 workflows and a dispatcher:

- `workflows/rna_underlier.wdl` computes RNA underlier calls and per-gene underlier prevalence.
- `workflows/q2_incidence.wdl` filters cohort VCFs and computes per-gene q² incidence from ClinVar P/LP small variants.
- `workflows/main.wdl` conditionally runs one or both branches.

All WDL tasks use the configurable `docker_image` input, which defaults to `underlier-prevalence:test`. Build that image before local execution:

```bash
docker build -f envs/Dockerfile -t underlier-prevalence:test .
```

## Main workflow flags and inputs

`main.run_rna` and `main.run_q2` select the RNA and q² branches respectively. Both default to `false`. Set one flag to run its branch independently, or set both to `true` to run both branches. When a branch is enabled, all of its required `main.rna_*` or `main.q2_*` inputs must be supplied; optional top-level inputs are selected by the enabled WDL call and are therefore not silently substituted.

The runnable synthetic example enables both flags in `tests/test_wdl_inputs.json`. A complete local smoke test is:

```bash
docker build -f envs/Dockerfile -t underlier-prevalence:test .
miniwdl run workflows/main.wdl -i tests/test_wdl_inputs.json
```

For a branch-only invocation, use the matching workflow and an input JSON whose keys are namespaced by the workflow name (for example, `rna_underlier.counts_gct` or `q2_incidence.genotype_vcfs`).

## RNA underlier prevalence

### Required input contract

`rna_counts_gct` must be a standard GCT v1.2 file:

1. Line 1 is exactly `#1.2`.
2. Line 2 contains tab-delimited integer counts of genes and samples.
3. Line 3 starts with `Name` and `Description`, followed by one unique sample identifier per count column.
4. Each subsequent row contains a unique gene identifier, a description, and finite, non-negative numeric counts.

The workflow converts this GCT to its internal `gene_id`-by-sample TSV. It matches genes to protein-coding `gene` records in `rna_gencode_gff`; the GFF or GFF3 must include `gene_id` and `gene_type` or `gene_biotype` attributes.

`rna_genotype_covariates_tsv` is a tab-delimited table with this schema:

```text
sample_id	Genotype_PC1	Genotype_PC2	...
SAMPLE_A	0.12	-1.03	...
```

`sample_id` is the first, unique column. Every following column must be named `Genotype_PC1`, `Genotype_PC2`, and so on, and contain finite numeric values. Its sample IDs must exactly match the GCT sample IDs. Set `rna_n_genotype_pcs` to the number of consecutive genotype-PC columns to use. These PCs are precomputed inputs: this workflow intentionally performs no genotype-PC computation, VCF processing, PLINK conversion, GDS conversion, or LD pruning.

### Analysis and interpretation

After protein-coding filtering, expression filtering, Freeman–Tukey normalization, and connectivity QC, the workflow runs `PCAtools::pca` on the RNA expression data. It selects phenotype PCs with `PCAtools::chooseGavishDonoho`. By default, the noise variance is the median of the lower half of the PCA variance spectrum. Supply `rna_phenotype_pc_noise` to override that default with a validated, non-negative noise variance; the output metadata records whether the source was `lower_half_median` or `override`.

The selected phenotype PCs and requested precomputed genotype PCs are regressed from the normalized expression used for z-scores. The selected phenotype PCs are also regressed from log2-CPM before calls are made. The `selected_phenotype_pcs.tsv` metadata records the Gavish–Donoho result, available rank, noise value/source, PC columns, design ranks, residual degrees of freedom, and post-QC sample/gene counts.

For each gene, an observation must first have PC-adjusted log2-CPM less than that gene's mean PC-adjusted log2-CPM minus `rna_logcpm_drop` (default `1`). This expression-drop definition is the **haplo** call. The **strict** call is the same haplo call with an additional residual expression-z-score threshold: `expression_zscore < -3` by default. The RNA WDL uses the script default strict cutoff, so its output tag is `z_-3`.

### RNA outputs

The RNA branch exposes these outputs (concrete in `rna_underlier.wdl`, optional from `main.wdl` when `run_rna` is false):

| Output | Meaning |
| --- | --- |
| `converted_counts_tsv` | Validated GCT converted to `gene_id` plus sample columns. |
| `selected_phenotype_pcs.tsv` | Gavish–Donoho PC-selection and residualization metadata. |
| `expr_z_join.tsv.gz` | One row per tested gene/sample with residual `expression_zscore` and PC-adjusted `expression_logcpm`. |
| `underliers_haplo.tsv.gz` | Calls meeting the expression-drop criterion. |
| `underliers_z_-3.tsv.gz` | Calls meeting both the haplo criterion and `expression_zscore < -3`. |
| `rna_outlier_prevalence_per_gene_haplo.tsv` | Per-gene prevalence for haplo calls. |
| `rna_outlier_prevalence_per_gene_z_-3.tsv` | Per-gene prevalence for strict calls. |

Each prevalence table contains `gene_id`, version-stripped `gene_nv`, `symbol`, `n_underlier_subjects`, `n_subjects_tested`, and `rna_outlier_prevalence = n_underlier_subjects / n_subjects_tested`. Subjects are counted once per gene even if a table could otherwise contain duplicate rows.

## q² incidence

### Required input contract

Supply equally sized `q2_genotype_vcfs` and `q2_genotype_vcf_indexes` arrays. Each index filename must be the matching VCF filename plus `.tbi`. Supply the same paired contract for `q2_clinvar_vcf` and `q2_clinvar_vcf_index`. Inputs may be VCF or VCF.GZ files. `q2_gene_whitelist` is optional and accepts one gene symbol per non-comment line.

Before aggregation, the workflow retains biallelic variants and recomputes `AC`, `AN`, `AF`, and `F_MISSING`; it then applies `q2_af_max` (default `0.01`), `q2_missing_max` (default `0.1`), and `q2_qual_min` (default `100`). It writes the resulting `all.filtered.vcf.gz` and tabix index.

### q² meaning and outputs

ClinVar contributes only pathogenic or likely-pathogenic, non-conflicting small variants assigned to a gene. For matching cohort alleles, duplicate observed records use the maximum allele frequency. Alleles with observed AF of zero are not counted as observed. For each gene:

```text
q            = sum of observed P/LP allele frequencies
incidence    = q²
carrier_freq = 1 - (1 - q)²
```

`carrier_freq` follows the supplied pipeline semantics: it is the probability of carrying at least one pathogenic allele under the model, including affected homozygotes. It is not a heterozygote-only carrier frequency.

`q2_incidence.tsv` has these columns:

| Column | Meaning |
| --- | --- |
| `gene` | ClinVar gene symbol. |
| `n_plp_alleles_clinvar` | Number of eligible unique ClinVar P/LP alleles for the gene. |
| `n_plp_alleles_observed` | Number of those alleles with observed AF greater than zero. |
| `q` | Sum of observed allele frequencies. |
| `incidence` | `q * q`. |
| `carrier_freq` | `1 - (1 - q) * (1 - q)`, including affected homozygotes. |

The branch outputs `filtered_vcf`, `filtered_vcf_index`, and `q2_incidence_tsv`; the dispatcher names them `q2_filtered_vcf`, `q2_filtered_vcf_index`, and `q2_incidence_tsv`.

## Validation, CI, and runtime limitations

GitHub Actions installs `pytest`, `bcftools`, `r-base`, and `miniwdl`; it runs the synthetic tests followed by validation of every WDL under `workflows/`. Run the equivalent local checks with:

```bash
python -m pip install pytest miniwdl
pytest -q
for wdl in workflows/*.wdl; do miniwdl check "$wdl"; done
git diff --check
```

The lightweight suite does not run the full RNA analysis unless an R runtime with `edgeR`, `limma`, `corral`, `PCAtools`, `scran`, `WGCNA`, and `rtracklayer` is available; that smoke test is automatically skipped otherwise. The full container smoke test above requires Docker, a container-capable local backend, and enough CPU, memory, and disk for the WDL task runtime requests. Production execution additionally requires real GCT, GENCODE, cohort VCF/index, and ClinVar VCF/index inputs and an execution backend configured to pull the selected image.

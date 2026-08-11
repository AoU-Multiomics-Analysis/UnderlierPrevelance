# Running the TOPMed workflows

The repository provides two WDL 1.0 workflows and a dispatcher:

- `workflows/rna_underlier.wdl` computes RNA underlier calls and per-gene underlier prevalence.
- `workflows/q2_incidence.wdl` filters cohort VCFs and computes per-gene q² incidence from ClinVar P/LP small variants.
- `workflows/main.wdl` conditionally runs one or both branches.

All WDL tasks use the configurable `docker_image` input, which defaults to `underlier-prevalence:test`. Build that image before local execution:

```bash
docker build -f envs/Dockerfile -t underlier-prevalence:test .
```

The production image is built for Linux/amd64, matching the current Bioconda
builds for the R/Bioconductor stack. On Apple Silicon, use Docker's amd64
emulation explicitly:

```bash
docker build --platform linux/amd64 -f envs/Dockerfile -t underlier-prevalence:test .
```

The image installs R, Bioconductor, Python, and bcftools from the conda-forge
and Bioconda channels. This avoids the previous Docker build step that
compiled the R/Bioconductor stack from source with BiocManager.

The runtime image owns the analysis scripts. The WDL inputs should not include
`convert_gct.py`, `run_rna_underlier.R`, `filter_variants.sh`, or
`compute_q2_incidence.py`; the tasks invoke the copies packaged at
`/opt/underlier-prevalence/scripts/`.

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

GCT conversion is streamed and validated one gene row at a time, then atomically renamed into place only after the declared dimensions and every row pass validation. Converter memory therefore scales with the sample header, one numeric row, and the set of gene IDs needed for duplicate detection rather than with the full count matrix. Every WDL task requests 32 GiB memory and 128 GiB of local disk by default; this includes the converter and downstream R analysis so large intermediate matrices and staged VCFs have consistent capacity.

`rna_genotype_covariates_tsv` is a tab-delimited table with this schema:

```text
sample_id	Genotype_PC1	Genotype_PC2	...
SAMPLE_A	0.12	-1.03	...
```

The table must contain one unique `sample_id` column, either first or last. The PC columns may be named `Genotype_PC1`, `Genotype_PC2`, ... or `GENETICPC1`, `GENETICPC2`, ...; they must be numbered consecutively and contain finite numeric values. Sample IDs must exactly match the GCT sample IDs. `rna_n_genotype_pcs` is optional: when supplied, it selects that many consecutive genotype-PC columns; when omitted, all genotype-PC columns are used. These PCs are precomputed inputs: this workflow intentionally performs no genotype-PC computation, VCF processing, PLINK conversion, GDS conversion, or LD pruning.

### Analysis and interpretation

After protein-coding filtering, expression filtering, Freeman–Tukey normalization, and connectivity QC, the workflow runs `PCAtools::pca` on the RNA expression data. It selects phenotype PCs with `PCAtools::chooseGavishDonoho`. By default, the noise variance is the median of the lower half of the PCA variance spectrum. Supply `rna_phenotype_pc_noise` to override that default with a validated, non-negative noise variance; the output metadata records whether the source was `lower_half_median` or `override`.

The selected phenotype PCs and requested precomputed genotype PCs are regressed from the normalized expression used for z-scores. The selected phenotype PCs are also regressed from log2-CPM before calls are made. The `selected_phenotype_pcs.tsv` metadata records the Gavish–Donoho result, available rank, noise value/source, PC columns, design ranks, residual degrees of freedom, and post-QC sample/gene counts.

For each gene, an observation must first have PC-adjusted log2-CPM less than that gene's mean PC-adjusted log2-CPM minus `rna_logcpm_drop` (default `1`). This expression-drop definition is the **haplo** call. The workflow then evaluates strict calls at configurable residual expression-z-score cutoffs. By default, `rna_z_cutoffs`/`z_cutoffs` contains `-1` through `-10` in increments of `1`, producing one underlier and prevalence output pair per cutoff. The thresholding loop reuses the same adjusted expression and z-score results; it does not rerun PCA or residualization.

### RNA outputs

At the `rna_underlier.wdl` interface, the RNA branch exposes named `converted_counts_tsv`, `selected_pc_metadata`, and `expr_z_join` outputs plus three `Array[File]` outputs: `underlier_artifacts`, `prevalence_artifacts`, and `all_rna_artifacts`. The named outputs map to `counts.tsv`, `selected_phenotype_pcs.tsv`, and `expr_z_join.tsv.gz`, respectively. `all_rna_artifacts` is the complete `glob("rna_outputs/*")` collection, so it includes all RNA output files rather than only the two underlier and two prevalence artifacts. The following files are members of the narrower artifact arrays; they are not individually named WDL output fields.

| Artifact array | File | Meaning |
| --- | --- |
| `underlier_artifacts` | `underliers_haplo.tsv.gz` | Calls meeting the expression-drop criterion. |
| `underlier_artifacts` | `underliers_z_<cutoff>.tsv.gz` | Calls meeting both the haplo criterion and `expression_zscore < cutoff`; defaults are `-1` through `-10`. |
| `prevalence_artifacts` | `rna_outlier_prevalence_per_gene_haplo.tsv` | Per-gene prevalence for haplo calls. |
| `prevalence_artifacts` | `rna_outlier_prevalence_per_gene_z_<cutoff>.tsv` | Per-gene prevalence for each configured strict cutoff. |

Each prevalence table contains `gene_id`, version-stripped `gene_nv`, `symbol`, `n_underlier_subjects`, `n_subjects_tested`, and `rna_outlier_prevalence = n_underlier_subjects / n_subjects_tested`. Subjects are counted once per gene even if a table could otherwise contain duplicate rows.

`main.wdl` exposes the same RNA results through optional aliases when `main.run_rna` is `true`; each alias is absent when that branch is disabled:

| `rna_underlier.wdl` output | `main.wdl` RNA output alias |
| --- | --- |
| `converted_counts_tsv` | `rna_converted_counts_tsv` |
| `selected_pc_metadata` | `rna_selected_pc_metadata` |
| `expr_z_join` | `rna_expr_z_join` |
| `underlier_artifacts` | `rna_underlier_artifacts` |
| `prevalence_artifacts` | `rna_prevalence_artifacts` |
| `all_rna_artifacts` (`rna_outputs/*`) | `rna_all_artifacts` |

## q² incidence

### Required input contract

Supply equally sized `q2_genotype_vcfs` and `q2_genotype_vcf_indexes` arrays. Each index filename must be the matching VCF filename plus `.tbi`. The genotype VCF array is a set of chromosome or region shards from one cohort, not a collection of independent cohorts: every genotype-bearing shard must contain identical sample IDs in identical order. Filtering fails before `bcftools concat` if a shard differs. Supply the same paired contract for `q2_clinvar_vcf` and `q2_clinvar_vcf_index`. Inputs may be VCF or VCF.GZ files. `q2_gene_whitelist` is optional and accepts one gene symbol per non-comment line.

Before aggregation, the workflow retains biallelic variants and recomputes `AC`, `AN`, `AF`, and `F_MISSING`; it then applies `q2_af_max` (default `0.01`), `q2_missing_max` (default `0.1`), and `q2_qual_min` (default `100`). Both AF and missingness thresholds must be within `[0, 1]`. It writes the resulting `all.filtered.vcf.gz` and tabix index. Filtering temporarily normalizes and indexes every shard before concatenation, so local-disk sizing must accommodate those staged copies plus the filtered output.

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

GitHub Actions installs `pytest`, `bcftools`, `r-base` (which provides `Rscript`), and `miniwdl`; it runs the synthetic tests followed by validation of every WDL under `workflows/`.

For local validation, `pytest` and `miniwdl` are required Python commands. The full lightweight `pytest` suite also requires native `bcftools` and `Rscript`: VCF-filter tests invoke `bcftools`, and dependency-free RNA interface tests invoke `Rscript`. Install the two Python commands, install `bcftools` and R through your operating system or environment manager, and confirm all four commands are available before running the complete suite:

```bash
python -m pip install pytest miniwdl
command -v pytest
command -v miniwdl
command -v bcftools
command -v Rscript
pytest -q
for wdl in workflows/*.wdl; do miniwdl check "$wdl"; done
git diff --check
```

The lightweight suite does not run the full RNA analysis unless an R runtime with `edgeR`, `limma`, `corral`, `PCAtools`, `scran`, `WGCNA`, and `rtracklayer` is available; that smoke test is automatically skipped otherwise. The full container smoke test above requires Docker, a container-capable local backend, and enough CPU, memory, and disk for the WDL task runtime requests. Production execution additionally requires real GCT, GENCODE, cohort VCF/index, and ClinVar VCF/index inputs and an execution backend configured to pull the selected image.

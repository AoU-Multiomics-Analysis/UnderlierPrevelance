version 1.0

task ConvertGct {
  input {
    File counts_gct
    String docker_image
  }

  command <<<
    set -euo pipefail

    mkdir -p converted
    python3 /opt/underlier-prevalence/scripts/convert_gct.py \
      --input "~{counts_gct}" \
      --output converted/counts.tsv
    test -s converted/counts.tsv
  >>>

  output {
    File counts_tsv = "converted/counts.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    maxRetries: 2
  }
}

task RunRnaUnderlier {
  input {
    File counts_tsv
    File genotype_covariates_tsv
    File gencode_gff
    Int n_genotype_pcs
    Float? phenotype_pc_noise
    Float conn_z
    Float logcpm_drop
    Int threads
    String docker_image
  }

  command <<<
    set -euo pipefail

    mkdir -p rna_outputs

    phenotype_pc_noise="~{default="" phenotype_pc_noise}"
    phenotype_pc_noise_args=()
    if [[ -n "$phenotype_pc_noise" ]]; then
      phenotype_pc_noise_args=(--phenotype-pc-noise "$phenotype_pc_noise")
    fi

    Rscript /opt/underlier-prevalence/scripts/run_rna_underlier.R \
      --counts "~{counts_tsv}" \
      --genotype-covariates "~{genotype_covariates_tsv}" \
      --gencode "~{gencode_gff}" \
      --out-dir rna_outputs \
      --n-geno-pcs "~{n_genotype_pcs}" \
      --connectivity-z "~{conn_z}" \
      --logcpm-drop "~{logcpm_drop}" \
      --threads "~{threads}" \
      "${phenotype_pc_noise_args[@]}"

    test -s rna_outputs/selected_phenotype_pcs.tsv
    test -s rna_outputs/expr_z_join.tsv.gz

    shopt -s nullglob
    underlier_artifacts=(rna_outputs/underliers_*.tsv.gz)
    prevalence_artifacts=(rna_outputs/rna_outlier_prevalence_per_gene_*.tsv)
    if [[ ${#underlier_artifacts[@]} -ne 2 ]]; then
      echo "error: RNA analysis did not emit both underlier artifacts" >&2
      exit 2
    fi
    if [[ ${#prevalence_artifacts[@]} -ne 2 ]]; then
      echo "error: RNA analysis did not emit both prevalence artifacts" >&2
      exit 2
    fi
  >>>

  output {
    File selected_pc_metadata = "rna_outputs/selected_phenotype_pcs.tsv"
    File expr_z_join = "rna_outputs/expr_z_join.tsv.gz"
    Array[File] underlier_artifacts = glob("rna_outputs/underliers_*.tsv.gz")
    Array[File] prevalence_artifacts = glob("rna_outputs/rna_outlier_prevalence_per_gene_*.tsv")
    Array[File] all_rna_artifacts = glob("rna_outputs/*")
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "16 GiB"
    disks: "local-disk 100 HDD"
    maxRetries: 2
  }
}

workflow rna_underlier {
  input {
    File counts_gct
    File genotype_covariates_tsv
    File gencode_gff
    Int n_genotype_pcs
    Float? phenotype_pc_noise
    Float conn_z = -3.0
    Float logcpm_drop = 1.0
    Int threads = 1
    String docker_image = "underlier-prevalence:test"
  }

  call ConvertGct {
    input:
      counts_gct = counts_gct,
      docker_image = docker_image
  }

  call RunRnaUnderlier {
    input:
      counts_tsv = ConvertGct.counts_tsv,
      genotype_covariates_tsv = genotype_covariates_tsv,
      gencode_gff = gencode_gff,
      n_genotype_pcs = n_genotype_pcs,
      phenotype_pc_noise = phenotype_pc_noise,
      conn_z = conn_z,
      logcpm_drop = logcpm_drop,
      threads = threads,
      docker_image = docker_image
  }

  output {
    File converted_counts_tsv = ConvertGct.counts_tsv
    File selected_pc_metadata = RunRnaUnderlier.selected_pc_metadata
    File expr_z_join = RunRnaUnderlier.expr_z_join
    Array[File] underlier_artifacts = RunRnaUnderlier.underlier_artifacts
    Array[File] prevalence_artifacts = RunRnaUnderlier.prevalence_artifacts
    Array[File] all_rna_artifacts = RunRnaUnderlier.all_rna_artifacts
  }
}

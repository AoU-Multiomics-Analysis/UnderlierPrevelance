version 1.0

import "rna_underlier.wdl" as RNA
import "q2_incidence.wdl" as Q2

workflow main {
  input {
    Boolean run_rna = false
    Boolean run_q2 = false
    String docker_image = "underlier-prevalence:test"

    File? rna_counts_gct
    File? rna_genotype_covariates_tsv
    File? rna_gencode_gff
    Int? rna_n_genotype_pcs
    Float? rna_phenotype_pc_noise
    Float rna_conn_z = -3.0
    Float rna_logcpm_drop = 1.0
    Int rna_threads = 1

    Array[File]? q2_genotype_vcfs
    Array[File]? q2_genotype_vcf_indexes
    File? q2_clinvar_vcf
    File? q2_clinvar_vcf_index
    Float q2_af_max = 0.01
    Float q2_missing_max = 0.1
    Float q2_qual_min = 100.0
    Int q2_threads = 1
    File? q2_gene_whitelist
  }

  if (run_rna) {
    call RNA.rna_underlier as RunRna {
      input:
        counts_gct = select_first([rna_counts_gct]),
        genotype_covariates_tsv = select_first([rna_genotype_covariates_tsv]),
        gencode_gff = select_first([rna_gencode_gff]),
        n_genotype_pcs = select_first([rna_n_genotype_pcs]),
        phenotype_pc_noise = rna_phenotype_pc_noise,
        conn_z = rna_conn_z,
        logcpm_drop = rna_logcpm_drop,
        threads = rna_threads,
        docker_image = docker_image
    }
  }

  if (run_q2) {
    call Q2.q2_incidence as RunQ2 {
      input:
        genotype_vcfs = select_first([q2_genotype_vcfs]),
        genotype_vcf_indexes = select_first([q2_genotype_vcf_indexes]),
        clinvar_vcf = select_first([q2_clinvar_vcf]),
        clinvar_vcf_index = select_first([q2_clinvar_vcf_index]),
        af_max = q2_af_max,
        missing_max = q2_missing_max,
        qual_min = q2_qual_min,
        threads = q2_threads,
        gene_whitelist = q2_gene_whitelist,
        docker_image = docker_image
    }
  }

  output {
    File? rna_converted_counts_tsv = RunRna.converted_counts_tsv
    File? rna_selected_pc_metadata = RunRna.selected_pc_metadata
    File? rna_expr_z_join = RunRna.expr_z_join
    Array[File]? rna_underlier_artifacts = RunRna.underlier_artifacts
    Array[File]? rna_prevalence_artifacts = RunRna.prevalence_artifacts
    Array[File]? rna_all_artifacts = RunRna.all_rna_artifacts

    File? q2_filtered_vcf = RunQ2.filtered_vcf
    File? q2_filtered_vcf_index = RunQ2.filtered_vcf_index
    File? q2_incidence_tsv = RunQ2.q2_incidence_tsv
  }
}

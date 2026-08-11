version 1.0

task StageVcfPairs {
  input {
    Array[File] genotype_vcfs
    Array[File] genotype_vcf_indexes
    File clinvar_vcf
    File clinvar_vcf_index
    String docker_image
  }

  command <<<
    set -euo pipefail

    if [[ ~{length(genotype_vcfs)} -eq 0 ]]; then
      echo "error: at least one genotype VCF/index pair is required" >&2
      exit 2
    fi
    if [[ ~{length(genotype_vcfs)} -ne ~{length(genotype_vcf_indexes)} ]]; then
      echo "error: genotype VCF and index arrays must have equal lengths" >&2
      exit 2
    fi

    mkdir -p staged/genotype staged/clinvar
    vcf_list="~{write_lines(genotype_vcfs)}"
    index_list="~{write_lines(genotype_vcf_indexes)}"

    ordinal=0
    while IFS=$'\t' read -r vcf index; do
      vcf_name=$(basename "$vcf")
      index_name=$(basename "$index")
      if [[ "$index_name" != "$vcf_name.tbi" ]]; then
        echo "error: genotype VCF/index pair does not match: $vcf_name and $index_name" >&2
        exit 2
      fi

      case "$vcf_name" in
        *.vcf.gz) suffix=.vcf.gz ;;
        *.vcf) suffix=.vcf ;;
        *)
          echo "error: genotype input is not a VCF or VCF.GZ: $vcf_name" >&2
          exit 2
          ;;
      esac

      prefix=$(printf '%04d' "$ordinal")
      cp -L "$vcf" "staged/genotype/${prefix}${suffix}"
      cp -L "$index" "staged/genotype/${prefix}${suffix}.tbi"
      ordinal=$((ordinal + 1))
    done < <(paste "$vcf_list" "$index_list")

    clinvar_name=$(basename "~{clinvar_vcf}")
    clinvar_index_name=$(basename "~{clinvar_vcf_index}")
    if [[ "$clinvar_index_name" != "$clinvar_name.tbi" ]]; then
      echo "error: ClinVar VCF/index pair does not match: $clinvar_name and $clinvar_index_name" >&2
      exit 2
    fi

    case "$clinvar_name" in
      *.vcf.gz) clinvar_suffix=.vcf.gz ;;
      *.vcf) clinvar_suffix=.vcf ;;
      *)
        echo "error: ClinVar input is not a VCF or VCF.GZ: $clinvar_name" >&2
        exit 2
        ;;
    esac

    staged_clinvar_vcf="staged/clinvar/clinvar${clinvar_suffix}"
    staged_clinvar_index="${staged_clinvar_vcf}.tbi"
    cp -L "~{clinvar_vcf}" "$staged_clinvar_vcf"
    cp -L "~{clinvar_vcf_index}" "$staged_clinvar_index"
    printf '%s\n' "$staged_clinvar_vcf" > staged_clinvar_vcf.path
    printf '%s\n' "$staged_clinvar_index" > staged_clinvar_index.path
  >>>

  output {
    Array[File] staged_genotype_files = glob("staged/genotype/*")
    File staged_clinvar_vcf = read_string("staged_clinvar_vcf.path")
    File staged_clinvar_index = read_string("staged_clinvar_index.path")
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "32 GiB"
    disks: "local-disk 128 HDD"
    maxRetries: 2
  }
}

task FilterVariants {
  input {
    Array[File] staged_genotype_files
    Float af_max
    Float missing_max
    Float qual_min
    Int threads
    String docker_image
  }

  command <<<
    set -euo pipefail

    mkdir -p genotype_inputs filtered
    staged_list="~{write_lines(staged_genotype_files)}"
    while IFS= read -r staged_file; do
      cp -L "$staged_file" "genotype_inputs/$(basename "$staged_file")"
    done < "$staged_list"

    bash /opt/underlier-prevalence/scripts/filter_variants.sh \
      --input-dir genotype_inputs \
      --output-dir filtered \
      --af-max "~{af_max}" \
      --missing-max "~{missing_max}" \
      --qual-min "~{qual_min}" \
      --threads "~{threads}"

    test -s filtered/all.filtered.vcf.gz
    test -s filtered/all.filtered.vcf.gz.tbi
  >>>

  output {
    File filtered_vcf = "filtered/all.filtered.vcf.gz"
    File filtered_vcf_index = "filtered/all.filtered.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "32 GiB"
    disks: "local-disk 128 HDD"
    maxRetries: 2
  }
}

task ComputeQ2Incidence {
  input {
    File filtered_vcf
    File filtered_vcf_index
    File clinvar_vcf
    File clinvar_vcf_index
    File? gene_whitelist
    String docker_image
  }

  command <<<
    set -euo pipefail

    if [[ "$(basename "~{filtered_vcf_index}")" != "$(basename "~{filtered_vcf}").tbi" ]]; then
      echo "error: filtered VCF/index pair does not match" >&2
      exit 2
    fi
    if [[ "$(basename "~{clinvar_vcf_index}")" != "$(basename "~{clinvar_vcf}").tbi" ]]; then
      echo "error: staged ClinVar VCF/index pair does not match" >&2
      exit 2
    fi

    mkdir -p cohort results
    cp -L "~{filtered_vcf}" cohort/0000.filtered.vcf.gz
    cp -L "~{filtered_vcf_index}" cohort/0000.filtered.vcf.gz.tbi

    gene_whitelist="~{default="" gene_whitelist}"
    gene_args=()
    if [[ -n "$gene_whitelist" ]]; then
      gene_args=(--genes "$gene_whitelist")
    fi

    python3 /opt/underlier-prevalence/scripts/compute_q2_incidence.py \
      --clinvar "~{clinvar_vcf}" \
      --vcf-glob 'cohort/*.filtered.vcf.gz' \
      --out results/q2_incidence.tsv \
      "${gene_args[@]}"

    test -s results/q2_incidence.tsv
  >>>

  output {
    File q2_incidence_tsv = "results/q2_incidence.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "32 GiB"
    disks: "local-disk 128 HDD"
    maxRetries: 2
  }
}

workflow q2_incidence {
  input {
    Array[File] genotype_vcfs
    Array[File] genotype_vcf_indexes
    File clinvar_vcf
    File clinvar_vcf_index
    Float af_max = 0.01
    Float missing_max = 0.1
    Float qual_min = 100.0
    Int threads = 1
    File? gene_whitelist
    String docker_image = "ghcr.io/aou-multiomics-analysis/underlierprevelance:main"
  }

  call StageVcfPairs {
    input:
      genotype_vcfs = genotype_vcfs,
      genotype_vcf_indexes = genotype_vcf_indexes,
      clinvar_vcf = clinvar_vcf,
      clinvar_vcf_index = clinvar_vcf_index,
      docker_image = docker_image
  }

  call FilterVariants {
    input:
      staged_genotype_files = StageVcfPairs.staged_genotype_files,
      af_max = af_max,
      missing_max = missing_max,
      qual_min = qual_min,
      threads = threads,
      docker_image = docker_image
  }

  call ComputeQ2Incidence {
    input:
      filtered_vcf = FilterVariants.filtered_vcf,
      filtered_vcf_index = FilterVariants.filtered_vcf_index,
      clinvar_vcf = StageVcfPairs.staged_clinvar_vcf,
      clinvar_vcf_index = StageVcfPairs.staged_clinvar_index,
      gene_whitelist = gene_whitelist,
      docker_image = docker_image
  }

  output {
    File filtered_vcf = FilterVariants.filtered_vcf
    File filtered_vcf_index = FilterVariants.filtered_vcf_index
    File q2_incidence_tsv = ComputeQ2Incidence.q2_incidence_tsv
  }
}

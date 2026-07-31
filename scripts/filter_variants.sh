#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: filter_variants.sh --input-dir DIR --output-dir DIR --af-max FLOAT \
  --missing-max FLOAT --qual-min FLOAT --threads INT

Stages genotype VCFs from DIR, computes AC/AN/AF/F_MISSING with bcftools
+fill-tags, and writes output-dir/all.filtered.vcf.gz plus its tabix index.
USAGE
}

input_dir=''
output_dir=''
af_max=''
missing_max=''
qual_min=''
threads=''

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) input_dir=${2:-}; shift 2 ;;
    --output-dir) output_dir=${2:-}; shift 2 ;;
    --af-max) af_max=${2:-}; shift 2 ;;
    --missing-max) missing_max=${2:-}; shift 2 ;;
    --qual-min) qual_min=${2:-}; shift 2 ;;
    --threads) threads=${2:-}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'error: unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

for required in input_dir output_dir af_max missing_max qual_min threads; do
  if [[ -z ${!required} ]]; then
    printf 'error: --%s is required\n' "${required//_/-}" >&2
    usage >&2
    exit 2
  fi
done

if ! command -v bcftools >/dev/null 2>&1; then
  printf 'error: bcftools is required but was not found on PATH\n' >&2
  exit 127
fi

if [[ ! -d $input_dir ]]; then
  printf 'error: input directory does not exist: %s\n' "$input_dir" >&2
  exit 2
fi

decimal='^([0-9]+([.][0-9]*)?|[.][0-9]+)$'
for value_name in af_max missing_max qual_min; do
  if [[ ! ${!value_name} =~ $decimal ]]; then
    printf 'error: --%s must be a non-negative decimal\n' "${value_name//_/-}" >&2
    exit 2
  fi
done
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
  printf 'error: --threads must be a positive integer\n' >&2
  exit 2
fi

mkdir -p "$output_dir"
input_dir=$(cd "$input_dir" && pwd)
output_dir=$(cd "$output_dir" && pwd)
output_vcf="$output_dir/all.filtered.vcf.gz"
staging_dir=$(mktemp -d "${TMPDIR:-/tmp}/filter_variants.XXXXXX")
trap 'rm -rf "$staging_dir"' EXIT

inputs=()
while IFS= read -r input; do
  [[ $input == "$output_dir"/* ]] && continue
  inputs+=("$input")
done < <(find "$input_dir" -maxdepth 1 -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) -print | LC_ALL=C sort)

if [[ ${#inputs[@]} -eq 0 ]]; then
  printf 'error: no VCF or VCF.GZ files found in %s\n' "$input_dir" >&2
  exit 2
fi

genotype_vcfs=()
for index in "${!inputs[@]}"; do
  staged="$staging_dir/$(printf '%04d' "$index").vcf.gz"
  bcftools view --threads "$threads" -Oz -o "$staged" "${inputs[$index]}"
  bcftools index --threads "$threads" -f -t "$staged"
  if [[ -n $(bcftools query -l "$staged") ]]; then
    genotype_vcfs+=("$staged")
  fi
done

if [[ ${#genotype_vcfs[@]} -eq 0 ]]; then
  printf 'error: no input VCF contains genotype samples\n' >&2
  exit 2
fi

bcftools concat --threads "$threads" --allow-overlaps -Ou "${genotype_vcfs[@]}" \
  | bcftools +fill-tags -Ou -- -t AC,AN,AF,F_MISSING \
  | bcftools view --threads "$threads" -m2 -M2 \
      -i "N_ALT=1 && INFO/AF[0]<=${af_max} && INFO/F_MISSING<=${missing_max} && QUAL>=${qual_min}" \
      -Oz -o "$output_vcf"
bcftools index --threads "$threads" -f -t "$output_vcf"

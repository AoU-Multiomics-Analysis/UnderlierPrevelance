import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def read_paired_vcf_inputs():
    manifest_lines = (FIXTURES / "vcf_inputs.tsv").read_text().splitlines()
    entries = [
        dict(zip(manifest_lines[0].split("\t"), line.split("\t")))
        for line in manifest_lines[1:]
    ]
    genotype = [entry for entry in entries if entry["role"] == "genotype"]
    clinvar = next(entry for entry in entries if entry["role"] == "clinvar")
    genotype_vcfs = [FIXTURES / entry["vcf"] for entry in genotype]
    genotype_indexes = [FIXTURES / entry["index"] for entry in genotype]
    clinvar_vcf = FIXTURES / clinvar["vcf"]
    clinvar_index = FIXTURES / clinvar["index"]
    assert len(genotype_vcfs) == len(genotype_indexes)
    assert all(path.exists() for path in genotype_vcfs + genotype_indexes)
    assert clinvar_vcf.exists() and clinvar_index.exists()
    return genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index


def run_q2(output):
    genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index = read_paired_vcf_inputs()
    subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "compute_q2_incidence.py"),
            "--clinvar",
            str(clinvar_vcf),
            "--vcf-glob",
            str(genotype_vcfs[0]),
            "--out",
            str(output),
        ],
        check=True,
        cwd=ROOT,
    )
    return output


def run_filtering(output_dir):
    read_paired_vcf_inputs()
    subprocess.run(
        [
            "bash",
            str(ROOT / "scripts" / "filter_variants.sh"),
            "--input-dir",
            str(FIXTURES),
            "--output-dir",
            str(output_dir),
            "--af-max",
            "0.25",
            "--missing-max",
            "0.1",
            "--qual-min",
            "125",
            "--threads",
            "1",
        ],
        check=True,
        cwd=ROOT,
    )
    return output_dir / "all.filtered.vcf.gz"


def read_variant_count(path):
    result = subprocess.run(
        ["bcftools", "view", "-H", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return len(result.stdout.splitlines())


def test_q2_output_contains_required_columns(tmp_path):
    output = tmp_path / "q2.tsv"
    run_q2(output)
    columns = output.read_text().splitlines()[0].split("\t")
    assert columns == [
        "gene",
        "n_plp_alleles_clinvar",
        "n_plp_alleles_observed",
        "q",
        "incidence",
        "carrier_freq",
    ]


def test_filtering_removes_nonbiallelic_and_high_af_sites(tmp_path):
    filtered = run_filtering(tmp_path)
    assert read_variant_count(filtered) == 3

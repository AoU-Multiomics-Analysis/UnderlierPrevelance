import subprocess
from pathlib import Path


def run_q2(output):
    subprocess.run(
        [
            "python3",
            "scripts/compute_q2_incidence.py",
            "--clinvar",
            "tests/fixtures/clinvar.vcf",
            "--vcf-glob",
            "tests/fixtures/cohort.vcf",
            "--out",
            str(output),
        ],
        check=True,
    )
    return output


def run_filtering(output_dir):
    subprocess.run(
        [
            "bash",
            "scripts/filter_variants.sh",
            "--input-dir",
            "tests/fixtures",
            "--output-dir",
            str(output_dir),
            "--af-max",
            "0.01",
            "--missing-max",
            "0.1",
            "--qual-min",
            "125",
            "--threads",
            "1",
        ],
        check=True,
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

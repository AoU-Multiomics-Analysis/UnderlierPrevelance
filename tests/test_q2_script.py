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
    if len(genotype_vcfs) != len(genotype_indexes):
        raise ValueError("genotype VCF and index arrays must have equal lengths")
    if not all(path.exists() for path in [*genotype_vcfs, *genotype_indexes, clinvar_vcf, clinvar_index]):
        raise ValueError("all VCF and index paths must exist")
    if any(vcf.name + ".tbi" != index.name for vcf, index in zip(genotype_vcfs, genotype_indexes)):
        raise ValueError("genotype VCF/index paths must be paired")
    if clinvar_vcf.name + ".tbi" != clinvar_index.name:
        raise ValueError("ClinVar VCF/index paths must be paired")
    return genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index


def run_q2(output, vcf_glob=None, genes=None):
    genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index = read_paired_vcf_inputs()
    assert len(genotype_vcfs) == len(genotype_indexes) == 2
    command = [
        "python3",
        str(ROOT / "scripts" / "compute_q2_incidence.py"),
        "--clinvar",
        str(clinvar_vcf),
        "--vcf-glob",
        vcf_glob or str(FIXTURES / "cohort*.vcf"),
        "--out",
        str(output),
    ]
    if genes is not None:
        command.extend(["--genes", str(genes)])
    subprocess.run(
        command,
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


def test_q2_uses_maximum_duplicate_allele_af_and_gene_level_q(tmp_path):
    clinvar = tmp_path / "clinvar.vcf"
    clinvar.write_text(
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">\n"
        "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;GENE=GENE1\n"
        "chr1\t20\t.\tC\tT\t.\tPASS\tCLNSIG=Likely_pathogenic;GENE=GENE1\n"
    )
    cohort_header = (
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    )
    (tmp_path / "cohort_low.vcf").write_text(
        cohort_header + "chr1\t10\t.\tA\tG\t200\tPASS\t.\tGT\t0/1\t0/0\n"
    )
    (tmp_path / "cohort_high.vcf").write_text(
        cohort_header
        + "chr1\t10\t.\tA\tG\t200\tPASS\t.\tGT\t1/1\t0/0\n"
        + "chr1\t20\t.\tC\tT\t200\tPASS\t.\tGT\t0/1\t0/0\n"
    )

    output = tmp_path / "q2.tsv"
    subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "compute_q2_incidence.py"),
            "--clinvar",
            str(clinvar),
            "--vcf-glob",
            str(tmp_path / "cohort_*.vcf"),
            "--out",
            str(output),
        ],
        check=True,
        cwd=ROOT,
    )

    row = dict(zip(output.read_text().splitlines()[0].split("\t"), output.read_text().splitlines()[1].split("\t")))
    assert row == {
        "gene": "GENE1",
        "n_plp_alleles_clinvar": "2",
        "n_plp_alleles_observed": "2",
        "q": "0.75",
        "incidence": "0.5625",
        "carrier_freq": "0.9375",
    }


def test_q2_rejects_a_vcf_glob_that_matches_no_files(tmp_path):
    result = subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "compute_q2_incidence.py"),
            "--clinvar",
            str(FIXTURES / "clinvar.vcf"),
            "--vcf-glob",
            str(tmp_path / "missing_*.vcf"),
            "--out",
            str(tmp_path / "q2.tsv"),
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert result.returncode != 0
    assert "matched no VCF files" in result.stderr


def test_q2_gene_whitelist_limits_output_genes(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("GENE2\n")

    output = run_q2(tmp_path / "q2.tsv", genes=genes)

    assert output.read_text().splitlines()[1].split("\t")[0] == "GENE2"
    assert len(output.read_text().splitlines()) == 2


def test_q2_helper_covers_all_genotype_vcf_index_pairs():
    genotype_vcfs, genotype_indexes, _, _ = read_paired_vcf_inputs()
    assert len(genotype_vcfs) == len(genotype_indexes) == 2
    assert {path.name for path in genotype_vcfs} == {"cohort.vcf", "cohort_2.vcf"}


def test_filtering_removes_nonbiallelic_and_high_af_sites(tmp_path):
    filtered = run_filtering(tmp_path)
    assert read_variant_count(filtered) == 3

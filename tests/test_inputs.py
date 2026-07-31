from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def parse_vcf_manifest(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    entries = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    genotype = [entry for entry in entries if entry["role"] == "genotype"]
    clinvar = next(entry for entry in entries if entry["role"] == "clinvar")
    return (
        [FIXTURES / entry["vcf"] for entry in genotype],
        [FIXTURES / entry["index"] for entry in genotype],
        FIXTURES / clinvar["vcf"],
        FIXTURES / clinvar["index"],
    )


def validate_vcf_index_pairs(genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index):
    if len(genotype_vcfs) != len(genotype_indexes):
        raise ValueError("genotype VCF and index arrays must have equal lengths")
    if not genotype_vcfs:
        raise ValueError("at least one genotype VCF/index pair is required")
    paths = [*genotype_vcfs, *genotype_indexes, clinvar_vcf, clinvar_index]
    if not all(path.exists() for path in paths):
        raise ValueError("all VCF and index paths must exist")
    if any(vcf.name + ".tbi" != index.name for vcf, index in zip(genotype_vcfs, genotype_indexes)):
        raise ValueError("genotype VCF/index paths must be paired")
    if clinvar_vcf.name + ".tbi" != clinvar_index.name:
        raise ValueError("ClinVar VCF/index paths must be paired")


def read_gct_contract(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3 or lines[0] != "#1.2":
        raise ValueError("GCT must start with #1.2")
    try:
        expected_genes, expected_samples = map(int, lines[1].split("\t"))
    except (ValueError, TypeError):
        raise ValueError("GCT dimensions must be two integers")
    header = lines[2].split("\t")
    if header[:2] != ["Name", "Description"]:
        raise ValueError("GCT header must start with Name and Description")
    if len(header) != expected_samples + 2 or len(lines[3:]) != expected_genes:
        raise ValueError("GCT dimensions do not match the data")
    for row_number, line in enumerate(lines[3:], start=4):
        fields = line.split("\t")
        if len(fields) != expected_samples + 2:
            raise ValueError(f"GCT row {row_number} has the wrong field count")
        if not fields[0] or not fields[1]:
            raise ValueError(f"GCT row {row_number} has missing Name or Description")
    return header[2:], lines[3:]


def validate_covariates(path, expected_samples, n_pcs=2):
    rows = [line.split("\t") for line in Path(path).read_text().splitlines()]
    header = rows[0]
    expected_columns = ["sample_id"] + [f"Genotype_PC{i}" for i in range(1, n_pcs + 1)]
    if header != expected_columns:
        raise ValueError("covariate columns do not match the required schema")
    for row_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise ValueError(f"covariate row {row_number} has the wrong field count")
    sample_ids = [row[0] for row in rows[1:]]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("covariate sample IDs must be unique")
    if sample_ids != list(expected_samples):
        raise ValueError("covariate sample IDs must exactly match the GCT samples")
    for row in rows[1:]:
        for value in row[1:]:
            try:
                float(value)
            except ValueError:
                raise ValueError("genotype PC values must be numeric")


def test_gct_fixture_has_standard_header():
    lines = (FIXTURES / "counts.gct").read_text().splitlines()
    assert lines[0] == "#1.2"
    assert lines[2].split("\t")[:2] == ["Name", "Description"]


def test_gct_fixture_declares_six_genes_and_four_samples():
    lines = (FIXTURES / "counts.gct").read_text().splitlines()
    assert lines[1].split("\t") == ["6", "4"]
    assert lines[2].split("\t")[2:] == ["S1", "S2", "S3", "S4"]
    assert len(lines[3:]) == 6


def test_covariates_use_sample_id_and_genotype_pc_columns():
    header = (FIXTURES / "genotype_covariates.tsv").read_text().splitlines()[0].split("\t")
    assert header[:2] == ["sample_id", "Genotype_PC1"]
    assert header == ["sample_id", "Genotype_PC1", "Genotype_PC2"]


def test_covariates_align_to_all_gct_samples():
    gct_header = (FIXTURES / "counts.gct").read_text().splitlines()[2].split("\t")
    covariate_header, *rows = (FIXTURES / "genotype_covariates.tsv").read_text().splitlines()
    assert [row.split("\t")[0] for row in rows] == gct_header[2:]


def test_gct_contract_rejects_malformed_header(tmp_path):
    malformed = tmp_path / "malformed.gct"
    malformed.write_text((FIXTURES / "counts.gct").read_text().replace("#1.2", "#1.1", 1))
    with pytest.raises(ValueError, match="#1.2"):
        read_gct_contract(malformed)


def test_gct_contract_rejects_wrong_dimensions(tmp_path):
    malformed = tmp_path / "wrong-dimensions.gct"
    malformed.write_text((FIXTURES / "counts.gct").read_text().replace("6\t4", "5\t4", 1))
    with pytest.raises(ValueError, match="dimensions"):
        read_gct_contract(malformed)


def test_gct_contract_rejects_truncated_data_row(tmp_path):
    malformed = tmp_path / "truncated-row.gct"
    malformed.write_text(
        (FIXTURES / "counts.gct").read_text().replace("GENE1\tGene 1\t10\t12\t14\t16", "GENE1\tGene 1\t10\t12\t14", 1)
    )
    with pytest.raises(ValueError, match="field count"):
        read_gct_contract(malformed)


@pytest.mark.parametrize("replacement", ["\tGene 1\t10", "GENE1\t\t10"])
def test_gct_contract_rejects_missing_metadata(tmp_path, replacement):
    malformed = tmp_path / "missing-metadata.gct"
    malformed.write_text(
        (FIXTURES / "counts.gct").read_text().replace("GENE1\tGene 1\t10", replacement, 1)
    )
    with pytest.raises(ValueError, match="missing Name or Description"):
        read_gct_contract(malformed)


def test_covariate_contract_rejects_duplicate_sample_ids(tmp_path):
    malformed = tmp_path / "duplicate-sample.tsv"
    malformed.write_text(
        (FIXTURES / "genotype_covariates.tsv").read_text().replace("S2\t-0.5", "S1\t-0.5", 1)
    )
    with pytest.raises(ValueError, match="unique"):
        validate_covariates(malformed, ["S1", "S2", "S3", "S4"])


def test_covariate_contract_rejects_missing_sample_id(tmp_path):
    malformed = tmp_path / "missing-sample.tsv"
    malformed.write_text(
        (FIXTURES / "genotype_covariates.tsv").read_text().replace("S4\t1.5", "S5\t1.5", 1)
    )
    with pytest.raises(ValueError, match="exactly match"):
        validate_covariates(malformed, ["S1", "S2", "S3", "S4"])


def test_covariate_contract_rejects_missing_genotype_pc_column(tmp_path):
    malformed = tmp_path / "missing-pc.tsv"
    malformed.write_text(
        (FIXTURES / "genotype_covariates.tsv").read_text().replace("\tGenotype_PC2", "", 1)
    )
    with pytest.raises(ValueError, match="columns"):
        validate_covariates(malformed, ["S1", "S2", "S3", "S4"])


def test_covariate_contract_rejects_non_numeric_genotype_pc(tmp_path):
    malformed = tmp_path / "nonnumeric-pc.tsv"
    malformed.write_text(
        (FIXTURES / "genotype_covariates.tsv").read_text().replace("S1\t-1.5", "S1\tnot-a-number", 1)
    )
    with pytest.raises(ValueError, match="numeric"):
        validate_covariates(malformed, ["S1", "S2", "S3", "S4"])


def test_covariate_contract_rejects_truncated_row(tmp_path):
    malformed = tmp_path / "truncated-row.tsv"
    malformed.write_text(
        (FIXTURES / "genotype_covariates.tsv").read_text().replace("S1\t-1.5\t0.5", "S1\t-1.5", 1)
    )
    with pytest.raises(ValueError, match="field count"):
        validate_covariates(malformed, ["S1", "S2", "S3", "S4"])


def test_gct_and_covariate_contracts_accept_fixture():
    samples, rows = read_gct_contract(FIXTURES / "counts.gct")
    assert len(rows) == 6
    validate_covariates(FIXTURES / "genotype_covariates.tsv", samples)


def test_q2_inputs_have_explicit_vcf_index_pairs():
    genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index = parse_vcf_manifest(
        FIXTURES / "vcf_inputs.tsv"
    )
    assert len(genotype_vcfs) == 2
    assert len(genotype_indexes) == 2
    assert len(genotype_vcfs) == len(genotype_indexes)
    validate_vcf_index_pairs(genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index)
    # These placeholders validate path pairing only; WDL/bcftools must validate real tabix indexes later.


def test_q2_inputs_reject_shortened_genotype_index_list():
    genotype_vcfs, genotype_indexes, clinvar_vcf, clinvar_index = parse_vcf_manifest(
        FIXTURES / "vcf_inputs.tsv"
    )
    with pytest.raises(ValueError, match="equal lengths"):
        validate_vcf_index_pairs(
            genotype_vcfs, genotype_indexes[:-1], clinvar_vcf, clinvar_index
        )


def test_annotation_contains_two_protein_coding_genes():
    records = [
        line.split("\t")
        for line in (FIXTURES / "annotation.gff3").read_text().splitlines()
        if line and not line.startswith("#")
    ]
    protein_coding = [record for record in records if "protein_coding" in record[8]]
    assert len(protein_coding) == 2


def test_clinvar_and_cohort_share_three_pathogenic_alleles():
    clinvar = {
        tuple(line.split("\t")[:5])
        for line in (FIXTURES / "clinvar.vcf").read_text().splitlines()
        if line and not line.startswith("#") and "CLNSIG=Pathogenic" in line
    }
    cohort = {
        tuple(line.split("\t")[:5])
        for line in (FIXTURES / "cohort.vcf").read_text().splitlines()
        if line and not line.startswith("#")
    }
    assert len(clinvar & cohort) >= 3

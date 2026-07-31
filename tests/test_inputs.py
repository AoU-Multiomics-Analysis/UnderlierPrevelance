import pytest
from pathlib import Path


FIXTURES = Path("tests/fixtures")


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

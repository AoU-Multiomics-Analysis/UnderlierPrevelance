import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def run_converter(input_path, output_path):
    return subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "convert_gct.py"),
            "--input",
            str(input_path),
            "--output",
            str(output_path),
        ],
        check=True,
    )


def test_convert_gct_writes_gene_id_and_sample_columns(tmp_path):
    output = tmp_path / "counts.tsv"
    run_converter(FIXTURES / "counts.gct", output)
    header = output.read_text().splitlines()[0].split("\t")
    assert header == ["gene_id", "S1", "S2", "S3", "S4"]


def test_converter_rejects_duplicate_gene_ids(tmp_path):
    bad = tmp_path / "bad.gct"
    bad.write_text(
        (FIXTURES / "counts.gct").read_text().replace("GENE2", "GENE1", 1)
    )
    with pytest.raises(subprocess.CalledProcessError):
        run_converter(bad, tmp_path / "out.tsv")


def test_converter_rejects_malformed_header(tmp_path):
    bad = tmp_path / "bad-header.gct"
    bad.write_text((FIXTURES / "counts.gct").read_text().replace("#1.2", "#1.1", 1))
    with pytest.raises(subprocess.CalledProcessError):
        run_converter(bad, tmp_path / "out.tsv")


def test_converter_rejects_wrong_dimensions(tmp_path):
    bad = tmp_path / "bad-dimensions.gct"
    bad.write_text((FIXTURES / "counts.gct").read_text().replace("6\t4", "5\t4", 1))
    with pytest.raises(subprocess.CalledProcessError):
        run_converter(bad, tmp_path / "out.tsv")


def test_converter_rejects_nonnumeric_counts(tmp_path):
    bad = tmp_path / "bad-count.gct"
    bad.write_text(
        (FIXTURES / "counts.gct").read_text().replace(
            "\t10\t12\t14\t16", "\tnot-a-count\t12\t14\t16", 1
        )
    )
    with pytest.raises(subprocess.CalledProcessError):
        run_converter(bad, tmp_path / "out.tsv")


def test_converter_rejects_negative_counts_with_row_and_sample_context(tmp_path):
    bad = tmp_path / "negative-count.gct"
    bad.write_text(
        (FIXTURES / "counts.gct").read_text().replace(
            "\t10\t12\t14\t16", "\t-1\t12\t14\t16", 1
        )
    )
    output = tmp_path / "out.tsv"

    result = subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "convert_gct.py"),
            "--input",
            str(bad),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )

    assert result.returncode == 2
    assert "GCT row 4, sample 'S1' must have a non-negative count" in result.stderr
    assert not output.exists()


def test_converter_does_not_replace_output_when_validation_fails(tmp_path):
    bad = tmp_path / "bad-count.gct"
    bad.write_text(
        (FIXTURES / "counts.gct").read_text().replace(
            "\t10\t12\t14\t16", "\tnot-a-count\t12\t14\t16", 1
        )
    )
    output = tmp_path / "counts.tsv"
    output.write_text("existing output\n")

    with pytest.raises(subprocess.CalledProcessError):
        run_converter(bad, output)

    assert output.read_text() == "existing output\n"


def test_rna_cli_accepts_documented_hyphenated_flags_without_dependencies():
    parser_program = (
        "source('scripts/run_rna_underlier.R'); "
        "options <- parse_cli(commandArgs(trailingOnly = TRUE)); "
        "cat(options$counts, options$genotype_covariates, options$gencode, "
        "options$out_dir, options$n_geno_pcs, options$phenotype_pc_noise, "
        "options$connectivity_z, options$logcpm_drop, options$z_cutoff, "
        "options$threads, sep = '\\t')"
    )
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        [
            "Rscript",
            "-e",
            parser_program,
            "--counts",
            "counts.tsv",
            "--genotype-covariates",
            "covariates.tsv",
            "--gencode",
            "annotation.gff3",
            "--out-dir",
            "out",
            "--n-geno-pcs",
            "2",
            "--phenotype-pc-noise",
            "0.25",
            "--connectivity-z",
            "-3",
            "--logcpm-drop",
            "1",
            "--z-cutoff",
            "-3",
            "--threads",
            "1",
        ],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.stdout == "counts.tsv\tcovariates.tsv\tannotation.gff3\tout\t2\t0.25\t-3\t1\t-3\t1"


@pytest.mark.parametrize(
    ("matrix_values", "n_samples", "expected_error"),
    [
        ("1, 2, 3, 4, 1, 2, 3, 4", 4, "is rank-deficient"),
        ("1, 2, 3, 2, 3, 5", 3, "leaves no residual degrees of freedom"),
    ],
)
def test_rna_covariate_design_rejects_invalid_rank_without_dependencies(
    matrix_values, n_samples, expected_error
):
    program = (
        "source('scripts/run_rna_underlier.R'); "
        f"covariates <- matrix(c({matrix_values}), nrow = {n_samples}, ncol = 2); "
        "rownames(covariates) <- paste0('S', seq_len(nrow(covariates))); "
        "colnames(covariates) <- c('PC1', 'Genotype_PC1'); "
        "validate_covariate_design(covariates, rownames(covariates), 'test')"
    )
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        ["Rscript", "-e", program],
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.returncode != 0
    assert expected_error in result.stderr


def test_rna_strict_output_tag_uses_the_requested_cutoff_without_dependencies():
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        [
            "Rscript",
            "-e",
            "source('scripts/run_rna_underlier.R'); cat(strict_output_tag(-2))",
        ],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.stdout == "z_-2"


def test_rna_gd_metadata_records_exact_covariate_design_without_dependencies():
    program = (
        "source('scripts/run_rna_underlier.R'); "
        "metadata <- selected_pc_metadata(1L, 1L, 4L, 0.25, 'override', 2L, "
        "list(columns = 'PC1', design_rank = 2L), "
        "c('Genotype_PC1', 'Genotype_PC2'), "
        "list(columns = c('PC1', 'Genotype_PC1', 'Genotype_PC2'), design_rank = 4L, "
        "residual_degrees_freedom = 6L), 10L, 16L); "
        "cat(metadata$phenotype_pc_method, metadata$phenotype_pc_columns, metadata$genotype_pc_columns, "
        "metadata$residualization_covariate_columns, metadata$residualization_design_rank, "
        "metadata$residual_degrees_freedom, sep = '\\t')"
    )
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        ["Rscript", "-e", program],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.stdout == "PCAtools::chooseGavishDonoho\tPC1\tGenotype_PC1,Genotype_PC2\tPC1,Genotype_PC1,Genotype_PC2\t4\t6"


def test_rna_counts_rejects_a_trailing_empty_count_without_dependencies(tmp_path):
    counts = tmp_path / "trailing-empty.tsv"
    counts.write_text("gene_id\tR1\tR2\nRNA_GENE01\t10\t\n")
    program = (
        "source('scripts/run_rna_underlier.R'); "
        f"read_counts('{counts}')"
    )
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        ["Rscript", "-e", program],
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.returncode != 0
    assert "Counts TSV contains nonnumeric or non-finite values" in result.stderr


def test_rna_counts_rejects_matching_trailing_delimiters_without_dependencies(tmp_path):
    counts = tmp_path / "matching-trailing-empty.tsv"
    counts.write_text("gene_id\tR1\t\nRNA_GENE01\t10\t\n")
    program = (
        "source('scripts/run_rna_underlier.R'); "
        f"read_counts('{counts}')"
    )
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        ["Rscript", "-e", program],
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.returncode != 0
    assert "Counts TSV sample IDs must be nonempty and unique" in result.stderr


def test_rna_covariates_reject_matching_trailing_delimiters_without_dependencies(tmp_path):
    covariates = tmp_path / "matching-trailing-empty.tsv"
    covariates.write_text("sample_id\tGenotype_PC1\t\nR1\t0.5\t\n")
    program = (
        "source('scripts/run_rna_underlier.R'); "
        f"read_genotype_covariates('{covariates}', 'R1', 1L)"
    )
    environment = os.environ | {"TOPMED_RNA_UNDERLIER_NO_MAIN": "1"}
    result = subprocess.run(
        ["Rscript", "-e", program],
        capture_output=True,
        env=environment,
        text=True,
    )
    assert result.returncode != 0
    assert "Genotype-covariate TSV columns after sample_id" in result.stderr


def test_rna_smoke_fixture_has_capacity_for_two_genotype_and_three_phenotype_pcs():
    samples, rows = read_gct_contract(FIXTURES / "rna_smoke_counts.gct")
    validate_covariates(FIXTURES / "rna_smoke_genotype_covariates.tsv", samples)
    records = [
        line
        for line in (FIXTURES / "rna_smoke_annotation.gff3").read_text().splitlines()
        if line and not line.startswith("#")
    ]
    assert len(samples) == 10
    assert len(rows) == 16
    assert len(records) == 16
    assert len(samples) - (1 + 2 + 3) > 0


def rna_runtime_available():
    if shutil.which("Rscript") is None:
        return False
    package_probe = (
        "required <- c('edgeR', 'limma', 'corral', 'PCAtools', 'scran', 'WGCNA', 'rtracklayer'); "
        "quit(status = as.integer(!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))))"
    )
    return subprocess.run(["Rscript", "-e", package_probe], check=False).returncode == 0


@pytest.mark.skipif(
    not rna_runtime_available(),
    reason="RNA smoke test requires edgeR, limma, corral, PCAtools, scran, WGCNA, and rtracklayer",
)
def test_rna_smoke_fixture_emits_core_outputs(tmp_path):
    counts_tsv = tmp_path / "counts.tsv"
    run_converter(FIXTURES / "rna_smoke_counts.gct", counts_tsv)
    out_dir = tmp_path / "rna-out"
    subprocess.run(
        [
            "Rscript",
            str(ROOT / "scripts" / "run_rna_underlier.R"),
            "--counts",
            str(counts_tsv),
            "--genotype-covariates",
            str(FIXTURES / "rna_smoke_genotype_covariates.tsv"),
            "--gencode",
            str(FIXTURES / "rna_smoke_annotation.gff3"),
            "--out-dir",
            str(out_dir),
            "--n-geno-pcs",
            "2",
            "--phenotype-pc-noise",
            "0.25",
        ],
        check=True,
    )
    expected_outputs = [
        "selected_phenotype_pcs.tsv",
        "expr_z_join.tsv.gz",
        "underliers_haplo.tsv.gz",
        "underliers_z_-3.tsv.gz",
        "rna_outlier_prevalence_per_gene_haplo.tsv",
        "rna_outlier_prevalence_per_gene_z_-3.tsv",
    ]
    assert all((out_dir / output).is_file() for output in expected_outputs)
    metadata = dict(
        zip(
            (out_dir / "selected_phenotype_pcs.tsv").read_text().splitlines()[0].split("\t"),
            (out_dir / "selected_phenotype_pcs.tsv").read_text().splitlines()[1].split("\t"),
        )
    )
    assert metadata["phenotype_pc_method"] == "PCAtools::chooseGavishDonoho"
    assert metadata["noise_source"] == "override"
    assert metadata["noise_variance"] == "0.25"
    assert metadata["selected_phenotype_pcs"] == "3"
    assert metadata["gavish_donoho_raw"] == "3"
    assert metadata["phenotype_pc_columns"] == "PC1,PC2,PC3"
    assert metadata["residualization_design_rank"] == "6"
    assert metadata["residual_degrees_freedom"] == "4"


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

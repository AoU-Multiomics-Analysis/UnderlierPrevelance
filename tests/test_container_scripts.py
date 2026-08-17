import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_IMAGE = "ghcr.io/aou-multiomics-analysis/underlierprevelance:main"
SCRIPT_NAMES = (
    "convert_gct.py",
    "run_rna_underlier.R",
    "filter_variants.sh",
    "compute_q2_incidence.py",
)
IMAGE_SCRIPT_DIR = "/opt/underlier-prevalence/scripts"


def test_runtime_scripts_are_packaged_and_not_wdl_file_inputs():
    dockerfile = (ROOT / "envs" / "Dockerfile").read_text()
    assert "COPY scripts/ /opt/underlier-prevalence/scripts/" in dockerfile
    for script_name in SCRIPT_NAMES:
        assert (ROOT / "scripts" / script_name).is_file()

    wdl_text = "\n".join(
        (ROOT / "workflows" / filename).read_text()
        for filename in ("main.wdl", "rna_underlier.wdl", "q2_incidence.wdl")
    )
    assert "File convert_gct_script" not in wdl_text
    assert "File rna_underlier_script" not in wdl_text
    assert "File filter_variants_script" not in wdl_text
    assert "File compute_q2_script" not in wdl_text
    for script_name in SCRIPT_NAMES:
        assert f"{IMAGE_SCRIPT_DIR}/{script_name}" in wdl_text

    test_inputs = json.loads((ROOT / "tests" / "test_wdl_inputs.json").read_text())
    assert not any("script" in input_name for input_name in test_inputs)


def test_wdl_tasks_use_expected_memory_and_128_gib_disk():
    for workflow_path in (ROOT / "workflows").glob("*.wdl"):
        workflow_text = workflow_path.read_text()
        runtime_count = workflow_text.count("runtime {")
        if workflow_path.name == "rna_underlier.wdl":
            assert workflow_text.count('memory: "32 GiB"') == runtime_count - 1
            assert workflow_text.count('memory: "256 GiB"') == 1
        else:
            assert workflow_text.count('memory: "32 GiB"') == runtime_count
        assert workflow_text.count('disks: "local-disk 128 HDD"') == workflow_text.count("runtime {")


def test_runtime_r_stack_is_declared_in_conda_environment():
    dockerfile = (ROOT / "envs" / "Dockerfile").read_text()
    environment = (ROOT / "envs" / "bioconda-environment.yml").read_text()

    assert "FROM mambaorg/micromamba:2.0.5" in dockerfile
    assert "BiocManager" not in dockerfile
    assert "install.packages" not in dockerfile
    for package in (
        "r-base",
        "bioconductor-edger",
        "bioconductor-limma",
        "bioconductor-corral",
        "bioconductor-pcatools",
        "bioconductor-scran",
        "bioconductor-scuttle",
        "bioconductor-rtracklayer",
        "r-wgcna",
        "r-vroom",
        "r-tidyverse",
        "r-magrittr",
    ):
        assert f"- {package}" in environment


def test_production_image_is_published_and_used_by_wdl_defaults():
    docker_workflow = (ROOT / ".github" / "workflows" / "docker-image.yml").read_text()
    assert "push: true" in docker_workflow

    for workflow_name in ("main.wdl", "rna_underlier.wdl", "q2_incidence.wdl"):
        workflow_text = (ROOT / "workflows" / workflow_name).read_text()
        assert f'docker_image = "{PRODUCTION_IMAGE}"' in workflow_text


def test_rna_workflow_exposes_default_z_cutoff_range():
    workflow_text = (ROOT / "workflows" / "rna_underlier.wdl").read_text()
    main_text = (ROOT / "workflows" / "main.wdl").read_text()
    r_script = (ROOT / "scripts" / "run_rna_underlier.R").read_text()
    expected_cutoffs = ", ".join(f"{value}.0" for value in range(-1, -11, -1))
    assert f"Array[Float] z_cutoffs = [{expected_cutoffs}]" in workflow_text
    assert "Array[Float] rna_z_cutoffs" in main_text
    assert "--z-cutoffs-file" in workflow_text
    assert "expected_artifacts=$((1 + 2 * ~{length(z_cutoffs)}))" in workflow_text
    assert "z_cutoffs" in r_script


def test_genotype_pc_count_is_optional_and_defaults_to_all_columns():
    workflow_text = (ROOT / "workflows" / "rna_underlier.wdl").read_text()
    main_text = (ROOT / "workflows" / "main.wdl").read_text()
    r_script = (ROOT / "scripts" / "run_rna_underlier.R").read_text()
    assert workflow_text.count("Int? n_genotype_pcs") == 2
    assert "Int? rna_n_genotype_pcs" in main_text
    assert "if (is.null(n_geno_pcs)) canonical_pc_columns" in r_script


def test_rna_pca_scores_are_coerced_to_numeric_matrix():
    r_script = (ROOT / "scripts" / "run_rna_underlier.R").read_text()
    assert "rotated_scores <- as.matrix(pca_result$rotated)" in r_script
    assert 'storage.mode(rotated_scores) <- "double"' in r_script

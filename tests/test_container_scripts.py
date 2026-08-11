import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
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


def test_all_wdl_tasks_use_default_32_gib_memory_and_128_gib_disk():
    for workflow_path in (ROOT / "workflows").glob("*.wdl"):
        workflow_text = workflow_path.read_text()
        assert workflow_text.count('memory: "32 GiB"') == workflow_text.count("runtime {")
        assert workflow_text.count('disks: "local-disk 128 HDD"') == workflow_text.count("runtime {")


def test_rna_workflow_exposes_default_z_cutoff_range():
    workflow_text = (ROOT / "workflows" / "rna_underlier.wdl").read_text()
    main_text = (ROOT / "workflows" / "main.wdl").read_text()
    r_script = (ROOT / "scripts" / "run_rna_underlier.R").read_text()
    expected_cutoffs = ", ".join(f"{value}.0" for value in range(-1, -11, -1))
    assert f"Array[Float] z_cutoffs = [{expected_cutoffs}]" in workflow_text
    assert "Array[Float] rna_z_cutoffs" in main_text
    assert "--z-cutoffs-file" in workflow_text
    assert "z_cutoffs" in r_script

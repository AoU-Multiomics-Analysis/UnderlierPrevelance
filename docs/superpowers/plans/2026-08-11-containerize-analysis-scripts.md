# Containerize Analysis Scripts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Docker image self-contained by packaging every runtime analysis script and removing script-file inputs from the WDL interfaces.

**Architecture:** `envs/Dockerfile` copies the four repository runtime scripts into `/opt/underlier-prevalence/scripts/`. WDL commands invoke those fixed paths, so Cromwell/WDL inputs contain only scientific data, configuration, and the Docker image reference.

**Tech Stack:** Docker, WDL 1.0, Python, R, Bash, pytest, miniwdl.

## Global Constraints

- Runtime script paths are fixed at `/opt/underlier-prevalence/scripts/` inside the image.
- The four packaged scripts are `convert_gct.py`, `run_rna_underlier.R`, `filter_variants.sh`, and `compute_q2_incidence.py`.
- No WDL workflow or task may require a runtime analysis script as a `File` input.
- Scientific input files and optional gene whitelist remain external WDL inputs.

---

### Task 1: Add the self-contained image/WDL contract test

**Files:**
- Create: `tests/test_container_scripts.py`

- [x] **Step 1: Write the failing test**

Assert that the Dockerfile copies all four scripts to the fixed image directory, that the WDL files contain no script `File` inputs, and that the main workflow/test JSON contain no script-file parameters.

- [x] **Step 2: Run the test to verify it fails**

Run: `pytest -q tests/test_container_scripts.py`

Expected: FAIL because the current Dockerfile does not copy scripts and the WDLs still accept script files.

### Task 2: Package scripts and update WDL interfaces

**Files:**
- Modify: `envs/Dockerfile`
- Modify: `workflows/rna_underlier.wdl`
- Modify: `workflows/q2_incidence.wdl`
- Modify: `workflows/main.wdl`
- Modify: `tests/test_wdl_inputs.json`

- [x] **Step 1: Add scripts to the image**

Copy `scripts/` into `/opt/underlier-prevalence/scripts/` in the final image.

- [x] **Step 2: Remove script-file WDL inputs**

Remove the four script `File` inputs and use the fixed paths in the command blocks and workflow call wiring.

- [x] **Step 3: Remove obsolete test inputs**

Delete the four script-file entries from `tests/test_wdl_inputs.json`.

- [x] **Step 4: Run focused validation**

Run: `pytest -q tests/test_container_scripts.py && miniwdl check workflows/main.wdl workflows/rna_underlier.wdl workflows/q2_incidence.wdl`

Expected: PASS.

### Task 3: Update execution documentation and run regression checks

**Files:**
- Modify: `docs/running.md`
- Modify: `README.md` if script-file inputs are documented there

- [x] **Step 1: Document image-owned scripts**

State that callers no longer provide runtime scripts and must use an image containing the packaged scripts.

- [x] **Step 2: Run the relevant test suite**

Run: `pytest -q`, `git diff --check`, and the three `miniwdl check` commands.

- [x] **Step 3: Inspect the final diff**

Confirm only the image, WDL interfaces, test contract, and affected documentation changed.

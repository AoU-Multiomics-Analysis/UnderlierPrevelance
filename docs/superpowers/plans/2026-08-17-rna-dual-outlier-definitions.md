# RNA Dual Outlier Definitions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Emit RNA prevalence results for z-score-only calls and for the haplo-plus-z-score intersection, with filenames that make the distinction explicit.

**Architecture:** Extend the existing R call helper with an explicit `require_haplo` control. Keep the haplo-only definition, then generate two definitions per configured z cutoff: `z_<cutoff>` from z-score filtering alone and `haplo_z_<cutoff>` from sequential haplo and z-score filtering. Reuse the existing prevalence aggregation and update the WDL artifact validation, tests, and documentation to reflect the expanded output set.

**Tech Stack:** R/Rscript, Python `pytest` tests, WDL/miniwdl validation, Markdown documentation.

## Global Constraints

- The haplo condition remains `expression_logcpm < gene_mean_logcpm - logcpm_drop`.
- Z-score calls remain strict: `expression_zscore < cutoff`.
- Prevalence counts unique `(gene_id, subjectid)` pairs and uses the same tested-subject denominator.
- The existing haplo-only output names remain unchanged.
- For every z cutoff, `z_<cutoff>` means z-only and `haplo_z_<cutoff>` means the haplo-plus-z intersection.
- Do not alter normalization, PCA, residualization, z-score calculation, or cutoff configuration.

---

## File map

- Modify `tests/test_inputs.py`: add dependency-free R behavior coverage and update RNA smoke output expectations.
- Modify `scripts/run_rna_underlier.R`: support z-only and haplo-plus-z call definitions and emit both pairs.
- Modify `workflows/rna_underlier.wdl`: validate two underlier/prevalence artifacts per z cutoff.
- Modify `tests/test_container_scripts.py`: update the expected WDL artifact-count expression.
- Modify `docs/running.md`: document z-only and haplo-plus-z artifact meanings.
- Modify `README.md`: summarize the expanded RNA output set.

### Task 1: Add the failing call-definition and prevalence tests

**Files:**
- Modify: `tests/test_inputs.py`

**Interfaces:**
- Consumes: `source('scripts/run_rna_underlier.R')` with `TOPMED_RNA_UNDERLIER_NO_MAIN=1`.
- Produces: a regression test requiring `call_underliers(expr_z_join, z_cutoff, logcpm_drop, require_haplo)` and checking the resulting prevalence counts.

- [ ] **Step 1: Write the failing test**

Add a test that runs R without loading the container-only packages and constructs one gene with these rows:

```text
subjectid  expression_logcpm  expression_zscore
S1         0.0                0.0
S2        -4.0               -4.0
S3        -1.5               -3.5
S4        -1.5               -3.0
S5        -1.5                0.0
```

With `logcpm_drop = 1`, only `S2` passes haplo; `S2` and `S3` pass z `< -3`; `S4` is excluded because the cutoff is strict; and only `S2` passes the intersection. Assert those subject sets using:

```r
z_only <- call_underliers(expr_z_join, -3, 1, require_haplo = FALSE)
haplo_z <- call_underliers(expr_z_join, -3, 1, require_haplo = TRUE)
stopifnot(identical(z_only$subjectid, c("S2", "S3")))
stopifnot(identical(haplo_z$subjectid, "S2"))
```

Also pass each result to `prevalence_by_gene()` with one `G1` metadata row and `n_subjects_tested = 5`; assert counts of `2` for z-only and `1` for the intersection.

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `python3 -m pytest tests/test_inputs.py -k 'z_only or dual or call_underliers' -q`

Expected: FAIL because the current `call_underliers()` does not accept `require_haplo`.

- [ ] **Step 3: Commit the failing test**

```bash
git add tests/test_inputs.py
git commit -m "test: specify dual RNA outlier definitions"
```

### Task 2: Implement the two R call definitions

**Files:**
- Modify: `scripts/run_rna_underlier.R:384-395,509-523`
- Test: `tests/test_inputs.py`

**Interfaces:**
- Consumes: `expr_z_join`, configured z cutoffs, and `logcpm_drop`.
- Produces: `call_underliers(expr_z_join, z_cutoff, logcpm_drop, require_haplo = TRUE)` and output tags `haplo`, `z_<cutoff>`, and `haplo_z_<cutoff>`.

- [ ] **Step 1: Add the minimal helper behavior**

Change the helper signature to include `require_haplo = TRUE`. Always add `gene_mean_logcpm` to the working data, apply the haplo predicate only when `require_haplo` is true, then apply the z predicate when `z_cutoff` is non-NULL:

```r
call_underliers <- function(expr_z_join, z_cutoff, logcpm_drop, require_haplo = TRUE) {
  gene_means <- stats::ave(expr_z_join$expression_logcpm, expr_z_join$gene_id, FUN = mean)
  expr_z_join$gene_mean_logcpm <- gene_means
  underliers <- expr_z_join
  if (require_haplo) {
    underliers <- underliers[
      underliers$expression_logcpm < underliers$gene_mean_logcpm - logcpm_drop,
      , drop = FALSE
    ]
  }
  if (!is.null(z_cutoff)) {
    underliers <- underliers[underliers$expression_zscore < z_cutoff, , drop = FALSE]
  }
  underliers
}
```

- [ ] **Step 2: Run the focused test to verify it passes**

Run: `python3 -m pytest tests/test_inputs.py -k 'z_only or dual or call_underliers' -q`

Expected: PASS for the z-only, intersection, strict-cutoff, and prevalence assertions.

- [ ] **Step 3: Generate both definitions for every cutoff**

Replace the current one-definition-per-cutoff construction with explicit entries:

```r
definitions <- list(list(tag = "haplo", z_cutoff = NULL, require_haplo = TRUE))
for (z_cutoff in options$z_cutoffs) {
  tag <- strict_output_tag(z_cutoff)
  definitions <- c(
    definitions,
    list(
      list(tag = tag, z_cutoff = z_cutoff, require_haplo = FALSE),
      list(tag = paste0("haplo_", tag), z_cutoff = z_cutoff, require_haplo = TRUE)
    )
  )
}
for (definition in definitions) {
  underliers <- call_underliers(
    expr_z_join,
    definition$z_cutoff,
    options$logcpm_drop,
    definition$require_haplo
  )
  # existing write_tsv and prevalence_by_gene calls remain
}
```

- [ ] **Step 4: Run the focused test again**

Run: `python3 -m pytest tests/test_inputs.py -k 'z_only or dual or call_underliers' -q`

Expected: PASS, with no change to haplo-only behavior.

- [ ] **Step 5: Commit the R implementation**

```bash
git add scripts/run_rna_underlier.R tests/test_inputs.py
git commit -m "feat: emit z-only and haplo-z RNA calls"
```

### Task 3: Update WDL artifact validation and smoke expectations

**Files:**
- Modify: `workflows/rna_underlier.wdl:79-89`
- Modify: `tests/test_container_scripts.py:90`
- Modify: `tests/test_inputs.py:405-413`

**Interfaces:**
- Consumes: the R output naming contract from Task 2.
- Produces: WDL validation expecting one haplo pair plus two pairs per z cutoff.

- [ ] **Step 1: Update the WDL expected artifact expression**

Change:

```bash
expected_artifacts=$((1 + ~{length(z_cutoffs)}))
```

to:

```bash
expected_artifacts=$((1 + 2 * ~{length(z_cutoffs)}))
```

- [ ] **Step 2: Update the dependency-free WDL text test**

Change its expected string to `expected_artifacts=$((1 + 2 * ~{length(z_cutoffs)}))`.

- [ ] **Step 3: Update the RNA smoke expected outputs**

For the smoke cutoff `-3`, require both:

```text
underliers_z_-3.tsv.gz
underliers_haplo_z_-3.tsv.gz
rna_outlier_prevalence_per_gene_z_-3.tsv
rna_outlier_prevalence_per_gene_haplo_z_-3.tsv
```

Keep the existing haplo-only filenames in the expected list.

- [ ] **Step 4: Run the focused WDL/input tests**

Run: `python3 -m pytest tests/test_container_scripts.py tests/test_inputs.py -k 'workflow or artifact or smoke' -q`

Expected: dependency-free WDL assertions pass; the RNA smoke test may require the container R dependencies.

- [ ] **Step 5: Commit the interface updates**

```bash
git add workflows/rna_underlier.wdl tests/test_container_scripts.py tests/test_inputs.py
git commit -m "test: cover dual RNA output artifacts"
```

### Task 4: Update user-facing documentation

**Files:**
- Modify: `docs/running.md:84,90-97`
- Modify: `README.md:50-52`

**Interfaces:**
- Consumes: final filenames and semantics from Tasks 2–3.
- Produces: documentation distinguishing z-only from haplo-plus-z prevalence.

- [ ] **Step 1: Update the analysis description**

State that each configured cutoff emits z-only calls and an intersection call requiring both haplo and z-score criteria.

- [ ] **Step 2: Update the artifact table**

Document:

```text
underliers_z_<cutoff>                  z-score-only calls
underliers_haplo_z_<cutoff>            haplo plus z-score calls
rna_outlier_prevalence_per_gene_z_*    z-score-only prevalence
rna_outlier_prevalence_per_gene_haplo_z_*  intersection prevalence
```

Retain the existing statement that prevalence counts each subject once per gene.

- [ ] **Step 3: Update the README summary**

Replace the “two underlier call tables, two per-gene prevalence tables” wording with a summary covering haplo-only, z-only, and haplo-plus-z outputs.

- [ ] **Step 4: Review documentation for stale meanings**

Run: `rg -n 'underliers_z_|prevalence_per_gene_z_|both the haplo|two underlier|two per-gene' README.md docs/running.md workflows tests`

Expected: no remaining statement says `z_<cutoff>` itself is the haplo intersection; intersection references use `haplo_z_<cutoff>`.

- [ ] **Step 5: Commit the documentation**

```bash
git add README.md docs/running.md
git commit -m "docs: describe dual RNA outlier outputs"
```

### Task 5: Run full verification

**Files:**
- Verify: all modified files and git diff

- [ ] **Step 1: Run dependency-free tests**

Run: `python3 -m pytest tests/test_inputs.py tests/test_container_scripts.py -q`

Expected: all tests that do not require the RNA container dependencies pass.

- [ ] **Step 2: Run the R function-level assertion directly if pytest is unavailable**

Run the same `Rscript` source/assertion used by Task 1 with `TOPMED_RNA_UNDERLIER_NO_MAIN=1`.

Expected: z-only returns both z `< -3` rows, the intersection returns only haplo-positive z `< -3` rows, and prevalence counts are `2` versus `1`.

- [ ] **Step 3: Validate the WDL syntax**

Run: `miniwdl check workflows/rna_underlier.wdl`

Expected: exit 0.

- [ ] **Step 4: Inspect the final diff and status**

Run: `git diff HEAD~4 --check && git status --short`

Expected: no whitespace errors; the pre-existing untracked `scripts/compute_full_int_pca.R` remains untouched.

- [ ] **Step 5: Commit any final verification-only corrections**

Only if verification identifies a documentation or assertion mismatch, make the smallest correction, rerun the relevant check, and commit it with:

```bash
git add <corrected-files>
git commit -m "fix: align RNA dual-output verification"
```

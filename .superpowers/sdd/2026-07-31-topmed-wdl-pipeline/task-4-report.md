# Task 4 report: container definitions

## Status

Implemented the Task 4 container/runtime definitions without changing analytical
scripts. Docker verification is incomplete because it was stopped at the user's
request after a source-package compatibility failure was diagnosed and addressed
in the Dockerfile.

## Changed files

- `envs/Dockerfile`
  - Builds a pinned q2 runtime with Bioconda (`bcftools 1.20`, Python 3.12.7,
    pip 24.2) and installs the pinned q2 Python requirements.
  - Uses `rocker/r-ver:4.4.2` for the final image and installs the full
    Bioconductor 3.20 RNA dependency set: `edgeR`, `limma`, `corral`,
    `PCAtools`, `scran`, `scuttle`, `WGCNA`, `rtracklayer`, `vroom`,
    `tidyverse`, and `magrittr`.
  - Pins `Rcpp 1.0.14` before the remaining RNA packages to avoid the observed
    incompatibility between the frozen snapshot's `Rcpp 1.0.13` and the Rocker
    R headers.
- `envs/bioconda-environment.yml`
  - Declares the q2 Bioconda/conda-forge environment.
- `envs/requirements-q2.txt`
  - Pins `numpy 2.1.3` and `pandas 2.2.3`.

## Validation

Passed local checks:

```text
git diff --check
container manifest structural checks: PASS
```

Docker checks attempted:

```text
docker build -f envs/Dockerfile -t underlier-prevalence:test .
```

- Baseline: failed as expected because `envs/Dockerfile` was empty.
- First implementation build: q2 environment resolution succeeded, but a
  combined R/Bioconda environment could not resolve the requested RNA package
  versions and `bcftools 1.20` together. The Dockerfile was revised to use
  Bioconda only for q2 and Bioconductor source installation for RNA.
- Second implementation build: successfully created the q2 environment and
  resolved Bioconductor 3.20 RNA packages, then failed compiling the frozen
  snapshot's `Rcpp 1.0.13` against the Rocker R headers. The Dockerfile now
  installs `Rcpp 1.0.14` explicitly before the RNA set.
- The post-fix full rebuild and the required `Rscript`/`bcftools` image probes
  were not run because the user asked to stop all long-running Docker work.

## Concerns / follow-up

- Run the three Task 4 Docker smoke commands after the `Rcpp 1.0.14` update to
  confirm the final image and runtime probes.
- No WDL files exist yet (Task 5 scope), so WDL image references cannot be
  centralized or made overrideable in this task without expanding scope.

## Fix round 1

- Preserved the q2 environment at `/opt/conda/envs/underlier-prevalence` in the
  final image, added its `bin` directory to `PATH`, and set
  `BCFTOOLS_PLUGINS` to its `libexec/bcftools` directory so `+fill-tags` is
  discoverable.
- Replaced the Bookworm-specific Posit URL with the distribution-neutral frozen
  Posit CRAN source snapshot URL.
- Added Docker build-time probes covering every required R package,
  `PCAtools::chooseGavishDonoho`, Python `numpy`/`pandas`, `bcftools --version`,
  and `bcftools +fill-tags -h`.
- Added Docker CI path triggers for both q2 dependency manifests.

### Validation status

- The fix-round structural check passed: preserved prefix and plugin path,
  distribution-neutral snapshot URL, required package probes, and CI triggers
  are all present.
- A full Docker build was started. It completed the Bioconda q2 environment,
  copied `/opt/conda` at the preserved prefix, installed system dependencies,
  resolved the neutral Posit snapshot, and compiled `Rcpp 1.0.14` past the
  prior failure. It was manually canceled during the longer remaining RNA
  package compilation at the user's request.
- Consequently, the final image was not emitted and the post-build
  `Rscript`, Python, `bcftools --version`, and `bcftools +fill-tags -h` probes
  remain unexecuted. The Dockerfile now runs all four categories during a
  successful image build.

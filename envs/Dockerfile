FROM mambaorg/micromamba:2.0.5 AS q2-runtime

COPY --chown=${MAMBA_USER}:${MAMBA_USER} envs/bioconda-environment.yml /tmp/bioconda-environment.yml
COPY --chown=${MAMBA_USER}:${MAMBA_USER} envs/requirements-q2.txt /tmp/requirements-q2.txt

RUN micromamba create --yes --name underlier-prevalence --file /tmp/bioconda-environment.yml \
    && micromamba run --name underlier-prevalence python -m pip install --no-cache-dir -r /tmp/requirements-q2.txt \
    && micromamba clean --all --yes

FROM rocker/r-ver:4.4.2

ARG BIOCONDUCTOR_VERSION=3.20
ARG CRAN_SNAPSHOT=2024-11-01
ENV PATH=/opt/conda/envs/underlier-prevalence/bin:${PATH}
ENV BCFTOOLS_PLUGINS=/opt/conda/envs/underlier-prevalence/libexec/bcftools
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

COPY --from=q2-runtime /opt/conda /opt/conda

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        gfortran \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgit2-dev \
        libglpk-dev \
        libgsl-dev \
        libharfbuzz-dev \
        libhdf5-dev \
        libjpeg-dev \
        liblzma-dev \
        libpcre2-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libudunits2-dev \
        libxml2-dev \
        libxt-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(repos = c(CRAN = sprintf("https://packagemanager.posit.co/cran/%s", Sys.getenv("CRAN_SNAPSHOT")))); install.packages("BiocManager"); install.packages("https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.14.tar.gz", repos = NULL, type = "source"); BiocManager::install(version = Sys.getenv("BIOCONDUCTOR_VERSION"), ask = FALSE, update = FALSE); BiocManager::install(c("edgeR", "limma", "corral", "PCAtools", "scran", "scuttle", "WGCNA", "rtracklayer", "vroom", "tidyverse", "magrittr"), ask = FALSE, update = FALSE)'

RUN Rscript -e 'packages <- c("edgeR", "limma", "corral", "PCAtools", "scran", "scuttle", "WGCNA", "rtracklayer", "vroom", "tidyverse", "magrittr"); missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", ")); stopifnot(is.function(PCAtools::chooseGavishDonoho))' \
    && python -c 'import numpy, pandas; print(numpy.__version__, pandas.__version__)' \
    && bcftools --version \
    && bcftools +fill-tags -h

COPY scripts/ /opt/underlier-prevalence/scripts/

WORKDIR /work

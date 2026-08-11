FROM mambaorg/micromamba:2.0.5 AS runtime-builder

COPY --chown=${MAMBA_USER}:${MAMBA_USER} envs/bioconda-environment.yml /tmp/bioconda-environment.yml
COPY --chown=${MAMBA_USER}:${MAMBA_USER} envs/requirements-q2.txt /tmp/requirements-q2.txt

RUN micromamba create --yes --name underlier-prevalence --file /tmp/bioconda-environment.yml \
    && micromamba run --name underlier-prevalence python -m pip install --no-cache-dir -r /tmp/requirements-q2.txt \
    && micromamba clean --all --yes

FROM mambaorg/micromamba:2.0.5

ENV PATH=/opt/conda/envs/underlier-prevalence/bin:${PATH}
ENV BCFTOOLS_PLUGINS=/opt/conda/envs/underlier-prevalence/libexec/bcftools
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

COPY --from=runtime-builder /opt/conda /opt/conda

RUN Rscript -e 'packages <- c("edgeR", "limma", "corral", "PCAtools", "scran", "scuttle", "WGCNA", "rtracklayer", "vroom", "tidyverse", "magrittr"); missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", ")); stopifnot(is.function(PCAtools::chooseGavishDonoho))' \
    && python -c 'import numpy, pandas; print(numpy.__version__, pandas.__version__)' \
    && bcftools --version \
    && bcftools +fill-tags -h

COPY scripts/ /opt/underlier-prevalence/scripts/

USER root
RUN mkdir -p /work \
    && chown -R ${MAMBA_USER}:${MAMBA_USER} /work /opt/underlier-prevalence
USER ${MAMBA_USER}

WORKDIR /work

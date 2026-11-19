FROM mambaorg/micromamba:latest

# ps is used by Nextflow inside containers
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# pinned to the versions used to produce the results in this repository
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      fastp=1.4.0 \
      bwa-mem2=2.3 \
      samtools=1.25 \
      bcftools=1.25 \
      snpeff=5.4.0c \
    && micromamba clean -a -y

# make the tools visible without relying on the image entrypoint (Nextflow overrides it)
ENV PATH=/opt/conda/bin:$PATH

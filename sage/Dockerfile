FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN mkdir /opt/sage
RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*
ADD target/sage-$VERSION-jar-with-dependencies.jar /opt/sage/sage.jar

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-genomicranges \
    bioconductor-variantannotation \
    font-ttf-dejavu \
    'openjdk >=8' \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    xorg-libxt \
    zlib && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"


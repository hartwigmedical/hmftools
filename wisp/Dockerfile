FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN mkdir -p /opt/wisp
RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

ADD target/wisp-$VERSION-jar-with-dependencies.jar /opt/wisp/wisp.jar

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'openjdk >=8' \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    r-gtable \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

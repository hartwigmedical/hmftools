FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN mkdir /opt/cobalt
RUN \
  apt-get update && \
  apt-get install -y procps && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*
ADD target/cobalt-$VERSION-jar-with-dependencies.jar /opt/cobalt/cobalt.jar

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    r-dplyr 'openjdk >=8' && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

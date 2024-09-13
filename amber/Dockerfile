FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN mkdir /opt/amber
RUN \
  apt-get update && \
  apt-get install -y procps && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*
ADD target/amber-$VERSION-jar-with-dependencies.jar /opt/amber/amber.jar

USER mambauser

RUN micromamba install -y -n base -c bioconda -c conda-forge r-dplyr 'openjdk >=8' && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"


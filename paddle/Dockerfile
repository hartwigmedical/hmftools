FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN mkdir -p /opt/paddle
RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

ADD target/paddle-$VERSION-jar-with-dependencies.jar /opt/paddle/paddle.jar

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-genomicranges \
    'openjdk >=8' \
    r-dplyr \
    r-devtools \
    r-tidyr && \
  micromamba clean --all --yes

# dndscv does not appear to be available via Conda
RUN /opt/conda/bin/R -e 'library(devtools); install_github("im3sanger/dndscv");'

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

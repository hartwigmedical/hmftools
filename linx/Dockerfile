FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN mkdir /opt/linx
RUN \
  apt-get update && \
  apt-get install -y procps wget && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

ADD target/linx-$VERSION-jar-with-dependencies.jar /opt/linx/linx.jar

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-gviz \
    'circos >=0.69.6' \
    'openjdk >=8' \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-magick \
    r-tidyr \
    xorg-libxtst && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

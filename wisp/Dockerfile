FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'openjdk >=8' \
    r-dplyr \
    r-ggplot2 \
    r-gtable \
    r-tidyr && \
  micromamba clean --all --yes

USER root
RUN mkdir -p /opt/wisp
ADD target/wisp-$VERSION-jar-with-dependencies.jar /opt/wisp/wisp.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

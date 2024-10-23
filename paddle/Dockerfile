FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-genomicranges \
    r-devtools \
    r-dplyr \
    r-tidyr && \
  micromamba clean --all --yes

# dndscv does not appear to be available via Conda
RUN /opt/conda/bin/R -e 'library(devtools); install_github("im3sanger/dndscv");'

USER root
RUN mkdir -p /opt/paddle
ADD target/paddle-$VERSION-jar-with-dependencies.jar /opt/paddle/paddle.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bcftools \
    bioconductor-variantannotation \
    'circos >=0.69.6' \
    r-cairo \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    xorg-libxtst \
    zlib && \
  micromamba clean --all --yes

USER root
RUN mkdir -p /opt/purple
ADD target/purple-$VERSION-jar-with-dependencies.jar /opt/purple/purple.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

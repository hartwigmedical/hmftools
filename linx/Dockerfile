FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-gviz \
    'circos >=0.69.6' \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-magick \
    r-tidyr \
    xorg-libxtst && \
  micromamba clean --all --yes

USER root
RUN mkdir /opt/linx
ADD target/linx-$VERSION-jar-with-dependencies.jar /opt/linx/linx.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

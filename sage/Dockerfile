FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-genomicranges \
    bioconductor-variantannotation \
    font-ttf-dejavu \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    xorg-libxt \
    zlib && \
  micromamba clean --all --yes

USER root
RUN mkdir /opt/sage
ADD target/sage-$VERSION-jar-with-dependencies.jar /opt/sage/sage.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"


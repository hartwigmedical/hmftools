FROM circos

ARG VERSION

USER mambauser
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bioconductor-gviz \
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

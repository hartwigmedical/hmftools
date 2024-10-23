FROM circos

ARG VERSION

USER mambauser
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    bcftools \
    bioconductor-variantannotation \
    r-cairo \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-tidyr \
    xorg-libxtst \
    zlib && \ 
  micromamba clean --all --yes --force-pkgs-dirs

USER root
RUN mkdir -p /opt/purple
ADD target/purple-$VERSION-jar-with-dependencies.jar /opt/purple/purple.jar

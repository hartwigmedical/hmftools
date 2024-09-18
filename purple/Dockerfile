FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER root
RUN apt update && \
  apt install -y libfindbin-libs-perl libconfig-general-perl libfont-ttf-perl \
    libclone-perl libgd-perl libmath-round-perl libmath-bezier-perl \
    liblist-moreutils-perl libparams-validate-perl libmath-vecstat-perl \
    libreadonly-perl libset-intspan-perl libsvg-perl libregexp-common-perl \
    libtext-format-perl libstatistics-basic-perl circos && \
  apt clean && rm -rf /var/lib/apt/lists/*
USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' \
  && micromamba clean --all --yes
#RUN micromamba install -y -n base -c bioconda -c conda-forge 'circos >=0.69.6' perl-gd && micromamba clean --all --yes
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
  micromamba clean --all --yes

USER root
RUN mkdir -p /opt/purple
ADD target/purple-$VERSION-jar-with-dependencies.jar /opt/purple/purple.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

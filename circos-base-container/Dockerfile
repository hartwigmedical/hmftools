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
RUN mkdir -p /usr/bin/etc
ADD etc/* /usr/bin/etc/
ADD data.tar.gz /usr/bin/
ADD fonts.tar.gz /usr/bin/

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && \
  micromamba clean --all --yes --force-pkgs-dirs

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

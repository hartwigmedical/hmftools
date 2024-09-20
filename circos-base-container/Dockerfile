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
RUN mkdir -p /usr/bin/data /usr/bin/etc /usr/bin/fonts
ADD data/* /usr/bin/data/
ADD etc/* /usr/bin/etc/
ADD fonts/* /usr/bin/fonts/

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && \
  micromamba clean --all --yes --force-pkgs-dirs

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

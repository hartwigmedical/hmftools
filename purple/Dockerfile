FROM mambaorg/micromamba:0.24.0

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    "openjdk>=8" \
    xorg-libxtst \
    r-tidyr \
    r-dplyr \
    r-ggplot2 \
    r-cowplot \
    r-magick \
    bioconductor-gviz \
    "circos>=0.69.9" \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

# Fix missing GD perl module error (despite GD already being installed with circos)
RUN apt update && apt install -y libconfig-general-perl
RUN perl -MCPAN -e 'install GD'

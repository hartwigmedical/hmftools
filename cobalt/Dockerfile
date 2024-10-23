FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser
RUN micromamba install -y -n base -c bioconda -c 'openjdk >=8' && micromamba clean --all --yes
RUN micromamba install -y -n base -c bioconda -c conda-forge r-dplyr && micromamba clean --all --yes

USER root
RUN mkdir /opt/cobalt
ADD target/cobalt-$VERSION-jar-with-dependencies.jar /opt/cobalt/cobalt.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

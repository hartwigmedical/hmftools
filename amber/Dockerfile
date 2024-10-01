FROM mambaorg/micromamba:0.24.0

ARG VERSION

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN micromamba install -y -n base -c bioconda -c conda-forge r-dplyr && micromamba clean --all --yes

USER root
RUN mkdir /opt/amber
ADD target/amber-$VERSION-jar-with-dependencies.jar /opt/amber/amber.jar

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"


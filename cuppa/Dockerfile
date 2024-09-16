FROM mambaorg/micromamba:0.24.0

ARG VERSION
USER root

RUN \
  apt-get update && \
  apt-get install -y python3 python3-pip && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

ADD target/cuppa-$VERSION-jar-with-dependencies.jar /opt/cuppa/cuppa.jar

USER mambauser

RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    'openjdk >=8' \
    r-ggplot2 \
    r-reshape2 \
    r-patchwork && \
  micromamba clean --all --yes

USER root

RUN cd /opt/cuppa && /opt/conda/bin/jar xvf /opt/cuppa/cuppa.jar pycuppa
RUN pip3 install --upgrade pip && pip3 install /opt/cuppa/pycuppa

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

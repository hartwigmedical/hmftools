FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG JAR=/usr/share/java/cuppa_v${VERSION}.jar

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN \
  micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    unzip \
    r-ggplot2 \
    r-patchwork \
    r-reshape2 \
    r-ggh4x && \
  micromamba clean --all --yes

USER root
ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

ADD target/cuppa-$VERSION-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/cuppa
RUN chmod +x /usr/local/bin/cuppa

RUN apt update && apt install -y python3 python3-pip && apt clean && rm -rf /var/lib/apt/lists/*
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN unzip -n ${JAR} -d /tmp/cuppa_jar
RUN pip3 install --upgrade pip && pip3 install /tmp/cuppa_jar/pycuppa/
RUN rm -r /tmp/cuppa_jar

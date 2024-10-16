FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG TOOL_NAME=cuppa
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar

USER root

ADD target/${TOOL_NAME}-$VERSION-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/${TOOL_NAME}
RUN chmod +x /usr/local/bin/${TOOL_NAME}

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    "openjdk >=8" \
    procps-ng \
    unzip \
    python=3.11 \
    r-ggplot2 \
    r-patchwork \
    r-reshape2 \
    r-ggh4x \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

RUN unzip -n ${JAR} -d /tmp/${TOOL_NAME}_jar
RUN pip3 install --upgrade pip && pip3 install /tmp/${TOOL_NAME}_jar/pycuppa/
RUN rm -r /tmp/${TOOL_NAME}_jar

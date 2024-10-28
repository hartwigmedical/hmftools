FROM mambaorg/micromamba:2.0.2

ARG VERSION
ARG TOOL_NAME=purple
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar
ARG RUN_SCRIPT=/usr/local/bin/${TOOL_NAME}

USER root

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh ${RUN_SCRIPT}
RUN chmod +x ${RUN_SCRIPT}

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    openjdk>=8 \
    xorg-libxtst \
    r-tidyr \
    r-dplyr \
    r-ggplot2 \
    r-cowplot \
    r-cairo \
    bioconductor-variantannotation \
    circos>=0.69.6 \
    perl-gd>=2.76 \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG TOOL_NAME=linx
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar

USER root

ADD target/${TOOL_NAME}-$VERSION-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/${TOOL_NAME}
RUN chmod +x /usr/local/bin/${TOOL_NAME}

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    openjdk>=8 \
    xorg-libxtst \
    r-tidyr \
    r-dplyr \
    r-ggplot2 \
    r-cowplot \
    r-magick \
    bioconductor-gviz \
    circos>=0.69.9 \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

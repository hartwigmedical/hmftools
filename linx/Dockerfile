FROM circos

ARG VERSION
ARG TOOL_NAME=linx
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar

USER root

ADD target/${TOOL_NAME}-$VERSION-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/${TOOL_NAME}
RUN chmod +x /usr/local/bin/${TOOL_NAME}

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    "openjdk >=8" \
    procps-ng \
    r-cowplot \
    r-dplyr \
    r-ggplot2 \
    r-magick \
    r-tidyr \
    bioconductor-gviz \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

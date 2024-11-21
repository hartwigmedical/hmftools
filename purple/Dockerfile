FROM mambaorg/micromamba:0.24.0

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
    "openjdk>=8,<=21" \
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

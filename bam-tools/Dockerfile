FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG TOOL_NAME=bam-tools
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar

ARG TOOL_ALIAS=bamtools
ARG RUN_SCRIPT=/usr/local/bin/${TOOL_ALIAS}

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    "openjdk >=8,<=21" \
    "sambamba>=1.0.1" \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh ${RUN_SCRIPT}
RUN chmod +x ${RUN_SCRIPT}

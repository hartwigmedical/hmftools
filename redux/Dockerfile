FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG TOOL_NAME=redux
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar
ARG RUN_SCRIPT=/usr/local/bin/${TOOL_NAME}

USER root

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh ${RUN_SCRIPT}
RUN chmod +x ${RUN_SCRIPT}

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    "openjdk=11" \
    "sambamba>=1.0.1" \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

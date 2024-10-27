FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG TOOL_NAME=chord
ARG JAR=/usr/share/java/${TOOL_NAME}_v${VERSION}.jar
ARG RUN_SCRIPT=/usr/local/bin/${TOOL_NAME}

USER root

ADD target/${TOOL_NAME}-${VERSION}-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh ${RUN_SCRIPT}
RUN chmod +x ${RUN_SCRIPT}

RUN micromamba install -y -n base -c bioconda -c conda-forge \
    "openjdk>=8" \
    procps-ng \
    unzip \
    r-base \
    r-randomforest \
    r-stringr \
    bioconductor-bsgenome \
    bioconductor-bsgenome.hsapiens.ucsc.hg19 \
    bioconductor-bsgenome.hsapiens.ucsc.hg38 \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

RUN unzip -n ${JAR} -d /tmp/${TOOL_NAME}_jar
RUN Rscript -e "install.packages(c('${TOOL_NAME}_jar/mutSigExtractor', '${TOOL_NAME}_jar/CHORD'), repos=NULL, type='source')"
RUN rm -r /tmp/${TOOL_NAME}_jar

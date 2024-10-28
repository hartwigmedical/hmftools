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
    "openjdk>=8" \
    xorg-libxtst \
    r-tidyr \
    r-dplyr \
    r-ggplot2 \
    r-cowplot \
    r-magick \
    bioconductor-gviz \
    wget \
    && micromamba clean --all --yes

RUN apt update && apt install -y \
    libconfig-general-perl \
    libclone-perl \
    libfindbin-libs-perl \
	libfont-ttf-perl \
    libgd-perl \
    liblist-moreutils-perl \
    libmath-round-perl \
    libmath-bezier-perl \
    libmath-vecstat-perl \
    libparams-validate-perl \
    libreadonly-perl \
    libregexp-common-perl \
    libset-intspan-perl \
    libstatistics-basic-perl \
    libsvg-perl \
    libtext-format-perl \
    circos \
    && apt clean && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

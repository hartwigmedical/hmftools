FROM mambaorg/micromamba:0.24.0

ARG VERSION
ARG JAR=/usr/share/java/amber_v${VERSION}.jar

USER mambauser
RUN micromamba install -y -n base -c bioconda -c conda-forge 'openjdk >=8' && micromamba clean --all --yes
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    procps-ng \
    r-dplyr  \
    bioconductor-copynumber \
    && micromamba clean --all --yes

USER root
ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"

ADD target/amber-$VERSION-jar-with-dependencies.jar ${JAR}
ADD target/entrypoint.sh /usr/local/bin/amber
RUN chmod +x /usr/local/bin/amber

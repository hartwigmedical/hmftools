FROM openjdk:17-slim

RUN mkdir /usr/share/hartwig
ARG VERSION
RUN apt update && apt install -y curl gnupg2 
RUN curl https://europe-west4-apt.pkg.dev/doc/repo-signing-key.gpg | apt-key add - && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN echo 'deb http://packages.cloud.google.com/apt apt-transport-artifact-registry-stable main' | tee -a /etc/apt/sources.list.d/artifact-registry.list
RUN apt update && apt install -y apt-transport-artifact-registry
RUN echo 'deb ar+https://europe-west4-apt.pkg.dev/projects/hmf-build hmf-apt main' | tee -a  /etc/apt/sources.list.d/artifact-registry.list
RUN apt update && apt install -y hmf-r-base hmf-circos hmf-r-ggplot2 hmf-r-dplyr hmf-r-bioconductor-variantannotation

ADD target/purple-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/purple.jar
#&& apt install -y circos liblist-moreutils-perl

#ENTRYPOINT ["java", "-jar", "/usr/share/hartwig/purple.jar"]
ENV PATH="$PATH:/opt/tools/R/bin"

ENTRYPOINT ["sleep", "3600"]

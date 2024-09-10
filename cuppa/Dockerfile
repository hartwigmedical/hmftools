FROM openjdk:17-slim

ARG VERSION
RUN mkdir -p /usr/share/hartwig /opt/tools
RUN apt update && apt install -y curl gnupg2 python3 python3-pip
RUN curl https://europe-west4-apt.pkg.dev/doc/repo-signing-key.gpg | apt-key add - && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN echo 'deb http://packages.cloud.google.com/apt apt-transport-artifact-registry-stable main' | tee -a /etc/apt/sources.list.d/artifact-registry.list
RUN apt update && apt install -y apt-transport-artifact-registry
RUN echo 'deb ar+https://europe-west4-apt.pkg.dev/projects/hmf-build hmf-apt main' | tee -a  /etc/apt/sources.list.d/artifact-registry.list
RUN apt update && apt install -y hmf-r-base 

ADD target/cuppa-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/cuppa.jar
RUN cd /opt/tools && jar xvf /usr/share/hartwig/cuppa.jar pycuppa 
RUN pip3 install --upgrade pip && pip3 install /opt/tools/pycuppa

ENV PATH="$PATH:/opt/tools/R/bin"
ENTRYPOINT ["java", "-jar", "/usr/share/hartwig/cuppa.jar"]

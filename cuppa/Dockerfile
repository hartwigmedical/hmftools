FROM openjdk:17-slim

ARG VERSION
RUN mkdir -p /usr/share/hartwig /opt/tools
RUN apt update && apt install -y python3 python3-pip

ADD target/cuppa-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/cuppa.jar
RUN cd /opt/tools && jar xvf /usr/share/hartwig/cuppa.jar pycuppa 
RUN pip3 install --upgrade pip && pip3 install /opt/tools/pycuppa

ENTRYPOINT ["java", "-cp", "/usr/share/hartwig/cuppa.jar", "com.hartwig.hmftools.cup.prep.CuppaDataPrep"]

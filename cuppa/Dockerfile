FROM openjdk:17-slim

ARG VERSION
RUN mkdir -p /usr/share/hartwig /opt/tools
ADD target/cuppa-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/cuppa.jar
RUN cd /opt/tools && jar xvf /usr/share/hartwig/cuppa.jar pycuppa 

ENTRYPOINT ["java", "-cp", "/usr/share/hartwig/cuppa.jar", "com.hartwig.hmftools.cup.prep.CuppaDataPrep"]

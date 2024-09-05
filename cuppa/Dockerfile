FROM openjdk:17-slim

ARG VERSION
RUN mkdir -p /usr/share/hartwig /opt/tools/pycuppa
ADD target/cuppa-$VERSION-jar-with-dependencies.jar /usr/share/hartwig/cuppa.jar
RUN jar xf /usr/share/hartwig/cuppa.jar pycuppa/* -d /opt/tools/pycuppa

ENTRYPOINT ["java", "-cp", "/usr/share/hartwig/cuppa.jar", "com.hartwig.hmftools.cup.prep.CuppaDataPrep"]

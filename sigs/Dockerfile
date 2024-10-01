FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/sigs
ADD target/sigs-$VERSION-jar-with-dependencies.jar /opt/sigs/sigs.jar

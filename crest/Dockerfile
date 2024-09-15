FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/crest
ADD target/crest-$VERSION-jar-with-dependencies.jar /opt/crest/crest.jar

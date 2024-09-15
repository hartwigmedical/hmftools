FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/gene-utils
ADD target/gene-utils-$VERSION-jar-with-dependencies.jar /opt/gene-utils/gene-utils.jar

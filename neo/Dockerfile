FROM openjdk:17-alpine

ARG VERSION
RUN mkdir /opt/neo
ADD target/neo-$VERSION-jar-with-dependencies.jar /opt/neo/neo.jar
